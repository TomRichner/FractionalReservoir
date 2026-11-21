function N = inm_run_network(seed, plan)
%INM_RUN_NETWORK One modular network of the network_mode stage: its baseline and stimulation blocks, and the full Jacobian's network mode along them and around discharges.
%
%   N = inm_run_network(seed, plan)
%
% plan from ied_network_mode_stage_plan (it documents what is computed).
% Fields of N:
%   seed, groups (n x 1 subnet), G, tau_d, seconds
%   block(b), b = 1 baseline, 2 stimulation: bundle_block (its number in
%     the bundle), amplitude, noise_seed, events (time, groups: the
%     detector's over the analysed part), n_events, repro (largest
%     difference of the mean E potential from the bundle's block file;
%     NaN without one), t (s, time in the run; the analysed part starts at
%     sim.settle), net (complex, 1/s), net_coh (coherence over the ten
%     populations), share (ns x G), n_net_unstable, abscissa, sim_s, eig_s
%   event: lags, t_event, t_ctrl, and for 'ev' and 'ctrl' the net, net_coh
%     and share at each centre (rows) and lag (columns)
%   summary: one table row (run_ied_network_mode_stage stacks them)
%
% Ported from the pilot's soied_run_network and soied_event_full
% (scripts/explorations/stim_outlier, 2026-11-21); changed: the bundle's
% own block numbering and whole block (so the run is the bundle's block,
% checked), the full Jacobian only (no frozen-gain series), and one
% simulation per block for both the series and the discharges.
%
% See also: ied_network_mode_stage_plan, run_ied_network_mode_stage

t_all = tic;
sim = plan.sim;
base = ied_base_config(); base.fs = sim.fs;
N = struct('seed', seed);
for b = 1:2
    t0 = tic;
    k_bundle = plan.blocks(b);
    A = sim.amplitude * (k_bundle > sim.n_no);
    cfg = base; cfg.seed = seed;
    cfg.noise_seed = base.noise_seed + (seed - base.seed) + 1000 * k_bundle;
    cfg.duration = sim.settle + sim.block_len;
    ic = struct('dc_levels', A, 'hold_dur', cfg.duration, 'ramp_dur', 0, ...
        'noise_intensity', 0, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
    [m, groups] = ied_build_network(cfg, ic, 'lya_method', 'none', ...
        'T_range', [-cfg.warmup cfg.duration], 'T_plot', [0 cfg.duration], 'store_full_state', true);
    m.run();
    sim_s = toc(t0);
    params = m.get_params();
    ws = ospec_weight_stats(params);
    if b == 1
        [U10, G] = population_basis(params, groups);
        N.groups = groups(:); N.G = G; N.tau_d = params.tau_d;
    end
    % The detector, exactly as ied_stim_blocks' run_block calls it.
    p = m.plot_data; t = p.t(:)';
    x = [p.x.E; p.x.I]; r = [p.r.E; p.r.I];
    sel = t >= sim.settle & t < cfg.duration;
    det = ied_detect_events(t(sel), x(:, sel), r(:, sel), groups, struct());
    ev = det.events;
    mean_x_E = mean(x(m.type_indices{1}, sel), 1);
    Bk = struct('bundle_block', k_bundle, 'amplitude', A, 'noise_seed', cfg.noise_seed, ...
        'events', struct('time', [ev.time], 'groups', {{ev.groups}}), 'n_events', numel(ev), ...
        'repro', repro(plan.run_dir, seed, k_bundle, mean_x_E), 'sim_s', sim_s);
    % The series on the first window_s s of the analysed part.
    tf = m.t_out(:);
    ts = (sim.settle:plan.full_dt:sim.settle + plan.window_s - 1e-9)';
    t1 = tic;
    X = modes_at(m.S_out, tf, ts, sim.fs, params, ws, plan.x_thr, U10, groups);
    Bk.t = ts; Bk.net = X.net; Bk.net_coh = X.net_coh; Bk.share = X.share;
    Bk.n_net_unstable = X.n_net_unstable; Bk.abscissa = X.abscissa;
    % Around discharges and control times, baseline block only.
    if b == 1
        te = [ev.time];
        use = te(1:min(plan.n_events, numel(te)));
        use = use(use + plan.lags(1) >= sim.settle & use + plan.lags(end) < cfg.duration);
        rng(seed, 'twister');
        ctrl = [];
        tries = 0;
        while numel(ctrl) < numel(use) && tries < 10000
            tries = tries + 1;
            c = sim.settle + plan.ctrl_gap + rand * (sim.block_len - 2 * plan.ctrl_gap);
            if all(abs(te - c) > plan.ctrl_gap); ctrl(end + 1) = c; end %#ok<AGROW>
        end
        E = struct('lags', plan.lags, 't_event', use, 't_ctrl', ctrl);
        kinds = {'ev', 'ctrl'}; centres = {use, ctrl};
        for kk = 1:2
            c = centres{kk}; nl = numel(plan.lags);
            Y = struct('net', complex(nan(numel(c), nl)), 'net_coh', nan(numel(c), nl), 'share', nan(numel(c), nl, G));
            for e = 1:numel(c)
                Z = modes_at(m.S_out, tf, c(e) + plan.lags(:), sim.fs, params, ws, plan.x_thr, U10, groups);
                Y.net(e, :) = Z.net.'; Y.net_coh(e, :) = Z.net_coh.'; Y.share(e, :, :) = Z.share;
            end
            E.(kinds{kk}) = Y;
        end
        N.event = E;
    end
    Bk.eig_s = toc(t1);
    N.block(b) = Bk;
    clear m
end
N.seconds = toc(t_all);
N.summary = summary_row(N, plan);
end

%% ------------------------------------------------------------------------
function X = modes_at(S, tf, ts, fs, params, ws, x_thr, U10, groups)
ns = numel(ts); G = max(groups);
X = struct('net', complex(nan(ns, 1)), 'net_coh', nan(ns, 1), 'share', nan(ns, G), ...
    'n_net_unstable', nan(ns, 1), 'abscissa', nan(ns, 1));
for q = 1:ns
    i = round((ts(q) - tf(1)) * fs) + 1;
    J = SRNNCellTypePairs.compute_Jacobian_fast(S(i, :)', params);
    f = full_jacobian_network_mode(J, ws.is_x, x_thr, U10, groups);
    X.net(q) = f.net; X.net_coh(q) = f.net_coh; X.share(q, :) = f.share;
    X.n_net_unstable(q) = f.n_net_unstable; X.abscissa(q) = f.abscissa;
end
end

function [U10, G] = population_basis(params, groups)
% The ten populations: subnet g E is population 2g - 1, subnet g I is 2g;
% orthonormal indicators.
n = params.n; G = max(groups);
type_of = zeros(n, 1);
for q = 1:params.n_cellTypes, type_of(params.type_indices{q}) = q; end
pop = 2 * (groups(:) - 1) + type_of;
U10 = zeros(n, 2 * G);
for k = 1:2 * G, U10(pop == k, k) = 1 / sqrt(nnz(pop == k)); end
end

function d = repro(run_dir, seed, k_bundle, mean_x_E)
% Largest difference of the mean E potential from the bundle's own block file.
d = NaN;
f = fullfile(run_dir, sprintf('seed%03d_block%d.mat', seed, k_bundle));
if ~isfile(f); return; end
L = load(f, 'signals');
a = double(L.signals.mean_x_E(:)); c = double(mean_x_E(:));
if numel(a) == numel(c); d = max(abs(a - c)); end
end

function T = summary_row(N, plan)
b1 = N.block(1); b2 = N.block(2); E = N.event;
pre = E.lags >= -0.3 & E.lags <= 0;
ev_pos = mean(any(real(E.ev.net(:, pre)) > 0, 2));
ctrl_pos = mean(any(real(E.ctrl.net(:, pre)) > 0, 2));
med = median(real(E.ev.net), 1);
[pk, ip] = max(med);
T = table(N.seed, b1.n_events, b2.n_events, b1.n_events / (plan.sim.block_len / 60), b2.n_events / (plan.sim.block_len / 60), ...
    mean(real(b1.net) > 0), mean(real(b2.net) > 0), median(real(b1.net)), median(real(b2.net)), ...
    median(b1.net_coh), median(b2.net_coh), numel(E.t_event), ev_pos, ctrl_pos, pk, E.lags(ip), ...
    max(b1.repro, b2.repro), b1.sim_s + b2.sim_s, b1.eig_s + b2.eig_s, N.seconds, ...
    'VariableNames', {'seed', 'events_base', 'events_stim', 'events_per_min_base', 'events_per_min_stim', ...
    'frac_unstable_base', 'frac_unstable_stim', 'median_re_base', 'median_re_stim', ...
    'median_coh_base', 'median_coh_stim', 'n_events_dense', 'events_unstable_pre', 'ctrl_unstable_pre', ...
    'peak_median_re', 'peak_lag', 'repro_max_diff', 'sim_seconds', 'eig_seconds', 'seconds'});
end
