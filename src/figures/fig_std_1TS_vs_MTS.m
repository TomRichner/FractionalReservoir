function out = fig_std_1TS_vs_MTS(cfg)
% FIG_STD_1TS_VS_MTS What single- and multiple-timescale depression share, and where they differ.
%
%   out = FIG_STD_1TS_VS_MTS()
%   out = FIG_STD_1TS_VS_MTS('preset_name', p, 'conditions', names)
%
% Figure tag Fig_STD_1TS_vs_MTS, with Fig_STD_1TS_vs_MTS_table.md beside it.
% FD-014 (Plan A, step 7). No run directory: every number comes from the
% route models of the conditions of cfg.preset_name.
%
% THE CLAIM IT SHOWS. The single-timescale condition (1TS, sfa1_std1) has one
% depression pool. The multiple-timescale condition (MTS, sfa3_stdc2) is the
% chained vesicle pools: the SAME releasable-pool equation plus a reserve
% that refills it (docs/EquationsParametersDocs/Equations_stability_paper.md,
% "Chained vesicle-pool depression"). The reserve starts full, so the two
% coincide at the onset of activity and come apart only as sustained activity
% draws the reserve down. Under a sinusoidal rate the peak synaptic output on
% the FIRST cycle is therefore close between the two, and at the periodic
% steady state the chain's is lower. They are not matched in steady-state
% depth, by decision. The product MTS (sfa3_std2, two independent pools), when
% selected, is drawn in green as the reference the reviewer objected to: it
% is deeper from the first cycle on.
%
% PANELS (one figure, 16 x 9 tiles):
%   (a)-(c)  a sinusoidal rate r(t) = (r_max/2)(1 - cos 2 pi f t), from rest,
%            over cfg.n_cycles cycles: the rate, the depression factor of each
%            condition (the chain's reserve b_2 dashed), and the synaptic
%            output theta = r F. Open circles mark each condition's first-cycle
%            peak output, filled diamonds its last-cycle peak.
%   (d)      the peak output of every cycle, per condition: the first-cycle
%            match and the separation to the periodic steady state at a glance.
%   (e)-(g)  rate steps from rest (cfg.step_rates, cfg.on_s on, cfg.off_s
%            off): the same three rows, with each condition's steady state at
%            the step's rate dotted. Coincidence at each onset, separation
%            during sustained activity, and the chain's slower recovery.
%   (h)      the steady-state depression factor against a constant rate (the
%            chain's reserve dashed), with the probe rates marked.
%   (i)      the sustained synaptic output against rate, with each monotone
%            curve's high-rate ceiling dotted and the product's peak marked.
%
% EVERY NUMBER IS THE ROUTE MODEL'S. The route is compiled by
% SRNNCellTypePairs.compile_route from the condition's synapse_config; time
% courses integrate its `rhs` with classical RK4 at cfg.fs
% (integrate_route; for the sinusoid the rate is evaluated at each step's
% start, middle and end), the depression factor is its `depression`, the
% output factor its `output`, and the steady state its closed-form
% `steady_state`. The figure contains no equation of its own.
%
% THE TABLE (written when cfg.save; statistics are there, not on the
% panels): per condition the first- and last-cycle peak outputs, their phase
% relative to the rate's peak, their ratios to 1TS, the periodic-state
% convergence, how long the factors stay within cfg.sep_tol of 1TS after
% onset, the recovery after each step, the steady-state factor and output at
% cfg.probe_rates, and the ceiling; then the first- and last-cycle ratios to
% 1TS over a grid of sinusoid frequencies (cfg.grid_freq) and peak rates
% (cfg.grid_rmax), so the robustness of the first-cycle statement can be
% judged. The 1TS reference is the selected condition whose route has ONE
% product pool; without one the ratios are NaN.
%
% Conditions with no depression on the route (no_adaptation) are skipped and
% listed in out.skipped. Titles and colours by condition name
% (srnn_condition_titles, manuscript_style). The 1TS trace is drawn wide and
% underneath, the others narrower on top, so a coincidence reads as a thin
% line inside a wide one rather than as one line hiding the other.
%
% See also: fig_STD_steady_state, integrate_route, route_pool, route_product,
%           scripts/explorations/adaptation_normalization/demo_chain_1TS_vs_MTS

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none'
    cfg.preset_name (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'
    % Condition names to draw, in order; empty = every depressing condition of
    % the preset (select_conditions).
    cfg.conditions                = {}
    cfg.route_pre   (1,:) char    = 'E'
    cfg.route_post  (1,:) char    = 'E'
    % The sinusoid Tom showed Brian Lundstrom (demo_chain_1TS_vs_MTS): 0 to
    % r_max at f. Ten cycles reach the periodic steady state to well under 0.1%
    % at the paper's parameters (the table reports the last cycle's change).
    cfg.sine_freq   (1,1) double  = 0.1
    cfg.sine_rmax   (1,1) double  = 0.5
    cfg.n_cycles    (1,1) double  = 10
    % Rate steps: settle, then on/off per rate. 15 s off lets the 1TS pool
    % recover almost fully (tau_rec 2 s) and shows the chain's reserve lagging.
    cfg.step_rates  (1,:) double  = [0.25 0.5 1]
    cfg.settle_s    (1,1) double  = 2
    cfg.on_s        (1,1) double  = 5
    cfg.off_s       (1,1) double  = 15
    cfg.fs          (1,1) double  = 1000
    % Rate at which the traces are DRAWN (see draw). 20 Hz is 200 points per
    % sinusoid cycle and about 4 per release time constant at r = 1; at 50 Hz
    % MATLAB's SVG export warns that the vector content is too large.
    cfg.plot_fs     (1,1) double  = 20
    cfg.probe_rates (1,:) double  = [0.1 0.25 0.5 1]
    % Robustness grid for the first-cycle ratio (table only).
    cfg.grid_freq   (1,:) double  = [0.05 0.1 0.2 0.5 1]
    cfg.grid_rmax   (1,:) double  = [0.25 0.5 1]
    % "Coincide": |F - F_1TS| <= sep_tol * F_1TS.
    cfg.sep_tol     (1,1) double  = 0.01
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
fig_tag = 'Fig_STD_1TS_vs_MTS';
label   = sprintf('%s->%s', cfg.route_pre, cfg.route_post);
dt      = 1 / cfg.fs;

[preset, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
names = select_conditions(cond_names, cfg.conditions);
zero_floor = isfield(preset, 'std_zero_floor') && preset.std_zero_floor;

%% ---- Routes ---------------------------------------------------------------
C = struct('name', {}, 'title', {}, 'color', {}, 'route_cfg', {}, 'model', {}, ...
    'P', {}, 'has_reserve', {}, 'is_single', {}, 'desc', {});
skipped = {};
for k = 1:numel(names)
    cnd = conditions{strcmp(cond_names, names{k})};
    std_cfg = route_std(cnd, cfg.route_pre, cfg.route_post);
    if isempty(std_cfg)
        skipped{end+1} = names{k}; %#ok<AGROW>
        continue
    end
    stf = [];
    e = cnd.synapse_config.(cfg.route_pre).(cfg.route_post);
    if isfield(e, 'stf'); stf = e.stf; end
    c = struct();
    c.name      = names{k};
    t_k         = titles_by_name(names(k));
    c.title     = t_k{1};
    c.color     = condition_color(st, names{k}, k);
    c.route_cfg = struct('std', std_cfg, 'stf', stf);
    R = compile_n(c.route_cfg, label, zero_floor, 1);
    c.model     = R.model;
    c.P         = R.P;
    c.has_reserve = strcmp(R.P.model, 'pool');
    c.is_single = strcmp(R.P.model, 'product') && R.P.n_b == 1 && R.P.n_g == 0;
    c.desc      = route_description(R.P);
    C(end+1) = c; %#ok<AGROW>
end
if isempty(C)
    error('fig_std_1TS_vs_MTS:NoDepression', ...
        'None of the conditions {%s} of preset ''%s'' carries depression on route %s.', ...
        strjoin(names, ', '), cfg.preset_name, label);
end
if ~isempty(skipped)
    vprintf(cfg.verbose, 'verbose', '[fig_std_1TS_vs_MTS] no depression on %s, not drawn: %s\n', ...
        label, strjoin(skipped, ', '));
end
nC = numel(C);
ref = find([C.is_single], 1);            % the 1TS reference for ratios
if isempty(ref); ref = NaN; end

%% ---- Sinusoid -------------------------------------------------------------
[sine, n_per] = sine_protocol(cfg.sine_freq, cfg.sine_rmax, cfg.n_cycles, cfg.fs);
for k = 1:nC
    S = integrate_route(C(k).model, C(k).P, sine.r0, dt, sine.rh, sine.r1);
    tr = trace_of(C(k), S, sine.r1);
    pk = cycle_peaks(tr.theta, sine.r, n_per);
    tr.peaks = pk;
    sine.traces(k) = tr;
end
sine.first = arrayfun(@(tr) tr.peaks.theta(1), sine.traces);
sine.last  = arrayfun(@(tr) tr.peaks.theta(end), sine.traces);
sine.first_phase = arrayfun(@(tr) tr.peaks.lead_deg(1), sine.traces);
sine.last_phase  = arrayfun(@(tr) tr.peaks.lead_deg(end), sine.traces);
sine.converged   = arrayfun(@(tr) abs(tr.peaks.theta(end) - tr.peaks.theta(end-1)) ...
    / tr.peaks.theta(end), sine.traces);
sine.sep_s = separation_time(sine.traces, ref, 1, cfg.sep_tol, dt);

%% ---- Steps ----------------------------------------------------------------
seg_rate = [0, reshape([cfg.step_rates; zeros(1, numel(cfg.step_rates))], 1, [])];
seg_dur  = [cfg.settle_s, repmat([cfg.on_s, cfg.off_s], 1, numel(cfg.step_rates))];
seg_n    = round(seg_dur * cfg.fs);
steps = struct();
steps.r_steps = repelem(seg_rate, seg_n);             % rate over step i
steps.t  = (0:sum(seg_n)) * dt;                        % sample 1 is t = 0, at rest
steps.r  = [0, steps.r_steps];                         % rate shown at each sample
seg_end  = [0, cumsum(seg_n)] + 1;                     % sample index at each segment boundary
steps.onset_idx = seg_end(2:2:end-1);                  % start of each "on" segment
steps.on_end_idx  = seg_end(3:2:end);
steps.off_end_idx = seg_end(4:2:end);
steps.seg_start_s = [0, cumsum(seg_dur)];
for k = 1:nC
    S = integrate_route(C(k).model, C(k).P, steps.r_steps, dt);
    steps.traces(k) = trace_of(C(k), S, steps.r_steps);
end
nst = numel(cfg.step_rates);
steps.sep_s = nan(nC, nst);
steps.recovery = nan(nC, nst);
steps.b2_after = nan(nC, nst);
for j = 1:nst
    steps.sep_s(:, j) = separation_time(steps.traces, ref, steps.onset_idx(j), cfg.sep_tol, dt);
    for k = 1:nC
        D = steps.traces(k).D;
        d_on = D(steps.on_end_idx(j)); d_off = D(steps.off_end_idx(j));
        steps.recovery(k, j) = (d_off - d_on) / (1 - d_on);
        if C(k).has_reserve
            steps.b2_after(k, j) = steps.traces(k).b2(steps.off_end_idx(j));
        end
    end
end
% Each condition's steady state at the step rates (dotted on the panels).
for k = 1:nC
    [Fp, b2p] = steady_factor(C(k), cfg.step_rates, label, zero_floor);
    steps.ss_D(k, :) = Fp.D; steps.ss_theta(k, :) = cfg.step_rates .* Fp.F;
    steps.ss_b2(k, :) = b2p;
end

%% ---- Steady state against rate -------------------------------------------
r_ss = linspace(0, 1, 2001);
rf   = linspace(0, 1, 4e5);
r_inf = 1e9;
ssd = struct('r', r_ss);
for k = 1:nC
    [Fs, b2s] = steady_factor(C(k), r_ss, label, zero_floor);
    ssd.D(k, :) = Fs.D; ssd.theta(k, :) = r_ss .* Fs.F; ssd.b2(k, :) = b2s;
    Ff = steady_factor(C(k), rf, label, zero_floor);
    [ssd.theta_peak(k), ip] = max(rf .* Ff.F);
    ssd.r_peak(k)   = rf(ip);
    ssd.has_peak(k) = ip < numel(rf);
    Fi = steady_factor(C(k), r_inf, label, zero_floor);
    ceil_k = r_inf * Fi.F;
    if ssd.has_peak(k) || ceil_k < 1e-6; ceil_k = NaN; end
    ssd.ceiling(k) = ceil_k;
    [Fp, b2p] = steady_factor(C(k), cfg.probe_rates, label, zero_floor);
    ssd.probe_D(k, :) = Fp.D; ssd.probe_theta(k, :) = cfg.probe_rates .* Fp.F;
    ssd.probe_b2(k, :) = b2p;
end

%% ---- Robustness grid (table only) -----------------------------------------
grid = grid_ratios(C, cfg, label, zero_floor, ref);

%% ---- Figure ---------------------------------------------------------------
fig = draw(C, ref, sine, steps, ssd, cfg, st);

%% ---- Output ---------------------------------------------------------------
if cfg.visible; set(fig, 'Visible', 'on'); end
summary = struct('name', {C.name}, 'title', {C.title}, 'route', {C.desc}, ...
    'first_peak', num2cell(sine.first), 'last_peak', num2cell(sine.last), ...
    'first_phase_deg', num2cell(sine.first_phase), 'last_phase_deg', num2cell(sine.last_phase), ...
    'ceiling', num2cell(ssd.ceiling), 'probe_D', num2cell(ssd.probe_D, 2)', ...
    'probe_theta', num2cell(ssd.probe_theta, 2)');
out = struct('figs', fig, 'files', {{}}, 'source', ['preset: ' cfg.preset_name], ...
    'conditions', {{C.name}}, 'skipped', {skipped}, 'reference', '', ...
    'summary', summary, 'sine', rmfield(sine, {'r0', 'rh', 'r1'}), 'steps', steps, ...
    'steady_state', ssd, 'grid', grid, 'table', '');
if ~isnan(ref); out.reference = C(ref).name; end
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    out.table = fullfile(out_dir, [fig_tag '_table.md']);
    write_table(out.table, C, ref, sine, steps, ssd, grid, cfg, label);
end
end

%% ========================================================================
%  Computation
%% ========================================================================
function [sine, n_per] = sine_protocol(f, rmax, n_cycles, fs)
% Rate at every sample (sine.t, sine.r) and at each RK4 step's start, middle
% and end (r0, rh, r1). One cycle must be a whole number of samples.
n_per = fs / f;
if abs(n_per - round(n_per)) > 1e-9
    error('fig_std_1TS_vs_MTS:Grid', 'fs/f = %g must be an integer number of samples.', n_per);
end
n_per = round(n_per);
nt = n_per * n_cycles;
dt = 1 / fs;
rfun = @(tt) (rmax / 2) * (1 - cos(2 * pi * f * tt));
ts = (0:nt - 1) * dt;                               % step starts
sine = struct('f', f, 'rmax', rmax, 'n_cycles', n_cycles);
sine.t  = (0:nt) * dt;
sine.r  = rfun(sine.t);
sine.r0 = rfun(ts);
sine.rh = rfun(ts + dt / 2);
sine.r1 = rfun(ts + dt);
end

function tr = trace_of(c, S, r_end)
% Factor, reserve and output at every sample, the initial state (rest)
% prepended as sample 1 (t = 0). S is nt x ns (one neuron).
s0 = c.model.initial(c.P);
S = [s0; S];
s3 = reshape(S.', 1, c.P.ns, []);                  % 1 x ns x (nt + 1)
tr = struct();
tr.D = reshape(c.model.depression(s3(:, 1:c.P.n_b, :), c.P), 1, []);
tr.b2 = [];
if c.has_reserve; tr.b2 = S(:, 2).'; end
F = reshape(c.model.output(s3, c.P), 1, []);
tr.theta = [0, r_end] .* F;                         % the rate is 0 at t = 0
tr.peaks = [];
end

function pk = cycle_peaks(theta, r, n_per)
% Peak output in each cycle [c-1, c]/f (both ends included) and its lead on
% the rate's peak, in degrees of the cycle (positive: output peaks first).
n_cyc = (numel(theta) - 1) / n_per;
pk = struct('theta', zeros(1, n_cyc), 'lead_deg', zeros(1, n_cyc));
for c = 1:n_cyc
    idx = (c - 1) * n_per + (1:n_per + 1);
    [pk.theta(c), it] = max(theta(idx));
    [~, ir] = max(r(idx));
    pk.lead_deg(c) = 360 * (ir - it) / n_per;
end
end

function sep = separation_time(traces, ref, i0, tol, dt)
% Time from sample i0 until each condition's factor first leaves
% tol * F_1TS of the 1TS factor (NaN for the reference itself or when there
% is none; Inf when it never leaves within the trace).
n = numel(traces);
sep = nan(n, 1);
if isnan(ref); return; end
Dr = traces(ref).D(i0:end);
for k = 1:n
    if k == ref; continue; end
    d = abs(traces(k).D(i0:end) - Dr) > tol * Dr;
    j = find(d, 1);
    if isempty(j); sep(k) = Inf; else; sep(k) = (j - 1) * dt; end
end
end

function [F, b2] = steady_factor(c, r, label, zero_floor)
% The route model's closed-form steady state at each rate in r (a row):
% F.F the output factor, F.D the depression factor; b2 the reserve (NaN
% without one).
r = r(:);
P = compile_n(c.route_cfg, label, zero_floor, numel(r)).P;
ss = c.model.steady_state(r, P);
F.F = reshape(ss.F, 1, []);
F.D = reshape(c.model.depression(ss.s(:, 1:P.n_b), P), 1, []);
b2 = nan(1, numel(r));
if c.has_reserve; b2 = reshape(ss.s(:, 2), 1, []); end
end

function g = grid_ratios(C, cfg, label, zero_floor, ref)
% First- and last-cycle peak output over (frequency x peak rate), every grid
% point a row of one integration per condition, all run for the same time
% (cfg.n_cycles cycles of the slowest frequency).
fq = cfg.grid_freq; am = cfg.grid_rmax;
[FF, AA] = ndgrid(fq, am);
rows = numel(FF);
T = cfg.n_cycles / min(fq);
dt = 1 / cfg.fs;
nt = round(T * cfg.fs);
n_per = cfg.fs ./ FF(:);
if any(abs(n_per - round(n_per)) > 1e-9) || any(abs(nt ./ n_per - round(nt ./ n_per)) > 1e-9)
    error('fig_std_1TS_vs_MTS:Grid', ...
        'Every grid frequency must give whole cycles of whole samples at fs = %g over %g s.', cfg.fs, T);
end
n_per = round(n_per);
ts = (0:nt - 1) * dt;
rfun = @(tt) (AA(:) / 2) .* (1 - cos(2 * pi * FF(:) .* tt));
r0 = rfun(ts); rh = rfun(ts + dt / 2); r1 = rfun(ts + dt);
first = nan(numel(fq), numel(am), numel(C));
last  = first;
for k = 1:numel(C)
    R = compile_n(C(k).route_cfg, label, zero_floor, rows);
    S = integrate_route(R.model, R.P, r0, dt, rh, r1);    % nt x ns x rows
    F = reshape(R.model.output(permute(S, [3 2 1]), R.P), rows, nt);
    th = [zeros(rows, 1), r1 .* F];
    for q = 1:rows
        n1 = n_per(q);
        first(q + numel(FF) * (k - 1)) = max(th(q, 1:n1 + 1));
        last(q + numel(FF) * (k - 1))  = max(th(q, end - n1:end));
    end
end
g = struct('freq', fq, 'rmax', am, 'first', first, 'last', last, 'T', T, ...
    'names', {{C.name}});
if isnan(ref)
    g.ratio_first = nan(size(first)); g.ratio_last = nan(size(last));
else
    g.ratio_first = first ./ first(:, :, ref);
    g.ratio_last  = last  ./ last(:, :, ref);
end
end

function R = compile_n(route_cfg, label, zero_floor, n)
opts = struct('npre', n, 'label', label, 'zero_floor', zero_floor, ...
    'err_id_prefix', 'fig_std_1TS_vs_MTS');
R = SRNNCellTypePairs.compile_route(route_cfg, opts);
end

function std_cfg = route_std(cnd, pre, post)
% The route's std configuration, or [] when the condition has no depression there.
std_cfg = [];
if ~isfield(cnd, 'synapse_config') || ~isstruct(cnd.synapse_config); return; end
sc = cnd.synapse_config;
if ~isfield(sc, pre) || ~isfield(sc.(pre), post) || ~isfield(sc.(pre).(post), 'std'); return; end
s = sc.(pre).(post).std;
if ~isempty(s) && ~(isstruct(s) && isempty(fieldnames(s)))
    std_cfg = s;
end
end

function s = route_description(P)
% The route's configured values, from the compiled route's nominal (plain text).
n = P.nominal;
v = @(x) strjoin(compose('%g', x), ' ');
if strcmp(P.model, 'pool')
    s = sprintf('chained pools, tau_rec = [%s] s, tau_rel = %s s, n2 = %g', ...
        v(n.tau_rec), v(n.tau_rel), n.n2);
elseif P.n_b == 1
    s = sprintf('one pool, tau_rec = %s s, tau_rel = %s s', v(n.tau_rec), v(n.tau_rel));
else
    s = sprintf('product of %d pools, tau_rec = [%s] s, tau_rel = [%s] s', ...
        P.n_b, v(n.tau_rec), v(n.tau_rel));
end
if ~isempty(n.tau_fac)
    s = sprintf('%s, facilitation tau_fac = [%s] s', s, v(n.tau_fac));
end
end

function col = condition_color(st, name, k)
if isKey(st.condition_color, name)
    col = st.condition_color(name);
else
    L = lines(max(k, 7));
    col = L(k, :);
end
end

%% ========================================================================
%  Drawing
%% ========================================================================
function fig = draw(C, ref, sine, steps, ssd, cfg, st)
nC = numel(C);
tick_fs = st.tick_fs - 2; label_fs = st.label_fs - 2; title_fs = st.title_fs;
rate_color = [0.25 0.25 0.25];
% Traces are DRAWN at cfg.plot_fs (every number is computed at cfg.fs): at
% 1 kHz each line held 1e5 points and the SVG export warned (it still does
% at 50 Hz). The steps keep
% the samples on either side of every rate change, so the edges stay sharp
% and each onset's output (the full rate) is drawn exactly.
dec = max(1, round(cfg.fs / cfg.plot_fs));
js = unique([1:dec:numel(sine.t), numel(sine.t)]);
jt = unique([1:dec:numel(steps.t), numel(steps.t), steps.onset_idx, steps.onset_idx + 1, ...
    steps.on_end_idx, steps.on_end_idx + 1]);
% Drawing order and widths: 1TS first, wide, underneath; then the others in
% reverse selection order, so the first after 1TS ends on top.
order = numel(C):-1:1;
if ~isnan(ref); order = [ref, order(order ~= ref)]; end
lw = st.line_lw * ones(1, nC);
if ~isnan(ref); lw(:) = 0.8 * st.line_lw; lw(ref) = 1.7 * st.line_lw; end

fig_size = [1400, 1150];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', 'Visible', 'off', ...
    'Position', [scr(1:2) + max((scr(3:4) - fig_size) / 2, 0), fig_size]);
nrow = 16; ncol = 9;
tl = tiledlayout(fig, nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
tile = @(r, c) (r - 1) * ncol + c;
L = panel_letters(9);
fmt = @(ax) set(ax, 'FontSize', tick_fs, 'Box', 'off', 'TickDir', 'out', 'LineWidth', st.axis_lw);

%% Sinusoid: rate, factor, output
ax_sr = nexttile(tl, tile(1, 1), [2 6]); hold(ax_sr, 'on');
plot(ax_sr, sine.t(js), sine.r(js), '-', 'Color', rate_color, 'LineWidth', 1.2);
fmt(ax_sr); ylim(ax_sr, [0 1.1 * sine.rmax]); set(ax_sr, 'XTickLabel', []);
ylabel(ax_sr, 'rate r', 'FontSize', label_fs);
panel_title(ax_sr, L{1}, sprintf('Sinusoidal rate: %g Hz, 0 to %g, from rest', ...
    sine.f, sine.rmax), title_fs);

ax_sd = nexttile(tl, tile(3, 1), [3 6]); hold(ax_sd, 'on');
for k = order
    if C(k).has_reserve
        plot(ax_sd, sine.t(js), sine.traces(k).b2(js), '--', 'Color', C(k).color, 'LineWidth', 1.2);
    end
    plot(ax_sd, sine.t(js), sine.traces(k).D(js), '-', 'Color', C(k).color, 'LineWidth', lw(k));
end
fmt(ax_sd); ylim(ax_sd, [0 1.05]); set(ax_sd, 'XTickLabel', [], 'YTick', [0 0.5 1]);
ylabel(ax_sd, 'depression factor', 'FontSize', label_fs);
panel_title(ax_sd, L{2}, 'Depression factor (chain: releasable pool b_1; reserve b_2 dashed)', title_fs);

ax_so = nexttile(tl, tile(6, 1), [3 6]); hold(ax_so, 'on');
n_per = (numel(sine.t) - 1) / sine.n_cycles;
for k = order
    plot(ax_so, sine.t(js), sine.traces(k).theta(js), '-', 'Color', C(k).color, 'LineWidth', lw(k));
end
for k = order
    th = sine.traces(k).theta;
    [~, i1] = max(th(1:n_per + 1));
    [~, iL] = max(th(end - n_per:end)); iL = iL + numel(th) - n_per - 1;
    plot(ax_so, sine.t(i1), th(i1), 'o', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', C(k).color, 'MarkerFaceColor', 'white');
    plot(ax_so, sine.t(iL), th(iL), 'd', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', C(k).color);
end
fmt(ax_so);
ymax_s = max(arrayfun(@(tr) max(tr.theta), sine.traces)) * 1.2;
ylim(ax_so, [0 ymax_s]);
xlabel(ax_so, 'time (s)', 'FontSize', label_fs);
ylabel(ax_so, 'synaptic output  \theta = r F', 'FontSize', label_fs);
panel_title(ax_so, L{3}, 'Synaptic output (open circle: first-cycle peak; filled diamond: last-cycle peak)', title_fs);
linkaxes([ax_sr, ax_sd, ax_so], 'x'); xlim(ax_so, [0 sine.t(end)]);
% One tick per two cycles. (With one per cycle the PNG export dropped the
% '0' of the '70' label, though the figure's XTickLabel held '70'.)
xt = 0:(max(1, round(sine.n_cycles / 5)) / sine.f):sine.t(end);
set([ax_sr, ax_sd, ax_so], 'XTick', xt);
set(ax_so, 'XTickLabel', compose('%g', xt));

%% Peak output per cycle + legend
ax_lg = nexttile(tl, tile(1, 7), [2 3]); axis(ax_lg, 'off'); hold(ax_lg, 'on');
hl = gobjects(0); tl_txt = {};
for k = 1:nC
    hl(end+1) = plot(ax_lg, NaN, NaN, '-', 'Color', C(k).color, 'LineWidth', 2.5); %#ok<AGROW>
    tl_txt{end+1} = C(k).title; %#ok<AGROW>
end
if any([C.has_reserve])
    hl(end+1) = plot(ax_lg, NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    tl_txt{end+1} = 'reserve b_2 (chain)';
end
hl(end+1) = plot(ax_lg, NaN, NaN, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);
tl_txt{end+1} = 'steady state';
legend(ax_lg, hl, tl_txt, 'Box', 'off', 'FontSize', tick_fs - 1, 'Location', 'west');

ax_pk = nexttile(tl, tile(3, 7), [6 3]); hold(ax_pk, 'on');
cyc = 1:sine.n_cycles;
for k = order
    p = sine.traces(k).peaks.theta;
    plot(ax_pk, cyc, p, '-', 'Color', C(k).color, 'LineWidth', max(1.2, lw(k) * 0.75));
    plot(ax_pk, cyc(1), p(1), 'o', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', C(k).color, 'MarkerFaceColor', 'white');
    plot(ax_pk, cyc(end), p(end), 'd', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', C(k).color);
end
fmt(ax_pk); xlim(ax_pk, [0.5, sine.n_cycles + 0.5]); ylim(ax_pk, [0 ymax_s]);
set(ax_pk, 'XTick', unique([1, 5:5:sine.n_cycles, sine.n_cycles]));
xlabel(ax_pk, 'cycle', 'FontSize', label_fs);
ylabel(ax_pk, 'peak output per cycle', 'FontSize', label_fs);
panel_title(ax_pk, L{4}, 'Peak output, cycle by cycle', title_fs);

%% Steps: rate, factor, output
ax_tr = nexttile(tl, tile(9, 1), [2 6]); hold(ax_tr, 'on');
plot(ax_tr, steps.t(jt), steps.r(jt), '-', 'Color', rate_color, 'LineWidth', 1.2);
fmt(ax_tr); ylim(ax_tr, [0 1.1 * max(steps.r)]); set(ax_tr, 'XTickLabel', []);
ylabel(ax_tr, 'rate r', 'FontSize', label_fs);
panel_title(ax_tr, L{5}, sprintf('Rate steps from rest: %s, %g s on, %g s off', ...
    strjoin(compose('%g', cfg.step_rates), ', '), cfg.on_s, cfg.off_s), title_fs);

ax_td = nexttile(tl, tile(11, 1), [3 6]); hold(ax_td, 'on');
ax_to = nexttile(tl, tile(14, 1), [3 6]); hold(ax_to, 'on');
on_s = steps.seg_start_s;
for k = order
    for j = 1:numel(cfg.step_rates)
        x = on_s([2 * j, 2 * j + 1]);
        plot(ax_td, x, steps.ss_D(k, [j j]), ':', 'Color', C(k).color, 'LineWidth', 1.4);
        plot(ax_to, x, steps.ss_theta(k, [j j]), ':', 'Color', C(k).color, 'LineWidth', 1.4);
    end
    if C(k).has_reserve
        plot(ax_td, steps.t(jt), steps.traces(k).b2(jt), '--', 'Color', C(k).color, 'LineWidth', 1.2);
    end
    plot(ax_td, steps.t(jt), steps.traces(k).D(jt), '-', 'Color', C(k).color, 'LineWidth', lw(k));
    plot(ax_to, steps.t(jt), steps.traces(k).theta(jt), '-', 'Color', C(k).color, 'LineWidth', lw(k));
end
fmt(ax_td); ylim(ax_td, [0 1.05]); set(ax_td, 'XTickLabel', [], 'YTick', [0 0.5 1]);
ylabel(ax_td, 'depression factor', 'FontSize', label_fs);
panel_title(ax_td, L{6}, 'Depression factor: together at each onset, apart while the rate is held', title_fs);
% Log scale: at onset the output is the full rate r (the pools are full), so on
% a linear axis the onset transients (up to r = 1) flatten the sustained
% levels, which is where the conditions differ, into the bottom tenth. The
% output is zero while the rate is (not drawn on a log axis).
fmt(ax_to); set(ax_to, 'YScale', 'log', 'YTick', [0.01 0.03 0.1 0.3 1]);
ylim(ax_to, [0.01, max(arrayfun(@(tr) max(tr.theta), steps.traces)) * 1.3]);
xlabel(ax_to, 'time (s)', 'FontSize', label_fs);
ylabel(ax_to, 'synaptic output  \theta = r F', 'FontSize', label_fs);
panel_title(ax_to, L{7}, 'Synaptic output (log scale): the full rate at onset, then each sustained level', title_fs);
linkaxes([ax_tr, ax_td, ax_to], 'x'); xlim(ax_to, [0 steps.t(end)]);

%% Steady state against rate
ax_ss = nexttile(tl, tile(9, 7), [4 3]); hold(ax_ss, 'on');
ax_st = nexttile(tl, tile(13, 7), [4 3]); hold(ax_st, 'on');
for k = order
    if C(k).has_reserve
        plot(ax_ss, ssd.r, ssd.b2(k, :), '--', 'Color', C(k).color, 'LineWidth', 1.2);
    end
    plot(ax_ss, ssd.r, ssd.D(k, :), '-', 'Color', C(k).color, 'LineWidth', lw(k));
    plot(ax_ss, cfg.probe_rates, ssd.probe_D(k, :), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', C(k).color, 'MarkerEdgeColor', 'none');
    plot(ax_st, ssd.r, ssd.theta(k, :), '-', 'Color', C(k).color, 'LineWidth', lw(k));
    if ssd.has_peak(k)
        plot(ax_st, ssd.r_peak(k), ssd.theta_peak(k), 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', C(k).color, 'MarkerEdgeColor', 'none');
    elseif ~isnan(ssd.ceiling(k))
        plot(ax_st, [0 1], ssd.ceiling(k) * [1 1], ':', 'Color', C(k).color, 'LineWidth', 1.4);
    end
end
fmt(ax_ss); xlim(ax_ss, [0 1]); ylim(ax_ss, [0 1.05]);
set(ax_ss, 'XTickLabel', [], 'YTick', [0 0.5 1]);
ylabel(ax_ss, 'steady-state factor', 'FontSize', label_fs);
panel_title(ax_ss, L{8}, 'Steady state against rate', title_fs);
fmt(ax_st); xlim(ax_st, [0 1]);
ylim(ax_st, [0 max([ssd.theta(:); ssd.ceiling(:)]) * 1.15]);
xlabel(ax_st, 'constant rate r', 'FontSize', label_fs);
ylabel(ax_st, 'sustained output', 'FontSize', label_fs);
panel_title(ax_st, L{9}, 'Sustained output (dotted: ceiling)', title_fs);
end

function panel_title(ax, letter, txt, fs)
title(ax, sprintf('%s  %s', letter, txt), 'FontWeight', 'normal', 'FontSize', fs, ...
    'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1 0]);
end

%% ========================================================================
%  Table
%% ========================================================================
function write_table(path, C, ref, sine, steps, ssd, grid, cfg, label)
fid = fopen(path, 'w');
if fid < 0
    warning('fig_std_1TS_vs_MTS:Table', 'Could not write %s.', path);
    return
end
closer = onCleanup(@() fclose(fid));
nC = numel(C);
if isnan(ref)
    ref_txt = 'none of the drawn conditions has a single pool, so the ratios are NaN';
else
    ref_txt = sprintf('`%s` (%s)', C(ref).name, C(ref).title);
end
fprintf(fid, '# Single- vs multiple-timescale depression (Fig_STD_1TS_vs_MTS)\n\n');
fprintf(fid, 'Preset `%s`, route %s. Every number is the route model''s (compile_route; time courses by RK4 on its rhs at %g Hz; steady states from its closed form). Ratios are to the 1TS reference, %s.\n\n', ...
    cfg.preset_name, label, cfg.fs, ref_txt);
fprintf(fid, '| Condition | Title | Route |\n|---|---|---|\n');
for k = 1:nC
    fprintf(fid, '| `%s` | %s | %s |\n', C(k).name, strrep(C(k).title, '\tau', 'tau'), C(k).desc);
end

fprintf(fid, '\n## Sinusoidal rate\n\n');
fprintf(fid, ['Input: r(t) = (r_max/2)(1 - cos(2 pi f t)) with f = %g Hz and r_max = %g, ' ...
    'starting at rest (every pool full, r(0) = 0), for %d cycles. The FIRST cycle is ' ...
    't in [0, 1/f] = [0, %g] s; the LAST cycle is the final one, t in [%g, %g] s, by which ' ...
    'the response is periodic (the last column: relative change of the peak from the cycle ' ...
    'before). Peak = the largest synaptic output theta = r F in that cycle. Phase = how far ' ...
    'the output peak leads the rate peak, in degrees of the cycle (positive: output first).\n\n'], ...
    sine.f, sine.rmax, sine.n_cycles, 1 / sine.f, (sine.n_cycles - 1) / sine.f, sine.n_cycles / sine.f);
fprintf(fid, '| Condition | first-cycle peak | phase (deg) | last-cycle peak | phase (deg) | first / 1TS first | last / 1TS last | last / first | last-cycle change |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
for k = 1:nC
    rf_ = ratio(sine.first, k, ref); rl_ = ratio(sine.last, k, ref);
    fprintf(fid, '| `%s` | %.4f | %+.1f | %.4f | %+.1f | %.4f | %.4f | %.4f | %.1e |\n', ...
        C(k).name, sine.first(k), sine.first_phase(k), sine.last(k), sine.last_phase(k), ...
        rf_, rl_, sine.last(k) / sine.first(k), sine.converged(k));
end
fprintf(fid, ['\nHow long the depression factor stays within %g%% of 1TS''s (|F - F_1TS| <= %g F_1TS), ' ...
    'measured from rest: under the sinusoid from t = 0, and from the onset of each rate step ' ...
    '(the second and third onsets follow %g s of recovery, not full rest).\n\n'], ...
    100 * cfg.sep_tol, cfg.sep_tol, cfg.off_s);
hdr = '| Condition | sinusoid (s) |'; sep = '|---|---|';
for j = 1:numel(cfg.step_rates)
    hdr = [hdr, sprintf(' step r = %g (s) |', cfg.step_rates(j))]; %#ok<AGROW>
    sep = [sep, '---|']; %#ok<AGROW>
end
fprintf(fid, '%s\n%s\n', hdr, sep);
for k = 1:nC
    if k == ref; continue; end
    fprintf(fid, '| `%s` | %s |', C(k).name, num_s(sine.sep_s(k)));
    fprintf(fid, ' %s |', strjoin(arrayfun(@num_s, steps.sep_s(k, :), 'UniformOutput', false), ' | '));
    fprintf(fid, '\n');
end

fprintf(fid, '\n## Rate steps: recovery\n\n');
fprintf(fid, ['Fraction of the depression recovered in the %g s off after each %g s step: ' ...
    '(F_end_off - F_end_on) / (1 - F_end_on). For the chain, its reserve b_2 at the end of the off period too.\n\n'], ...
    cfg.off_s, cfg.on_s);
hdr = '| Condition |'; sep = '|---|';
for j = 1:numel(cfg.step_rates)
    hdr = [hdr, sprintf(' after r = %g |', cfg.step_rates(j))]; %#ok<AGROW>
    sep = [sep, '---|']; %#ok<AGROW>
end
fprintf(fid, '%s\n%s\n', hdr, sep);
for k = 1:nC
    cells = arrayfun(@(j) rec_cell(steps, k, j, C(k).has_reserve), 1:numel(cfg.step_rates), ...
        'UniformOutput', false);
    fprintf(fid, '| `%s` | %s |\n', C(k).name, strjoin(cells, ' | '));
end

fprintf(fid, '\n## Steady state at constant rate\n\n');
fprintf(fid, ['Depression factor F (the chain: its releasable pool b_1), sustained output r F, and the ' ...
    'high-rate limit of r F (the ceiling; a product of pools has none, its output peaks and falls, ' ...
    'so its peak is given instead).\n\n']);
hdr = '| Condition |'; sep = '|---|';
for j = 1:numel(cfg.probe_rates)
    hdr = [hdr, sprintf(' F at r = %g |', cfg.probe_rates(j))]; %#ok<AGROW>
    sep = [sep, '---|']; %#ok<AGROW>
end
for j = 1:numel(cfg.probe_rates)
    hdr = [hdr, sprintf(' r F at r = %g |', cfg.probe_rates(j))]; %#ok<AGROW>
    sep = [sep, '---|']; %#ok<AGROW>
end
hdr = [hdr, ' ceiling of r F |']; sep = [sep, '---|'];
fprintf(fid, '%s\n%s\n', hdr, sep);
for k = 1:nC
    if ~isnan(ssd.ceiling(k))
        ceil_txt = sprintf('%.4f', ssd.ceiling(k));
    elseif ssd.has_peak(k)
        ceil_txt = sprintf('none (peak %.4f at r = %.3f)', ssd.theta_peak(k), ssd.r_peak(k));
    else
        ceil_txt = 'none';
    end
    fprintf(fid, '| `%s` | %s | %s | %s |\n', C(k).name, ...
        strjoin(compose('%.4f', ssd.probe_D(k, :)), ' | '), ...
        strjoin(compose('%.4f', ssd.probe_theta(k, :)), ' | '), ceil_txt);
end

fprintf(fid, '\n## How robust is the first-cycle match?\n\n');
fprintf(fid, ['The same sinusoid at other frequencies f (rows) and peak rates r_max (columns), each ' ...
    'from rest and run for %g s (%d cycles of the slowest). Entries are the first-cycle peak output ' ...
    'over 1TS''s first-cycle peak, then the same for the last cycle.\n'], grid.T, cfg.n_cycles);
for k = 1:nC
    if k == ref; continue; end
    for what = {'first', 'last'}
        M = grid.(['ratio_' what{1}])(:, :, k);
        fprintf(fid, '\n`%s`, %s-cycle peak / 1TS:\n\n', C(k).name, what{1});
        fprintf(fid, '| f (Hz) | %s |\n', strjoin(compose('r_max = %g', grid.rmax), ' | '));
        fprintf(fid, '|---|%s\n', repmat('---|', 1, numel(grid.rmax)));
        for i = 1:numel(grid.freq)
            fprintf(fid, '| %g | %s |\n', grid.freq(i), strjoin(compose('%.3f', M(i, :)), ' | '));
        end
    end
end
end

function v = ratio(x, k, ref)
if isnan(ref); v = NaN; else; v = x(k) / x(ref); end
end

function s = num_s(x)
if isinf(x); s = 'never'; elseif isnan(x); s = '-'; else; s = sprintf('%.2f', x); end
end

function s = rec_cell(steps, k, j, has_reserve)
s = sprintf('%.3f', steps.recovery(k, j));
if has_reserve
    s = sprintf('%s (b_2 = %.3f)', s, steps.b2_after(k, j));
end
end
