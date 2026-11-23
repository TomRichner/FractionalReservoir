function E = soied_event_full(seed, events, S, C)
% SOIED_EVENT_FULL The full Jacobian's rightmost network mode densely around baseline discharges and around quiet control times.
%
%   E = soied_event_full(seed, events, S, C)
%
% seed    network seed; events: R.block(1).events from soied_run_network
% S       the settings soied_run_network was run with (same block 1:
%         the run is deterministic given the noise seed, so this
%         re-simulation reproduces it; checked by comparing mean_x_E)
% C       this call's settings: full_event_dt, n_event_full (number of
%         discharges, the first ones, and of control times), coh_thr
%
% Around each of the first n_use discharges and n_use control times (drawn
% with rng(seed), at least 1 s from every discharge), the full Jacobian is
% decomposed (dense eig with vectors) every S.full_event_dt s over
% [t + dense_win(1), t + dense_win(2)]. Fields per sample: the rightmost
% network mode (>= 10% of norm^2 in the dendritic rows), its coherence over
% the ten populations, and its subnet shares; and the frozen-gain matrix's
% rightmost mode with coh10 >= S.coh_thr at the same state.
%
% See also: soied_run_network, soied_run

n_use = C.n_event_full; S.full_event_dt = C.full_event_dt; S.coh_thr = C.coh_thr;
base = ied_base_config(); base.fs = S.fs;
cfg = base; cfg.seed = seed;
cfg.noise_seed = base.noise_seed + (seed - base.seed) + 1000 * 1;
cfg.duration = S.settle + S.block_len;
ic = struct('dc_levels', 0, 'hold_dur', cfg.duration, 'ramp_dur', 0, ...
    'noise_intensity', 0, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
[m, groups] = ied_build_network(cfg, ic, 'lya_method', 'none', ...
    'T_range', [-cfg.warmup cfg.duration], 'T_plot', [0 cfg.duration], 'store_full_state', true);
m.run();
params = m.get_params(); ws = ospec_weight_stats(params);
n = params.n; G = max(groups);
type_of = zeros(n, 1);
for q = 1:2, type_of(params.type_indices{q}) = q; end
pop = 2 * (groups(:) - 1) + type_of;
U10 = zeros(n, 2 * G);
for k = 1:2 * G, U10(pop == k, k) = 1 / sqrt(nnz(pop == k)); end

p = m.plot_data;
selp = p.t >= S.settle;
E = struct('seed', seed, 'repro_mean_x_E', mean(p.x.E(:, selp & p.t < cfg.duration), 'all'), ...
    'full_event_dt', S.full_event_dt, 'coh_thr', S.coh_thr, 'n_use', n_use);

ev = events.time(:);
use = ev(1:min(n_use, numel(ev)));
rng(seed, 'twister');
ctrl = [];
while numel(ctrl) < numel(use)
    c = S.settle + 1 + rand * (S.block_len - 2);
    if all(abs(ev - c) > 1), ctrl(end + 1) = c; end %#ok<AGROW>
end
lags = S.dense_win(1):S.full_event_dt:S.dense_win(2);
E.lags = lags; E.t_event = use'; E.t_ctrl = ctrl;
tf = m.t_out(:);
kinds = {'event', 'ctrl'}; centres = {use', ctrl};
for kk = 1:2
    c = centres{kk}; ne = numel(c); nl = numel(lags);
    X = struct('net', complex(nan(ne, nl)), 'net_coh10', nan(ne, nl), 'net_share', nan(ne, nl, G), ...
        'netc', complex(nan(ne, nl)), 'M_coh', complex(nan(ne, nl)), 'M_coh_share', nan(ne, nl, G), ...
        'M_lam1', complex(nan(ne, nl)), 'n_net_unstable', nan(ne, nl));
    for e = 1:ne
        for l = 1:nl
            i = round((c(e) + lags(l) - tf(1)) * S.fs) + 1;
            st = m.S_out(i, :)';
            J = SRNNCellTypePairs.compute_Jacobian_fast(st, params);
            [V, D] = eig(full(J), 'vector'); V = V ./ vecnorm(V);
            en = abs(V).^2; xf = sum(en(ws.is_x, :), 1)';
            vx = V(ws.is_x, :);
            ch = (sum(abs(U10' * vx).^2, 1) ./ sum(abs(vx).^2, 1))';
            [~, o] = sort(real(D), 'descend');
            net = o(xf(o) >= 0.1);
            X.n_net_unstable(e, l) = nnz(real(D(net)) > 0);
            k = net(1);
            X.net(e, l) = complex(real(D(k)), abs(imag(D(k)))); X.net_coh10(e, l) = ch(k);
            X.net_share(e, l, :) = accumarray(groups(:), en(ws.is_x, k), [G 1]) / sum(en(ws.is_x, k));
            kc = net(find(ch(net) >= S.coh_thr, 1));
            if ~isempty(kc), X.netc(e, l) = complex(real(D(kc)), abs(imag(D(kc)))); end
            Mf = full(J(ws.is_x, ws.is_x)) * params.tau_d + eye(n);
            [Vm, Dm] = eig(Mf, 'vector'); Vm = Vm ./ vecnorm(Vm);
            cm = sum(abs(U10' * Vm).^2, 1)';
            [~, om] = sort(real(Dm), 'descend');
            X.M_lam1(e, l) = Dm(om(1));
            km = om(find(cm(om) >= S.coh_thr, 1));
            if ~isempty(km)            % no mode that coherent: left NaN
                X.M_coh(e, l) = complex(real(Dm(km)), abs(imag(Dm(km))));
                X.M_coh_share(e, l, :) = accumarray(groups(:), abs(Vm(:, km)).^2, [G 1]);
            end
        end
    end
    E.(kinds{kk}) = X;
end
end
