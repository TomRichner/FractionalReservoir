function out = p2_run_arms(spec)
% P2_RUN_ARMS The intact, frozen and control arms of block B pilot 2 on one network.
%
%   out = p2_run_arms(spec)
%
% spec fields:
%   label        text for the figure and the report
%   build        @(gain) a BUILT SRNNCellTypePairs with store_full_state on
%                and lya_method 'none'; gain 1 = the stimulus, 0 = none.
%                Both must draw the same network, initial state and noise.
%   t_on         nominal onset time of the stimulus (s)
%   deci         keep every deci-th sample for the traces
%
% Arms, all on one seed, input and noise path:
%   intact        the model's own run() with the stimulus
%   frozen        the intact run up to the freeze sample k (the last sample
%                 before the input changes), then x integrated on with a, b,
%                 g held at their values at k (p2_frozen_rhs), through the
%                 class's own integrate() and the class's noise increments
%   nostim        the model's own run() without the stimulus
%   nostim_frozen the same freezing applied to the no-stimulus run
% and one arm that is only a check:
%   resumed       the second segment integrated the same way with NOTHING
%                 frozen; it must reproduce the intact run.
%
% Returns the checks as numbers (out.check) and decimated traces
% (out.arm.<name>) for p2_frozen_figure.

m = spec.build(1);
params = m.cached_params;
params.u_interpolant = m.u_interpolant;
layout = params.state_layout;
x_idx = layout.x;
adapt_idx = setdiff(1:layout.N, x_idx);
t = m.t_ex(:);

% Freeze sample: the last one before the input differs from the no-stimulus
% input. The model interpolates u linearly between samples, so the stimulus
% starts to act in the step AFTER this sample.
m0 = spec.build(0);
changed = find(any(m.u_ex ~= m0.u_ex, 1), 1, 'first');
k = changed - 1;

m.run();
S_intact = m.S_out;
m0.run();
S_nostim = m0.S_out;

chk = struct();
chk.t_freeze = t(k);
chk.solver = m.ode_solver;
chk.n_states = layout.N;
chk.n_frozen_states = numel(adapt_idx);
% The two runs share network, initial state and noise: identical up to k.
chk.stim_vs_nostim_before_onset = max(abs(S_intact(1:k, :) - S_nostim(1:k, :)), [], 'all');

% Second segments. run() drops the noise increments when it is done;
% build_noise regenerates them from the seed, indexed by absolute time, so
% a segment starting at t(k) draws the columns the run drew there.
tspan = t(k:end);
m.build_noise();
p = params;                                         % the stimulus run's input
[~, S2_resumed] = m.integrate(@(tt, S) SRNNCellTypePairs.dynamics_fast(tt, S, p), tspan, S_intact(k, :)');
[~, S2_frozen]  = m.integrate(@(tt, S) p2_frozen_rhs(tt, S, p, x_idx), tspan, S_intact(k, :)');

% Candidate 1 as a cross-check (it has to know the route models' parameter
% names, which is why it is not the method): every adaptation time constant
% infinite on the params struct.
chk.infinite_tau_vs_wrapper = NaN;
chk.infinite_tau_drift = NaN;
chk.infinite_tau_fields = '';
chk.infinite_tau_error = '';
try
    pinf = p;
    for q = 1:numel(pinf.tau_a_matrix)
        pinf.tau_a_matrix{q} = inf(size(pinf.tau_a_matrix{q}));
    end
    known = {'tau_rec', 'tau_rel', 'tau_rec1', 'tau_rec2', 'tau_dec', 'tau_fac', 'tau_h'};
    hit = {};
    for r = 1:numel(pinf.route)
        for f = 1:numel(known)
            if isfield(pinf.route{r}.P, known{f})
                pinf.route{r}.P.(known{f}) = inf(size(pinf.route{r}.P.(known{f})));
                hit{end + 1} = known{f}; %#ok<AGROW>
            end
        end
    end
    chk.infinite_tau_fields = strjoin(unique(hit), ', ');
    [~, S2_inf] = m.integrate(@(tt, S) SRNNCellTypePairs.dynamics_fast(tt, S, pinf), tspan, S_intact(k, :)');
    chk.infinite_tau_vs_wrapper = max(abs(S2_inf - S2_frozen), [], 'all');
    chk.infinite_tau_drift = max(abs(S2_inf(:, adapt_idx) - S2_inf(1, adapt_idx)), [], 'all');
    clear S2_inf
catch err
    chk.infinite_tau_error = err.message;
end

m0.build_noise();
p0 = m0.cached_params;
p0.u_interpolant = m0.u_interpolant;
[~, S2_nostim_frozen] = m0.integrate(@(tt, S) p2_frozen_rhs(tt, S, p0, x_idx), tspan, S_nostim(k, :)');
chk.noise_sigma = m.sigma_u_noise;

% The checks
S_frozen = [S_intact(1:k - 1, :); S2_frozen];
S_nostim_frozen = [S_nostim(1:k - 1, :); S2_nostim_frozen];
chk.before_onset_frozen_vs_intact = max(abs(S_frozen(1:k, :) - S_intact(1:k, :)), [], 'all');
chk.adaptation_drift_frozen = max(abs(S2_frozen(:, adapt_idx) - S2_frozen(1, adapt_idx)), [], 'all');
chk.adaptation_drift_nostim_frozen = max(abs(S2_nostim_frozen(:, adapt_idx) - S2_nostim_frozen(1, adapt_idx)), [], 'all');
chk.adaptation_change_intact = max(abs(S_intact(k:end, adapt_idx) - S_intact(k, adapt_idx)), [], 'all');
chk.restart_state_mismatch = max(abs(S2_frozen(1, :) - S_intact(k, :)));
dx_frozen = abs(diff(S_frozen(:, x_idx), 1, 1));
dx_intact = abs(diff(S_intact(:, x_idx), 1, 1));
chk.x_step_at_onset_frozen = max(dx_frozen(k, :));          % the step k -> k+1
chk.x_step_at_onset_intact = max(dx_intact(k, :));
chk.x_first_step_frozen_vs_intact = max(abs(S_frozen(k + 1, x_idx) - S_intact(k + 1, x_idx)));
chk.x_step_typical_intact = median(max(dx_intact, [], 2));
chk.resumed_vs_intact_x = max(abs(S2_resumed(:, x_idx) - S_intact(k:end, x_idx)), [], 'all');
chk.resumed_vs_intact_all = max(abs(S2_resumed - S_intact(k:end, :)), [], 'all');
clear S2_resumed S2_frozen S2_nostim_frozen dx_frozen dx_intact

% Traces and rates
sel = 1:spec.deci:numel(t);
names = {'intact', 'frozen', 'nostim', 'nostim_frozen'};
states = {S_intact, S_frozen, S_nostim, S_nostim_frozen};
inputs = {m.u_ex, m.u_ex, m0.u_ex, m0.u_ex};
pars = {p, p, p0, p0};
arm = struct();
after = t(sel) > t(k);
for i = 1:numel(names)
    arm.(names{i}) = traces(states{i}(sel, :), inputs{i}(:, sel), pars{i}, after);
end
stim_rows = any(m.u_ex ~= m0.u_ex, 2);
stim_per_type = cellfun(@(idx) nnz(stim_rows(idx)), params.type_indices);

out = struct('label', spec.label, 't', t(sel), 't_on', spec.t_on, 'k', k, ...
    'check', chk, 'arm', arm, 'cell_type_names', {params.cell_type_names}, ...
    'n_per_type', params.n_per_type, 'stim_rows', stim_rows, ...
    'stim_per_type', stim_per_type, 'u_stim', m.u_ex(:, end));
end

function T = traces(S, u, params, after)
[~, a, b, r] = SRNNCellTypePairs.unpack_and_compute_states(S, params);
names = params.cell_type_names;
C = numel(names);
T = struct('rate', [], 'adapt', [], 'dep', [], 'dep_names', {{}}, 'u', []);
nt = size(S, 1);
T.rate = zeros(C, nt);  T.adapt = nan(C, nt);  T.u = zeros(C, nt);
for q = 1:C
    T.rate(q, :) = mean(r.(names{q}), 1);
    if ~isempty(a.(names{q}))
        % the adaptation term in the rate's argument, (c/K) sum_k a_k
        T.adapt(q, :) = mean(params.c_eff(q) .* reshape(sum(a.(names{q}), 2), [], nt), 1);
    end
    T.u(q, :) = mean(u(params.type_indices{q}, :), 1);
end
for pre = 1:C
    for post = 1:C
        bb = b.(names{pre}).(names{post});
        if isempty(bb), continue; end
        D = SRNNCellTypePairs.std_route_factor(bb, params, pre, post);   % the model's factor
        T.dep(end + 1, :) = mean(reshape(D, [], nt), 1);
        T.dep_names{end + 1} = [names{pre} ' to ' names{post}];
    end
end
T.rate_mean_after = mean(T.rate(:, after), 2)';
T.rate_sd_after = std(T.rate(:, after), 0, 2)';      % over time, of the population mean
T.rate_mean_before = mean(T.rate(:, ~after), 2)';
T.rate_max_after = max(T.rate(:, after), [], 2)';
end
