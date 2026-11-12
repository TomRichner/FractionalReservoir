function X = p4_simulate_extra(out_file)
% P4_SIMULATE_EXTRA The per-neuron rates pilot 4 needs beyond the saved traces of pilots 0a and 0b.
%
%   X = p4_simulate_extra(out_file)
%
% Seven short simulations, serial, about 3 minutes in all. X is a struct
% array of cases, each with: name, network, t (s, 100 Hz), r (n x nt
% per-neuron rate), is_E (n x 1), window ([t0 t1] to analyse), and where it
% applies x, groups (run13 modules) for the IED detector. Saved to out_file
% (-v7.3) when given.
%
%   bursting arms (4)   pilot 2's network (ii): bursting_pairs, sfa3_std1,
%                       seed [42 42], 0..60 s, +0.4 to a random half of E
%                       from 20 s; intact / frozen x stimulus / none, made
%                       exactly as p2_run_arms makes them (the model's run,
%                       then p2_frozen_rhs through the class's integrate).
%                       Analysed over 25..60 s (5 s after onset dropped).
%   run13 stimulated (2) run13 'original' and 'chain_mts' with the protocol
%                       of ied_stim_experiment: uniform input 0 for 60 s,
%                       then 0.20 for 60 s (dc_staircase_stimulus), seed
%                       42. Two windows each: 10..60 s (none) and 70..120 s
%                       (stimulated). 0.20 is one of that experiment's four
%                       amplitudes, chosen once here, not tuned.
%   paper network (2)   the mu7chain preset, conditions no_adaptation and
%                       sfa3_stdc2 (chained MTS), no input, seed [7 8]
%                       (pilot 2's), 10 s settling + 100 s: a reference for
%                       what the measures read on the manuscript's network.
%
% Block B pilot 4 (2026-11-12). See also: run_p4_aberrant_measure,
% p2_run_arms, run13_condition_model, bursting_condition_model

X = struct('name', {}, 'network', {}, 't', {}, 'r', {}, 'is_E', {}, 'window', {}, ...
    'x', {}, 'groups', {});

%% Bursting arms (pilot 2, network ii)
build = @(gain) build_from_preset('bursting_pairs', '', ...
    'rng_seeds', [42 42], 'T_range', [0 60], 'fs', 400, 'lya_method', 'none', ...
    'store_full_state', true, 'store_decimated_state', false, ...
    'input_config', tonic_config(gain, 0.4, 0.5, 20));
m = build(1);  m0 = build(0);
p = m.cached_params;   p.u_interpolant = m.u_interpolant;
p0 = m0.cached_params; p0.u_interpolant = m0.u_interpolant;
x_idx = p.state_layout.x;
t = m.t_ex(:);
k = find(any(m.u_ex ~= m0.u_ex, 1), 1, 'first') - 1;
m.run();  S_intact = m.S_out;
m0.run(); S_nostim = m0.S_out;
tspan = t(k:end);
m.build_noise();
[~, S2] = m.integrate(@(tt, S) p2_frozen_rhs(tt, S, p, x_idx), tspan, S_intact(k, :)');
S_frozen = [S_intact(1:k - 1, :); S2];
m0.build_noise();
[~, S2] = m0.integrate(@(tt, S) p2_frozen_rhs(tt, S, p0, x_idx), tspan, S_nostim(k, :)');
S_nostim_frozen = [S_nostim(1:k - 1, :); S2];
names = {'intact_nostim', 'intact_stim', 'frozen_stim', 'frozen_nostim'};
states = {S_nostim, S_intact, S_frozen, S_nostim_frozen};
pars = {p0, p, p, p0};
sel = 1:4:numel(t);
for i = 1:4
    [~, ~, ~, r] = SRNNCellTypePairs.unpack_and_compute_states(states{i}(sel, :), pars{i});
    X(end + 1) = struct('name', ['bursting_' names{i}], 'network', 'bursting arms', ...
        't', t(sel)', 'r', [r.E; r.I], ...
        'is_E', [true(size(r.E, 1), 1); false(size(r.I, 1), 1)], ...
        'window', [25 60], 'x', [], 'groups', []); %#ok<AGROW>
end
clear S_intact S_nostim S_frozen S_nostim_frozen S2 states

%% run13, original and chain MTS, 0 then 0.20 uniform input
hold_s = 60; A = 0.20;
ic = struct('dc_levels', [0 A], 'hold_dur', hold_s, 'ramp_dur', 0, ...
    'noise_intensity', 0, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
conds = {'original', 'chain_mts'};
for c = 1:2
    [m, groups] = run13_condition_model(conds{c}, 'duration', 2 * hold_s, 'input_config', ic);
    m.run();
    pd = m.plot_data;
    tt = pd.t(:)';
    step = max(1, round(0.01 / median(diff(tt))));
    s = 1:step:numel(tt);
    rr = [pd.r.E; pd.r.I];  xx = [pd.x.E; pd.x.I];
    isE = [true(size(pd.r.E, 1), 1); false(size(pd.r.I, 1), 1)];
    wins = {[10 hold_s], [hold_s + 10, 2 * hold_s]};
    tags = {'nostim', 'stim'};
    for w = 1:2
        X(end + 1) = struct('name', sprintf('run13_%s_%s', conds{c}, tags{w}), ...
            'network', 'run13 stimulated', 't', tt(s), 'r', single(rr(:, s)), 'is_E', isE, ...
            'window', wins{w}, 'x', single(xx(:, s)), 'groups', groups); %#ok<AGROW>
    end
end

%% The paper's network (mu7chain), no input
preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
conds = {'no_adaptation', 'sfa3_stdc2'};
for c = 1:2
    m = build_from_preset(preset, conds{c}, 'rng_seeds', [7 8], 'T_range', [-10 100], ...
        'T_plot', [0 100], 'fs', 400, 'plot_deci', 4, 'lya_method', 'none');
    m.run();
    pd = m.plot_data;
    tt = pd.t(:)';
    step = max(1, round(0.01 / median(diff(tt))));
    s = 1:step:numel(tt);
    rr = [pd.r.E; pd.r.I];
    isE = [true(size(pd.r.E, 1), 1); false(size(pd.r.I, 1), 1)];
    X(end + 1) = struct('name', ['mu7chain_' conds{c}], 'network', 'paper network', ...
        't', tt(s), 'r', single(rr(:, s)), 'is_E', isE, 'window', [0 100], ...
        'x', [], 'groups', []); %#ok<AGROW>
end

if nargin >= 1 && ~isempty(out_file)
    save(out_file, 'X', '-v7.3');
end
end

%% Pilot 2's tonic stimulus (a local function of run_p2_frozen_state, copied)
function ic = tonic_config(gain, amp, frac, t_on)
ic = struct('intrinsic_drive', 0, ...
    'generator', @(params, T, fs, seed, config) tonic_input(params, T, fs, seed, gain, amp, frac, t_on));
end

function [u_ex, t_ex] = tonic_input(params, T, fs, rng_seed, gain, amp, frac, t_on)
% +amp to a random fraction of the first cell type (E) from t_on to the end.
rng(rng_seed);
t_ex = (0:1/fs:T)';
idx = params.type_indices{1};
chosen = idx(rand(numel(idx), 1) < frac);
u_ex = zeros(params.n, numel(t_ex));
u_ex(chosen, t_ex >= t_on) = gain * amp;
end
