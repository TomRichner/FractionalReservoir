% RUN_P2_FROZEN_STATE Block B pilot 2: can adaptation be frozen at its pre-stimulation value?
%
% A feasibility pilot (Tom approved the block B pilots on 2026-11-12). The
% planned causal experiment has an arm "stimulation on, adaptation clamped
% at its pre-stimulation values". This script shows that the arm can be run
% from OUTSIDE the model class, with no change to the class or the route
% models: run the model to stimulus onset, then integrate on with the
% class's own integrate() and noise increments and a right-hand side that
% is the class's dynamics_fast with the derivative of every state but x set
% to zero (p2_frozen_rhs). See p2_run_arms for the arms and the checks.
%
% Two networks:
%   (i)  the paper's network, chained MTS (mu7chain, sfa3_stdc2, seed [7 8]),
%        30 s, with the signed step of example_mu7chain_step_signed
%        (+ to the stimulated E neurons, - to the stimulated I neurons; the
%        same neurons and magnitudes) on from 10 s to the end;
%   (ii) bursting_pairs under its most-adapted condition (sfa3_std1),
%        seed [42 42] (the seed of fig_stim_engages_adaptation), 60 s, with a
%        tonic input of +0.4 (that figure's DC level) to a random half of the
%        E neurons from 20 s on. Amplitude and fraction were chosen once,
%        before any run, and not tuned.
%
% Run setup_paths first (cwd at the repo root). Serial, no pool.
% Output: figs/explore_B/p2_frozen_state/ of the MAIN checkout (out_dir
% below): two figures as .png and .fig, and p2_checks.mat / .txt.
%
% See also: p2_run_arms, p2_frozen_rhs, p2_frozen_figure

out_dir = 'C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p2_frozen_state';
if ~isfolder(out_dir), mkdir(out_dir); end
t_start = tic;

%% (i) the paper's network, chained MTS
preset_i = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
spec_i = struct('label', 'mu7chain, chained MTS (sfa3_stdc2), seed [7 8], step +E / -I from 10 s', ...
    'build', @(gain) build_from_preset(preset_i, 'sfa3_stdc2', ...
        'rng_seeds', [7 8], 'T_range', [0 30], 'fs', 400, 'lya_method', 'none', ...
        'store_full_state', true, 'store_decimated_state', false, ...
        'input_config', held_step_config(gain)), ...
    't_on', 10, 'deci', 4);
res_i = p2_run_arms(spec_i);

%% (ii) bursting_pairs
spec_ii = struct('label', 'bursting_pairs (sfa3_std1), seed [42 42], +0.4 to half of E from 20 s', ...
    'build', @(gain) build_from_preset('bursting_pairs', '', ...
        'rng_seeds', [42 42], 'T_range', [0 60], 'fs', 400, 'lya_method', 'none', ...
        'store_full_state', true, 'store_decimated_state', false, ...
        'input_config', tonic_config(gain, 0.4, 0.5, 20)), ...
    't_on', 20, 'deci', 2);
res_ii = p2_run_arms(spec_ii);

%% Figures and the numbers
results = {res_i, res_ii};
tags = {'p2_frozen_mu7chain', 'p2_frozen_bursting'};
fid = fopen(fullfile(out_dir, 'p2_checks.txt'), 'w');
for i = 1:2
    R = results{i};
    fig = p2_frozen_figure(R, 'off');
    exportgraphics(fig, fullfile(out_dir, [tags{i} '.png']), 'Resolution', 110);
    set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');   % the saved .fig opens visible
    savefig(fig, fullfile(out_dir, [tags{i} '.fig']));
    close(fig);
    txt = evalc('disp(R.check)');
    for f = [1 fid]
        fprintf(f, '\n== %s ==\n%s\n', R.label, txt);
        fprintf(f, 'stimulated neurons per type: %s of %s\n', mat2str(R.stim_per_type), mat2str(R.n_per_type));
        fprintf(f, 'stimulus amplitude, min / max over stimulated: %.4f / %.4f\n', ...
            min(R.u_stim(R.stim_rows)), max(R.u_stim(R.stim_rows)));
        arms = fieldnames(R.arm);
        for a = 1:numel(arms)
            A = R.arm.(arms{a});
            fprintf(f, '%-14s rate after onset, mean [%s]  SD over time [%s]  max [%s]   (before: [%s])\n', arms{a}, ...
                num2str(A.rate_mean_after, '%.4f '), num2str(A.rate_sd_after, '%.4f '), ...
                num2str(A.rate_max_after, '%.4f '), num2str(A.rate_mean_before, '%.4f '));
        end
    end
end
run_time_s = toc(t_start);
fprintf(fid, '\nrun time %.0f s\n', run_time_s);
fclose(fid);
save(fullfile(out_dir, 'p2_checks.mat'), 'results', 'run_time_s');
fprintf('run time %.0f s\n', run_time_s);

%% The stimuli
function ic = held_step_config(gain)
% The stimulus of example_mu7chain_step_signed (about 15 % of each type,
% magnitude |0.5 randn|, + to E and - to I), held from the start of the
% middle one of three equal steps (10 s of 30 s) to the end.
ic = struct('n_steps', 3, 'step_density', struct('E', 0.15, 'I', 0.15), 'amp', 0.5, ...
    'no_stim_pattern', [true false true], 'intrinsic_drive', 0, 'positive_only', false, ...
    'generator', @(params, T, fs, seed, config) held_step_input(params, T, fs, seed, config, gain));
end

function [u_ex, t_ex] = held_step_input(params, T, fs, rng_seed, config, gain)
% The stock generator draws the neurons and magnitudes (and leaves the
% random stream where it always does, so x0 is the same for every gain);
% the middle step's amplitudes are then held to the end and signed.
stock = rmfield(config, 'generator');
[u_ex, t_ex] = SRNNCellTypePairs.generate_external_input(params, T, fs, rng_seed, stock);
first = find(any(u_ex ~= 0, 1), 1, 'first');
amp = abs(u_ex(:, first));
signs = [+1, -1];                                   % E, I
for q = 1:params.n_cellTypes
    idx = params.type_indices{q};
    amp(idx) = signs(q) .* amp(idx);
end
u_ex(:) = 0;
u_ex(:, first:end) = repmat(gain .* amp, 1, size(u_ex, 2) - first + 1);
end

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
