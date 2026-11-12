% EXAMPLE_MU7CHAIN_STEP The mu7chain network under a step stimulus: time series, stimulus, eigenvalues of W.
%
% Reproduces the example Tom looked at on 2026-11-12 before the essential
% run (FD-023): one network of the mu7chain bundle (n = 500, per-neuron
% SFA and depression spreads), seed [7 8], simulated for 15 s under all
% four adaptation conditions, with a step applied to about 15 % of the
% neurons of each type for the middle 5 s. Three figures open:
%
%   1. the example time series, one column per condition (x, rate, synaptic
%      output, SFA, depression factor, local and accumulated lambda_1);
%   2. the stimulus: which neurons receive input, when, and how much;
%   3. the eigenvalues of the connectivity matrix W that is simulated, with
%      the bulk radius R, the spectral abscissa and the spectral radius.
%
% Run setup_paths once in the session first (cwd at the repository root);
% like every script in scripts/examples this file has no path code. It takes
% about 75 s. Nothing is saved unless save_figs is set.
%
% WHAT TO EXPECT with the settings below, on the code as of 2026-11-12
% (main at 9c9c1f1 or later; preset with std_spread = [0.25 0.25]):
%   lambda_1:  no_adaptation +2.6269,  sfa1_std1 +0.0488,
%              sfa3_stdc2 -0.1402,     sfa3_std2 -0.1089
%   stimulus:  73 of 500 neurons (34 E, 39 I), on from 5 to 10 s, zero
%              elsewhere, amplitudes about -1 to +1.44
%   W:         spectral abscissa 4.8802, spectral radius 5.2761, R = 4.799
% One seed pair fixes W, the setpoints, the SFA and depression draws, the
% stimulus and the noise path, so a rerun gives these numbers exactly. A
% different lambda_1 means the model, the preset or this script changed.
%
% TO TRY SOMETHING ELSE, edit the settings block. The stimulus is the
% model's own step input (input_config): the run is cut into n_steps equal
% steps, each stimulated neuron gets its own constant amplitude amp * randn
% for a step, and no_stim_pattern says which steps are off. So "the middle
% 5 s" is the middle one of three steps of a 15 s run. Further name-value
% overrides of the preset go in extra_args, for example another
% connectivity:  extra_args = {'mu_tilde_relative', [6 -6; 6 -6]};
%
% See also: fig_example_timeseries, srnn_param_preset, mu7chain_config,
%           SRNNCellTypePairs

%% Settings
preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
rng_seeds   = [7 8];        % network seed pair
T_range     = [0 15];       % s; three steps of 5 s
stim_frac   = 0.15;         % share of each cell type that receives the step
stim_amp    = 0.5;          % each stimulated neuron's amplitude is stim_amp * randn
stim_on     = [false true false];   % which of the equal steps carry input (the middle one)
extra_args  = {};           % further name-value overrides of the preset
save_figs   = false;        % true: also write the three figures to out_dir
out_dir     = '';           % '' = the figure function's default folder under figs/

%% The stimulus, in the model's own terms
input_config = struct( ...
    'n_steps',         numel(stim_on), ...
    'step_density',    struct('E', stim_frac, 'I', stim_frac), ...
    'amp',             stim_amp, ...
    'no_stim_pattern', ~stim_on, ...
    'intrinsic_drive', 0, ...
    'positive_only',   false);

%% Simulate the four conditions on one network and draw the three figures
out = fig_example_timeseries( ...
    'preset_name',     preset_name, ...
    'rng_seeds',       rng_seeds, ...
    'T_range',         T_range, ...
    'model_args',      [{'input_config', input_config}, extra_args], ...
    'show_stimulus',   true, ...
    'show_W_spectrum', true, ...
    'save',            save_figs, ...
    'out_dir',         out_dir, ...
    'verbose',         'verbose');      % prints lambda_1 per condition

%% Name the windows and say what was drawn
window_names = {'time series', 'stimulus', 'eigenvalues of W'};
for k = 1:min(numel(out.figs), numel(window_names))
    set(out.figs(k), 'NumberTitle', 'off', 'Name', ...
        sprintf('mu7chain, seed %s: %s', mat2str(rng_seeds), window_names{k}));
end
fprintf('Stimulated neurons: %d.  Spectral abscissa of W: %.4f, spectral radius: %.4f.\n', ...
    out.stimulus.n_stimulated, out.W_spectrum.abscissa, out.W_spectrum.radius);
