% EXAMPLE_MU7CHAIN_STEPS_TOP30 The mu7chain network: rest, 20 s of back-to-back random steps, rest; top-30 Lyapunov exponents.
%
% The same network and seed as example_mu7chain_step
% (scripts/examples/mu7chain_step_example): the mu7chain bundle's network
% (n = 500, per-neuron SFA and depression spreads), seed [7 8], all four
% adaptation conditions. What differs:
%
%   - 40 s in all: NO stimulus for the first 10 s, then a NEW random step
%     every 5 s for 20 s (four steps back to back, 10 to 30 s), then no
%     stimulus for the last 10 s. Every step is a fresh draw: about 15 % of
%     each cell type is chosen anew and each chosen neuron gets its own
%     constant amplitude stim_amp * randn (either sign) for that step.
%   - The lambda row shows the TOP 30 Lyapunov exponents, accumulating over
%     5 to 40 s: lambda_1 thick, lambda_2..30 thin beneath it, with the
%     final lambda_1 and how many of the 30 ended positive in the corner.
%   - In that row columns 2 to 4 (the three adapting conditions) share
%     their y axis; column 1 (no adaptation) keeps its own, because its
%     exponents are an order of magnitude larger.
%
% Three figures open: the time series, the stimulus, and the eigenvalues of
% W (the same W as in the other mu7chain examples).
%
% Run setup_paths once in the session first (cwd at the repository root);
% like every script in scripts/examples this file has no path code. It takes
% about three and a half minutes (40 s of simulation and 30 exponents per
% condition). Nothing is saved unless save_figs is set.
%
% WHAT TO EXPECT with the settings below, on the code as of 2026-11-12
% (fig_example_timeseries with lambda_K and lambda_link; preset with
% std_spread = [0.25 0.25]):
%   lambda_1:  no_adaptation +2.5514,  sfa1_std1 -0.0991,
%              sfa3_stdc2 -0.1015,     sfa3_std2 -0.0933
%   positive exponents among the top 30:  5, 0, 0, 0
%   stimulus:  zero before 10 s and after 30 s. The four steps stimulate
%              66 (34 E, 32 I), 87 (41 E, 46 I), 74 (37 E, 37 I) and
%              89 (46 E, 43 I) neurons; 253 of 500 are stimulated in at
%              least one step; max |u| = 1.508
%   W:         spectral abscissa 4.8802, spectral radius 5.2761, R = 4.799
% A rerun gives these numbers exactly. A different lambda_1 means the
% model, the preset or this script changed.
%
% READING THE LAMBDA ROW. The accumulation window is the figure function's:
% 5 s to the end, after a 5 s warm-up. So each exponent is an average over
% the rest before the steps, the 20 s of steps and the rest after them
% together; they are not separate values for the stimulated and
% unstimulated stretches. The thin lines are the estimator's columns in the
% order it carries them, not sorted at each instant, so early in the
% accumulation a thin line can lie above the thick lambda_1.
%
% TO TRY SOMETHING ELSE, edit the settings block. The run is cut into equal
% steps of step_len seconds and stim_on says which carry input, so the
% pattern below is [off off on on on on off off] at 5 s a step.
%
% See also: example_mu7chain_step, example_mu7chain_step_signed,
%           fig_example_timeseries, lyapunov_topk

%% Settings
preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
rng_seeds   = [7 8];        % network seed pair, as in example_mu7chain_step
step_len    = 5;            % s; a new random step every step_len seconds
stim_on     = logical([0 0 1 1 1 1 0 0]);   % 10 s off, 20 s of steps, 10 s off
stim_frac   = 0.15;         % share of each cell type chosen anew for every step
stim_amp    = 0.5;          % a chosen neuron's amplitude for a step is stim_amp * randn
lambda_K    = 30;           % how many Lyapunov exponents the lambda row shows
lambda_link = [2 3 4];      % columns of the lambda row that share a y axis
extra_args  = {};           % further name-value overrides of the preset
save_figs   = false;        % true: also write the three figures to out_dir
out_dir     = '';           % '' = the figure function's default folder under figs/

T_range = [0, step_len * numel(stim_on)];   % 40 s

%% The stimulus, in the model's own terms
input_config = struct( ...
    'n_steps',         numel(stim_on), ...
    'step_density',    struct('E', stim_frac, 'I', stim_frac), ...
    'amp',             stim_amp, ...
    'no_stim_pattern', ~stim_on, ...
    'intrinsic_drive', 0, ...
    'positive_only',   false);

%% Simulate the four conditions on one network and draw the three figures
t_start = tic;
out = fig_example_timeseries( ...
    'preset_name',     preset_name, ...
    'rng_seeds',       rng_seeds, ...
    'T_range',         T_range, ...
    'model_args',      [{'input_config', input_config}, extra_args], ...
    'lambda_K',        lambda_K, ...
    'lambda_link',     lambda_link, ...
    'show_stimulus',   true, ...
    'show_W_spectrum', true, ...
    'save',            save_figs, ...
    'out_dir',         out_dir, ...
    'verbose',         'verbose');      % prints lambda_1 per condition

%% Name the windows and say what was drawn
window_names = {'time series, top-30 exponents', 'stimulus', 'eigenvalues of W'};
for k = 1:min(numel(out.figs), numel(window_names))
    set(out.figs(k), 'NumberTitle', 'off', 'Name', ...
        sprintf('mu7chain, seed %s, steps 10-30 s: %s', mat2str(rng_seeds), window_names{k}));
end
fprintf('Neurons stimulated in at least one step: %d.  Spectral abscissa of W: %.4f, spectral radius: %.4f.  %.0f s.\n', ...
    out.stimulus.n_stimulated, out.W_spectrum.abscissa, out.W_spectrum.radius, toc(t_start));
