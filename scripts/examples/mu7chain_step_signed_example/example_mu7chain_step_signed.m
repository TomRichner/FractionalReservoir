% EXAMPLE_MU7CHAIN_STEP_SIGNED The mu7chain step example with the step positive to E and negative to I.
%
% The same network, seed and stimulated neurons as example_mu7chain_step
% (scripts/examples/mu7chain_step_example), with ONE change: the sign of
% the step. There each stimulated neuron's amplitude is stim_amp * randn,
% of either sign. Here a stimulated EXCITATORY neuron receives
% +|stim_amp * randn| and a stimulated INHIBITORY neuron -|stim_amp * randn|:
% excitation is pushed up and inhibition down for the middle 5 s.
%
% Everything else is identical neuron for neuron: the mu7chain bundle's
% network (n = 500, per-neuron SFA and depression spreads), seed [7 8],
% 15 s under all four adaptation conditions, about 15 % of each cell type
% stimulated from 5 to 10 s, the same neurons with the same magnitudes, and
% the same initial state and noise path (the generator below draws exactly
% what the stock one draws and only then sets the signs, so the random
% stream is left where the stock generator leaves it). Three figures open:
% the time series, the stimulus, and the eigenvalues of W (the same W as in
% the unsigned example).
%
% Run setup_paths once in the session first (cwd at the repository root);
% like every script in scripts/examples this file has no path code. It takes
% about 75 s. Nothing is saved unless save_figs is set.
%
% WHAT TO EXPECT with the settings below, on the code as of 2026-11-12
% (main at 9c9c1f1 or later; preset with std_spread = [0.25 0.25]):
%   lambda_1:  no_adaptation +1.8595,  sfa1_std1 -0.0391,
%              sfa3_stdc2 -0.1502,     sfa3_std2 -0.1073
%              (the unsigned example: +2.6269, +0.0488, -0.1402, -0.1089)
%   stimulus:  73 of 500 neurons (34 E, 39 I), on from 5 to 10 s, zero
%              elsewhere; every E amplitude positive (0.013 to 0.749),
%              every I amplitude negative (-1.440 to -0.002); |u| equals
%              the unsigned example's input at every neuron and sample
%   W:         spectral abscissa 4.8802, spectral radius 5.2761, R = 4.799
% A rerun gives these numbers exactly. A different lambda_1 means the
% model, the preset or this script changed.
%
% HOW THE SIGN IS SET. The model's step input has no per-type sign, but
% input_config may carry a `generator`, a function the model calls in place
% of its own (SRNNCellTypePairs.generate_external_input). signed_step_input,
% at the end of this file, calls the stock generator and then makes the E
% rows positive and the I rows negative. Change e_sign / i_sign in the
% settings for another combination (both +1 is the model's positive_only).
%
% See also: example_mu7chain_step, fig_example_timeseries,
%           SRNNCellTypePairs.generate_external_input

%% Settings
preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
rng_seeds   = [7 8];        % network seed pair, as in example_mu7chain_step
T_range     = [0 15];       % s; three steps of 5 s
stim_frac   = 0.15;         % share of each cell type that receives the step
stim_amp    = 0.5;          % a stimulated neuron's magnitude is |stim_amp * randn|
stim_on     = [false true false];   % which of the equal steps carry input (the middle one)
e_sign      = +1;           % sign of the step on stimulated E neurons
i_sign      = -1;           % sign of the step on stimulated I neurons
extra_args  = {};           % further name-value overrides of the preset
save_figs   = false;        % true: also write the three figures to out_dir
out_dir     = '';           % '' = the figure function's default folder under figs/

%% The stimulus: the model's step input, then the sign per cell type
input_config = struct( ...
    'n_steps',         numel(stim_on), ...
    'step_density',    struct('E', stim_frac, 'I', stim_frac), ...
    'amp',             stim_amp, ...
    'no_stim_pattern', ~stim_on, ...
    'intrinsic_drive', 0, ...
    'positive_only',   false, ...
    'generator',       @(params, T, fs, seed, config) ...
                           signed_step_input(params, T, fs, seed, config, e_sign, i_sign));

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
        sprintf('mu7chain, seed %s, step +E / -I: %s', mat2str(rng_seeds), window_names{k}));
end
fprintf('Stimulated neurons: %d.  Spectral abscissa of W: %.4f, spectral radius: %.4f.\n', ...
    out.stimulus.n_stimulated, out.W_spectrum.abscissa, out.W_spectrum.radius);

%% The generator
function [u_ex, t_ex] = signed_step_input(params, T, fs, rng_seed, config, e_sign, i_sign)
% The model's own step input with a fixed sign per cell type.
%
% Calls the stock generator on the same config (without the `generator`
% field, or it would call this function again), so the neurons chosen, the
% magnitudes and the state of the random stream afterwards are exactly the
% stock ones; then every E row becomes e_sign * |u| and every I row
% i_sign * |u|. Refuses a nonzero intrinsic_drive, which |u| would fold in.
if any(config.intrinsic_drive ~= 0, 'all')
    error('example_mu7chain_step_signed:IntrinsicDrive', ...
        'signed_step_input takes |u| per neuron and cannot carry an intrinsic drive.');
end
stock = rmfield(config, 'generator');
[u_ex, t_ex] = SRNNCellTypePairs.generate_external_input(params, T, fs, rng_seed, stock);
signs = struct('E', e_sign, 'I', i_sign);
for q = 1:params.n_cellTypes
    name = params.cell_type_names{q};
    if ~isfield(signs, name)
        error('example_mu7chain_step_signed:CellType', ...
            'No sign is given for cell type ''%s'' (this example has E and I).', name);
    end
    idx = params.type_indices{q};
    u_ex(idx, :) = signs.(name) .* abs(u_ex(idx, :));
end
end
