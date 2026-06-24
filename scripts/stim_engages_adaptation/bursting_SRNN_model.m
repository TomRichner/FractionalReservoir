% bursting_SRNN_model.m - SRNNModel2 scratch driver with ALL major params exposed
%
% Goal: tune defaults until the network produces BURSTS of activity (stimulus
% engages adaptation/depression). Every knob that meaningfully shapes the
% dynamics is pulled up to the top of this script and grouped by category, so
% you can edit in place instead of digging into SRNNModel2.m.
%
% Lifecycle is unchanged: construct -> build() -> run() -> plot().
%
% Each parameter below mirrors a property of src/SRNNModel2.m; the comment
% after it gives the class default. Anything you leave commented falls back to
% that default. See CLAUDE.md and SRNNModel2.m for the governing equations.

close all; clc;

% --- RNG seeds: persist & advance across successive runs --------------------
% On the first run, initialize rng_seeds = [1 2]. On every subsequent run,
% increment both seeds so each execution explores a new network / stimulus
% realization. clearvars wipes everything EXCEPT rng_seeds, so the value
% survives from one run to the next (it lives in the base workspace).
if exist('rng_seeds', 'var')
    rng_seeds = rng_seeds + 1
else
    rng_seeds = [1 2]
end
clearvars -except rng_seeds;

% Add src/ and scripts/ to the MATLAB path
setup_paths();

%% ======================================================================
%  1. NETWORK ARCHITECTURE
%  ======================================================================
n          = 50;      % total number of neurons              (default 300; old bursting net = 10)
f          = 0.5;     % fraction excitatory, n_E = round(f*n) (default 0.5; old EI = 0.7)
indegree   = 4;       % expected in-degree -> alpha = indegree/n (default 100; old mean degree ~4)

%% ======================================================================
%  2. CONNECTIVITY  (Random Matrix Theory, Harris 2023 tilde-notation)
%  ======================================================================
%  W is built by RMTMatrix from the tilde means/stds, then scaled by
%  level_of_chaos. The tilde params are conventionally expressed as multiples
%  of the normalization factor F = 1/sqrt(N*alpha*(2-alpha)) (so that R=1 when
%  all four are equal to F). F is computed below from n & indegree so you can
%  keep editing in "F units" exactly like the class defaults do.
alpha = indegree / n;
F     = 1 / sqrt(n * alpha * (2 - alpha));   % RMT normalization factor

mu_E_tilde    =  3*F;   % normalized E mean      (class default  3*F)
mu_I_tilde    = -3*F;   % normalized I mean      (class default -4*F)
sigma_E_tilde =  1*F;   % normalized E std dev   (class default  1*F)
sigma_I_tilde =  1*F;   % normalized I std dev   (class default  1*F)
E_W           =  0;     % common mean offset added to mu_E_tilde & mu_I_tilde (default 0)
zrs_mode      = 'none'; % 'none' | 'ZRS' | 'SZRS' | 'Partial_SZRS'           (default 'none')

level_of_chaos      = 1.0;    % multiplicative scale on W; >1 past edge of chaos (default 1.0)
rescale_by_abscissa = false;  % if true, rescale W by 1/abscissa_0 before chaos scaling (default false)

%% ======================================================================
%  3. SPIKE-FREQUENCY ADAPTATION (SFA)   a-variables
%  ======================================================================
%  da_{i,k}/dt = (-a_{i,k} + r_i)/tau_k ;  firing rate sees  x - c*sum_k a_k
n_a_E   = 3;     % # adaptation timescales for E neurons (default 0; set >0 for SFA)
n_a_I   = 0;     % # adaptation timescales for I neurons (default 0)
c_E     = 0.5/3;   % adaptation strength for E   (default 0.15/3 ~ 0.05; old bursting = 0.5/3 ~ 0.167)
c_I     = 0.15/3;  % adaptation strength for I   (default 0.15/3 ~ 0.05)

% Adaptation time constants (row vectors). Leave [] to auto-fill with
% logspace(log10(0.25), log10(10), n_a) inside build_network().
tau_a_E = logspace(log10(0.3), log10(15), n_a_E);    % old bursting: [0.3 ... 15] s (slower tail)
tau_a_I = [];                                        % auto if n_a_I>0

%% ======================================================================
%  4. SHORT-TERM SYNAPTIC DEPRESSION (STD)   b-variables
%  ======================================================================
%  db_i/dt = (1-b_i)/tau_rec - (b_i*r_i)/tau_rel ;  synaptic output is b*r
%  NOTE: n_b_E, n_b_I must be 0 or 1.
n_b_E       = 1;     % enable STD for E (0/1)            (default 0)
n_b_I       = 0;     % enable STD for I (0/1)            (default 0)
tau_b_E_rec = 1;     % E recovery time constant (s)      (default 1; old bursting = 1)
tau_b_E_rel = 1;     % E release  time constant (s)      (default 0.25; old bursting = 1 = tau_STD)
tau_b_I_rec = 1;     % I recovery time constant (s)      (default 1)
tau_b_I_rel = 0.25;  % I release  time constant (s)      (default 0.25)

%% ======================================================================
%  5. INTRINSIC DYNAMICS & NONLINEARITY
%  ======================================================================
%  dx_i/dt = (-x_i + sum_j w_ij r_j + u_i)/tau_d ;  r_i = b_i * phi(x_i - c*sum_k a_k)
tau_d = 0.025;   % dendritic time constant (s)   (default 0.1; old bursting = 0.025)

% --- Nonlinearity phi ---------------------------------------------------
% 'piecewise' : piecewiseSigmoid(x, S_a, S_c)  -- bounded [0,1], parameterized below
% 'logistic'  : logisticSigmoid(x)             -- bounded [0,1]
% 'tanh'      : tanhActivation(x)              -- bounded [-1,1]
% 'relu'      : max(0, x)  -- UNBOUNDED, matches old SRNN model
nonlinearity = 'relu';
S_a = 0.9;       % piecewiseSigmoid slope param a (default 0.9)
S_c = 0.35;      % piecewiseSigmoid center  param c (default 0.35)

switch lower(nonlinearity)
    case 'piecewise'
        phi      = @(x) SRNNModel2.piecewiseSigmoid(x, S_a, S_c);
        phi_deriv = @(x) SRNNModel2.piecewiseSigmoidDerivative(x, S_a, S_c);
    case 'logistic'
        phi      = @(x) SRNNModel2.logisticSigmoid(x);
        phi_deriv = @(x) SRNNModel2.logisticSigmoidDerivative(x);
    case 'tanh'
        phi      = @(x) SRNNModel2.tanhActivation(x);
        phi_deriv = @(x) SRNNModel2.tanhActivationDerivative(x);
    case 'relu'
        phi      = @(x) max(0, x);
        phi_deriv = @(x) double(x > 0);
    otherwise
        error('Unknown nonlinearity: %s', nonlinearity);
end

%% ======================================================================
%  6. EXTERNAL STIMULUS  (old IED-example drive: tonic DC + elevated epochs)
%  ======================================================================
%  Reproduces the stimulus from SRNN_caller_example_tseries.m: a small tonic
%  DC applied UNIFORMLY to every neuron (ramped on over the first ramp_dur
%  seconds), with two elevated-drive epochs added on top. Built by the custom
%  generator dc_step_stimulus() at the BOTTOM of this file, wired in via
%  input_config.generator (which overrides the built-in sparse-step generator).
%
%  NOTE: SRNNModel2 forces u = 0 during the negative warmup (t < 0), so this
%  profile applies to t in [0, T_range(2)]. The ramp therefore runs over the
%  first ramp_dur seconds of the positive window, not during the warmup.
input_config = struct();
input_config.DC       = 0.02;   % tonic baseline drive, all neurons (old DC = 0.02)
input_config.ramp_dur = 10;     % seconds to ramp 0 -> DC at start  (old ramp_duration = 10)
input_config.epochs   = [20 35 0.08;    % [t_start  t_end  amp_added] elevated epochs
                         35 45 0.48];    % old: +0.08 (20-35s) then +0.48 (35-45s), all neurons
input_config.intrinsic_drive = [];      % unused by this generator; kept for class compatibility
input_config.generator = @dc_step_stimulus;  % custom generator (see bottom of file)

u_ex_scale = 1.0;     % global scale on the external input (default 1.0)

%% ======================================================================
%  7. SIMULATION / INTEGRATION
%  ======================================================================
fs         = 400;          % sampling frequency (Hz)            (default 400)
T_range    = [-30, 100];   % integration interval [start end] s (default [0 50]); negative start = warmup
T_plot     = [];           % plotting window; [] -> T_range
ode_solver = @ode45;       % ODE solver handle                  (default @ode45)
% rng_seeds is set & auto-incremented at the TOP of this script (see header).
% [network seed, stimulus seed]. Do not reassign it here, or the per-run
% increment will be clobbered.

%% ======================================================================
%  8. LYAPUNOV / ANALYSIS
%  ======================================================================
%  'none' is fastest while tuning; 'benettin' adds LLE; 'qr' = full spectrum.
lya_method = 'none';       % 'none' | 'benettin' | 'qr'         (default 'benettin')

%% ======================================================================
%  9. PLOTTING
%  ======================================================================
%  plot_deci sets the plot decimation: effective plot rate = fs / plot_deci.
%  Set explicitly here so it overrides the class default (round(fs/plot_freq)
%  = 40 -> 10 Hz). plot_deci = 10 with fs = 400 -> 40 Hz plotting.
plot_deci = 10;            % decimation factor for plotting (fs/plot_deci = 40 Hz)
store_full_state = true;   % keep full-res S_out (needed for the PSD of x)

%% ======================================================================
%  10. PSD ANALYSIS (pwelch of mean dendritic potential x)
%  ======================================================================
%  Single overall PSD of the mean dendritic potential, full-resolution (fs Hz).
%  Mirrors the old template (SRNN_caller_example_tseries_raster_adapt_psd.m).
psd_t_start      = 0;     % analysis window start (s); use t > this (drop warmup/transient)
psd_win_len_s    = 15;    % Hamming window length (s)   [template used 15]
psd_overlap_frac = 0.75;  % segment overlap fraction    [template used 0.75]
psd_f            = logspace(log10(0.05), log10(100), 50);  % requested freqs (Hz)

%% ======================================================================
%  11. CONSTRUCT, BUILD, RUN, PLOT
%  ======================================================================
model = SRNNModel2( ...
    'n', n, 'f', f, 'indegree', indegree, ...
    'mu_E_tilde', mu_E_tilde, 'mu_I_tilde', mu_I_tilde, ...
    'sigma_E_tilde', sigma_E_tilde, 'sigma_I_tilde', sigma_I_tilde, ...
    'E_W', E_W, 'zrs_mode', zrs_mode, ...
    'level_of_chaos', level_of_chaos, 'rescale_by_abscissa', rescale_by_abscissa, ...
    'n_a_E', n_a_E, 'n_a_I', n_a_I, 'c_E', c_E, 'c_I', c_I, ...
    'tau_a_E', tau_a_E, 'tau_a_I', tau_a_I, ...
    'n_b_E', n_b_E, 'n_b_I', n_b_I, ...
    'tau_b_E_rec', tau_b_E_rec, 'tau_b_E_rel', tau_b_E_rel, ...
    'tau_b_I_rec', tau_b_I_rec, 'tau_b_I_rel', tau_b_I_rel, ...
    'tau_d', tau_d, 'S_a', S_a, 'S_c', S_c, ...
    'activation_function', phi, 'activation_function_derivative', phi_deriv, ...
    'input_config', input_config, 'u_ex_scale', u_ex_scale, ...
    'fs', fs, 'T_range', T_range, 'T_plot', T_plot, ...
    'ode_solver', ode_solver, 'rng_seeds', rng_seeds, ...
    'lya_method', lya_method, 'plot_deci', plot_deci, ...
    'store_full_state', store_full_state);

model.build();
model.run();
[fig_handle, ax_handles] = model.plot();

% --- Un-cap the firing-rate & synaptic-output panels -----------------------
% SRNNModel2.plot_firing_rate / plot_synaptic_output hardcode ylim([0 1]),
% which clips the unbounded ReLU output. Re-autoscale those axes here (matched
% by y-label, so panel order doesn't matter) rather than editing the class.
for k = 1:numel(ax_handles)
    ylab = get(get(ax_handles(k), 'YLabel'), 'String');
    if any(strcmp(ylab, {'firing rate', 'synaptic output'}))
        ylim(ax_handles(k), 'auto');
        yticks(ax_handles(k), 'auto');
    end
end

%% ======================================================================
%  PSD of mean dendritic potential (x)  --  pwelch, full-resolution
%  ======================================================================
% Extract full-resolution x (last n columns of S_out) and average over neurons.
nE = model.n_E; nI = model.n_I;
len_a = nE*model.n_a_E + nI*model.n_a_I;     % a_E, a_I block lengths
len_b = nE*model.n_b_E + nI*model.n_b_I;     % b_E, b_I block lengths
x_cols = (len_a + len_b) + (1:model.n);      % x is the last n state columns
x_mean = mean(model.S_out(:, x_cols), 2);    % nt x 1, mean dendritic potential
t_full = model.t_out;

% Post-transient window; remove DC so we get the PSD of fluctuations
sel   = t_full > psd_t_start;
x_seg = x_mean(sel) - mean(x_mean(sel));

% Clamp Hamming window to available signal length (guards short signals)
win_len  = min(round(psd_win_len_s * model.fs), numel(x_seg));
win      = hamming(win_len);
noverlap = floor(psd_overlap_frac * win_len);

[pxx, f] = pwelch(x_seg, win, noverlap, psd_f, model.fs);

figure('Name', 'PSD of Mean Dendritic Potential');
loglog(f, pxx, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Dendritic Potential^2/Hz');
title('PSD of Mean Dendritic Potential (x)');
grid on;

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================
function [u_ex, t_ex] = dc_step_stimulus(params, T, fs, ~, input_config)
    % DC_STEP_STIMULUS Reproduces the old SRNN IED-example drive.
    %
    % Uniform tonic DC (ramped on over the first ramp_dur seconds) plus
    % elevated-drive epochs, applied identically to every neuron, over [0, T].
    % Signature matches the SRNNModel2 generator hook:
    %   [u_ex, t_ex] = generator(params, T, fs, rng_seed, input_config)
    dt   = 1 / fs;
    t_ex = (0:dt:T)';          % nt x 1, matches built-in generator
    nt   = numel(t_ex);

    DC       = input_config.DC;
    ramp_dur = input_config.ramp_dur;
    epochs   = input_config.epochs;   % rows: [t_start  t_end  amp_added]

    % Tonic DC with a linear ramp 0 -> DC over the first ramp_dur seconds
    dc_profile = DC * ones(nt, 1);
    ramp_idx   = t_ex <= ramp_dur;
    dc_profile(ramp_idx) = linspace(0, DC, nnz(ramp_idx))';

    % Add elevated-drive epochs on top of the baseline
    for k = 1:size(epochs, 1)
        in_epoch = t_ex > epochs(k, 1) & t_ex < epochs(k, 2);
        dc_profile(in_epoch) = dc_profile(in_epoch) + epochs(k, 3);
    end

    % Same drive to every neuron: n x nt
    u_ex = repmat(dc_profile', params.n, 1);
end
