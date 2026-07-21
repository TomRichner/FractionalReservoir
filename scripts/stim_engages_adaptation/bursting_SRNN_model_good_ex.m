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
% if exist('rng_seeds', 'var')
%     rng_seeds = rng_seeds + 1
% else
%     rng_seeds = [1 2]
% end

rng_seeds = [19 20]

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
%  6. EXTERNAL STIMULUS  (DC staircase + white-noise probe)
%  ======================================================================
%  A tonic DC applied UNIFORMLY to every neuron, stepped through a sequence of
%  levels, each HELD for hold_dur seconds, PLUS independent white noise added
%  per neuron. The staircase sets the operating point; the noise probes how the
%  network filters its input, revealed by the per-level PSD (Section 10).
%  Built by the custom generator dc_staircase_stimulus() at the BOTTOM of this
%  file, wired in via input_config.generator (overrides the built-in generator).
%
%  NOTE: SRNNModel2 forces u = 0 during the negative warmup (t < 0), so this
%  profile applies to t in [0, T_range(2)]. The ramp into the first level runs
%  over the first ramp_dur seconds of the positive window, not the warmup.
dc_levels = [0.0 0.025 0.05 0.1 0.2];   % absolute DC per level (all neurons); edit this sweep
hold_dur  = 90;                       % seconds each level is held
% White-noise INTENSITY (fs-invariant): the generator adds noise_intensity*sqrt(fs)*randn
% per neuron, so the continuous-time noise PSD (~noise_intensity^2) is independent of fs.
% Effective per-sample std = noise_intensity*sqrt(fs); 0.001*sqrt(400) = 0.02 at fs=400.
noise_intensity = 0.001;              % white-noise intensity (input units * sqrt(s)); 0 = off

input_config = struct();
input_config.dc_levels = dc_levels;   % staircase levels
input_config.hold_dur  = hold_dur;    % seconds per level
input_config.ramp_dur  = 10;          % ramp 0 -> dc_levels(1) over first ramp_dur s
input_config.noise_intensity = noise_intensity;  % fs-invariant noise intensity (to the generator)
input_config.intrinsic_drive = [];    % unused by this generator; required by the class
input_config.generator = @dc_staircase_stimulus;  % custom generator (see bottom of file)

u_ex_scale = 1.0;     % global scale on the external input (default 1.0)

%% ======================================================================
%  7. SIMULATION / INTEGRATION
%  ======================================================================
fs         = 400;          % sampling frequency (Hz)            (default 400)
% Positive window must span the whole staircase: numel(dc_levels)*hold_dur.
T_range    = [-0, numel(dc_levels)*hold_dur];   % warmup + staircase (e.g. [-30 300])
T_plot     = [];           % plotting window; [] -> T_range
% ODE solver:
%   @ode45    - adaptive RK4(5); accurate, but VERY slow with noisy forcing.
%   @ode_rk4  - fixed-step classic RK4 at the fs grid (local fn, bottom of file);
%               fast and well-behaved when white noise is added. Recommended here.
ode_solver = @ode_rk4;     % @ode45 (adaptive) | @ode_rk4 (fixed-step, fast with noise)
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

% Figure saving: when true, save all open figures (time-series + PSD) to data/
% via save_some_figs_to_folder_2 (writes .fig/.png/.pdf). Each run goes to a
% timestamped subfolder so nothing is overwritten.
save_figs  = trueor;        % set true to save figures for sharing
save_types = {'fig', 'png', 'pdf'};   % formats; pdf bundles all figs into one _report.pdf

%% ======================================================================
%  10. PSD ANALYSIS (pwelch of mean dendritic potential x, per DC level)
%  ======================================================================
%  One PSD per DC staircase level (full-resolution, fs Hz), overlaid to show
%  how the drive level shapes the burst spectrum. For each level we skip the
%  first psd_settle seconds after the DC step (settling) before the PSD window.
%  Mirrors the per-period overlay in SRNN_caller_example_tseries_raster_adapt_psd.m.
%  NOTE: the slowest adaptation timescale is ~15 s, so psd_settle ~ 20 s covers
%  ~1.3 tau; raise hold_dur / psd_settle for cleaner level separation, and raise
%  hold_dur (and psd_win_len_s) for finer low-frequency resolution (~1/win).
psd_settle       = 15;    % seconds to skip after each DC step before the PSD window
psd_win_len_s    = 10;    % Hamming window length (s)   [template used 15]
psd_overlap_frac = 0.5;  % segment overlap fraction    [template used 0.75]
psd_f            = logspace(log10(0.1), log10(100), 100);  % requested freqs (Hz)

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
%  PSD of mean dendritic potential (x), per DC level  --  pwelch, full-res
%  ======================================================================
% Extract full-resolution x (last n columns of S_out) and average over neurons.
nE = model.n_E; nI = model.n_I;
len_a = nE*model.n_a_E + nI*model.n_a_I;     % a_E, a_I block lengths
len_b = nE*model.n_b_E + nI*model.n_b_I;     % b_E, b_I block lengths
x_cols = (len_a + len_b) + (1:model.n);      % x is the last n state columns
x_mean = mean(model.S_out(:, x_cols), 2);    % nt x 1, mean dendritic potential
t_full = model.t_out;

% One PSD per DC staircase level, overlaid on a log-log axis.
nL = numel(dc_levels);
figure('Name', 'PSD of Mean Dendritic Potential vs DC level');
hold on;
cmap = parula(nL+1);
labels = cell(nL, 1);
for k = 1:nL
    % Steady window for level k: skip the first psd_settle s after the step.
    lo  = (k-1)*hold_dur + psd_settle;
    hi  = k*hold_dur;
    sel = t_full > lo & t_full <= hi;
    
    % Remove DC so we get the PSD of fluctuations
    x_seg = x_mean(sel) - mean(x_mean(sel));
    
    % Clamp Hamming window to available signal length (guards short windows)
    win_len  = min(round(psd_win_len_s * model.fs), numel(x_seg));
    win      = hamming(win_len);
    noverlap = floor(psd_overlap_frac * win_len);
    
    [pxx, fpx] = pwelch(x_seg, win, noverlap, psd_f, model.fs);
    plot(fpx, pxx, 'LineWidth', 1.5, 'Color', cmap(k, :));
    labels{k} = sprintf('DC = %.3g', dc_levels(k));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Frequency (Hz)');
ylabel('Dendritic Potential^2/Hz');
title(sprintf('PSD of Mean Dendritic Potential (x) vs DC level  (noise intensity = %.3g)', noise_intensity));
legend(labels, 'Location', 'southwest');
grid on;

%% ======================================================================
%  SHORT STAIRCASE FOR THE PAPER  (3 DC levels, 10 s holds)
%  ======================================================================
%  The staircase above is long and busy -- great for the per-level PSD, too
%  cluttered for a figure. Re-run the SAME network/adaptation/nonlinearity with
%  a compact 3-level staircase so the time-series reads cleanly in a paper.
%  Only dc_levels, hold_dur and the derived T_range change; everything else is
%  inherited from the parameters set at the top of this script.
dc_levels_short = [0.0 0.05 0.2];   % compact staircase for the paper figure
hold_dur_short  = 20;               % seconds each level is held

input_config_short = input_config;                 % inherit noise, generator, etc.
input_config_short.dc_levels = dc_levels_short;
input_config_short.hold_dur  = hold_dur_short;

T_range_short = [-0, numel(dc_levels_short)*hold_dur_short];

model_short = SRNNModel2( ...
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
    'input_config', input_config_short, 'u_ex_scale', u_ex_scale, ...
    'fs', fs, 'T_range', T_range_short, 'T_plot', T_plot, ...
    'ode_solver', ode_solver, 'rng_seeds', rng_seeds, ...
    'lya_method', lya_method, 'plot_deci', plot_deci, ...
    'store_full_state', store_full_state);

model_short.build();
model_short.run();
[fig_handle_short, ax_handles_short] = model_short.plot();
set(fig_handle_short, 'Name', 'Short staircase (paper figure)');

% Un-cap the firing-rate & synaptic-output panels (same ReLU fix as above).
for k = 1:numel(ax_handles_short)
    ylab = get(get(ax_handles_short(k), 'YLabel'), 'String');
    if any(strcmp(ylab, {'firing rate', 'synaptic output'}))
        ylim(ax_handles_short(k), 'auto');
        yticks(ax_handles_short(k), 'auto');
    end
end

%% ======================================================================
%  SAVE FIGURES  (time-series + PSD) -> data/<timestamped subfolder>
%  ======================================================================
if save_figs
    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));  % .../FractionalResevoir
    stamp       = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    save_folder = fullfile(project_root, 'data', ['bursting_SRNN_model_' stamp]);
    save_name   = sprintf('bursting_seed%d_%d', rng_seeds(1), rng_seeds(2));
    % fig_vec = [] -> save every open figure (close all at top guarantees these
    % are just the time-series and PSD figures from this run).
    save_some_figs_to_folder_2(save_folder, save_name, [], save_types);
    fprintf('Saved figures to %s\n', save_folder);
end

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================
function [u_ex, t_ex] = dc_staircase_stimulus(params, T, fs, rng_seed, input_config)
% DC_STAIRCASE_STIMULUS Uniform tonic DC stepped through a sequence of levels,
% plus independent per-neuron white noise.
%
% Each level dc_levels(k) is held for hold_dur seconds, applied identically
% to every neuron over [0, T]. The first hold ramps in linearly from 0 over
% the first ramp_dur seconds. fs-invariant white noise
% (input_config.noise_intensity) is added per neuron on top. Signature
% matches the SRNNModel2 generator hook:
%   [u_ex, t_ex] = generator(params, T, fs, rng_seed, input_config)
dt   = 1 / fs;
t_ex = (0:dt:T)';          % nt x 1, matches built-in generator
nt   = numel(t_ex);

dc_levels = input_config.dc_levels;
hold_dur  = input_config.hold_dur;
ramp_dur  = input_config.ramp_dur;
nL        = numel(dc_levels);

% Staircase: level k over [(k-1)*hold_dur, k*hold_dur)
dc_profile = zeros(nt, 1);
for k = 1:nL
    seg = t_ex >= (k-1)*hold_dur & t_ex < k*hold_dur;
    dc_profile(seg) = dc_levels(k);
end
dc_profile(t_ex >= nL*hold_dur) = dc_levels(nL);   % final boundary sample

% Ramp the first hold in linearly: 0 -> dc_levels(1) over the first ramp_dur s
ramp_idx = t_ex <= ramp_dur;
dc_profile(ramp_idx) = linspace(0, dc_levels(1), nnz(ramp_idx))';

% Same drive to every neuron: n x nt
u_ex = repmat(dc_profile', params.n, 1);

% Add independent white noise per neuron over [0, T] (probes the network's
% filtering). The sqrt(fs) factor makes this an fs-invariant white noise: the
% continuous-time PSD ~ noise_intensity^2 is independent of fs (standard
% Euler-Maruyama 1/sqrt(dt) scaling). The model's linear interpolant
% band-limits it to ~Nyquist (fs/2), flat over our <100 Hz band. Seeded for
% reproducibility.
if isfield(input_config, 'noise_intensity') && input_config.noise_intensity > 0
    rng(rng_seed);
    u_ex = u_ex + input_config.noise_intensity * sqrt(fs) * randn(params.n, numel(t_ex));
end
end

function [t, Y] = ode_rk4(odefun, tspan, y0, ~)
% ODE_RK4 Fixed-step classic RK4, matching the @ode45 call signature used by
% SRNNModel2: solver(rhs, t_ex, S0, opts). Steps at the native spacing of the
% supplied time vector tspan (uniform fs grid) and returns the solution at
% exactly those times, so the class's output-time check passes. opts ignored.
%
% Much faster than adaptive ode45 when the forcing is noisy (no step-size
% control thrashing). rhs is evaluated 4x per step.
t  = tspan(:);
nt = numel(t);
y  = y0(:);
Y  = zeros(nt, numel(y));
Y(1, :) = y.';
for k = 1:nt-1
    h  = t(k+1) - t(k);
    tk = t(k);
    k1 = odefun(tk,       y);
    k2 = odefun(tk + h/2, y + (h/2)*k1);
    k3 = odefun(tk + h/2, y + (h/2)*k2);
    k4 = odefun(tk + h,   y + h*k3);
    y  = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    Y(k+1, :) = y.';
end
end
