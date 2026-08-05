% bursting_SRNN_model_good_ex_piecewise.m - piecewise-sigmoid variant of
% bursting_SRNN_model_good_ex.m
%
% Identical to bursting_SRNN_model_good_ex.m except the nonlinearity phi is the
% bounded piecewiseSigmoid(x, S_a, S_c) instead of the unbounded ReLU, matching
% the class defaults' family of bounded activations (see test_SRNN2_defaults.m).
% Because phi is now bounded in [0,1], the firing-rate / synaptic-output panel
% ylim([0 1]) set by SRNNModel2.plot is correct and is left alone.
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

% rng_seeds = [19 20]
rng_seeds = [42 42];

clearvars -except rng_seeds;

% Add src/ and scripts/ to the MATLAB path
setup_paths();

%% ======================================================================
%  1. NETWORK ARCHITECTURE
%  ======================================================================
n          = 50;      % total number of neurons              (default 300; old bursting net = 10)
f          = 0.7;     % fraction excitatory, n_E = round(f*n) (default 0.5; old EI = 0.7)
indegree   = 10;       % expected in-degree -> alpha = indegree/n (default 100; old mean degree ~4)
check_connectivity = true;

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
mu_I_tilde    = -2*F;   % normalized I mean      (class default -4*F)
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
tau_d = 0.1;   % dendritic time constant (s)   (default 0.1; old bursting = 0.025)

% --- Nonlinearity phi ---------------------------------------------------
% 'piecewise' : piecewiseSigmoid(x, S_a, S_c)  -- bounded [0,1], parameterized below
% 'logistic'  : logisticSigmoid(x, S_c)        -- bounded [0,1], centered on S_c
% 'tanh'      : tanhActivation(x)              -- bounded [-1,1]
% 'relu'      : max(0, x)  -- UNBOUNDED, matches old SRNN model
nonlinearity = 'piecewise';
% nonlinearity = 'logistic';

S_a = 0.9;       % piecewiseSigmoid slope param a (default 0.9)
% S_c is the center param c, shared by BOTH piecewise and logistic (class default 0.4).
% For the logistic it sets the operating point: phi(S_c)=0.5 with unit slope there, so
% the resting rate at x=0 is 1/(1+exp(4*S_c)) -- 0.48 -> ~0.13, 0.35 -> ~0.20, 1.0 -> ~0.02.
% Larger S_c => lower baseline firing and lower gain near rest.
S_c = 0.45;      % activation center param c (default 0.35)

switch lower(nonlinearity)
    case 'piecewise'
        phi      = @(x) SRNNModel2.piecewiseSigmoid(x, S_a, S_c);
        phi_deriv = @(x) SRNNModel2.piecewiseSigmoidDerivative(x, S_a, S_c);
    case 'logistic'
        phi      = @(x) SRNNModel2.logisticSigmoid(x, S_c);
        phi_deriv = @(x) SRNNModel2.logisticSigmoidDerivative(x, S_c);
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
dc_levels = [0.0 0.4];   % absolute DC per level (all neurons); edit this sweep
hold_dur  = 50; % 45                      % seconds each level is held
% White-noise INTENSITY (fs-invariant): the generator adds noise_intensity*sqrt(fs)*randn
% per neuron, so the continuous-time noise PSD (~noise_intensity^2) is independent of fs.
% Effective per-sample std = noise_intensity*sqrt(fs); 0.001*sqrt(400) = 0.02 at fs=400.
noise_intensity = 0.001;     % 0.001         % white-noise intensity (input units * sqrt(s)); 0 = off

input_config = struct();
input_config.dc_levels = dc_levels;   % staircase levels
input_config.hold_dur  = hold_dur;    % seconds per level
input_config.ramp_dur  = 0;          % ramp 0 -> dc_levels(1) over first ramp_dur s
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
save_figs  = true;        % set true to save figures for sharing
save_types = {'png'};   % formats; pdf bundles all figs into one _report.pdf

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
psd_f            = logspace(log10(0.3), log10(100), 100);  % requested freqs (Hz)

%% ======================================================================
%  11. POPULATION SYNCHRONY  (chi^2, Golomb-Rinzel), per DC level
%  ======================================================================
%  For a population matrix M (n_neurons x nt) over a settled window, with
%  population mean R(t) = mean_i M(i,t):
%
%      chi2 = var_t(R) / mean_i( var_t(M_i) )
%
%  i.e. the variance of the average divided by the average of the variances.
%  The only difference between numerator and denominator is WHEN you average
%  across neurons -- before or after taking the variance over time.
%
%    chi2 = 1   -> every neuron follows the same trace; averaging across
%                  neurons cancels nothing (perfect synchrony).
%    chi2 = 1/N -> neurons fluctuate independently; averaging N independent
%                  signals shrinks variance by 1/N (perfect asynchrony).
%
%  chi2 IS the mean pairwise correlation plus a 1/N floor:
%      chi2 = 1/N + (1 - 1/N)*rho_bar   =>   rho_bar = (chi2 - 1/N)/(1 - 1/N)
%  so rho_bar is the floor-corrected, N-independent form. Both are reported.
%
%  Computed on the E population only: with f = 0.7 and adaptation on E alone
%  (n_a_E = 3, n_a_I = 0), E and I cells do structurally different things and
%  blending them muddies the measure.
%
%  Computed on BOTH x (dendritic potential) and r (firing rate). x is primary:
%  in the burst regime r clips at phi's ceiling of 1, and saturation distorts
%  variance. x is unbounded and shows the burst cleanly.
%
%  Computed BOTH raw and with a per-neuron linear detrend. tau_a_E reaches 15 s,
%  so the network can still be drifting inside a hold; that drift is shared
%  across neurons (common-mode) and spuriously inflates chi2 in the suppressed
%  regime. Comparing the two lines shows how big that confound is. NOTE this is
%  an affine trend removal, NOT a bandpass -- no filtering is applied here.
chi2_settle    = psd_settle;  % seconds to skip after each DC step; matched to the PSD
                              % window so both analyses describe the same epoch
chi2_var_floor = 1e-24;   % denominator guard; below this chi2 -> NaN instead of 0/0
                          % (matters only if noise_intensity = 0, where a suppressed
                          %  fixed point has genuinely zero variance)

%% ======================================================================
%  12. CONSTRUCT, BUILD, RUN, PLOT
%  ======================================================================
model = SRNNModel2( ...
    'n', n, 'f', f, 'indegree', indegree, 'check_connectivity',check_connectivity, ...
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
% phi is bounded in [0,1] here, so the class's ylim([0 1]) on the firing-rate
% and synaptic-output panels is left as-is (the ReLU version re-autoscaled it).

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
figure('Name', 'PSD of Mean Dendritic Potential vs DC level', ...
    'Position', [20   297   505   418]);
hold on;

% Only show a subset of DC levels (the first = no stim, and the last = the
% strongest stim), but keep the parula(nL+1) colors indexed by the original
% level number so the traces match the full-staircase color scheme.
% Derived from dc_levels so editing the staircase can't ask for a level that
% was never simulated (a stale index yields an empty window and a pwelch error).
cmap1 = parula(5+1);
cmap = [cmap1(1,:); cmap1(5,:)];
psd_show_levels = [1 nL];
psd_show_labels = {'no-stim', 'stim'};   % legend text for each shown level
% Label each shown trace by its DC level, so the legend stays in sync with
% whatever staircase / subset of levels is being plotted.
% cmap = parula(nL+1);
% psd_show_levels = [1:nL];
% psd_show_labels = arrayfun(@(v) sprintf('DC = %g', v), ...
%     dc_levels(psd_show_levels), 'UniformOutput', false);

labels = cell(numel(psd_show_levels), 1);
for ii = 1:numel(psd_show_levels)
    k = psd_show_levels(ii);
    % Steady window for level k: skip the first psd_settle s after the step.
    lo  = (k-1)*hold_dur + psd_settle;
    hi  = k*hold_dur;
    sel = t_full > lo & t_full <= hi;
    
    % Remove DC so we get the PSD of fluctuations
    x_seg = x_mean(sel) - mean(x_mean(sel));
    
    % Clamp Hamming window to available signal length (guards short windows)
    win_len  = min(round(psd_win_len_s * model.fs), numel(x_seg));
    win      = hann(win_len);
    noverlap = floor(psd_overlap_frac * win_len);
    
    [pxx, fpx] = pwelch(x_seg, win, noverlap, psd_f, model.fs);
    plot(fpx, pxx, 'LineWidth', 1.5, 'Color', cmap(k, :));
    labels{ii} = psd_show_labels{ii};
end
hold off;
psd_ax = gca;
tick_fs = 14;
label_fs = 15.4;
legend_fs = 12.6;
xlim([0.2 90])
set(psd_ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', tick_fs);
xlim(psd_ax, [0.2 90]);
xlabel(psd_ax, 'Frequency (Hz)', 'FontSize', label_fs);
ylabel(psd_ax, 'PSD $x^2$/Hz', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
legend(psd_ax, labels, 'Location', 'northeast', 'Box', 'off', ...
    'FontSize', legend_fs);
grid(psd_ax, 'off');

%% ======================================================================
%  POPULATION SYNCHRONY chi^2 vs DC level
%  ======================================================================
% See Section 11 for the definition and the reasoning behind the four variants.
% Unpack the full-resolution state once. unpack_and_compute_states returns
% x/r as structs with .E (n_E x nt) and .I (n_I x nt) -- neurons in ROWS, time
% in COLUMNS -- so variance over time is var(M,0,2) and the population mean is
% mean(M,1).
[x_st, ~, ~, r_st, ~] = SRNNModel2.unpack_and_compute_states(model.S_out, model.get_params());
xE = x_st.E;    % n_E x nt, dendritic potential
rE = r_st.E;    % n_E x nt, firing rate

chi2_x    = nan(nL, 1);  rho_x    = nan(nL, 1);
chi2_xd   = nan(nL, 1);  rho_xd   = nan(nL, 1);
chi2_r    = nan(nL, 1);  rho_r    = nan(nL, 1);
chi2_rd   = nan(nL, 1);  rho_rd   = nan(nL, 1);
varR_x    = nan(nL, 1);  meanvar_x = nan(nL, 1);
varR_r    = nan(nL, 1);  meanvar_r = nan(nL, 1);
mean_rate = nan(nL, 1);

for k = 1:nL
    % Steady window for level k: skip the first chi2_settle s after the step.
    lo  = (k-1)*hold_dur + chi2_settle;
    hi  = k*hold_dur;
    sel = t_full > lo & t_full <= hi;

    xk = xE(:, sel);
    rk = rE(:, sel);
    % Per-neuron linear detrend (detrend works column-wise -> transpose).
    xkd = detrend(xk', 'linear')';
    rkd = detrend(rk', 'linear')';

    [chi2_x(k),  rho_x(k),  varR_x(k), meanvar_x(k)] = population_chi2(xk,  chi2_var_floor);
    [chi2_xd(k), rho_xd(k), ~, ~]                    = population_chi2(xkd, chi2_var_floor);
    [chi2_r(k),  rho_r(k),  varR_r(k), meanvar_r(k)] = population_chi2(rk,  chi2_var_floor);
    [chi2_rd(k), rho_rd(k), ~, ~]                    = population_chi2(rkd, chi2_var_floor);

    % Mean E firing rate: distinguishes suppression-to-quiescence from
    % suppression-by-saturation (both read as "no bursts", very different states).
    mean_rate(k) = mean(rk(:));
end

chi2_floor = 1 / model.n_E;   % asynchronous floor for the E population

% --- Figure: chi^2 vs DC, plus the numerator/denominator that produce it ----
figure('Name', 'Population synchrony (chi^2) vs DC level', ...
    'Position', [540 297 505 560]);

ax_chi = subplot(2, 1, 1);
hold(ax_chi, 'on');
plot(ax_chi, dc_levels, chi2_x,  '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(ax_chi, dc_levels, chi2_xd, '--o', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(ax_chi, dc_levels, chi2_r,  '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(ax_chi, dc_levels, chi2_rd, '--s', 'LineWidth', 1.5, 'MarkerSize', 6);
yline(ax_chi, chi2_floor, ':k', 'LineWidth', 1.5);
hold(ax_chi, 'off');
ylim(ax_chi, [0 1]);
set(ax_chi, 'FontSize', tick_fs);
xlabel(ax_chi, 'DC level', 'FontSize', label_fs);
ylabel(ax_chi, '$\chi^2$', 'Interpreter', 'latex', 'FontSize', label_fs);
legend(ax_chi, {'x (raw)', 'x (detrended)', 'r (raw)', 'r (detrended)', ...
    sprintf('async floor 1/n_E = %.3f', chi2_floor)}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', legend_fs);
grid(ax_chi, 'off');

% Diagnostic panel on the FIRING RATE r: the numerator and denominator that
% produce the chi2_r curve above. If the two collapse in lockstep the ratio
% cannot move, which is exactly why chi2 is blind to burst suppression here --
% the transition is sync -> low-variance, not sync -> async.
ax_var = subplot(2, 1, 2);
hold(ax_var, 'on');
plot(ax_var, dc_levels, varR_r,    '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(ax_var, dc_levels, meanvar_r, '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
hold(ax_var, 'off');
set(ax_var, 'YScale', 'log', 'FontSize', tick_fs);
xlabel(ax_var, 'DC level', 'FontSize', label_fs);
ylabel(ax_var, 'variance (r)', 'FontSize', label_fs);
legend(ax_var, {'var_t(R) [numerator]', '\langle var_t(r_i) \rangle_i [denominator]'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', legend_fs);
grid(ax_var, 'off');

% --- Console table ---------------------------------------------------------
fprintf('\nPopulation synchrony (E cells, n_E = %d, window = last %g s of each %g s hold)\n', ...
    model.n_E, hold_dur - chi2_settle, hold_dur);
fprintf('%8s %9s %11s %9s %11s %11s %11s %11s %11s %8s\n', ...
    'DC', 'chi2_x', 'chi2_x_dtr', 'chi2_r', 'chi2_r_dtr', ...
    'varR_x', 'mvar_x', 'varR_r', 'mvar_r', 'mean_r');
for k = 1:nL
    fprintf('%8.4g %9.4f %11.4f %9.4f %11.4f %11.3e %11.3e %11.3e %11.3e %8.4f\n', ...
        dc_levels(k), chi2_x(k), chi2_xd(k), chi2_r(k), chi2_rd(k), ...
        varR_x(k), meanvar_x(k), varR_r(k), meanvar_r(k), mean_rate(k));
end
fprintf('\n');

% %% ======================================================================
% %  SHORT STAIRCASE FOR THE PAPER  (3 DC levels, 20 s holds)
% %  ======================================================================
% % The staircase above is long and busy -- useful for the per-level PSD, but
% % too cluttered for a manuscript time-series figure. Re-run the same network
% % configuration with a compact three-level staircase. Only the DC levels,
% % hold duration, and derived simulation range change.
% dc_levels_short = [0.0 0.05 0.2];
% hold_dur_short  = 20;
% 
% input_config_short = input_config;
% input_config_short.dc_levels = dc_levels_short;
% input_config_short.hold_dur  = hold_dur_short;
% 
% T_range_short = [0, numel(dc_levels_short) * hold_dur_short];
% 
% model_short = SRNNModel2( ...
%     'n', n, 'f', f, 'indegree', indegree, ...
%     'mu_E_tilde', mu_E_tilde, 'mu_I_tilde', mu_I_tilde, ...
%     'sigma_E_tilde', sigma_E_tilde, 'sigma_I_tilde', sigma_I_tilde, ...
%     'E_W', E_W, 'zrs_mode', zrs_mode, ...
%     'level_of_chaos', level_of_chaos, 'rescale_by_abscissa', rescale_by_abscissa, ...
%     'n_a_E', n_a_E, 'n_a_I', n_a_I, 'c_E', c_E, 'c_I', c_I, ...
%     'tau_a_E', tau_a_E, 'tau_a_I', tau_a_I, ...
%     'n_b_E', n_b_E, 'n_b_I', n_b_I, ...
%     'tau_b_E_rec', tau_b_E_rec, 'tau_b_E_rel', tau_b_E_rel, ...
%     'tau_b_I_rec', tau_b_I_rec, 'tau_b_I_rel', tau_b_I_rel, ...
%     'tau_d', tau_d, 'S_a', S_a, 'S_c', S_c, ...
%     'activation_function', phi, 'activation_function_derivative', phi_deriv, ...
%     'input_config', input_config_short, 'u_ex_scale', u_ex_scale, ...
%     'fs', fs, 'T_range', T_range_short, 'T_plot', T_plot, ...
%     'ode_solver', ode_solver, 'rng_seeds', rng_seeds, ...
%     'lya_method', lya_method, 'plot_deci', plot_deci, ...
%     'store_full_state', store_full_state);
% 
% model_short.build();
% model_short.run();
% [fig_handle_short, ax_handles_short] = model_short.plot();
% set(fig_handle_short, 'Name', 'Short staircase (paper figure)');

%% ======================================================================
%  SAVE FIGURES  (time-series + PSD) -> this manuscript figure folder
%  ======================================================================
if save_figs
    save_folder = fileparts(mfilename('fullpath'));
    save_name   = sprintf('bursting_piecewise_seed%d_%d', rng_seeds(1), rng_seeds(2));
    % fig_vec = [] -> save every open figure (close all at top guarantees these
    % are just the time-series and PSD figures from this run).
    save_some_figs_to_folder_2(save_folder, save_name, [], save_types);
    fprintf('Saved figures to %s\n', save_folder);
end

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================
function [chi2, rho, varR, meanvar] = population_chi2(M, var_floor)
    % POPULATION_CHI2 Golomb-Rinzel population synchrony measure.
    %
    %   chi2 = var_t( mean_i M ) / mean_i( var_t(M_i) )
    %
    % Variance of the population average over the average of the per-neuron
    % variances. Equals 1 for perfectly synchronous neurons (averaging cancels
    % nothing) and 1/N for independent neurons (averaging N independent signals
    % shrinks variance by 1/N).
    %
    % Inputs:
    %   M         - n_neurons x nt population matrix (neurons in ROWS)
    %   var_floor - if the denominator falls below this, return NaN rather than
    %               a meaningless 0/0 (a noiseless fixed point has zero variance)
    %
    % Outputs:
    %   chi2    - synchrony measure, in [1/N, 1]
    %   rho     - floor-corrected form (chi2 - 1/N)/(1 - 1/N), the mean pairwise
    %             correlation; independent of N, so comparable across net sizes
    %   varR    - numerator,   var_t(mean_i M)
    %   meanvar - denominator, mean_i var_t(M_i)
    N       = size(M, 1);
    varR    = var(mean(M, 1), 0, 2);
    meanvar = mean(var(M, 0, 2));

    if meanvar < var_floor
        chi2 = NaN;
    else
        chi2 = varR / meanvar;
    end
    rho = (chi2 - 1/N) / (1 - 1/N);
end

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
