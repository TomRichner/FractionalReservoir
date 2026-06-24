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

close all; clear; clc;

% Add src/ and scripts/ to the MATLAB path
setup_paths();

%% ======================================================================
%  1. NETWORK ARCHITECTURE
%  ======================================================================
n          = 10;      % total number of neurons              (default 300; old bursting net = 10)
f          = 0.7;     % fraction excitatory, n_E = round(f*n) (default 0.5; old EI = 0.7)
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
mu_I_tilde    = -4*F;   % normalized I mean      (class default -4*F)
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
tau_b_E_rec = 1;     % E recovery time constant (s)      (default 1)
tau_b_E_rel = 0.25;  % E release  time constant (s)      (default 0.25)
tau_b_I_rec = 1;     % I recovery time constant (s)      (default 1)
tau_b_I_rel = 0.25;  % I release  time constant (s)      (default 0.25)

%% ======================================================================
%  5. INTRINSIC DYNAMICS & NONLINEARITY
%  ======================================================================
%  dx_i/dt = (-x_i + sum_j w_ij r_j + u_i)/tau_d ;  r_i = b_i * phi(x_i - c*sum_k a_k)
tau_d = 0.1;     % dendritic time constant (s)   (default 0.1)

% --- Nonlinearity phi ---------------------------------------------------
% 'piecewise' : piecewiseSigmoid(x, S_a, S_c)  -- default, parameterized below
% 'logistic'  : logisticSigmoid(x)
% 'tanh'      : tanhActivation(x)
nonlinearity = 'piecewise';
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
    otherwise
        error('Unknown nonlinearity: %s', nonlinearity);
end

%% ======================================================================
%  6. EXTERNAL STIMULUS  (sparse step input that "engages" adaptation)
%  ======================================================================
%  Built by generate_external_input from this struct. The default is a few
%  steps where a fraction of E neurons get a step amplitude, alternating with
%  silent windows (no_stim_pattern) so you can watch the network burst & relax.
input_config = struct();
input_config.n_steps        = 3;            % # stimulus steps                 (default 3)
input_config.step_density   = 0.2;          % legacy overall density           (default 0.2)
input_config.amp            = 0.5;          % step amplitude                   (default 0.5)
input_config.no_stim_pattern = [true false true];  % which steps are silent    (default alternating)
input_config.positive_only  = false;        % allow +/- amplitudes             (default false)
input_config.step_density_E = 0.15;         % fraction of E neurons driven     (default 0.15)
input_config.step_density_I = 0;            % fraction of I neurons driven      (default 0)
input_config.intrinsic_drive = [];          % per-neuron constant drive (n x 1); [] -> zeros

u_ex_scale = 1.0;     % global scale on the external input (default 1.0)

%% ======================================================================
%  7. SIMULATION / INTEGRATION
%  ======================================================================
fs         = 400;          % sampling frequency (Hz)            (default 400)
T_range    = [0, 50];      % integration interval [start end] s (default [0 50])
T_plot     = [];           % plotting window; [] -> T_range
ode_solver = @ode45;       % ODE solver handle                  (default @ode45)
rng_seeds  = [1 2];        % [network seed, stimulus seed]      (default [1 2])

%% ======================================================================
%  8. LYAPUNOV / ANALYSIS
%  ======================================================================
%  'none' is fastest while tuning; 'benettin' adds LLE; 'qr' = full spectrum.
lya_method = 'none';       % 'none' | 'benettin' | 'qr'         (default 'benettin')

%% ======================================================================
%  9. CONSTRUCT, BUILD, RUN, PLOT
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
    'lya_method', lya_method);

model.build();
model.run();
model.plot();
