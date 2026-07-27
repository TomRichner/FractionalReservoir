%% panelA_bottom_traces.m
% Figure 1, Panel A (bottom): example time series of a random recurrent
% network across three levels of dynamical gain, reproducing the
% Sompolinsky-Crisanti-Sommers 1988 result (chaos in random neural networks)
% using the FractionalReservoir SRNNModel2 class.
%
% Three networks are simulated with an identical underlying random weight
% matrix (same rng seed) scaled by three gain values (level_of_chaos = gamma).
% Only the gain differs -> isolates the effect of stability/regulation:
% "same brain, different regulation".
%
% Setup:  purely random (non-E/I) Gaussian connectivity, tanh nonlinearity,
%         no SFA, no STD, N = 200, no external input.
% Spectral radius R = gamma exactly, so chaos onset is at gamma = 1.
% Gammas 0.8 / 1.5 / 2.5 span the three regimes: gamma = 0.8 (< 1) is
% subcritical and relaxes to the trivial fixed point (LLE < 0); gamma = 1.5
% is rich edge-of-chaos activity; gamma = 2.5 is wild chaos (both LLE > 0).
%
% Reference drivers (older SRNNModel class):
%   ConnectivityAdaptation/StabilityAnalysis/scripts/Sompolinsky_N_200_g_{1p8,2p1}.m
%
% Reference: Sompolinsky, H., Crisanti, A. & Sommers, H. J. Chaos in Random
% Neural Networks. Phys. Rev. Lett. 61, 259-262 (1988).

%% Setup
close all; clear; clc;
setup_paths();

this_dir = fileparts(mfilename('fullpath'));

%% ---- Editable parameters -------------------------------------------------
N            = 200;                 % network size
tau_d        = 1;                  % dendritic time constant (s), matches reference
gammas       = [0.9, 1.6, 2.5];    % level_of_chaos per panel (left -> right)
                                   % 0.8 < 1: subcritical -> relaxes to fixed point
                                   % 1.5, 2.5 > 1: chaotic (rich -> wild)
n_traces     = 15;                 % number of membrane-potential traces to plot
fs           = 200;                % sampling frequency (Hz), dt = 5 ms
T_range      = [0, 60];            % simulation interval (s)
lya_interval = [30, 60];            % Lyapunov window (skip ~5 s transient)
T_plot       = [0, 60];            % time window to display
x0_std       = 1.0;                % std of random initial membrane potential x0
                                   % (larger -> more visible relaxation in the stable panel)
net_seed     = 0;                 % shared network rng seed (same W_raw across gammas)
% --------------------------------------------------------------------------

n_cases   = numel(gammas);
results   = struct('gamma', {}, 't', {}, 'x', {}, 'LLE', {}, 'R', {}, 'W', {});

%% Simulate the three networks
for k = 1:n_cases
    gamma = gammas(k);
    fprintf('\n=== Case %d/%d : gamma = %.2f ===\n', k, n_cases, gamma);

    model = SRNNModel2();

    % Network architecture
    model.n        = N;
    model.f        = 1;                 % single population (E/I irrelevant once stats equal)
    model.indegree = N;                 % fully connected, alpha = 1
    model.tau_d    = tau_d;             % dendritic time constant (matches reference)

    % Purely random Gaussian weights: zero mean, std 1/sqrt(N), no E/I structure
    model.mu_E_tilde    = 0;
    model.mu_I_tilde    = 0;
    model.sigma_E_tilde = 1/sqrt(N);
    model.sigma_I_tilde = 1/sqrt(N);
    model.E_W           = 0;
    model.zrs_mode      = 'none';
    model.rescale_by_abscissa = false;
    model.level_of_chaos      = gamma;  % scales W by gamma -> R = gamma

    % No SFA, no STD
    model.n_a_E = 0; model.n_a_I = 0;
    model.n_b_E = 0; model.n_b_I = 0;

    % tanh activation
    model.activation_function            = @SRNNModel2.tanhActivation;
    model.activation_function_derivative = @SRNNModel2.tanhActivationDerivative;

    % Simulation settings
    model.T_range   = T_range;
    model.fs        = fs;
    model.x0_std    = x0_std;           % larger random IC -> visible relaxation
    model.rng_seeds = [net_seed, 1];    % same network seed across all gammas

    % No external input
    model.input_config.n_steps         = 1;
    model.input_config.step_density    = 0;
    model.input_config.step_density_E  = 0;
    model.input_config.step_density_I  = 0;
    model.input_config.amp             = 0;
    model.input_config.no_stim_pattern = true;
    model.input_config.intrinsic_drive = zeros(N, 1);
    model.u_ex_scale = 0;

    % Lyapunov (largest exponent via Benettin)
    model.lya_method      = 'benettin';
    model.lya_T_interval  = lya_interval;
    model.filter_local_lya = false;

    % Storage
    model.store_full_state      = false;   % benettin runs before S_out is cleared
    model.store_decimated_state = true;
    model.plot_freq             = fs;      % keep full temporal detail for plotting

    % Build + run
    model.build();
    model.run();

    LLE = model.lya_results.LLE;
    fprintf('  R (spectral radius) = %.3f | LLE = %.4f -> %s\n', ...
        model.R, LLE, ternary(LLE > 0, 'CHAOTIC', 'non-chaotic'));

    results(k).gamma = gamma;
    results(k).t     = model.plot_data.t;
    results(k).x     = model.plot_data.x.E;   % 200 x nt membrane potentials
    results(k).LLE   = LLE;
    results(k).R     = model.R;
    results(k).W     = model.W;                % gamma already folded in (W = gamma*W_raw)
end

%% Plot: 1 x 3 linked time-series panels
fig = figure('Color', 'white', 'Position', [100, 300, 1200, 320]);
tl  = tiledlayout(1, n_cases, 'TileSpacing', 'compact', 'Padding', 'compact');

idx  = round(linspace(1, N, n_traces));   % which neurons to show
cols = lines(n_traces);
ax   = gobjects(1, n_cases);

for k = 1:n_cases
    ax(k) = nexttile; hold on;
    t = results(k).t;
    x = results(k).x;
    for j = 1:n_traces
        plot(t, x(idx(j), :), 'LineWidth', 0.6, 'Color', cols(j, :));
    end
    box off;
    set(ax(k), 'Color', 'none');
    ax(k).XAxis.Visible = 'off';   % hide time axis (scale bar added instead)
    if k == 1
        ylabel('x (AU)');
    end
end

linkaxes(ax, 'xy');       % shared x and y scale across the three panels
xlim(ax(1), T_plot);
set(ax, 'YTick', [-3, 0, 3]);

% 10 s time scale bar in the lower-right of the 1st subplot
bar_len = 10;                        % seconds
xr = xlim(ax(1)); yr = ylim(ax(1));
wx = xr(2) - xr(1); wy = yr(2) - yr(1);
x_end   = xr(2) - 0.03*wx;
x_start = x_end - bar_len;
y_bar   = yr(1) + 0.10*wy;
plot(ax(1), [x_start, x_end], [y_bar, y_bar], 'k', 'LineWidth', 4);
text(ax(1), mean([x_start, x_end]), y_bar - 0.03*wy, sprintf('%d s', bar_len), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);
hold(ax(1), 'off');

%% Save (fig, png, svg via repo helper)
save_some_figs_to_folder_2(this_dir, 'panelA_bottom_traces', fig.Number, {'fig', 'png', 'svg'});

%% Eigenspectrum of the Jacobian at x = 0:  J = (-I + gamma*W)/tau_d
% (tanh'(0) = 1, so the linearization uses W directly; gamma is already in W.)
% Eigenvalues fill a disk centered at -1/tau_d with radius gamma/tau_d.
% The vertical (Im) axis at Re = 0 is the stability boundary: once the disk
% crosses it (gamma > 1) the fixed point is unstable.
%
% Styling matches ConnectivityAdaptation/RandomMatrixTheory (RMT.plot_spectrum /
% Fig_1_RMT_examples.m): black open-circle eigenvalues, solid black theoretical
% radius circle, axis off with hand-drawn Re/Im axes through the origin, equal
% aspect, common scaling, and (a)/(b)/(c) letters.
fig2  = figure('Color', 'white', 'Position', [100, 300, 977, 380]);
tl2   = tiledlayout(1, n_cases, 'TileSpacing', 'tight', 'Padding', 'compact'); %#ok<NASGU>
axe   = gobjects(1, n_cases);
theta = linspace(0, 2*pi, 200);
xc    = -1 / tau_d;                 % common disk center (shift from -I/tau_d)
mSize = 4;

% Plot eigenvalues + theoretical radius circle (RMT style)
for k = 1:n_cases
    axe(k) = nexttile; hold(axe(k), 'on');
    ev = eig((-eye(N) + results(k).W) / tau_d);
    Rk = results(k).gamma / tau_d;                 % theoretical radius
    unstable = real(ev) > 0;                        % right of the imaginary axis
    plot(axe(k), real(ev(~unstable)), imag(ev(~unstable)), 'ko', ...   % stable modes
        'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
    plot(axe(k), real(ev(unstable)), imag(ev(unstable)), 'o', ...      % unstable modes (Re>0)
        'MarkerSize', mSize, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', [0.7 0 0], 'LineWidth', 0.5);
    plot(axe(k), xc + Rk*cos(theta), Rk*sin(theta), 'k-', 'LineWidth', 2);
end

% Common scaling centered on the disk center (RMT-style equal aspect)
% margin 1.15, then zoomed out a further 20% so the Im axis/label clears the disk
max_R = max([results.gamma]) / tau_d;
common_radius = max_R * 1.15 * 1.20;
for k = 1:n_cases
    xlim(axe(k), [xc - common_radius, xc + common_radius]);
    ylim(axe(k), [-common_radius, common_radius]);
    daspect(axe(k), [1 1 1]);
end

% Format axes: hide box, draw Re/Im axes through the origin (RMT style)
for k = 1:n_cases
    axes(axe(k)); %#ok<LAXES>
    x_lim = xlim; y_lim = ylim;
    y_lim_axis = min(0.75*y_lim, 1.1*max_R);
    axis off;
    hold on;
    h_x = plot(x_lim, [0, 0], 'k', 'LineWidth', 1.5);      % Re axis
    h_y = plot([0, 0], y_lim_axis, 'k', 'LineWidth', 1.5); % Im axis = Re=0 stability boundary
    uistack([h_x, h_y], 'bottom');
    text(x_lim(2), 0, ' Re', 'Interpreter', 'latex', ...
        'VerticalAlignment', 'middle', 'FontSize', 16);
    text(0, y_lim_axis(2), 'Im', 'Interpreter', 'latex', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 16);
    xlim(x_lim); ylim(y_lim);
    hold off;
end

% (a), (b), (c) letters omitted for the eigenspectrum panels
drawnow;

save_some_figs_to_folder_2(this_dir, 'panelA_eigenspectrum', fig2.Number, {'fig', 'png', 'svg'});

%% Local helper
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
