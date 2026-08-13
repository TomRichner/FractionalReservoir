% example_SRNNCellTypePairs_noise.m - The same network as
% example_SRNNCellTypePairs_S_c_by_type, run as an SDE with additive Wiener
% noise on the dendritic state and integrated with Roessler SRA1.
%
% The network, seeds, stimulus and per-type setpoints are identical to that
% example. The only changes are
%
%   'ode_solver',    'sra1', ...
%   'sigma_u_noise', 0.03, ...
%
% so everything below is attributable to the noise alone.
%
% WHAT sigma_u_noise MEANS. It is INPUT-REFERRED: the same units as u, so it
% reads directly against the intrinsic_drive of 0.1 used here -- 0.03 is noise
% at 30% of the DC drive. Two Dependent views come off it:
%
%   sigma_x_raw = sigma_u_noise / tau_d          the coefficient on dW
%   x_noise_std = sigma_u_noise / sqrt(2*tau_d)  the stationary std of x it
%                                                would produce with W = 0
%
% x_noise_std is the number to compare against the sigmoid width. With
% 'piecewise' at S_a = 0.8 the nonlinearity spans S_c +/- 0.6, so a noise std of
% ~0.07 moves a neuron across roughly a tenth of its input range.
%
% WHY SRA1. Noise requires a stochastic integrator ('euler', 'heun' or 'sra1');
% asking for it with rk4/ode45 is an error rather than a silent drop. Prefer
% SRA1: it costs the same two drift evaluations per step as stochastic Heun but
% is strong order 1.5 instead of 1.0 -- measured ~85x more accurate at the same
% step size (see test_sde_integrators).
%
% WHAT THE NOISE DOES. Additive noise does not simply add jitter on top of the
% same dynamics. It pushes x_eff away from the setpoint into the flatter part of
% the sigmoid, which LOWERS the mean effective gain phi', and a lower gain means
% a lower Lyapunov exponent. That is noise-induced synchronisation. The sweep at
% the bottom measures the whole chain (3 noise realisations per level):
%
%   sigma_u   LLE (mean +/- sd)   <phi'>   std(x)   E rate
%   0.000     +0.1178 +/- 0.0000  0.2836   0.4940   0.0102
%   0.010     +0.3062 +/- 0.2781  0.2819   0.4928   0.0103
%   0.020     +0.1521 +/- 0.1397  0.2765   0.4959   0.0106
%   0.030     +0.0082 +/- 0.1205  0.2703   0.4979   0.0109
%   0.050     -0.0904 +/- 0.0130  0.2552   0.5116   0.0118
%   0.080     -0.1011 +/- 0.0026  0.2346   0.5450   0.0128
%
% Three things in that table are worth more than the headline:
%
%   * The network crosses from chaotic to contracting at sigma_u ~ 0.03. The
%     gain falls monotonically the whole way, which is the mechanism; the LLE
%     follows it.
%   * std(x) RISES while the LLE falls. A network whose LLE falls under noise
%     really has become less sensitive to its initial conditions -- two
%     trajectories driven by the SAME noise path converge -- but it has NOT
%     become quieter. "More stable" here means "not chaotic", not "better
%     behaved", and reporting the LLE without the variance would mislead.
%   * The LLE scatter COLLAPSES as sigma rises, from +/-0.28 to +/-0.003.
%     Finite-time exponents of a chaotic trajectory scatter; those of a
%     contracting one do not. This is also why the sweep averages over
%     realisations at all -- at one realisation per level the sigma = 0.01
%     point alone would read as a spurious non-monotonicity.
%
% The Benettin estimate stays valid at any noise level because the perturbed
% re-integration replays the same pre-generated Wiener increments as the
% fiducial run, so the additive noise cancels in their difference and what is
% left is the linearised stretching. Nothing about the LLE here is a
% small-noise approximation.
%
% The W plots from the S_c example are omitted: noise moves the operating
% point, not the connectivity, and W is bit-identical to the noise-free twin.
%
% See also: example_SRNNCellTypePairs_S_c_by_type, sde_fixed_step,
%           test_SRNNCellTypePairs_noise, test_sde_integrators

close all; clear; clc;

compare_to_deterministic = true;   % also run the sigma = 0 twin
run_sigma_sweep          = true;   % LLE / gain / var(x) against sigma

sigma_u_noise = 0.03;              % input-referred: 30% of the 0.1 DC drive

%% Short-term depression on E->E and I->I
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
synapse_config.I.I.std = struct('tau_rec', 4, 'tau_rel', 0.5);

%% Create model
% Identical to example_SRNNCellTypePairs_S_c_by_type except for the integrator.
% 'sra1' is used for the sigma = 0 twin as well, so the comparison isolates the
% NOISE rather than confounding it with a change of integrator (at sigma = 0
% SRA1 degenerates exactly to a deterministic two-stage RK with node 3/4).
Tend = 25;
common = { ...
    'rng_seeds', [19 20] + 7, ...
    'activation', 'piecewise', ...
    'S_a', 0.8, ...
    'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], ...
    'mu_tilde_relative',    [4.5 -3; 4 -3], ...   % multiples of F, (post <- pre)
    'sigma_tilde_relative', [1 1; 1 1], ...
    'level_of_chaos', 1.5, ...
    'n_a', [3 0], ...                     % SFA on E only, 3 timescales
    'c', [0.5/3, 0], ...
    'mu_S_c',    [0.0 0.25], ...          % E centred at 0, I at 0.25
    'sigma_S_c', [0.0 0.0], ...
    'T_range', [-5 Tend], ...
    'ode_solver', 'sra1', ...             % Roessler SRA1, strong order 1.5
    'fs', 200, ...
    'T_plot', [0 Tend], ...
    'lya_T_interval', [max(0, Tend-10), Tend], ...
    'synapse_config', synapse_config};

model = SRNNCellTypePairs(common{:}, 'sigma_u_noise', sigma_u_noise);
model.input_config.intrinsic_drive = 0.1 * ones(model.n, 1);

% plot_eigenvalues reads the Jacobian off the full trajectory, and the gain
% diagnostic below reads x from it, so the full state has to be kept.
model.store_full_state = true;

%% Build and run
model.build();
model.run();

%% Plots
model.plot();                       % compact summary: one panel per quantity
model.plot_celltypes();             % one column per type, every neuron

% Effective-Jacobian spectrum at three times. The noise moves the operating
% point sample by sample, so these spectra wander more than in the noise-free
% twin even though W is untouched.
model.plot_eigenvalues([2 10 18]);

%% Report
fprintf('\n=== Additive noise, integrated with SRA1 ===\n');
fprintf('sigma_u_noise (input-referred) = %.4f\n', model.sigma_u_noise);
fprintf('sigma_x_raw   (coefficient on dW) = %.4f\n', model.sigma_x_raw);
fprintf('x_noise_std   (nominal, W = 0)    = %.4f\n', model.x_noise_std);
fprintf('  ...as a fraction of the sigmoid half-width (S_a = %.2f): %.1f%%\n', ...
    model.S_a, 100 * model.x_noise_std / 0.6);
fprintf('intrinsic_drive                = %.4f  (noise is %.0f%% of it)\n', ...
    0.1, 100 * model.sigma_u_noise / 0.1);
fprintf('state dimension                = %d\n', model.N_sys_eqs);
fprintf('mean rate E / I                = %.4f / %.4f\n', ...
    mean(model.plot_data.r.E, 'all'), mean(model.plot_data.r.I, 'all'));
fprintf('measured std(x)                = %.4f  (nominal %.4f; they differ\n', ...
    std_of_x(model), model.x_noise_std);
fprintf('                                  because recurrence amplifies it)\n');
fprintf('mean effective gain <phi''>     = %.4f\n', mean_gain(model));
fprintf('LLE                            = %.6f\n', model.lya_results.LLE);

%% The noise-free twin, for comparison
if compare_to_deterministic
    twin = SRNNCellTypePairs(common{:}, 'sigma_u_noise', 0);
    twin.input_config.intrinsic_drive = 0.1 * ones(twin.n, 1);
    twin.store_full_state = true;
    twin.build();
    twin.run();

    fprintf('\n=== sigma_u_noise = 0 (same W, same stimulus, same integrator) ===\n');
    fprintf('W identical to above           = %d\n', ...
        max(abs(full(model.W) - full(twin.W)), [], 'all') == 0);
    fprintf('mean rate E / I                = %.4f / %.4f\n', ...
        mean(twin.plot_data.r.E, 'all'), mean(twin.plot_data.r.I, 'all'));
    fprintf('measured std(x)                = %.4f\n', std_of_x(twin));
    fprintf('mean effective gain <phi''>     = %.4f\n', mean_gain(twin));
    fprintf('LLE                            = %.6f\n', twin.lya_results.LLE);

    fprintf('\nEffect of adding sigma_u_noise = %.3f:\n', sigma_u_noise);
    fprintf('  E rate     %.4f -> %.4f\n', ...
        mean(twin.plot_data.r.E, 'all'), mean(model.plot_data.r.E, 'all'));
    fprintf('  I rate     %.4f -> %.4f\n', ...
        mean(twin.plot_data.r.I, 'all'), mean(model.plot_data.r.I, 'all'));
    fprintf('  std(x)     %.4f -> %.4f\n', std_of_x(twin), std_of_x(model));
    fprintf('  <phi''>     %.4f -> %.4f\n', mean_gain(twin), mean_gain(model));
    fprintf('  LLE        %.6f -> %.6f\n', ...
        twin.lya_results.LLE, model.lya_results.LLE);

    % One neuron, with and without noise, on the same axes. Same W, same
    % stimulus, same integrator -- the two traces separate purely because one
    % is being driven by a Wiener process.
    figure('Name', 'One E neuron with and without noise');
    plot(twin.plot_data.t, twin.plot_data.x.E(1, :), 'LineWidth', 1.2, ...
        'DisplayName', '\sigma_u = 0');
    hold on;
    plot(model.plot_data.t, model.plot_data.x.E(1, :), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('\\sigma_u = %.3f', sigma_u_noise));
    hold off; box on; grid on;
    xlabel('time (s)'); ylabel('x (E neuron 1)');
    title('Additive noise on the dendritic state');
    legend('Location', 'best');
end

%% How the exponent moves with the noise level
% The headline of this example. Every run reuses the same W, stimulus and
% setpoints; only sigma_u_noise changes.
% Several noise realisations per level, because a finite-time LLE measured on a
% single noise path over a 10 s window genuinely scatters -- more so as sigma
% rises. One realisation per level produces visible non-monotonicity that is
% sampling noise rather than physics, so the trend is only readable once it is
% averaged. This is the same reason a production sweep wants a reps axis.
if run_sigma_sweep
    sigmas = [0 0.01 0.02 0.03 0.05 0.08];
    noise_seeds = [1 2 3];
    n_s = numel(sigmas);
    lle = nan(n_s, numel(noise_seeds));
    gain = nan(n_s, numel(noise_seeds));
    xstd = nan(n_s, numel(noise_seeds));
    rateE = nan(n_s, numel(noise_seeds));

    fprintf('\n=== sigma sweep (%d noise realisations per level) ===\n', ...
        numel(noise_seeds));
    fprintf('  %-8s %-20s %-10s %-10s %-10s\n', ...
        'sigma_u', 'LLE (mean +/- sd)', '<phi''>', 'std(x)', 'E rate');
    for k = 1:n_s
        % At sigma = 0 the realisation is irrelevant, so run it once.
        seeds_here = noise_seeds;
        if sigmas(k) == 0, seeds_here = noise_seeds(1); end
        for j = 1:numel(seeds_here)
            s = SRNNCellTypePairs(common{:}, 'sigma_u_noise', sigmas(k), ...
                'noise_seed', seeds_here(j));
            s.input_config.intrinsic_drive = 0.1 * ones(s.n, 1);
            s.store_full_state = true;
            s.build();
            evalc('s.run();');          % keep the sweep output readable
            lle(k, j) = s.lya_results.LLE;
            gain(k, j) = mean_gain(s);
            xstd(k, j) = std_of_x(s);
            rateE(k, j) = mean(s.plot_data.r.E, 'all');
        end
        fprintf('  %-8.3f %+7.4f +/- %-8.4f %-10.4f %-10.4f %-10.4f\n', ...
            sigmas(k), mean(lle(k, :), 'omitnan'), std(lle(k, :), 'omitnan'), ...
            mean(gain(k, :), 'omitnan'), mean(xstd(k, :), 'omitnan'), ...
            mean(rateE(k, :), 'omitnan'));
    end

    figure('Name', 'Effect of noise level');
    tiledlayout(3, 1, 'TileSpacing', 'compact');

    nexttile;
    plot(sigmas, lle, 'o', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    hold on;
    plot(sigmas, mean(lle, 2, 'omitnan'), '-o', 'LineWidth', 1.6, ...
        'DisplayName', 'mean over realisations');
    yline(0, 'k:', 'edge of chaos', 'HandleVisibility', 'off');
    hold off; box on; grid on; ylabel('LLE');
    legend('Location', 'best');
    title('Largest Lyapunov exponent falls as noise rises (grey: individual paths)');

    nexttile;
    plot(sigmas, mean(gain, 2, 'omitnan'), '-o', 'LineWidth', 1.6);
    box on; grid on; ylabel('mean \phi''');
    title('...because the mean effective gain falls with it');

    nexttile;
    plot(sigmas, mean(xstd, 2, 'omitnan'), '-o', 'LineWidth', 1.6);
    box on; grid on; ylabel('std(x)');
    xlabel('\sigma_{u,noise}   (input-referred)');
    title('while the state itself gets NOISIER, not quieter');
end

%% Local helpers
function v = std_of_x(m)
% Std of the dendritic state over every neuron and plotted sample.
x = [];
for q = 1:m.n_cellTypes
    xq = m.plot_data.x.(m.cell_type_names{q});
    x = [x; xq(:)]; %#ok<AGROW>
end
v = std(x);
end

function g = mean_gain(m)
% Mean phi'(x_eff) over every neuron and plotted sample -- the effective gain
% the linearisation sees, and the quantity noise actually acts on.
%
% x_eff is reassembled as a full n x nt matrix in NEURON order rather than
% evaluated per cell type. With per-type setpoints S_c_vec is an n x 1 vector,
% so the derivative handle lines its centre up with x elementwise and is only
% valid for length-n input; it broadcasts correctly across n x nt, but a
% flattened per-type slice would silently pair the wrong centres.
nt = numel(m.plot_data.t);
x_eff = zeros(m.n, nt);
for q = 1:m.n_cellTypes
    name = m.cell_type_names{q};
    idx = m.type_indices{q};
    xq = m.plot_data.x.(name);
    if m.n_a(q) > 0
        xq = xq - m.c(q) * reshape(sum(m.plot_data.a.(name), 2), numel(idx), nt);
    end
    x_eff(idx, :) = xq;
end
g = mean(m.activation_function_derivative(x_eff), 'all');
end
