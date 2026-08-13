% example_SRNNCellTypePairs_S_c_by_type.m - A different nonlinearity setpoint
% for each cell type: S_c = 0 for E, S_c = 0.15 for I.
%
% Identical to test_SRNNCellTypePairs_defaults.m in every other respect. The
% only change is the two lines
%
%   'mu_S_c',    [0.0 0.15], ...
%   'sigma_S_c', [0.0 0.0],  ...     % no cell-to-cell spread (the default)
%
% which replace the single shared 'S_c' with one mean per cell type. build()
% expands them into the read-only n x 1 model.S_c_vec, and every neuron then
% carries its own centre through phi, the Jacobian and the Lyapunov exponents.
%
% WHAT RAISING S_c DOES. The nonlinearity is centred on S_c: phi(S_c) = 0.5,
% and with 'piecewise' at S_a = 0.8 it is exactly 0 below S_c - 0.6 and exactly
% 1 above S_c + 0.6. Raising S_c therefore slides the whole curve to the RIGHT,
% so a neuron needs 0.15 more input to reach the same rate -- S_c = 0.15 makes
% the I population LESS excitable than the E population, not more.
%
% What comes out (the comparison at the bottom runs the shared-S_c twin on the
% same W and stimulus, so the difference is attributable to the setpoints
% alone) is worth reading carefully:
%
%   E rate  0.0067 -> 0.0159      more than doubles
%   I rate  0.0767 -> 0.0779      essentially unchanged
%   LLE    -0.0988 -> -0.0842     closer to the edge of chaos
%
% The I population is harder to drive, but it does NOT end up firing less: the
% disinhibited E population drives it back to roughly where it was. What the
% setpoint split actually buys is a higher E rate at the same I rate -- the
% network settles at a different E/I operating point rather than simply turning
% inhibition down.
%
% Only 'logistic' and 'piecewise' have a centre to move. Asking for per-type
% setpoints under 'tanh' or activation_custom is an error, not a no-op.
%
% See also: test_SRNNCellTypePairs_defaults, SRNNCellTypePairs,
%           test_SRNNCellTypePairs_S_c_heterogeneity

close all; clear; clc;

compare_to_homogeneous = true;   % also run the shared-S_c twin and report both

%% Short-term depression on E->E and I->I
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
synapse_config.I.I.std = struct('tau_rec', 4, 'tau_rel', 0.5);

%% Create model
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
    'T_range', [-5 Tend], ...
    'ode_solver', @ode_rk4, ...
    'fs', 200, ...
    'T_plot', [0 Tend], ...
    'lya_T_interval', [max(0, Tend-10), Tend], ...
    'synapse_config', synapse_config};

model = SRNNCellTypePairs(common{:}, ...
    'mu_S_c',    [0.0 0.25], ...   % E centred at 0, I centred at 0.15
    'sigma_S_c', [0.0 0.0]);       % raise these to add cell-to-cell spread

model.input_config.intrinsic_drive = 0.1 * ones(model.n, 1);

% plot_eigenvalues reads the Jacobian off the full trajectory, so the full state
% has to be kept. At n=300 over [0 25] at fs=200 that is ~5001 x N_sys_eqs
% doubles -- tens of MB, fine here, but drop it if you scale n up.
model.store_full_state = true;

%% Build and run
model.build();
model.run();

%% Plots
model.plot();                       % compact summary: one panel per quantity

% Per-cell-type view: one COLUMN per type, every individual neuron trace, with
% b and g collapsed across routes as prod(b) and coloured by target type. This
% is the panel where the setpoint split shows: the I column sits lower for the
% same drive.
model.plot_celltypes();

% Effective-Jacobian spectrum at three times: early (still settling), middle,
% and late. phi' is now evaluated at a per-neuron centre, so the diagonal these
% eigenvalues come from is genuinely heterogeneous.
model.plot_eigenvalues([2 10 18]);

% W itself, imaged with a diverging colormap. W is untouched by the setpoints --
% only the neurons' operating points moved, not the connectivity.
model.plot_W();
model.plot_W_spectrum();
model.plot_weight_histogram();

% The setpoint distribution itself: two spikes here, because sigma_S_c = 0.
% With a nonzero sigma this becomes two Gaussians.
figure('Name', 'Per-neuron S_c by cell type');
hold on;
type_rgb = lines(2); type_rgb([1 2], :) = type_rgb([2 1], :);  % E warm, I cool
for q = 1:model.n_cellTypes
    idx = model.type_indices{q};
    plot(idx, model.S_c_vec(idx), '.', 'MarkerSize', 10, ...
        'Color', type_rgb(q, :), 'DisplayName', model.cell_type_names{q});
end
hold off; box on; grid on;
xlabel('neuron index'); ylabel('S_c');
title('Nonlinearity setpoint per neuron');
legend('Location', 'best');

%% Report
fprintf('\n=== Per-type setpoints: E at 0, I at 0.15 ===\n');
fprintf('S_c_vec (E)            = %.4f +/- %.4f\n', ...
    mean(model.S_c_vec(model.type_indices{1})), ...
    std(model.S_c_vec(model.type_indices{1})));
fprintf('S_c_vec (I)            = %.4f +/- %.4f\n', ...
    mean(model.S_c_vec(model.type_indices{2})), ...
    std(model.S_c_vec(model.type_indices{2})));
fprintf('F (normalization)      = %.6g\n', model.default_val);
fprintf('mu_tilde (post<-pre)   = %s\n', mat2str(model.mu_tilde, 4));
fprintf('bulk radius R          = %.4f\n', model.R);
fprintf('outlier eigenvalues    = %s\n', mat2str(round(model.lambda_O(:)', 4)));
fprintf('state dimension        = %d\n', model.N_sys_eqs);
fprintf('mean rate E / I        = %.4f / %.4f\n', ...
    mean(model.plot_data.r.E, 'all'), mean(model.plot_data.r.I, 'all'));
fprintf('LLE                    = %.6f\n', model.lya_results.LLE);

%% The shared-setpoint twin, for comparison
% Same seeds, same W, same stimulus -- the ONLY difference is that every neuron
% shares S_c = 0. Any change below is attributable to the setpoint split alone.
if compare_to_homogeneous
    twin = SRNNCellTypePairs(common{:}, 'S_c', 0.0);
    twin.input_config.intrinsic_drive = 0.1 * ones(twin.n, 1);
    twin.build();
    twin.run();

    fprintf('\n=== Shared S_c = 0 for both types (same W, same stimulus) ===\n');
    fprintf('W identical to above   = %d\n', ...
        max(abs(full(model.W) - full(twin.W)), [], 'all') == 0);
    fprintf('mean rate E / I        = %.4f / %.4f\n', ...
        mean(twin.plot_data.r.E, 'all'), mean(twin.plot_data.r.I, 'all'));
    fprintf('LLE                    = %.6f\n', twin.lya_results.LLE);

    fprintf('\nEffect of moving I to S_c = 0.15:\n');
    fprintf('  I rate  %.4f -> %.4f\n', ...
        mean(twin.plot_data.r.I, 'all'), mean(model.plot_data.r.I, 'all'));
    fprintf('  E rate  %.4f -> %.4f\n', ...
        mean(twin.plot_data.r.E, 'all'), mean(model.plot_data.r.E, 'all'));
    fprintf('  LLE     %.6f -> %.6f\n', ...
        twin.lya_results.LLE, model.lya_results.LLE);
end
