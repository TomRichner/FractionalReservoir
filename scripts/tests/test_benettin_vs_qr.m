% test_benettin_vs_qr.m - Compare the Benettin LLE and the QR largest exponent
% on the SAME SRNNModel2 network (SFA + STD).
%
% The two Lyapunov methods are independent algorithms:
%   'benettin' - shadow-trajectory rescaling  -> lya_results.LLE
%   'qr'       - variational flow + QR (uses the analytic Jacobian)
%                -> lya_results.LE_spectrum (sorted descending; (1) is largest)
% On the same fiducial trajectory their largest exponents should agree closely.
% This is an end-to-end cross-check of the corrected QR Jacobian.
%
% QR cost scales ~ N_states^2 and is infeasible at the default n=300, so this
% test uses a small network at the default connection regime (alpha=1/3):
% n=40, indegree=20. Everything else is at SRNNModel2 defaults (fs=400,
% ode_solver=@ode45, f=0.5, rng_seeds=[1 2]). ode45 is required: the QR
% variational step integrates on a 2-point span that ode_rk4 rejects.
%
% T_range starts negative ([-10, 15]) so the settling transient happens at
% t<0. Both methods accumulate over lya_T_interval, which defaults here to
% [T_range(1)+15, T_range(2)] = [5, 15], and both iterate for lya_warmup
% seconds before that so the perturbation / basis can align first.
%
% lya_warmup is raised to 10 s here rather than left at the class default of
% 5. This test cross-validates two independent algorithms, so both must be run
% at CONVERGED settings or it measures warmup bias instead of method
% agreement. On this network QR is still ~12% from its plateau at 5 s of
% warmup while Benettin is within 0.2% (see the table on SRNNModel2.lya_warmup);
% at 10 s both have plateaued.

close all; clear; clc;

%% Shared configuration
% Same rng_seeds (default [1 2]) => identical W, stimulus, and fiducial
% trajectory for both models; only the Lyapunov computation differs.
common = {'n', 40, 'indegree', 20, 'n_a_E', 3, 'n_b_E', 1, ...
    'T_range', [-10, 15], 'lya_warmup', 10};

%% Benettin
fprintf('\n=== Benettin (largest Lyapunov exponent) ===\n');
model_b = SRNNModel2(common{:}, 'lya_method', 'benettin');
model_b.build();
model_b.run();

%% QR (full spectrum)
fprintf('\n=== QR (full Lyapunov spectrum) ===\n');
model_q = SRNNModel2(common{:}, 'lya_method', 'qr');
model_q.build();
model_q.run();

%% Compare largest exponents
LLE_benettin = model_b.lya_results.LLE;
LE_spectrum  = model_q.lya_results.LE_spectrum;   % sorted descending
LLE_qr       = LE_spectrum(1);

abs_diff = abs(LLE_benettin - LLE_qr);
rel_diff = abs_diff / max(abs(LLE_benettin), eps);
tol = 0.05;   % documented, generous tolerance for two independent methods

n_show = min(20, numel(LE_spectrum));
fprintf('\n========================================\n');
fprintf('Benettin vs QR largest Lyapunov exponent\n');
fprintf('========================================\n');
fprintf('  Benettin LLE          : %+.4f\n', LLE_benettin);
fprintf('  QR largest exponent   : %+.4f\n', LLE_qr);
fprintf('  |difference|          : %.4f\n', abs_diff);
fprintf('  relative difference   : %.1f%%\n', 100 * rel_diff);
fprintf('  QR spectrum (top %d)   : %s\n', n_show, mat2str(LE_spectrum(1:n_show)', 4));
if abs_diff < tol
    fprintf('  RESULT: PASS (|diff| < %.3f)\n', tol);
else
    fprintf('  RESULT: FAIL (|diff| >= %.3f)\n', tol);
end
fprintf('========================================\n');

%% Full time-series plots (same panels as test_SRNN2_defaults)
% Input, dendritic potential, firing rates, adaptation, depression, synaptic
% output, and the local Lyapunov exponent. The fiducial trajectory is identical
% for both models (same seeds), so the Benettin model's panels represent both.
model_b.plot();

%% Overlay: Benettin local LLE vs the QR top-K local exponents
% The two methods use different rescaling intervals (lya_dt 0.02 vs 0.1), so
% each local series is plotted against its own time base. The QR columns are
% sorted descending, so columns 1:K are the K largest exponents.
K = min(20, size(model_q.lya_results.local_LE_spectrum_t, 2));
t_qr = model_q.lya_results.t_lya;
qr_local = model_q.lya_results.local_LE_spectrum_t;

figure('Name', 'Benettin vs QR: local Lyapunov exponents');
hold on;
% QR top-K local exponents in a gradient (dark = largest lambda_1, light =
% smaller). Only the first and last are labelled to keep the legend readable.
cmap = parula(K);
for j = 1:K
    if j == 1
        name = 'QR local \lambda_{1} (largest)';
    elseif j == K
        name = sprintf('QR local \\lambda_{%d}', K);
    else
        name = '';
    end
    if isempty(name)
        plot(t_qr, qr_local(:, j), 'Color', cmap(j, :), 'HandleVisibility', 'off');
    else
        plot(t_qr, qr_local(:, j), 'Color', cmap(j, :), 'DisplayName', name);
    end
end
% Benettin local LLE on top, in black, so it stands out against the QR band.
plot(model_b.lya_results.t_lya, model_b.lya_results.local_lya, ...
    'k', 'LineWidth', 1.5, 'DisplayName', 'Benettin (local LLE)');
yline(LLE_benettin, 'k--', 'DisplayName', 'Benettin LLE (mean)');
hold off;
xlabel('time (s)');
ylabel('local Lyapunov exponent');
title(sprintf('Benettin LLE = %+.4f   |   QR largest = %+.4f   (top-%d QR local exponents shown)', ...
    LLE_benettin, LE_spectrum(1), K));
legend('Location', 'best');
