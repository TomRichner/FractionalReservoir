function out_file = amp_outlier_compute(varargin)
%AMP_OUTLIER_COMPUTE Non-normal amplification against the effective outlier, from the saved outlier_spectra runs (no integration).
%
%   out_file = amp_outlier_compute()
%   out_file = amp_outlier_compute('run_dir', d, 'out_dir', o, 't_grid', 0:0.05:2)
%
% For every run of data/narrative_essential_spectral/outlier_spectra and
% every calibration state, with the frozen-gain effective connectivity
%   M(i, j) = W(i, j) gain_route(type of i, j)      (the class's own M)
% and J = (-I + M)/tau_d (the paper's J_xx, run_transient_gain), computes
%   1. alpha = max Re eig J, omega = max eig (J + J')/2, margin = omega - alpha
%   2. G(t) = ||expm(J t)||_2 on t_grid; peak over the grid, G(tau_d), G(0.3)
%   3. the same for the uniform-gain control M_u = W rms(d) (rms over every
%      neuron's route gains, the rule of ospec_null / R.cal.ev_uniform)
%   4. the mean part P(i, j) = mu(type i, type j) gain_route(type i, j)
%      (mu: realized block means of W, zeros included): its two nonzero
%      eigenvalues (= eig K_eff), ||P||_F, and w_FF = sqrt(||P||_F^2 -
%      sum |lambda_k|^2); the same at uniform gain
%   5. Henrici's departure of M and of M_u
% Checks (on the first state of the EEx1p5 sfa3_stdc2 run, and the max over
% all states): eig M against R.cal.ev_M, eig M_u against R.cal.ev_uniform,
% eig P against R.K_eig at the state's column, M against
% ospec_effective_matrix from the saved state vector.
%
% See also: amp_outlier_figs, ospec_load_run, run_transient_gain

p = inputParser;
p.addParameter('run_dir', fullfile(repo_root(), 'data', 'narrative_essential_spectral'));
p.addParameter('out_dir', fullfile(repo_root(), 'data', 'explorations', 'amplification_outlier'));
p.addParameter('t_grid', 0:0.05:2);
p.parse(varargin{:});
o = p.Results;
if ~isfolder(o.out_dir); mkdir(o.out_dir); end

runs = ospec_read_index(o.run_dir);
rows = {};
checks = struct('run', {}, 'state', {}, 'ev_M', {}, 'ev_u', {}, 'K_eig', {}, 'M_vs_class', {}, 'M_vs_Wdiag', {});
tg = o.t_grid(:)';
for r = 1:height(runs)
    if ~runs.ok(r); continue; end
    R = ospec_load_run(o.run_dir, r);
    [~, params, ws] = ospec_rebuild_model(R);
    n = R.n; C = numel(R.n_per_type); I = eye(n); tau = R.tau_d;
    Wf = ws.Wf; mu = ws.mu; tof = ws.type_of;
    Pw = mu(tof, tof);                        % block means as an n x n matrix
    evW = eig(Wf); hW2 = norm(Wf, 'fro')^2 - sum(abs(evW).^2);
    Ws = (Wf + Wf') / 2; evWs = eig(Ws);
    n_cal = numel(R.cal.t);
    fprintf('%-10s %-14s', R.network, R.condition);
    for k = 1:n_cal
        G = R.cal.gain_route(:, :, k); G(isnan(G)) = 0;
        M = Wf .* G(tof, :);
        evM = eig(M);
        d_rms = sqrt(mean(G(:).^2));
        % -- checks
        e1 = max(abs(sort(evM, 'ComparisonMethod', 'real') - sort(R.cal.ev_M(:, k), 'ComparisonMethod', 'real')));
        evU = d_rms * evW;
        e2 = max(abs(sort(evU, 'ComparisonMethod', 'real') - sort(R.cal.ev_uniform(:, k), 'ComparisonMethod', 'real')));
        P = Pw .* G(tof, :);
        evP = eig(P); [~, ip] = sort(abs(evP), 'descend'); lamP = evP(ip(1:C));
        Ksv = R.K_eig(:, R.cal.idx(k));
        e3 = max(abs(sort(lamP, 'ComparisonMethod', 'real') - sort(Ksv, 'ComparisonMethod', 'real')));
        e4 = NaN; e5 = norm(M - Wf * diag(R.cal.gain(:, k)), 'fro') / max(norm(M, 'fro'), eps);
        if k == 1
            Mc = ospec_effective_matrix(R.cal.S(k, :)', params, ws);
            e4 = norm(M - Mc, 'fro') / max(norm(M, 'fro'), eps);
        end
        checks(end + 1) = struct('run', sprintf('%s/%s', R.network, R.condition), 'state', k, ...
            'ev_M', e1, 'ev_u', e2, 'K_eig', e3, 'M_vs_class', e4, 'M_vs_Wdiag', e5); %#ok<AGROW>
        % -- 1, 2 measured
        J = (-I + M) / tau;
        a_m = max(real(evM) - 1) / tau;
        w_m = max(eig((J + J') / 2));
        Gt_m = gain_curve(J, tg);
        % -- 1, 2 uniform
        Ju = (-I + d_rms * Wf) / tau;
        a_u = (max(real(evU)) - 1) / tau;
        w_u = (d_rms * max(evWs) - 1) / tau;
        Gt_u = gain_curve(Ju, tg);
        % -- 4 mean part, measured and uniform
        PF = norm(P, 'fro');
        wff = sqrt(max(0, PF^2 - sum(abs(lamP).^2)));
        Ku = ws.K_struct * d_rms; lamU = eig(Ku);
        PFu = norm(Pw, 'fro') * d_rms;
        wffu = sqrt(max(0, PFu^2 - sum(abs(lamU).^2)));
        [~, i1] = max(real(lamP)); [~, i1u] = max(real(lamU));
        % -- 5 Henrici
        hM = sqrt(max(0, norm(M, 'fro')^2 - sum(abs(evM).^2)));
        hU = d_rms * sqrt(max(0, hW2));
        dE = mean(G(:, ws.ti{1}), 'all'); dI = mean(G(:, ws.ti{2}), 'all');
        rows(end + 1, :) = {R.network, R.condition, k, R.cal.t(k), d_rms, dE, dI, ...
            R.frac_on(1, R.cal.idx(k)), R.frac_on(2, R.cal.idx(k)), R.cal.R_eff(k), ...
            a_m, w_m, w_m - a_m, max(Gt_m), Gt_m(abs(tg - tau) < 1e-12), interp1(tg, Gt_m, 0.3), ...
            a_u, w_u, w_u - a_u, max(Gt_u), Gt_u(abs(tg - tau) < 1e-12), interp1(tg, Gt_u, 0.3), ...
            lamP(i1), lamP(3 - i1), PF, wff, lamU(i1u), lamU(3 - i1u), PFu, wffu, hM, hU, ...
            {Gt_m}, {Gt_u}}; %#ok<AGROW>
        fprintf('.');
    end
    fprintf('\n');
end
names = {'network', 'condition', 'state', 't', 'd_rms', 'dE', 'dI', 'frac_on_E', 'frac_on_I', 'R_eff', ...
    'alpha', 'omega', 'margin', 'Gpeak', 'G_tau', 'G_0p3', ...
    'alpha_u', 'omega_u', 'margin_u', 'Gpeak_u', 'G_tau_u', 'G_0p3_u', ...
    'lamO', 'lam2', 'PF', 'wFF', 'lamO_u', 'lam2_u', 'PF_u', 'wFF_u', 'henrici', 'henrici_u', ...
    'Gt', 'Gt_u'};
T = cell2table(rows, 'VariableNames', names);
CK = struct2table(checks);
t_grid = tg; %#ok<NASGU>
out_file = fullfile(o.out_dir, 'amp_outlier_results.mat');
save(out_file, 'T', 'CK', 't_grid', '-v7');
fprintf('checks, max over all states: eig M %.2g, eig M_u %.2g, eig P vs K_eig %.2g, M vs class %.2g, M vs W diag(gain) %.2g\n', ...
    max(CK.ev_M), max(CK.ev_u), max(CK.K_eig), max(CK.M_vs_class, [], 'omitnan'), max(CK.M_vs_Wdiag));
fprintf('saved %s\n', out_file);
end

function g = gain_curve(J, tg)
% ||expm(J t)||_2 on an even grid starting at 0: one expm, then powers.
dt = tg(2) - tg(1);
E = expm(J * dt); Phi = eye(size(J));
g = zeros(size(tg)); g(1) = 1;
for q = 2:numel(tg)
    Phi = E * Phi;
    g(q) = norm(Phi, 2);
end
end

function r = repo_root()
r = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
