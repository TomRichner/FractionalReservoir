% test_RMTBlocks_predictions.m
%
% Empirical validation of the RMTBlocks spectral predictions.
%
% No published result covers sparse + block mean + block variance simultaneously:
%
%   Harris et al. (2023)      sparse + Dale, but COLUMN-ONLY mu and sigma
%   Aljadeff et al. (2015)    full D x D block variances, but DENSE and ZERO-MEAN
%   Ahmadian et al. (2015)    arbitrary low-rank mean, but SEPARABLE variance only
%
% RMTBlocks composes them: Harris's blockwise sparse moment matching feeds Aljadeff's
% Perron-root radius, with the block mean contributing only outliers (Ahmadian Sec. V C).
% That composition is a CONJECTURE. This script tests it numerically, which is the same
% standard Harris meets for his own Eq. (18) (validated against numerical spectra in his
% Fig. 2(b), not proven).
%
% Measurement notes:
%
%  * The bulk edge is estimated with a PERCENTILE of |lambda|, not the max. Harris
%    Sec. III D 2 documents a small number of "local eigenvalue outliers" escaping the
%    circular support - which is exactly why his Eq. (18) is introduced as the
%    "approximate" radius - so max(|lambda|) would reject a correct formula. For a
%    uniform disk the q-quantile of |lambda| is R*sqrt(q), so p99 should land near
%    0.995*R; both statistics are reported.
%
%  * Harris also notes the escape distance GROWS as alpha -> 1, so both a dense and a
%    sparse case are included.
%
%  * Section 3 is the discriminating test: a configuration where the block Perron root
%    and the naive population-averaged variance disagree by 2x. This is Aljadeff's
%    Fig. 1(a,b) point - the mean synaptic gain gives the wrong answer - and it is the
%    evidence that the block structure genuinely matters rather than averaging out.
%
% Run with the matlab MCP (run_matlab_file), not matlab -batch.

%% Bootstrap
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(repo_root);
setup_paths();

rng(20260810);

n_pass = 0;
n_fail = 0;

fprintf('\n=========================================================================\n');
fprintf('  RMTBlocks prediction validation\n');
fprintf('=========================================================================\n');

%% ---------------------------------------------------------------
%  Configuration grid
%  ---------------------------------------------------------------
% All statistics are given in units of 1/sqrt(N) so that R stays O(1) as N varies.
% mu and sigma are specified as (postsynaptic <- presynaptic) blocks.

cfgs = {};

cfgs{end+1} = struct('name', 'column-uniform (Harris)', ...
    'f', [0.5 0.5], ...
    'mu',  [ 3.0 -4.0 ;  3.0 -4.0 ], ...
    'sig', [ 1.0  1.0 ;  1.0  1.0 ]);

cfgs{end+1} = struct('name', 'mu_EE elevated', ...
    'f', [0.5 0.5], ...
    'mu',  [ 4.5 -4.0 ;  3.0 -4.0 ], ...
    'sig', [ 1.0  1.0 ;  1.0  1.0 ]);

cfgs{end+1} = struct('name', 'sigma_II reduced', ...
    'f', [0.5 0.5], ...
    'mu',  [ 3.0 -4.0 ;  3.0 -4.0 ], ...
    'sig', [ 1.0  1.0 ;  1.0  0.3 ]);

cfgs{end+1} = struct('name', 'non-separable sigma, zero mean', ...
    'f', [0.5 0.5], ...
    'mu',  [ 0.0  0.0 ;  0.0  0.0 ], ...
    'sig', [ 1.0  0.4 ;  0.4  1.6 ]);

cfgs{end+1} = struct('name', 'D=3 cell types', ...
    'f', [0.6 0.2 0.2], ...
    'mu',  [ 2.0 -3.0 -1.0 ;  2.0 -3.0 -1.0 ;  1.0 -2.0 -1.0 ], ...
    'sig', [ 1.0  0.8  0.5 ;  0.6  1.2  0.4 ;  0.9  0.3  1.1 ]);

alphas = [1.0, 1/3];
Ns     = [300, 1000];
reps   = [ 20,   10];      % fewer realizations at the larger N

%% ---------------------------------------------------------------
%  Section 1 & 2: predicted vs measured bulk edge and outliers
%  ---------------------------------------------------------------
fprintf('\n%-32s %6s %6s %8s %8s %8s %8s %9s\n', ...
    'configuration', 'alpha', 'N', 'R_pred', 'p99/R', 'max/R', 'edge/R', 'outl.err');
fprintf('%s\n', repmat('-', 1, 96));

% edge_err(ci, ai, ni) = |edge_estimate/R_pred - 1|,  outl_err(ci, ai, ni) likewise.
% Thresholds are applied at the LARGEST N and convergence is checked across N: at
% N = 300 the finite-size bias is several percent even for Harris's own already-
% validated column-uniform case, so a tight bound there would reject a correct
% formula. What distinguishes "correct but finite-size limited" from "wrong" is
% whether the error shrinks as N grows.
edge_err = nan(numel(cfgs), numel(alphas), numel(Ns));
outl_err = nan(numel(cfgs), numel(alphas), numel(Ns));
max_ratios = [];

for ci = 1:numel(cfgs)
    c = cfgs{ci};
    for ai = 1:numel(alphas)
        for ni = 1:numel(Ns)
            N = Ns(ni);
            n_rep = reps(ni);
            sq = 1/sqrt(N);

            r = RMTBlocks(N);
            r.alpha = alphas(ai);
            r.set_types(c.f, c.mu * sq, c.sig * sq);

            R_pred = r.R;
            lam_pred = r.lambda_O;
            % Only outliers that actually escape the disk are identifiable
            outlier_pred = lam_pred(abs(lam_pred) > 1.10 * R_pred);

            p99_acc = zeros(n_rep, 1);
            max_acc = zeros(n_rep, 1);
            oerr_acc = nan(n_rep, 1);

            for rep = 1:n_rep
                % Fresh disorder for each realization
                r.A = randn(N, N);
                r.update_sparsity();

                lam = eig(r.W);

                % Peel off the eigenvalues matched to predicted outliers before
                % estimating the bulk edge.
                lam_bulk = lam;
                if ~isempty(outlier_pred)
                    d_rel = zeros(numel(outlier_pred), 1);
                    for k = 1:numel(outlier_pred)
                        [dmin, imin] = min(abs(lam_bulk - outlier_pred(k)));
                        d_rel(k) = dmin / abs(outlier_pred(k));
                        lam_bulk(imin) = [];
                    end
                    oerr_acc(rep) = max(d_rel);
                end

                p99_acc(rep) = RMTBlocks.percentile(abs(lam_bulk), 99);
                max_acc(rep) = max(abs(lam_bulk));
            end

            p99_r = mean(p99_acc) / R_pred;
            max_r = mean(max_acc) / R_pred;
            % p99 of a uniform disk sits at sqrt(0.99)*R, so undo that to get an
            % edge estimate directly comparable to R_pred.
            edge_r = p99_r / sqrt(0.99);
            oerr = mean(oerr_acc, 'omitnan');

            edge_err(ci, ai, ni) = abs(edge_r - 1);
            outl_err(ci, ai, ni) = oerr;
            max_ratios(end+1) = max_r; %#ok<SAGROW>

            if isnan(oerr)
                oerr_str = '     n/a';
            else
                oerr_str = sprintf('%8.4f', oerr);
            end
            fprintf('%-32s %6.3f %6d %8.4f %8.4f %8.4f %8.4f %s\n', ...
                c.name, alphas(ai), N, R_pred, p99_r, max_r, edge_r, oerr_str);
        end
    end
end

fprintf('%s\n', repmat('-', 1, 96));

edge_big = edge_err(:, :, end);
edge_small = edge_err(:, :, 1);
outl_big = outl_err(:, :, end);

[n_pass, n_fail] = report(all(edge_big(:) < 0.04), ...
    sprintf('bulk edge within 4%% of R_pred at N = %d, every configuration', Ns(end)), ...
    sprintf('worst %.4f', max(edge_big(:))), n_pass, n_fail);

[n_pass, n_fail] = report(all(edge_big(:) < edge_small(:)), ...
    sprintf('bulk-edge error shrinks from N = %d to N = %d in every configuration', Ns(1), Ns(end)), ...
    sprintf('median %.4f -> %.4f', median(edge_small(:)), median(edge_big(:))), n_pass, n_fail);

[n_pass, n_fail] = report(all(max_ratios > 0.95 & max_ratios < 1.35), ...
    'max(|lambda|)/R_pred within [0.95, 1.35] (local outliers escape, per Harris)', ...
    sprintf('range [%.4f, %.4f]', min(max_ratios), max(max_ratios)), n_pass, n_fail);

ob = outl_big(~isnan(outl_big));
[n_pass, n_fail] = report(all(ob < 0.12), ...
    sprintf('predicted outliers matched within 12%% of their modulus at N = %d', Ns(end)), ...
    sprintf('worst %.4f', max(ob)), n_pass, n_fail);

%% ---------------------------------------------------------------
%  Section 3: discrimination against the population-averaged variance
%  ---------------------------------------------------------------
% Aljadeff Fig. 1(a,b): the mean synaptic gain predicts the wrong phase for
% block-structured networks. Here f = [0.9 0.1] with a hyper-variable minority whose
% recurrent block dominates the Perron root; the naive average is ~2x too small.
fprintf('\n--- 3. Block Perron root vs naive averaged variance ---\n');

N3 = 1200;
sq3 = 1/sqrt(N3);
f3 = [0.9 0.1];
sig3 = [0.5 0.5; 0.5 4.0] * sq3;

r3 = RMTBlocks(N3);
r3.alpha = 1.0;
r3.set_types(f3, zeros(2), sig3);

R_block = r3.R;

% Naive alternative: treat the variance as homogeneous at its population average,
% i.e. sum_{a,b} f_a f_b sigma_s_sq(a,b). This is what a single-cell-type reading of
% Harris Eq. (18) would give.
ss = r3.sigma_s_sq;
R_avg = sqrt(N3 * (f3(:)' * ss * f3(:)));

n_rep3 = 12;
edge_acc = zeros(n_rep3, 1);
for rep = 1:n_rep3
    r3.A = randn(N3, N3);
    r3.update_sparsity();
    lam = eig(r3.W);
    edge_acc(rep) = RMTBlocks.percentile(abs(lam), 99) / sqrt(0.99);
end
edge_meas = mean(edge_acc);

err_block = abs(edge_meas - R_block) / R_block;
err_avg   = abs(edge_meas - R_avg)   / R_avg;

fprintf('  measured bulk edge      : %.4f\n', edge_meas);
fprintf('  block Perron root (ours): %.4f   (relative error %.4f)\n', R_block, err_block);
fprintf('  naive averaged variance : %.4f   (relative error %.4f)\n', R_avg, err_avg);

[n_pass, n_fail] = report(err_block < 0.05, ...
    'block Perron root predicts the measured edge within 5%', ...
    sprintf('%.4f', err_block), n_pass, n_fail);
[n_pass, n_fail] = report(err_avg > 5 * err_block, ...
    'naive averaged variance is decisively worse (>5x the error)', ...
    sprintf('%.4f vs %.4f', err_avg, err_block), n_pass, n_fail);

%% ---------------------------------------------------------------
%  Section 4: finite-size trend
%  ---------------------------------------------------------------
% If the formula is right, agreement should IMPROVE with N. If instead the error is
% flat or growing, the composition is wrong rather than merely finite-size limited.
fprintf('\n--- 4. Finite-size trend (non-separable sigma, alpha = 1/3) ---\n');

sig4 = [1.0 0.4; 0.4 1.6];
mu4  = [2.0 -3.0; 1.0 -2.5];
Ns4 = [200, 400, 800, 1600];
errs4 = zeros(size(Ns4));

for k = 1:numel(Ns4)
    N4 = Ns4(k);
    sq4 = 1/sqrt(N4);
    r4 = RMTBlocks(N4);
    r4.alpha = 1/3;
    r4.set_types([0.5 0.5], mu4 * sq4, sig4 * sq4);

    R4 = r4.R;
    lam_pred4 = r4.lambda_O;
    outl4 = lam_pred4(abs(lam_pred4) > 1.10 * R4);

    n_rep4 = max(3, round(2400 / N4));
    acc = zeros(n_rep4, 1);
    for rep = 1:n_rep4
        r4.A = randn(N4, N4);
        r4.update_sparsity();
        lam = eig(r4.W);
        for kk = 1:numel(outl4)
            [~, imin] = min(abs(lam - outl4(kk)));
            lam(imin) = [];
        end
        acc(rep) = RMTBlocks.percentile(abs(lam), 99) / sqrt(0.99);
    end
    errs4(k) = abs(mean(acc) - R4) / R4;
    fprintf('  N = %5d   R_pred = %.4f   measured = %.4f   rel err = %.4f\n', ...
        N4, R4, mean(acc), errs4(k));
end

[n_pass, n_fail] = report(errs4(end) < errs4(1), ...
    'relative error decreases from the smallest to the largest N', ...
    sprintf('%.4f -> %.4f', errs4(1), errs4(end)), n_pass, n_fail);
[n_pass, n_fail] = report(errs4(end) < 0.03, ...
    'relative error below 3% at the largest N', ...
    sprintf('%.4f', errs4(end)), n_pass, n_fail);

%% ---------------------------------------------------------------
%  Summary
%  ---------------------------------------------------------------
fprintf('\n=========================================================================\n');
fprintf('  passed: %d   failed: %d\n', n_pass, n_fail);
fprintf('=========================================================================\n\n');

if n_fail > 0
    warning('test_RMTBlocks_predictions:Failed', ...
        ['%d check(s) failed. These are empirical thresholds on a conjectured ', ...
         'formula, so inspect the printed table before concluding the formula is ', ...
         'wrong - a borderline ratio may just be finite-size scatter.'], n_fail);
end

%% ---------------------------------------------------------------
%  Local functions
%  ---------------------------------------------------------------
function [n_pass, n_fail] = report(ok, name, detail, n_pass, n_fail)
if ok
    n_pass = n_pass + 1;
    fprintf('  PASS  %s  [%s]\n', name, detail);
else
    n_fail = n_fail + 1;
    fprintf('  FAIL  %s  [%s]\n', name, detail);
end
end
