% test_RMTBlocks_equivalence.m
%
% Backward-compatibility and reduction tests for RMTBlocks.
%
% The central claim of RMTBlocks is that it is a strict generalization: with
% column-uniform blocks it must reproduce the original Harris/RMT parameterization
% EXACTLY, both in the matrix it builds and in its theoretical predictions.
%
% Sections:
%   1. Bit-for-bit W against RMT               (all 4 zrs_modes, both alphas, shift, f=1)
%   2. Bit-for-bit W against RMTMatrix         (the SRNNModel2 configuration)
%   3. Harris Eq. (17)/(18) reduction          (no reference class needed)
%   4. Alias round-trip and ambiguity errors
%   5. Empty populations: f = [1 0] vs a true D = 1 object
%   6. Block functionality: non-separable sigma, g_mu / g_sigma
%
% Sections 1 and 2 are skipped with a printed note when the reference classes are
% not on the path, so this repository stays standalone-clonable.
%
% Run with the matlab MCP (run_matlab_file), not matlab -batch.

%% Bootstrap
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(repo_root);
setup_paths();

% Probe for the reference classes in sibling repositories. All of the author's
% repos live side by side under github_repos/, so a relative probe is reliable;
% if they are not there the relevant sections simply skip.
siblings = fileparts(repo_root);
ref_paths = { ...
    fullfile(siblings, 'ConnectivityAdaptation', 'RandomMatrixTheory'), ...
    fullfile(siblings, 'FractionalReservoir', 'src', 'connectivity') ...
    };
for k = 1:numel(ref_paths)
    if exist(ref_paths{k}, 'dir')
        addpath(ref_paths{k});
    end
end

has_RMT       = exist('RMT', 'class') == 8;
has_RMTMatrix = exist('RMTMatrix', 'class') == 8;

n_pass = 0;
n_fail = 0;
n_skip = 0;

fprintf('\n=================================================\n');
fprintf('  RMTBlocks equivalence tests\n');
fprintf('=================================================\n');

%% ---------------------------------------------------------------
%  Section 1: bit-for-bit W against RMT
%  ---------------------------------------------------------------
% These are the six configurations from Fig_1_RMT_examples.m plus an explicit
% SZRS case, which Fig_1 never exercises.

N = 300;
sq = 1 / sqrt(N);
E_W = 0.05 / sqrt(N);

% Example (e) solves Harris Eq. (14) for sigma_tilde_i at fixed total variance.
sig_e_e = 0.35 * sq;
sig_i_e = sqrt((1/N - 0.5 * sig_e_e^2) / 0.5);

cfgs = {};
cfgs{end+1} = struct('name', '(a) dense unbalanced, no ZRS', ...
    'alpha', 1.0, 'mu_e', E_W, 'mu_i', E_W, 'sig_e', sq,      'sig_i', sq, ...
    'f', 1,   'zrs', 'none',          'shift', 0);
cfgs{end+1} = struct('name', '(b) dense balanced, shifted', ...
    'alpha', 1.0, 'mu_e', 0,   'mu_i', 0,   'sig_e', sq,      'sig_i', sq, ...
    'f', 1,   'zrs', 'none',          'shift', -1);
cfgs{end+1} = struct('name', '(c) dense balanced Dale', ...
    'alpha', 1.0, 'mu_e', sq,  'mu_i', -sq, 'sig_e', sq,      'sig_i', sq, ...
    'f', 0.5, 'zrs', 'none',          'shift', -1);
cfgs{end+1} = struct('name', '(d) dense Dale ZRS', ...
    'alpha', 1.0, 'mu_e', sq,  'mu_i', -sq, 'sig_e', sq,      'sig_i', sq, ...
    'f', 0.5, 'zrs', 'ZRS',           'shift', -1);
cfgs{end+1} = struct('name', '(e) dense Dale different sigmas, ZRS', ...
    'alpha', 1.0, 'mu_e', sq,  'mu_i', -sq, 'sig_e', sig_e_e, 'sig_i', sig_i_e, ...
    'f', 0.5, 'zrs', 'ZRS',           'shift', -1);
cfgs{end+1} = struct('name', '(f) sparse Dale Partial_SZRS', ...
    'alpha', 0.5, 'mu_e', sq,  'mu_i', -sq, 'sig_e', sig_e_e, 'sig_i', sig_i_e, ...
    'f', 0.5, 'zrs', 'Partial_SZRS',  'shift', -1);
cfgs{end+1} = struct('name', '(g) sparse Dale SZRS', ...
    'alpha', 0.5, 'mu_e', sq,  'mu_i', -sq, 'sig_e', sq,      'sig_i', sq, ...
    'f', 0.5, 'zrs', 'SZRS',          'shift', 0);

fprintf('\n--- 1. Bit-for-bit W vs RMT ---\n');
if ~has_RMT
    fprintf('  SKIP: class RMT not found on the path.\n');
    n_skip = n_skip + numel(cfgs);
else
    for k = 1:numel(cfgs)
        c = cfgs{k};

        rng(1000 + k);
        ref = RMT(N);
        ref = apply_cfg(ref, c, true);

        rng(1000 + k);
        new = RMTBlocks(N);
        new = apply_cfg(new, c, true);

        % Inject the reference realization explicitly (alpha's setter redraws S,
        % so this must come after apply_cfg).
        new.A = ref.A;
        new.S = ref.S;

        W_ref = ref.W;
        W_new = new.W;

        if isequal(W_ref, W_new)
            [n_pass, n_fail] = report(true, sprintf('%s: W bit-for-bit', c.name), '', n_pass, n_fail);
        else
            d = max(abs(W_ref(:) - W_new(:)));
            scale = max(abs(W_ref(:)));
            if d <= 1e-14 * max(scale, 1)
                [n_pass, n_fail] = report(true, ...
                    sprintf('%s: W equal to %.1e (not bitwise)', c.name, d), '', n_pass, n_fail);
            else
                [n_pass, n_fail] = report(false, sprintf('%s: W', c.name), ...
                    sprintf('max abs diff %.3e (scale %.3e)', d, scale), n_pass, n_fail);
            end
        end

        % R must agree too (it is read by Fig_1 under every zrs_mode).
        [n_pass, n_fail] = report(abs(ref.R - new.R) <= 1e-12 * max(abs(ref.R), 1), ...
            sprintf('%s: R matches RMT', c.name), ...
            sprintf('ref %.10g vs new %.10g', ref.R, new.R), n_pass, n_fail);
    end
end

%% ---------------------------------------------------------------
%  Section 2: bit-for-bit W against RMTMatrix (SRNNModel2 configuration)
%  ---------------------------------------------------------------
fprintf('\n--- 2. Bit-for-bit W vs RMTMatrix (SRNNModel2 config) ---\n');
if ~has_RMTMatrix
    fprintf('  SKIP: class RMTMatrix not found on the path.\n');
    n_skip = n_skip + 1;
else
    n_srnn = 300;
    indegree = 100;
    alpha_srnn = indegree / n_srnn;
    F = 1 / sqrt(n_srnn * alpha_srnn * (2 - alpha_srnn));

    c = struct('alpha', alpha_srnn, 'mu_e', 3*F, 'mu_i', -4*F, ...
        'sig_e', F, 'sig_i', F, 'f', 0.5, 'zrs', 'none', 'shift', 0);

    rng(2024);
    ref = RMTMatrix(n_srnn);
    ref = apply_cfg(ref, c, false);          % RMTMatrix has no shift property

    rng(2024);
    new = RMTBlocks(n_srnn);
    new = apply_cfg(new, c, true);
    new.A = ref.A;
    new.S = ref.S;

    W_ref = ref.W;
    W_new = new.W;
    [n_pass, n_fail] = report(isequal(W_ref, W_new), ...
        'SRNNModel2 config: W bit-for-bit vs RMTMatrix', ...
        sprintf('max abs diff %.3e', max(abs(W_ref(:) - W_new(:)))), n_pass, n_fail);

    [n_pass, n_fail] = report(abs(ref.R - new.R) <= 1e-12 * abs(ref.R), ...
        'SRNNModel2 config: R matches RMTMatrix', ...
        sprintf('ref %.10g vs new %.10g', ref.R, new.R), n_pass, n_fail);

    [n_pass, n_fail] = report(abs(ref.lambda_O - new.lambda_O(1)) <= 1e-12 * abs(ref.lambda_O), ...
        'SRNNModel2 config: lambda_O(1) matches RMTMatrix', ...
        sprintf('ref %.10g vs new %.10g', ref.lambda_O, new.lambda_O(1)), n_pass, n_fail);
end

%% ---------------------------------------------------------------
%  Section 3: Harris Eq. (17)/(18) reduction, self-contained
%  ---------------------------------------------------------------
fprintf('\n--- 3. Harris Eq. (17)/(18) reduction ---\n');

n3 = 300;
alpha3 = 1/3;
F3 = 1 / sqrt(n3 * alpha3 * (2 - alpha3));
f3 = 0.5;

r3 = RMTBlocks(n3);
r3.alpha = alpha3;
r3.mu_tilde_e = 3*F3;
r3.mu_tilde_i = -4*F3;
r3.sigma_tilde_e = F3;
r3.sigma_tilde_i = F3;
r3.f = f3;

% Harris Eqs. (15), (16)
mu_se_h  = alpha3 * (3*F3);
mu_si_h  = alpha3 * (-4*F3);
sig_se_h = alpha3*(1-alpha3)*(3*F3)^2  + alpha3*F3^2;
sig_si_h = alpha3*(1-alpha3)*(-4*F3)^2 + alpha3*F3^2;

lambda_O_h = n3 * (f3*mu_se_h + (1-f3)*mu_si_h);              % Eq. (17)
R_h        = sqrt(n3 * (f3*sig_se_h + (1-f3)*sig_si_h));      % Eq. (18)

lam = r3.lambda_O;
[n_pass, n_fail] = report(abs(lam(1) - lambda_O_h) <= 1e-12*abs(lambda_O_h), ...
    'lambda_O(1) == Harris Eq. (17)', ...
    sprintf('%.10g vs %.10g', lam(1), lambda_O_h), n_pass, n_fail);
[n_pass, n_fail] = report(abs(lam(2)) <= 1e-12*max(abs(lambda_O_h), 1), ...
    'lambda_O(2) == 0 for column-uniform blocks', ...
    sprintf('%.3e', abs(lam(2))), n_pass, n_fail);
[n_pass, n_fail] = report(abs(r3.R - R_h) <= 1e-12*R_h, ...
    'R == Harris Eq. (18)', ...
    sprintf('%.10g vs %.10g', r3.R, R_h), n_pass, n_fail);
[n_pass, n_fail] = report(r3.is_column_uniform(), ...
    'is_column_uniform() true for scalar-alias configuration', '', n_pass, n_fail);

fprintf('  [info] for this configuration the sparsity-induced term alpha(1-alpha)mu^2\n');
fprintf('         supplies %.1f%% (E) and %.1f%% (I) of the block variance.\n', ...
    100 * alpha3*(1-alpha3)*(3*F3)^2 / sig_se_h, ...
    100 * alpha3*(1-alpha3)*(4*F3)^2 / sig_si_h);

%% ---------------------------------------------------------------
%  Section 4: alias round-trip and ambiguity errors
%  ---------------------------------------------------------------
fprintf('\n--- 4. Alias behavior ---\n');

r4 = RMTBlocks(100);
r4.mu_tilde_e = 0.3;
[n_pass, n_fail] = report(r4.mu_tilde_e == 0.3 && all(r4.mu_tilde(:,1) == 0.3), ...
    'set/get mu_tilde_e round-trips and writes the whole column', '', n_pass, n_fail);

r4.mu_tilde(1, 1) = 0.9;                 % break column uniformity
ok = false;
try
    v = r4.mu_tilde_e;
catch ME
    ok = strcmp(ME.identifier, 'RMTBlocks:AmbiguousAlias');
end
[n_pass, n_fail] = report(ok, 'non-uniform column raises RMTBlocks:AmbiguousAlias', '', n_pass, n_fail);

r4b = RMTBlocks(100);
r4b.set_types([1/3 1/3 1/3], zeros(3), 0.1*ones(3));
ok = false;
try
    v = r4b.mu_tilde_e;
catch ME
    ok = strcmp(ME.identifier, 'RMTBlocks:NotTwoTypes');
end
[n_pass, n_fail] = report(ok, 'D = 3 raises RMTBlocks:NotTwoTypes on a D=2 alias', '', n_pass, n_fail);
[n_pass, n_fail] = report(r4b.n_types == 3 && sum(r4b.n_per_type) == 100, ...
    'set_types(D=3) allocates counts summing to N', ...
    mat2str(r4b.n_per_type), n_pass, n_fail);

% Inconsistent sizes must be caught at use time, not silently
r4c = RMTBlocks(50);
r4c.mu_tilde = zeros(3);
ok = false;
try
    W = r4c.W;
catch ME
    ok = strcmp(ME.identifier, 'RMTBlocks:InconsistentTypes');
end
[n_pass, n_fail] = report(ok, 'size mismatch raises RMTBlocks:InconsistentTypes at use time', '', n_pass, n_fail);

%% ---------------------------------------------------------------
%  Section 5: empty populations
%  ---------------------------------------------------------------
fprintf('\n--- 5. Empty populations (f = [1 0]) ---\n');

n5 = 200;
mu11 = 0.02;
sig11 = 0.05;

r5_two = RMTBlocks(n5);
r5_two.alpha = 0.4;
r5_two.mu_tilde = [mu11, -0.07; mu11, -0.07];
r5_two.sigma_tilde = [sig11, 0.09; sig11, 0.09];
r5_two.f = 1;                                    % -> [1 0], type 2 empty

r5_one = RMTBlocks(n5);
r5_one.alpha = 0.4;
r5_one.set_types(1, mu11, sig11);

[n_pass, n_fail] = report(isequal(r5_two.n_per_type, [n5 0]), ...
    'f = 1 gives counts [N 0] without error', mat2str(r5_two.n_per_type), n_pass, n_fail);
[n_pass, n_fail] = report(abs(r5_two.R - r5_one.R) <= 1e-12*max(r5_one.R, 1), ...
    'R with an empty population equals the true D = 1 answer', ...
    sprintf('%.10g vs %.10g', r5_two.R, r5_one.R), n_pass, n_fail);
[n_pass, n_fail] = report(abs(r5_two.lambda_O(1) - r5_one.lambda_O(1)) <= 1e-12*max(abs(r5_one.lambda_O(1)), 1), ...
    'lambda_O(1) with an empty population equals the D = 1 answer', ...
    sprintf('%.10g vs %.10g', r5_two.lambda_O(1), r5_one.lambda_O(1)), n_pass, n_fail);

W5 = r5_two.W;
[n_pass, n_fail] = report(isequal(size(W5), [n5 n5]) && all(isfinite(W5(:))), ...
    'W builds cleanly with an empty population', '', n_pass, n_fail);

%% ---------------------------------------------------------------
%  Section 6: block functionality
%  ---------------------------------------------------------------
fprintf('\n--- 6. Block functionality ---\n');

n6 = 400;
rng(77);
r6 = RMTBlocks(n6);
r6.alpha = 1.0;                                  % dense, so measured stats are clean
% Deliberately NON-SEPARABLE sigma: sig_EE*sig_II ~= sig_EI*sig_IE, which the
% Ahmadian LJR ensemble cannot represent.
r6.mu_tilde    = [ 0.30, -0.50 ;  0.10, -0.20 ];
r6.sigma_tilde = [ 0.05,  0.02 ;  0.03,  0.11 ];
r6.f = 0.5;

[n_pass, n_fail] = report(~r6.is_column_uniform(), ...
    'is_column_uniform() false for genuine block structure', '', n_pass, n_fail);

sep_lhs = r6.sigma_tilde(1,1) * r6.sigma_tilde(2,2);
sep_rhs = r6.sigma_tilde(1,2) * r6.sigma_tilde(2,1);
[n_pass, n_fail] = report(abs(sep_lhs - sep_rhs) > 1e-6, ...
    'test sigma is non-separable (outside the LJR ensemble)', ...
    sprintf('%.4g vs %.4g', sep_lhs, sep_rhs), n_pass, n_fail);

W6 = r6.W;
idx6 = r6.type_indices;
max_mu_err = 0;
max_sd_err = 0;
for a = 1:2
    for b = 1:2
        blk = W6(idx6{a}, idx6{b});
        max_mu_err = max(max_mu_err, abs(mean(blk(:)) - r6.mu_tilde(a,b)));
        max_sd_err = max(max_sd_err, abs(std(blk(:), 1) - r6.sigma_tilde(a,b)));
    end
end
% ~40000 samples per block, so the standard error on the mean is ~sigma/200.
[n_pass, n_fail] = report(max_mu_err < 0.005, ...
    'each (a<-b) block of W has the requested mean', ...
    sprintf('max err %.4g', max_mu_err), n_pass, n_fail);
[n_pass, n_fail] = report(max_sd_err < 0.005, ...
    'each (a<-b) block of W has the requested std dev', ...
    sprintf('max err %.4g', max_sd_err), n_pass, n_fail);

% g_mu / g_sigma must be equivalent to scaling the block matrices directly
r6b = r6.copy();
r6b.g_mu = 2.0;
r6b.g_sigma = 3.0;

r6c = r6.copy();
r6c.mu_tilde = 2.0 * r6.mu_tilde;
r6c.sigma_tilde = 3.0 * r6.sigma_tilde;

[n_pass, n_fail] = report(isequal(r6b.W, r6c.W), ...
    'g_mu / g_sigma equivalent to scaling mu_tilde / sigma_tilde', ...
    sprintf('max abs diff %.3e', max(abs(r6b.W(:) - r6c.W(:)))), n_pass, n_fail);
[n_pass, n_fail] = report(abs(r6b.R - r6c.R) <= 1e-12*r6c.R && ...
                          max(abs(sort(r6b.lambda_O) - sort(r6c.lambda_O))) <= 1e-12*max(abs(r6c.lambda_O(1)), 1), ...
    'g_mu / g_sigma flow through R and lambda_O', '', n_pass, n_fail);

% copy() must be a genuine deep copy of the parameters
r6d = r6.copy();
r6d.mu_tilde(1,1) = 999;
[n_pass, n_fail] = report(r6.mu_tilde(1,1) ~= 999, 'copy() is a deep copy', '', n_pass, n_fail);

% Two outliers where the column-uniform case has one
lam6 = r6.lambda_O;
[n_pass, n_fail] = report(numel(lam6) == 2 && abs(lam6(2)) > 1e-8, ...
    'block means produce a second nonzero outlier', ...
    sprintf('lambda_O = [%.4g, %.4g]', lam6(1), lam6(2)), n_pass, n_fail);

%% ---------------------------------------------------------------
%  Section 7: RNG stream alignment with RMT
%  ---------------------------------------------------------------
% Sections 1 and 2 inject A and S from the reference object, which proves the
% CONSTRUCTION is equivalent but says nothing about how much randomness each class
% consumes. examples/Fig_1_RMT_examples.m seeds once with rng(100) and then chains
% six copy() calls, so reproducing those figures additionally requires RMTBlocks to
% draw from the stream in exactly the same order as RMT:
%
%   constructor : randn(N)  then  rand(N) via update_sparsity
%   copy()      : the constructor's two draws, then another rand(N) because
%                 set.alpha re-runs update_sparsity when A is already populated
%
% That was designed for deliberately; this section is what makes it a guarantee
% rather than an assumption.

fprintf('\n--- 7. RNG stream alignment with RMT ---\n');
if ~has_RMT
    fprintf('  SKIP: class RMT not found on the path.\n');
    n_skip = n_skip + 1;
else
    N7 = 120;

    rng(100); a1 = RMT(N7);
    rng(100); b1 = RMTBlocks(N7);
    [n_pass, n_fail] = report(isequal(a1.A, b1.A) && isequal(a1.S, b1.S), ...
        'constructor consumes the RNG identically', '', n_pass, n_fail);

    a2 = a1.copy();  b2 = b1.copy();
    [n_pass, n_fail] = report(isequal(a2.A, b2.A) && isequal(a2.S, b2.S), ...
        'copy() consumes the RNG identically', '', n_pass, n_fail);

    % Fig_1 chains six copies, so confirm the streams have not drifted apart by the
    % time a later panel draws fresh disorder.
    a3 = a2.copy();  b3 = b2.copy();
    rng(7); a3.update_sparsity();
    rng(7); b3.update_sparsity();
    [n_pass, n_fail] = report(isequal(a3.S, b3.S), ...
        'streams still aligned after a second copy()', '', n_pass, n_fail);

    % The executable form of "the ported figure is the same figure": run Fig_1's
    % exact call sequence on both classes from one seed, with NO injection of A/S.
    N_fig = 600;
    rng(100); W_ref = fig1_sequence(@RMT, N_fig);
    rng(100); W_new = fig1_sequence(@RMTBlocks, N_fig);

    all_ok = true;
    worst = 0;
    for k = 1:numel(W_ref)
        if ~isequal(W_ref{k}, W_new{k})
            all_ok = false;
            worst = max(worst, max(abs(W_ref{k}(:) - W_new{k}(:))));
        end
    end
    [n_pass, n_fail] = report(all_ok, ...
        'Fig_1 call sequence: all 6 panels bit-for-bit from one seed, no injection', ...
        sprintf('worst diff %.3e', worst), n_pass, n_fail);
end

%% ---------------------------------------------------------------
%  Summary
%  ---------------------------------------------------------------
fprintf('\n=================================================\n');
fprintf('  passed: %d   failed: %d   skipped: %d\n', n_pass, n_fail, n_skip);
fprintf('=================================================\n\n');

if n_fail > 0
    error('test_RMTBlocks_equivalence:Failed', '%d assertion(s) failed.', n_fail);
end

%% ---------------------------------------------------------------
%  Local functions
%  ---------------------------------------------------------------
function obj = apply_cfg(obj, c, supports_shift)
% Configure a connectivity object through the SHARED D=2 alias API. The fact that
% this same function drives RMT, RMTMatrix and RMTBlocks is the backward-compat
% claim in executable form.
%
% alpha goes first because its setter redraws the sparsity mask.
obj.alpha = c.alpha;
obj.mu_tilde_e = c.mu_e;
obj.mu_tilde_i = c.mu_i;
obj.sigma_tilde_e = c.sig_e;
obj.sigma_tilde_i = c.sig_i;
obj.f = c.f;
obj.zrs_mode = c.zrs;
if supports_shift
    obj.shift = c.shift;
end
end

function Ws = fig1_sequence(ctor, N)
% Replay the exact construction sequence of Fig_1_RMT_examples.m, in the same
% property-assignment ORDER, since alpha's setter redraws the sparsity mask and any
% reordering would desynchronize the RNG stream. Returns the six panels' W matrices.
Ws = cell(1, 6);
G = cell(1, 6);
E_W = 0.05 / sqrt(N);

% (a) Dense random matrix, unbalanced with global outlier
G{1} = ctor(N);
G{1}.mu_tilde_e = 0 + E_W;
G{1}.mu_tilde_i = 0 + E_W;
G{1}.sigma_tilde_e = 1/sqrt(N);
G{1}.sigma_tilde_i = 1/sqrt(N);
G{1}.f = 1;
G{1}.alpha = 1.0;
G{1}.set_zrs_mode('none');

% (b) Dense balanced shifted
G{2} = G{1}.copy();
G{2}.mu_tilde_e = 0;
G{2}.mu_tilde_i = 0;
R_b = G{2}.R;
G{2}.shift = -R_b;

% (c) Dense balanced Dale's law
G{3} = G{2}.copy();
G{3}.mu_tilde_e = 1/sqrt(N);
G{3}.mu_tilde_i = -1/sqrt(N);
G{3}.sigma_tilde_e = 1/sqrt(N);
G{3}.sigma_tilde_i = 1/sqrt(N);
G{3}.f = 0.5;

% (d) Dense balanced Dale ZRS
G{4} = G{3}.copy();
G{4}.set_zrs_mode('ZRS');

% (e) Dense unbalanced Dale, different sigmas, ZRS
G{5} = G{4}.copy();
G{5}.sigma_tilde_e = 0.35/sqrt(N);
G{5}.sigma_tilde_i = G{5}.compute_sigma_tilde_i_for_target_variance(1/N);

% (f) Sparse unbalanced with Partial SZRS
G{6} = G{5}.copy();
G{6}.set_alpha(0.5);
G{6}.set_zrs_mode('Partial_SZRS');

for k = 1:6
    Ws{k} = G{k}.W;
end
end

function [n_pass, n_fail] = report(ok, name, detail, n_pass, n_fail)
if ok
    n_pass = n_pass + 1;
    fprintf('  PASS  %s\n', name);
else
    n_fail = n_fail + 1;
    if isempty(detail)
        fprintf('  FAIL  %s\n', name);
    else
        fprintf('  FAIL  %s  [%s]\n', name, detail);
    end
end
end
