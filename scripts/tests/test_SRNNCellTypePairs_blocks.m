% test_SRNNCellTypePairs_blocks.m
% Verify the RMT modernization of SRNNCellTypePairs: block statistics via
% RMTBlocks, tildes stored as multiples of F, a selectable F convention, and a
% named activation -- matching what SRNNModel2 gained.
%
% See also: SRNNCellTypePairs, SRNNModel2, RMTBlocks

fprintf('=== Testing SRNNCellTypePairs block RMT ===\n\n');
all_passed = true;
base = {'n', 120, 'indegree', 40, 'n_cellTypes', 2, ...
    'cell_type_names', {'E','I'}, 'f', [0.5 0.5], ...
    'tau_a', {[], []}, 'T_range', [0 0.5], 'fs', 200, ...
    'ode_solver', 'rk4', 'lya_method', 'none', 'rng_seeds', [41 42]};

%% F, and the tildes as multiples of it
m = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4], ...
    'sigma_tilde_relative', [1 1]);
F = m.default_val;
all_passed = check('F = 1/sqrt(n*alpha*(2-alpha))', ...
    abs(F - 1/sqrt(120 * (40/120) * (2 - 40/120))) < 1e-15) && all_passed;
all_passed = check('a 1xC row broadcasts down the columns', ...
    isequal(size(m.mu_tilde), [2 2]) && isequal(m.mu_tilde(1,:), m.mu_tilde(2,:))) && all_passed;
all_passed = check('absolute mu_tilde = relative * F', ...
    abs(m.mu_tilde(1,1) - 3*F) < 1e-15 && abs(m.mu_tilde(1,2) + 4*F) < 1e-15) && all_passed;

%% A full CxC block sets the postsynaptic dependence
mb = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4; 9 -4], ...
    'sigma_tilde_relative', [1 1]);
all_passed = check('a CxC block is accepted, (post, pre)', ...
    abs(mb.mu_tilde(1,1) - 3*F) < 1e-15 && abs(mb.mu_tilde(2,1) - 9*F) < 1e-15) && all_passed;

m.build(); mb.build();
nE = m.n_per_type(1);
dW = abs(full(mb.W) - full(m.W));
blockIE = dW(nE+1:end, 1:nE);            % rows = post (I), cols = pre (E)
rest = dW; rest(nE+1:end, 1:nE) = 0;
all_passed = check('the block changes only the I<-E submatrix of W', ...
    max(blockIE(:)) > 0 && max(rest(:)) == 0) && all_passed;
all_passed = check('W stays sparse', issparse(m.W)) && all_passed;

%% R and lambda_O exist and behave
all_passed = check('R is finite and positive', isfinite(m.R) && m.R > 0) && all_passed;
all_passed = check('lambda_O has one entry per type', numel(m.lambda_O) == 2) && all_passed;
fprintf('\n  mu_EE_rel   R        max Re(lambda_O)\n');
prev = -Inf; mono = true;
for rel = [3 6 12 24]
    mm = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [rel -4; 3 -4], ...
        'sigma_tilde_relative', [1 1]);
    lam = max(real(mm.lambda_O));
    fprintf('  %-11g %-8.4f %.4f\n', rel, mm.R, lam);
    mono = mono && lam > prev; prev = lam;
end
all_passed = check('the outlier grows with mu_EE', mono) && all_passed;

%% F convention
mf = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4], ...
    'sigma_tilde_relative', [1 1], 'F_tracks_network', false, ...
    'F_ref_n', 300, 'F_ref_indegree', 100);
all_passed = check('F_tracks_network=false pins F to the reference network', ...
    abs(mf.default_val - 1/sqrt(300*(1/3)*(2-1/3))) < 1e-15) && all_passed;
all_passed = check('...but alpha still tracks the real network', ...
    abs(mf.alpha - 40/120) < 1e-15) && all_passed;

%% Named activation
for nm = {'logistic', 'piecewise', 'tanh'}
    ma = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4], ...
        'sigma_tilde_relative', [1 1], 'activation', nm{1});
    ok = true;
    try
        ma.build(); ma.run();
    catch
        ok = false;
    end
    all_passed = check(sprintf('a %s model builds and runs', nm{1}), ok) && all_passed;
end

mt = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4], ...
    'sigma_tilde_relative', [1 1], 'activation', 'tanh');
mt.S_c = 0.9;
all_passed = check('tanh ignores S_c', ...
    abs(mt.activation_function(0.5) - tanh(0.5)) < 1e-15) && all_passed;

threw = false;
try
    SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4], ...
        'sigma_tilde_relative', [1 1], 'activation_function', @sin);
catch err
    threw = strcmp(err.identifier, 'SRNNCellTypePairs:RenamedProperty');
end
all_passed = check('setting activation_function errors, naming activation', threw) && all_passed;

%% Both model classes now build the SAME network from matched parameters
m2 = SRNNModel2('n', 120, 'indegree', 40, 'f', 0.5, 'rng_seeds', [41 42], ...
    'T_range', [0 0.1], 'fs', 100, 'lya_method', 'none');
m2.build();
mp = SRNNCellTypePairs(base{:}, 'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1]);
mp.build();
all_passed = check('SRNNModel2 and Pairs build identical W from matched params', ...
    max(abs(full(m2.W) - full(mp.W)), [], 'all') == 0) && all_passed;
all_passed = check('...and agree on R', abs(m2.R - mp.R) < 1e-12) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
