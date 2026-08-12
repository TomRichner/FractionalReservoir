% test_SRNN2_blocks.m
% Verify SRNNModel2's block connectivity: the RMT statistics are indexed by
% BOTH postsynaptic and presynaptic type, so E->E can differ from E->I.
%
%   mu_EE = E receives from E    mu_EI = E receives from I
%   mu_IE = I receives from E    mu_II = I receives from I
%
% Block overrides are empty by default, in which case the column shorthands
% (mu_E_tilde_relative / mu_I_tilde_relative) apply to every postsynaptic type
% and the model reduces exactly to the classical Harris (2023) setting.
%
% See also: SRNNModel2, RMTBlocks

fprintf('=== Testing SRNNModel2 block connectivity ===\n\n');
all_passed = true;
args = {'n', 200, 'indegree', 60, 'f', 0.5, 'rng_seeds', [11 12], ...
    'T_range', [0 0.1], 'fs', 100, 'lya_method', 'none'};

%% With no overrides the blocks are column-uniform
m = SRNNModel2(args{:});
F = m.default_val;
mu = m.mu_tilde_block;
all_passed = check('no overrides -> columns uniform', ...
    mu(1,1) == mu(2,1) && mu(1,2) == mu(2,2)) && all_passed;
all_passed = check('column 1 is mu_E_tilde_relative * F', ...
    abs(mu(1,1) - 3*F) < 1e-15) && all_passed;
all_passed = check('column 2 is mu_I_tilde_relative * F', ...
    abs(mu(1,2) + 4*F) < 1e-15) && all_passed;

%% A block override wins over the column shorthand, and only for its block
m2 = SRNNModel2(args{:}, 'mu_EE_relative', 9);
mu2 = m2.mu_tilde_block;
all_passed = check('mu_EE_relative sets only the (E<-E) block', ...
    abs(mu2(1,1) - 9*F) < 1e-15 && abs(mu2(2,1) - 3*F) < 1e-15 && ...
    isequal(mu2(:,2), mu(:,2))) && all_passed;

%% ...and it changes only the corresponding block of W
m.build(); m2.build();
nE = m.n_E;
dW = abs(full(m2.W) - full(m.W));
blockEE = dW(1:nE, 1:nE);              % rows = post, cols = pre
rest = dW; rest(1:nE, 1:nE) = 0;
all_passed = check('W changes in the E<-E block only', ...
    max(blockEE(:)) > 0 && max(rest(:)) == 0) && all_passed;

%% Raising mu_EE drives a positive outlier eigenvalue out past the bulk
fprintf('\n  mu_EE_rel   R        max Re(lambda_O)   outlier beyond R?\n');
prev = -Inf;
mono = true;
for rel = [3 6 12 24]
    mm = SRNNModel2(args{:}, 'mu_EE_relative', rel);
    lam = max(real(mm.lambda_O));
    fprintf('  %-11g %-8.4f %-18.4f %s\n', rel, mm.R, lam, ...
        tern(lam > mm.R, 'yes', 'no'));
    mono = mono && lam > prev;
    prev = lam;
end
all_passed = check('the dominant outlier grows with mu_EE', mono) && all_passed;

%% Sigma blocks work the same way
m3 = SRNNModel2(args{:}, 'sigma_II_relative', 5);
sg = m3.sigma_tilde_block;
all_passed = check('sigma_II_relative sets only the (I<-I) block', ...
    abs(sg(2,2) - 5*F) < 1e-15 && abs(sg(1,2) - F) < 1e-15) && all_passed;
all_passed = check('a sigma block raises R', m3.R > m.R) && all_passed;

%% R still reduces to the classical formula when column-uniform
m4 = SRNNModel2(args{:});
sse = m4.sigma_se; ssi = m4.sigma_si;
R_classic = sqrt(m4.n * (m4.f * sse^2 + (1 - m4.f) * ssi^2)) * m4.level_of_chaos;
all_passed = check('R matches the column-uniform Harris Eq. (18)', ...
    abs(m4.R - R_classic) < 1e-12) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

function s = tern(c, a, b), if c, s = a; else, s = b; end, end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
