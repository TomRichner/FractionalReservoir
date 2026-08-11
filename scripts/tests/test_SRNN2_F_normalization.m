% test_SRNN2_F_normalization.m
% Verify the relative RMT parameters and the F_tracks_network convention.
%
% The RMT tildes are stored as multipliers of F = 1/sqrt(n*alpha*(2-alpha)):
%   mu_E_tilde = mu_E_tilde_relative * F,  and likewise for the others.
%
% F_tracks_network = true  (default) computes F from the CURRENT n/indegree,
%                          which makes the theoretical spectral radius R exactly
%                          independent of n.
% F_tracks_network = false pins F to (F_ref_n, F_ref_indegree), so R varies with
%                          network size.
%
% No simulations are run -- everything here is closed-form via the dependent
% properties.
%
% See also: SRNNModel2

fprintf('=== Testing SRNNModel2 RMT normalization ===\n\n');
all_passed = true;

%% Relative multipliers produce the documented absolute tildes
m = SRNNModel2('n', 300, 'indegree', 100);
F = m.default_val;
all_passed = check('F = 1/sqrt(n*alpha*(2-alpha))', ...
    abs(F - 1/sqrt(300 * (100/300) * (2 - 100/300))) < 1e-15) && all_passed;
all_passed = check('class defaults are 3F / -4F / 1F / 1F / 0', ...
    abs(m.mu_E_tilde - 3*F) < 1e-15 && abs(m.mu_I_tilde + 4*F) < 1e-15 && ...
    abs(m.sigma_E_tilde - F) < 1e-15 && abs(m.sigma_I_tilde - F) < 1e-15 && ...
    m.E_W == 0) && all_passed;

m.mu_E_tilde_relative = 3.5;
m.E_W_relative = -0.5;
all_passed = check('setting a relative updates the absolute', ...
    abs(m.mu_E_tilde - 3.5*F) < 1e-15 && abs(m.E_W + 0.5*F) < 1e-15) && all_passed;

%% The absolute names are read-only, and name their replacement
threw = false;
try
    SRNNModel2('mu_E_tilde', 0.2);
catch err
    threw = strcmp(err.identifier, 'SRNNModel:RenamedProperty') && ...
        contains(err.message, 'mu_E_tilde_relative');
end
all_passed = check('constructing with mu_E_tilde errors and names the replacement', threw) && all_passed;

%% Tildes track a changed n (they are recomputed, never cached)
m = SRNNModel2('n', 300, 'indegree', 100);
mu_before = m.mu_E_tilde;
m.n = 600;
expected = 3 / sqrt(600 * (100/600) * (2 - 100/600));
all_passed = check('tildes recompute when n changes', ...
    abs(m.mu_E_tilde - expected) < 1e-15 && m.mu_E_tilde ~= mu_before) && all_passed;

%% A disconnected network (alpha = 0 -> F = Inf) keeps zero multipliers at zero
m = SRNNModel2('n', 1, 'f', 1, 'indegree', 0, ...
    'mu_E_tilde_relative', 0, 'mu_I_tilde_relative', 0, ...
    'sigma_E_tilde_relative', 0, 'sigma_I_tilde_relative', 0);
vals = [m.mu_E_tilde, m.mu_I_tilde, m.sigma_E_tilde, m.sigma_I_tilde, m.E_W];
all_passed = check('indegree=0 gives zeros, not 0*Inf = NaN', ...
    all(vals == 0) && ~any(isnan(vals))) && all_passed;

%% F tracking: R is EXACTLY invariant in n when alpha is held fixed
ns = [100 200 400 800];
Rs = arrayfun(@(k) SRNNModel2('n', k, 'indegree', k/2, 'f', 0.5).R, ns);
all_passed = check(sprintf('R invariant in n at fixed alpha (spread %.2g)', ...
    max(abs(Rs - Rs(1)))), max(abs(Rs - Rs(1))) < 1e-13) && all_passed;

%% Frozen F reproduces tracking F exactly at the reference point
m_track = SRNNModel2('n', 300, 'indegree', 100);
m_froz = SRNNModel2('n', 300, 'indegree', 100, ...
    'F_tracks_network', false, 'F_ref_n', 300, 'F_ref_indegree', 100);
all_passed = check('frozen F == tracking F at the reference network', ...
    abs(m_track.default_val - m_froz.default_val) < 1e-16 && ...
    abs(m_track.R - m_froz.R) < 1e-13) && all_passed;

%% Frozen F stays put across n, and R then does vary with network size
m_small = SRNNModel2('n', 100,  'indegree', 100, 'f', 0.5, 'F_tracks_network', false);
m_large = SRNNModel2('n', 1000, 'indegree', 100, 'f', 0.5, 'F_tracks_network', false);
all_passed = check('frozen F identical across n', ...
    abs(m_small.default_val - m_large.default_val) < 1e-16) && all_passed;
all_passed = check('frozen F lets R vary with n (the opt-in has an effect)', ...
    abs(m_small.R - m_large.R) > 1e-6) && all_passed;

% ... while the network itself still tracks the grid point: freezing F pins the
% weight SCALE, not the connectivity.
all_passed = check('freezing F does not freeze alpha', ...
    abs(m_small.alpha - 1) < 1e-15 && abs(m_large.alpha - 0.1) < 1e-15) && all_passed;

%% An invalid reference network is rejected by build()
threw = false;
try
    m = SRNNModel2('F_tracks_network', false, 'F_ref_n', 100, 'F_ref_indegree', 200, ...
        'T_range', [0 0.1], 'fs', 100, 'lya_method', 'none');
    m.build();
catch err
    threw = contains(err.message, 'F_ref_indegree');
end
all_passed = check('F_ref_indegree > F_ref_n errors from build()', threw) && all_passed;

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
