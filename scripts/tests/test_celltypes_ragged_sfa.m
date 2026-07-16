% test_celltypes_ragged_sfa.m - Ragged per-cell-type multi-timescale SFA
%
% Verifies the ragged SFA layout in SRNNModelCellTypes: per-type timescale counts
% n_a (K-vector, 0-3) with distinct logspaced taus, DC-gain-split strength (c/n_a), no
% orphan SFA state, backward compatibility at the single-timescale default, and the
% deferred-Jacobian guard for multi-timescale. Prints CELLTYPES_RAGGED_SFA_PASS/FAIL.

close all; clear; clc;
setup_paths();

pass = true;
act  = @(x) SRNNModelBase.piecewiseSigmoid(x, 0.9, 0.35, 1);
actd = @(x) SRNNModelBase.piecewiseSigmoidDerivative(x, 0.9, 0.35, 1);
mk = @(varargin) SRNNModelCellTypes('n',120,'fs',500,'T_range',[0 3], ...
        'lya_method','benettin','rng_seeds',[1 2], ...
        'activation_function',act,'activation_function_derivative',actd, varargin{:});

%% 1. Regression: default (n_a=1) -> one pool per adapting neuron, no orphan state
m = mk(); m.build();
if m.cached_params.sfa_len ~= m.n_ad
    fprintf('FAIL: default sfa_len=%d ~= n_ad=%d\n', m.cached_params.sfa_len, m.n_ad); pass = false;
end
m.run();
if ~isfinite(m.lya_results.LLE)
    fprintf('FAIL: default LLE not finite\n'); pass = false;
end
fprintf('1. default n_a=[%s] sfa_len=%d n_ad=%d LLE=%.4f\n', ...
    num2str(m.cached_params.n_a'), m.cached_params.sfa_len, m.n_ad, m.lya_results.LLE);

%% 2. c_gain=0 -> SFA drive is exactly zero (x_eff == x, r == phi(x))
m0 = mk('c_gain',0); m0.build(); m0.run();
P0 = m0.cached_params;
drive0 = 0;
for s = P0.sfa, drive0 = max(drive0, abs(s.c)); end
if drive0 ~= 0
    fprintf('FAIL: c_gain=0 but a block has nonzero c (%.3g)\n', drive0); pass = false;
end
fprintf('2. c_gain=0 -> max block c = %.3g (expect 0)\n', drive0);

%% 3. Ragged layout: n_a=[3 0 2 1], distinct logspaced taus, no orphan state
m2 = mk('n_a',[3 0 2 1],'tau_a_range',[0.05 2]); m2.build();
P = m2.cached_params; sfa = P.sfa;
cnt = arrayfun(@(s) s.count, sfa); na = arrayfun(@(s) s.n_a, sfa); typ = arrayfun(@(s) s.type, sfa);
if P.sfa_len ~= sum(cnt.*na)
    fprintf('FAIL: sfa_len=%d ~= sum(count*n_a)=%d\n', P.sfa_len, sum(cnt.*na)); pass = false;
end
if ~isequal(typ(:)', [1 3 4])   % type 2 (Pvalb, n_a=0) must carry NO block
    fprintf('FAIL: block types=[%s] (expected [1 3 4], Pvalb absent)\n', num2str(typ)); pass = false;
end
for s = sfa
    if numel(unique(s.tau)) ~= s.n_a
        fprintf('FAIL: type %d taus not distinct: [%s]\n', s.type, num2str(s.tau)); pass = false;
    end
end
m2.run();
if ~isfinite(m2.lya_results.LLE)
    fprintf('FAIL: ragged LLE not finite\n'); pass = false;
end
f = m2.plot_by_celltype(); close(f);   % plotting must survive ragged SFA
fprintf('3. ragged types=[%s] counts=[%s] n_a=[%s] sfa_len=%d LLE=%.4f\n', ...
    num2str(typ), num2str(cnt), num2str(na), P.sfa_len, m2.lya_results.LLE);

%% 4. DC-gain split: bl.c * n_a == c_gain * adapt_index(type) (invariant to n_a)
for s = m2.cached_params.sfa
    lhs = s.c * s.n_a; rhs = m2.c_gain * m2.adapt_index(s.type);
    if abs(lhs - rhs) > 1e-12
        fprintf('FAIL: DC-split type %d: c*n_a=%.6g ~= c_gain*ai=%.6g\n', s.type, lhs, rhs); pass = false;
    end
end
fprintf('4. DC-split c*n_a == c_gain*adapt_index for all blocks\n');

%% 5. Jacobian guard: multi-timescale errors; single-timescale works (via plot_eigenvalues)
mg = mk('n_a',[2 0 0 0],'store_full_state',true); mg.build(); mg.run();
guard_ok = false;
try
    mg.plot_eigenvalues(1.0);
catch e
    guard_ok = strcmp(e.identifier, 'SRNNModelCellTypes:JacobianDeferredMultiSFA');
end
if ~guard_ok
    fprintf('FAIL: multi-timescale Jacobian guard did not fire\n'); pass = false;
end
close all force;
ms = mk('n_a',[1 0 1 1],'store_full_state',true); ms.build(); ms.run();
single_ok = false;
try
    fe = ms.plot_eigenvalues(1.0); single_ok = isgraphics(fe); close(fe);
catch e
    fprintf('FAIL: single-timescale Jacobian errored: %s\n', e.message);
end
if ~single_ok, pass = false; end
fprintf('5. Jacobian guard fired for multi-SFA=%d, single-SFA works=%d\n', guard_ok, single_ok);

%% Result
close all force;
if pass
    fprintf('\nCELLTYPES_RAGGED_SFA_PASS\n');
else
    fprintf('\nCELLTYPES_RAGGED_SFA_FAIL\n');
end
