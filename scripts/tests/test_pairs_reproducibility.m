% test_pairs_reproducibility.m - the same SRNNCellTypePairs model built and
% run twice from the same seeds gives bit-identical output, with noise off
% and on, in every condition of the manuscript's network.
%
% Written 2026-11-10 as part of the safety net for changes to the model
% class (docs/ISSUES.md, FD-001). A regression reference is only usable if
% a rerun of unchanged code reproduces it, so this is the property the
% regression tests (FD-003) stand on. It also separates the two things a
% changed number can mean: if this passes and a frozen checksum moves, the
% model changed; if this fails, the run is not a function of its seeds.
%
% The network is the manuscript's preset at n = 40, so the run carries
% everything the paper's runs do: per-neuron setpoints (sigma_S_c), a
% per-neuron SFA ladder (tau_a_spread 0.25), depression on all four routes
% with one component (sfa1_std1) or two (sfa3_std2), and the top-K
% Lyapunov spectrum.
%
% Checks, per condition (no_adaptation, sfa1_std1, sfa3_std2) and per
% integrator ('rk4' at sigma_u_noise = 0; 'sra1' at the preset's 0.025):
%   1. Two fresh objects from the same arguments: W, S0, the input, the
%      setpoint draw, the SFA ladders, t_out, S_out, every field of
%      plot_data and the whole Lyapunov spectrum are identical (isequaln).
%   2. The global random stream at call time does not matter: the second
%      object is built and run after rng('shuffle').
%   3. Calling run() a second time on the same object reproduces its first
%      run.
% And once, on sfa3_std2:
%   4. The seeds do what they say: another network seed changes W; another
%      noise_seed changes the trajectory and leaves W, S0 and the input
%      alone, and reproduces.
%   5. Noise off under 'sra1' is a deterministic run (reproducible, and no
%      noise path is kept), and differs from the noisy run.
%
% lya_results.seconds is wall-clock time and is the one field left out of
% the comparison.
%
% NOT checked, because it is not true and not claimed: build() does not
% leave the caller's global random stream alone. It ends reseeded to
% rng_seeds(2) (generate_external_input calls rng(seed) and nothing
% restores it), so anything the caller draws after build() is the same
% every time. The setpoint, ladder and noise draws do save and restore.
%
% Prints PASS/FAIL per check and a final banner. About 5 s. Assumes
% setup_paths has run.
%
% See also: test_SRNNCellTypePairs_noise (the noise path itself, on a
%           network without depression), test_std_product_model

fprintf('=== Testing that SRNNCellTypePairs runs reproduce bit-for-bit ===\n\n');
all_passed = true;
rng_at_entry = rng;

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium';
[preset, ~, conditions] = srnn_param_preset(P);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
small = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'rng_seeds', [5 6], ...
    'T_range', [0 4], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, 'lya_dt', 0.1, ...
    'lya_T_interval', [2 4], 'lya_warmup', 1};
integrators = { ...
    'rk4, noise off',              {'sigma_u_noise', 0, 'ode_solver', 'rk4'}
    sprintf('sra1, noise %g', preset.sigma_u_noise), {'ode_solver', 'sra1'}};

%% 1-3. Per condition and integrator
for i = 1:numel(names)
    for s = 1:size(integrators, 1)
        args = [small, integrators{s, 2}];
        label = sprintf('%-13s %s', names{i}, integrators{s, 1});

        rng(1, 'twister');
        A = run_model(P, names{i}, args);
        rng('shuffle');                         % check 2: the global stream is irrelevant
        B = run_model(P, names{i}, args);

        [same, what] = same_run(A, B);
        all_passed = check(sprintf('%s: two fresh objects identical%s', label, what), same) && all_passed;
        all_passed = check(sprintf('%s: the run did something (finite, states moved, lambda_1 = %+.3f)', ...
            label, A.lya_results.LLE), ...
            all(isfinite(A.S_out), 'all') && any(A.S_out(end, :) ~= A.S_out(1, :)) && ...
            isfinite(A.lya_results.LLE)) && all_passed;

        S_first = A.S_out; lya_first = untimed(A.lya_results); pd_first = A.plot_data;
        evalc('A.run();');
        all_passed = check(sprintf('%s: run() twice on one object reproduces itself', label), ...
            isequal(A.S_out, S_first) && isequaln(untimed(A.lya_results), lya_first) && ...
            isequaln(A.plot_data, pd_first)) && all_passed;
    end
end

%% 4. The seeds do what they say
args = [small, {'ode_solver', 'sra1'}];
C = run_model(P, 'sfa3_std2', args);
D = run_model(P, 'sfa3_std2', [args, {'rng_seeds', [7 6]}]);
all_passed = check('another network seed gives another W', ~isequal(C.W, D.W)) && all_passed;

E = run_model(P, 'sfa3_std2', [args, {'noise_seed', 99}]);
F = run_model(P, 'sfa3_std2', [args, {'noise_seed', 99}]);
all_passed = check('another noise_seed: same W, S0 and input, different trajectory', ...
    isequal(E.W, C.W) && isequal(E.S0, C.S0) && isequal(E.u_ex, C.u_ex) && ...
    ~isequal(E.S_out, C.S_out)) && all_passed;
all_passed = check('an explicit noise_seed reproduces', isequal(E.S_out, F.S_out)) && all_passed;

%% 5. Noise off under sra1
G1 = run_model(P, 'sfa3_std2', [args, {'sigma_u_noise', 0}]);
G2 = run_model(P, 'sfa3_std2', [args, {'sigma_u_noise', 0}]);
all_passed = check('sra1 at sigma_u_noise = 0 reproduces and differs from the noisy run', ...
    isequal(G1.S_out, G2.S_out) && ~isequal(G1.S_out, C.S_out)) && all_passed;
all_passed = check('no noise path is left on the object after run()', ...
    isempty(C.noise_increments) && isempty(G1.noise_increments)) && all_passed;

rng(rng_at_entry);
fprintf('\n');
if all_passed
    fprintf('=== ALL reproducibility TESTS PASSED ===\n');
else
    fprintf('=== SOME reproducibility TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function m = run_model(P, cond, args) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, args{:}); m.run();');
end

function [same, what] = same_run(A, B)
% Compare everything a run produces. `what` names the first field that
% differs, so a failure says where to look.
fields = {'W', 'S0', 'u_ex', 't_ex', 'S_c_vec', 'tau_a_matrix', 't_out', 'S_out', ...
    'plot_data', 'lya_results'};
same = true; what = '';
for k = 1:numel(fields)
    a = A.(fields{k}); b = B.(fields{k});
    if strcmp(fields{k}, 'lya_results'); a = untimed(a); b = untimed(b); end
    if ~isequaln(a, b)
        same = false;
        what = sprintf(' (first difference: %s)', fields{k});
        return
    end
end
end

function r = untimed(r)
% lya_results without its wall-clock field.
if isfield(r, 'seconds'); r = rmfield(r, 'seconds'); end
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
