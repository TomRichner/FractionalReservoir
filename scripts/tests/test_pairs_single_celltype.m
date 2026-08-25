% TEST_PAIRS_SINGLE_CELLTYPE One-cell-type support in SRNNCellTypePairs.
%
% Covers the 2026-08-23 change that let build_network configure its RMTBlocks
% generator with set_types instead of three separate assignments, which is what
% makes n_cellTypes = 1 buildable.
%
% THE FIRST SECTION IS THE LOAD-BEARING ONE. The change had to be a no-op at
% D = 2, because every sweep preset, every committed figure and every archived
% run directory depends on the weight matrices being unchanged. The checksums
% below were captured by running the PRE-CHANGE code; they are not derived from
% anything, so they cannot drift with it. If one fails, the connectivity moved.
%
% sompolinsky_pairs is IN the table even though it moved from two cell types to
% one in the same commit. That was expected to redraw W -- changing D changes how
% the generator partitions its blocks -- and it does not: its two types had
% identical zero-mean statistics, so the per-block scaling was uniform and the
% underlying draw never depended on the partition. Measured identical to twelve
% significant figures over 4856 nonzeros. Figure 1 panel A is therefore
% unchanged, and this row is what keeps that true.
%
% single_neuron_stf and single_neuron_dualStd are exempt: W == 0 either way, so
% a checksum proves nothing. Their shape is checked further down instead.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs, RMTBlocks/set_types, srnn_adaptation_conditions,
%           singleCellTypeRefactor.md

all_passed = true;

%% D = 2 invariance: frozen pre-change checksums
% n and indegree are overridden to keep this fast; deterministic either way, and
% the frozen numbers were captured with the same overrides.
fprintf('\n-- D = 2 invariance (frozen pre-change checksums) --\n');
frozen = { ...
    'celltype_pairs',                                                 190.602117992,  38.2392507317; ...
    'celltype_pairs_S_c_by_type',                                     518.014669734,  45.8233392159; ...
    'celltype_pairs_S_c_by_type_n500',                                518.014669734,  45.8233392159; ...
    'celltype_pairs_S_c_by_type_n500_fixedF',                         315.253688591,  27.8871961653; ...
    'celltype_pairs_all_std_n500',                                    315.253688591,  27.8871961653; ...
    'celltype_pairs_uniform_std_n500',                               -1.13502310520,  30.1299675811; ...
    'celltype_pairs_uniform_std_n500_mu5p5',                         -1.03433758718,  29.1467059839; ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive',                 -1.03433758718,  29.1467059839; ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5',          -1.14156058489,  29.7665779012; ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p02',-1.14156058489,  29.7665779012; ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01',-1.14156058489,  29.7665779012; ...
    'celltype_pairs_Sc0p2_noise0p025',                               -1.14156058489,  29.7665779012; ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd',                       -1.14156058489,  29.7665779012; ...
    'bursting_pairs',                                                 891.629319107,  24.9210947937; ...
    'sompolinsky_pairs',                                            -0.563793163494, 13.5057856038};

for k = 1:size(frozen, 1)
    name = frozen{k, 1};
    W = preset_W(name, {'n', 120, 'indegree', 40});
    ok = abs(sum(W(:)) - frozen{k, 2}) < 1e-9 && ...
         abs(norm(W, 'fro') - frozen{k, 3}) < 1e-9 && nnz(W) == 4856;
    all_passed = check(sprintf('%-62s W unchanged', name), ok) && all_passed;
end

%% C = 1 builds, and its RMT maths is right
fprintf('\n-- C = 1 builds --\n');
gamma = 1.6;
s = one_type_model(gamma);
s.build();
s.run();          % the plotting checks below need states, not just a network
all_passed = check('C = 1 model builds and runs', true) && all_passed;
all_passed = check(sprintf('R == level_of_chaos exactly (%.4f)', s.R), ...
    abs(s.R - gamma) < 1e-12) && all_passed;
frac_neg = mean(nonzeros(s.W) < 0);
all_passed = check(sprintf('zero-mean -> ~50%% negative weights (%.1f%%)', 100*frac_neg), ...
    abs(frac_neg - 0.5) < 0.05) && all_passed;
all_passed = check('n_per_type holds every neuron', isequal(s.n_per_type, s.n)) && all_passed;

%% C = 1 at n = 1: a genuinely single neuron with all three mechanisms
fprintf('\n-- C = 1, n = 1, SFA + dual STD + STF --\n');
q = single_neuron_model();
q.build();
q.run();
p = q.plot_data;
% 3 SFA + 2 STD + 1 STF + 1 x
all_passed = check('N_sys_eqs == 7', q.N_sys_eqs == 7) && all_passed;
all_passed = check('W == 0 (1x1)', isequal(full(q.W), 0)) && all_passed;
all_passed = check('x.E is 1 x nt',   size(p.x.E, 1) == 1) && all_passed;
all_passed = check('a.E is 1 x 3 x nt', isequal(size(p.a.E, 1), 1) && size(p.a.E, 2) == 3) && all_passed;
all_passed = check('b.E.E is 1 x 2 x nt', size(p.b.E.E, 2) == 2) && all_passed;
all_passed = check('g.E.E is 1 x 1 x nt', size(p.g.E.E, 2) == 1) && all_passed;
% The step must actually engage the mechanisms, or the figure it backs is empty.
b_min = min(squeeze(prod(p.b.E.E(1, :, :), 2)));
g_max = max(squeeze(prod(p.g.E.E(1, :, :), 2)));
all_passed = check(sprintf('dual STD depresses (min prod(b) = %.3f < 0.5)', b_min), ...
    b_min < 0.5) && all_passed;
all_passed = check(sprintf('STF facilitates (max prod(g) = %.3f > 1.5)', g_max), ...
    g_max > 1.5) && all_passed;

%% Plotting at C = 1
% These were all unreachable while build() failed, so they were never exercised.
fprintf('\n-- plotting at C = 1 --\n');
for m = {'plot', 'plot_celltypes', 'plot_W', 'plot_W_spectrum', 'plot_weight_histogram'}
    ok = true;
    try
        s.(m{1})();
    catch ME
        ok = false;
        fprintf('    %s\n', ME.message);
    end
    all_passed = check(sprintf('%s runs', m{1}), ok) && all_passed;
end
close all force;

%% The two-type aliases must still refuse
% They are what enforces "C = 1 is for figures, C >= 2 is for sweeps".
fprintf('\n-- scalar aliases at C = 1 --\n');
for a = {'f_E', 'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'}
    all_passed = check(sprintf('%s errors NotTwoTypes', a{1}), ...
        throws_id(@() s.(a{1}), 'SRNNCellTypePairs:NotTwoTypes')) && all_passed;
end
% tau_a_E is different: it asserts only that type 1 is named E, so it works at
% any number of types. q is named E; s is named 'all'.
all_passed = check('tau_a_E works at C = 1 when type 1 is E', ...
    isequal(q.tau_a_E, q.tau_a{1})) && all_passed;

%% Adaptation conditions carry the right row length
fprintf('\n-- srnn_adaptation_conditions --\n');
c1 = srnn_adaptation_conditions('SRNNCellTypePairs', [], srnn_sfa_timescales(3), 1);
all_passed = check('C = 1 gives 1-element n_a rows', ...
    all(cellfun(@(c) isscalar(c.n_a),c1))) && all_passed;
all_passed = check('C = 1 sfa_only is n_a = 3', isequal(c1{2}.n_a, 3)) && all_passed;

c2 = srnn_adaptation_conditions('SRNNCellTypePairs', [], srnn_sfa_timescales(3), 2);
c2_default = srnn_adaptation_conditions('SRNNCellTypePairs', [], srnn_sfa_timescales(3));
all_passed = check('C = 2 unchanged from the default', isequal(c2, c2_default)) && all_passed;
all_passed = check('C = 2 sfa_only is still [3 0]', isequal(c2{2}.n_a, [3 0])) && all_passed;

all_passed = check('SRNNModel2 refuses a non-2 type count', ...
    throws_id(@() srnn_adaptation_conditions('SRNNModel2', [], srnn_sfa_timescales(3), 1), ...
              'srnn_adaptation_conditions:NotTwoTypes')) && all_passed;

%% The C = 1 presets
fprintf('\n-- presets that are now C = 1 --\n');
for nm = {'sompolinsky_pairs', 'single_neuron_stf', 'single_neuron_dualStd'}
    [d, cls, cond] = srnn_param_preset(nm{1});
    ok = strcmp(cls, 'SRNNCellTypePairs') && d.n_cellTypes == 1 && ...
         all(cellfun(@(c) isscalar(c.n_a),cond));
    all_passed = check(sprintf('%-22s is C = 1 with 1-element conditions', nm{1}), ok) && all_passed;
end
% The single-neuron cartoon must be ONE neuron. This is the check that would
% have caught fig_adaptation_methods building the 500-neuron network.
for nm = {'single_neuron_stf', 'single_neuron_dualStd'}
    d = srnn_param_preset(nm{1});
    all_passed = check(sprintf('%-22s is n = 1 with zero connectivity', nm{1}), ...
        d.n == 1 && all(d.mu_tilde_relative(:) == 0) && ...
        all(d.sigma_tilde_relative(:) == 0)) && all_passed;
end
% And it must carry the paper's adaptation, not its own invented values.
paper = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd');
sn    = srnn_param_preset('single_neuron_dualStd');
all_passed = check('single_neuron_dualStd inherits the paper''s c', ...
    abs(sn.c - paper.c(1)) < 1e-12) && all_passed;
[~, ~, paper_cond] = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd');
[~, ~, sn_cond]    = srnn_param_preset('single_neuron_dualStd');
all_passed = check('single_neuron_dualStd inherits the dual STD timescales', ...
    isequal(sn_cond{4}.synapse_config.E.E.std, ...
            paper_cond{4}.synapse_config.E.E.std)) && all_passed;

%% ------------------------------------------------------------------------
if all_passed
    fprintf('\ntest_pairs_single_celltype: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_pairs_single_celltype: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function W = preset_W(preset_name, extra)
% build() is chatty; the checksum is all we want here, so evalc swallows it.
% The argument list is assembled OUTSIDE the evalc string -- the alternative
% interpolates preset_name into code, which the Code Analyzer cannot see through.
args = [{preset_name, 'sfa_and_std'}, extra, ...
        {'lya_method', 'none', 'T_range', [0 1], 'rng_seeds', [42 42]}]; %#ok<NASGU>
evalc('m = build_from_preset(args{:});');
W = full(m.W);
end

function m = one_type_model(gamma)
% A Sompolinsky-shaped single-population network: Dale-free, zero mean, tanh.
m = SRNNCellTypePairs('n_cellTypes', 1, 'cell_type_names', {'all'}, 'f', 1, ...
    'mu_tilde_relative', 0, 'sigma_tilde_relative', 1, ...
    'n', 200, 'indegree', 200, 'level_of_chaos', gamma, ...
    'activation', 'tanh', 'mu_S_c', [], 'sigma_S_c', [], 'tau_d', 1.0, ...
    'c', 0, 'n_a', 0, 'synapse_config', struct(), 'x0_std', 1.0, ...
    'ode_solver', 'rk4', 'fs', 100, 'lya_method', 'none', 'T_range', [0 2], ...
    'u_ex_scale', 0, 'rng_seeds', [7 1], ...
    'store_decimated_state', true, 'plot_freq', 100);   % plot() needs the states
end

function m = single_neuron_model()
% One neuron, no recurrence, every mechanism on, driven by a step.
sc = struct('E', struct('E', struct( ...
    'std', struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]), ...
    'stf', struct('tau_dec', 6, 'tau_fac', 2.5, 'G', 1/0.35))));
ic = struct('intrinsic_drive', 0, 'generator', @step_gen);
m = SRNNCellTypePairs('n_cellTypes', 1, 'cell_type_names', {'E'}, 'f', 1, ...
    'mu_tilde_relative', 0, 'sigma_tilde_relative', 0, 'n', 1, 'indegree', 1, ...
    'activation', 'piecewise', 'S_a', 0.8, 'S_c', 0.2, 'tau_d', 0.1, ...
    'c', 0.5/3, 'n_a', 3, 'synapse_config', sc, 'x0_std', 0, ...
    'input_config', ic, 'ode_solver', 'rk4', 'fs', 200, ...
    'lya_method', 'none', 'T_range', [0 25], 'plot_deci', 1);
end

function [u_ex, t_ex] = step_gen(params, T, fs, ~, ~)
% Signature fixed by SRNNCellTypePairs: generator(params, T, fs, seed, config).
t_ex = (0:1/fs:T)';
u_ex = zeros(params.n, numel(t_ex));
u_ex(:, (t_ex >= 5) & (t_ex < 15)) = 0.5;
end

function tf = throws_id(f, id)
tf = false;
try
    f();
catch ME
    tf = strcmp(ME.identifier, id);
end
end
