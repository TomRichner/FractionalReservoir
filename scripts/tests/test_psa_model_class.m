% test_psa_model_class.m
% Verify ParamSpaceAnalysis2 can drive a model class other than SRNNModel2.
%
% PSA needs a constructor and a property list, not a common base class, so the
% model classes are duck-typed through the model_class property rather than
% sharing a hierarchy. This test pins the three things that had to generalize:
%   1. construction/introspection by class name,
%   2. conditions that carry ARBITRARY fields (a whole synapse_config struct,
%      not just SRNNModel2's four adaptation counts),
%   3. metric extraction over any number of cell types, and under either name
%      for the synaptic output ('br' vs 'synaptic_output').
%
% See also: ParamSpaceAnalysis2, SRNNModel2, SRNNCellTypePairs

ok = true;

%% 1. The default path still works and records the class
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 10, 'note', 'p3a', 'verbose', false);
psa.use_parallel = false;
psa.output_dir = fullfile(tempdir, 'psa_p3');
psa.set_conditions({struct('name','sfa_and_std','n_a_E',3,'n_b_E',1)});
psa.add_grid_parameter('level_of_chaos', [1, 1.4]);
psa.model_defaults.n = 24; psa.model_defaults.indegree = 8;
psa.model_defaults.T_range = [0 2]; psa.model_defaults.fs = 200;
psa.model_defaults.ode_solver = 'rk4';
psa.run();
r = psa.results.sfa_and_std;
lle_default = cellfun(@(x) x.LLE, r(~cellfun(@isempty, r)));
ok = rep('SRNNModel2 sweep runs', numel(lle_default) == 2) && ok;
ok = rep('resolved_defaults holds SRNNModel2-shaped params', ...
    isfield(psa.resolved_defaults, 'mu_E_tilde_relative')) && ok;
ok = rep('mean_rate finite', all(isfinite(cellfun(@(x) x.mean_rate, r(~cellfun(@isempty,r)))))) && ok;
ok = rep('mean_synaptic_output finite (br)', ...
    all(isfinite(cellfun(@(x) x.mean_synaptic_output, r(~cellfun(@isempty,r)))))) && ok;

%% 2. A SRNNCellTypePairs sweep: mu_EE swept, STD on E->E only vs none
sc_std = struct(); sc_std.E.E.std = struct('tau_rec', 2.0, 'tau_rel', 0.15);
p2 = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 10, 'note', 'p3b', 'verbose', false);
p2.model_class = 'SRNNCellTypePairs';
p2.use_parallel = false;
p2.output_dir = fullfile(tempdir, 'psa_p3b');
p2.integer_params = {'n', 'indegree'};
% Conditions differ by a whole synapse_config STRUCT -- impossible before.
p2.set_conditions({ ...
    struct('name', 'no_std',   'synapse_config', struct()), ...
    struct('name', 'std_EE',   'synapse_config', sc_std)});
p2.add_grid_parameter('mu_EE_relative', [3, 9]);
p2.model_defaults.n = 24;
p2.model_defaults.indegree = 8;
p2.model_defaults.n_cellTypes = 2;
p2.model_defaults.cell_type_names = {'E','I'};
p2.model_defaults.f = [0.5 0.5];
p2.model_defaults.mu_tilde_relative = [3 -4];
p2.model_defaults.sigma_tilde_relative = [1 1];
p2.model_defaults.n_a = [0 0];
p2.model_defaults.T_range = [0 2];
p2.model_defaults.fs = 200;
p2.model_defaults.ode_solver = 'rk4';
try
    p2.run();
    ran = true;
catch err
    ran = false;
    fprintf('   Pairs sweep error: %s\n', err.message);
end
ok = rep('Pairs sweep runs', ran) && ok;
if ran
    for cn = {'no_std','std_EE'}
        rr = p2.results.(cn{1});
        v = rr(~cellfun(@isempty, rr));
        fprintf('   %-8s LLE = %s   mean_rate = %s\n', cn{1}, ...
            mat2str(round(cellfun(@(x) x.LLE, v)', 4)), ...
            mat2str(round(cellfun(@(x) x.mean_rate, v)', 4)));
    end
    ok = rep('resolved_defaults holds Pairs-shaped params', ...
        isfield(p2.resolved_defaults, 'mu_tilde_relative') && ...
        ~isfield(p2.resolved_defaults, 'mu_E_tilde_relative')) && ok;
    ok = rep('synaptic output found under either name', ...
        all(isfinite(cellfun(@(x) x.mean_synaptic_output, ...
            p2.results.std_EE(~cellfun(@isempty, p2.results.std_EE)))))) && ok;
end

%% 3. Different classes must not pool
[tf, why] = psa.same_config(p2);
ok = rep(sprintf('same_config refuses across classes (%s)', why), ...
    ~tf && contains(why, 'model_class')) && ok;

fprintf('\n%s\n', tern(ok, 'PHASE 3 CHECKS PASS', 'FAILURES ABOVE'));

function p = rep(name, cond)
p = cond; fprintf('  [%s] %s\n', tern(cond,'PASS','FAIL'), name);
end
function s = tern(c, a, b), if c, s = a; else, s = b; end, end
