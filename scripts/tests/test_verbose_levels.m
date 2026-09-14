% test_verbose_levels.m - the three-level verbose setting.
%
% Checks:
%   1. verbose_level ranks the three names and the two logicals; a bad name
%      errors verbose_level:BadLevel; verbose_name canonicalises.
%   2. A SRNNCellTypePairs (n = 20) built and run under 'minimal' and
%      'near-none' prints NOTHING; under 'verbose' it prints the old lines.
%      The property accepts a logical too.
%   3. A serial 2 x 2 ParamSpaceAnalysis2 sweep under 'minimal' prints at
%      most 6 lines (one start, one per batch, one end, plus consolidation
%      counts); the models inside it inherit the level (the transcript has
%      no "Integration complete"); resolved_defaults and same_config ignore
%      verbose; an old-style logical verbose loads.
%   4. paper_config().verbose is 'minimal' and a bad level errors at config
%      time.
%
% ~30 s. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: verbose_level, vprintf, paper_config, ParamSpaceAnalysis2

fprintf('=== Testing verbose levels ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';

%% 1. helper
all_passed = check('ranks: verbose 2, minimal 1, near-none 0, true 2, false 1', ...
    verbose_level('verbose') == 2 && verbose_level('minimal') == 1 && verbose_level('near-none') == 0 && ...
    verbose_level(true) == 2 && verbose_level(false) == 1) && all_passed;
all_passed = check('bad level errors verbose_level:BadLevel', ...
    throws_id(@() verbose_level('loud'), 'verbose_level:BadLevel')) && all_passed;
all_passed = check('verbose_name canonicalises', strcmp(verbose_name(true), 'verbose') && ...
    strcmp(verbose_name(false), 'minimal') && strcmp(verbose_name('near-none'), 'near-none')) && all_passed;
out = evalc('vprintf(''minimal'', ''verbose'', ''x\n''); vprintf(''minimal'', ''minimal'', ''y\n'');');
all_passed = check('vprintf gates', strcmp(strtrim(out), 'y')) && all_passed;

%% 2. model
outs = struct();
for lvl = {'minimal', 'near-none', 'verbose'}
    outs.(strrep(lvl{1}, '-', '_')) = evalc(['m = build_from_preset(P, ''sfa1_std1'', ''n'', 20, ''indegree'', 5, ' ...
        '''T_range'', [0 2], ''lya_method'', ''none'', ''verbose'', ''' lvl{1} '''); m.run();']);
end
all_passed = check('model prints nothing at minimal and near-none', ...
    isempty(strtrim(outs.minimal)) && isempty(strtrim(outs.near_none))) && all_passed;
all_passed = check('model prints the build/run lines at verbose', ...
    contains(outs.verbose, 'Integration complete') && contains(outs.verbose, 'built successfully')) && all_passed;
m.verbose = true;
all_passed = check('model property accepts a logical', strcmp(m.verbose, 'verbose')) && all_passed;
out = evalc('m2 = build_from_preset(P, ''sfa1_std1'', ''n'', 20, ''indegree'', 5, ''T_range'', [0 4], ''lya_method'', ''topk'', ''lya_K'', 2, ''lya_T_interval'', [2 4], ''lya_warmup'', 1); m2.run();');
all_passed = check('default level is minimal: a top-K run prints nothing', isempty(strtrim(out))) && all_passed;

%% 3. sweep
[d, cls, conds] = srnn_param_preset(P);
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 2, 'model_class', cls);
psa.use_parallel = false;
psa.output_dir = fullfile(tempdir, 'verbose_levels_test');
if exist(psa.output_dir, 'dir'); rmdir(psa.output_dir, 's'); end
psa.model_defaults = d;
psa.model_defaults.n = 20; psa.model_defaults.indegree = 5;
psa.model_defaults.T_range = [0 3]; psa.model_defaults.fs = 200;
psa.model_defaults.lya_method = 'none'; psa.model_defaults.ode_solver = 'sra1';
psa.add_grid_parameter('level_of_chaos', [0.8, 1.2]);
psa.set_conditions(conds(1:2));
all_passed = check('psa.verbose defaults to minimal', strcmp(psa.verbose, 'minimal')) && all_passed;
out = evalc('psa.run();');
lines = strsplit(strtrim(out), newline);
fprintf('  sweep transcript at minimal (%d lines):\n', numel(lines));
fprintf('    | %s\n', lines{:});
all_passed = check('sweep at minimal prints <= 6 lines, none per model', ...
    numel(lines) <= 6 && ~contains(out, 'Integration complete')) && all_passed;
all_passed = check('verbose is not in resolved_defaults', ~isfield(psa.resolved_defaults, 'verbose')) && all_passed;
psa2 = ParamSpaceAnalysis2.from_dir(psa.output_dir);
psa2.verbose = false;
all_passed = check('same_config ignores verbose; logical setter maps to minimal', ...
    psa.same_config(psa2) && strcmp(psa2.verbose, 'minimal')) && all_passed;
out = evalc('psa2.verbose = ''near-none''; psa3 = ParamSpaceAnalysis2(''n_levels'', 2, ''batch_size'', 2, ''model_class'', cls); psa3.verbose = ''near-none''; psa3.use_parallel = false; psa3.output_dir = fullfile(tempdir, ''verbose_levels_test2''); psa3.model_defaults = psa.model_defaults; psa3.add_grid_parameter(''level_of_chaos'', [0.8 1.2]); psa3.set_conditions(conds(1:2)); psa3.run();');
all_passed = check('sweep at near-none prints <= 2 lines', numel(strsplit(strtrim(out), newline)) <= 2) && all_passed;
rmdir(psa.output_dir, 's'); rmdir(psa3.output_dir, 's');

%% 4. config
cfg = paper_config();
all_passed = check('paper_config().verbose is minimal', strcmp(cfg.verbose, 'minimal')) && all_passed;
all_passed = check('paper_config rejects a bad level', ...
    throws_id(@() paper_config('verbose', 'loud'), 'verbose_level:BadLevel')) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL verbose level TESTS PASSED ===\n');
else
    fprintf('=== SOME verbose level TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
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
