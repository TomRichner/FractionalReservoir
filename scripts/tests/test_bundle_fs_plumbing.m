function test_bundle_fs_plumbing()
% TEST_BUNDLE_FS_PLUMBING The bundle config's cfg.fs and cfg.ode_solver reach every stage (FD-031 item 14).
%
% Tom (2026-11-17): every simulation of every paper analysis uses SRA1 at
% 500 Hz, stated once in the bundle config (mu7chain_config's cfg.fs and
% cfg.ode_solver, inherited by mu7chain_AC_config and
% mu7chain_gainbias_config). The one exception is the numerics stage, which
% compares solvers and rates on purpose.
%
% PART 1 is a DRY RUN of run_all_paper_analyses: the real entry point, the
% real run_all_analyses and resolve_run_context, the real opt-in stage plans
% (mc_gain_sweep_stage_plan, noise_switch_stage_plan), with every stage that
% simulates replaced by a stub (a folder put at the front of the path) that
% records what it was called with. Nothing is simulated.
%
% Checked, at 'fast' and 'production', for mu7chain_config,
% mu7chain_AC_config and mu7chain_gainbias_config:
%   * the three sweeps' resolved model_defaults are fs 500 and ode_solver
%     'sra1', i.e. the analysis_run_config mode table's are replaced;
%   * eig_heatmap, lyapunov_spectrum, transient_gain, local_lyapunov,
%     dc_lle, memory_capacity and mc_example are called with 'fs', 500 and
%     'ode_solver', 'sra1';
%   * paper_illustrations gets the bundle cfg carrying both;
%   * every run of block C's mc_gain_sweep stage is planned with
%     'ode_solver', 'sra1' and a protocol with fs = 500;
%   * the noise_switch stage is planned with 'fs', 500, 'ode_solver',
%     'sra1' and storage rates that divide 500;
%   * the opt-in stages built after the plumbing (the merge of
%     2026-11-21): lyapunov_spectrum_n80, outlier_spectra and
%     outlier_structures are planned at fs 500 with 'sra1' (their plans
%     read cfg.fs and cfg.ode_solver); frequency_response and
%     eigenvalue_density simulate nothing and are only checked to be
%     reached. All five are stubbed like the others: the dry run of a
%     mu7chain bundle would otherwise run them for real, at 500 neurons;
%     and Plan 2's concept_examples and uncoupled_reference, stubbed the
%     same way and checked at fs 500 with 'sra1' (until they were, the
%     dry runs ran both for real, 29 minutes);
%   * numerics, the one stage that keeps its own rates and solvers, gets
%     neither;
%   * run_manifest records fs = 500 and ode_solver = 'sra1';
%   * the two registry entries that simulate while drawing
%     (fig_introductory_concepts, fig_stim_engages_adaptation) carry both.
% And with cfg.fs and cfg.ode_solver removed, no stage gets either (memory
% capacity keeps cfg.mc_ode_solver) and the sweeps keep the mode table's
% rate and solver, as before the fields existed; the three opt-in stages
% above fall back to their own 500 Hz and name no integrator. A stage's own
% knob (cfg.outlier_spectra_fs, cfg.outlier_structures_fs) overrides the
% bundle's rate for that stage.
%
% PART 2 checks that the stages USE what they are handed: build_from_preset
% is stubbed to record its arguments and stop, and the real
% run_eig_heatmap, run_local_lyapunov, run_paper_illustrations,
% fig_introductory_concepts and fig_stim_engages_adaptation are called with
% fs 500 and 'sra1'; each must hand both to the model it builds. (Those
% build in the client; the stages that build inside a parfor are covered by
% part 1 and the pilots; memory capacity, block C and the noise-switch
% rollouts construct their models directly and are covered the same way.)
%
% Finally stage_fs and stage_solver themselves.
%
% Quick: a few seconds.

global FS_PLUMB_LOG %#ok<GVMIS>   the stubs append here
global FS_PLUMB_STAGES_REAL %#ok<GVMIS>   false: the stage stubs record and return (part 1)
FS_PLUMB_STAGES_REAL = false;   % never inherit a true left by an interrupted run

stub_dir = tempname;
mkdir(stub_dir);
write_stubs(stub_dir);
addpath(stub_dir, '-begin');
% The stubs' folder is also the current folder for the test: the current
% folder outranks the path, and the session's is usually the project root,
% whose own setup_paths.m would otherwise win over the stub.
old_dir = cd(stub_dir);
cleaner = onCleanup(@() cleanup(stub_dir, old_dir)); %#ok<NASGU>
rehash;

FS = 500;
SOLVER = 'sra1';
PLUMBED = {'eig_heatmap', 'lyapunov_spectrum', 'transient_gain', 'local_lyapunov', 'dc_lle', ...
           'memory_capacity', 'mc_example'};
OWN_RATE = {'numerics'};
OPT_IN_PLUMBED = {'lyapunov_spectrum_n80', 'outlier_spectra', 'outlier_structures'};   % plans read cfg.fs, cfg.ode_solver
OPT_IN_PLAN2 = {'concept_examples', 'uncoupled_reference'};   % Plan 2: plans read cfg.fs, cfg.ode_solver (the second the joint sample's way)
OPT_IN_NO_SIM = {'frequency_response', 'eigenvalue_density'};                         % read outlier_spectra's output only
SWEEPS = {'sensitivity', 'tau_sensitivity', 'param_space'};

%% PART 1: the dry runs
for bundle = {'mu7chain_config', 'mu7chain_AC_config', 'mu7chain_gainbias_config'}
    for mode = {'fast', 'production'}
        cfg = feval(bundle{1}, mode{1}, '3cond');
        label = sprintf('%s(''%s'')', bundle{1}, mode{1});
        assert(isequal(cfg.fs, FS), '%s: cfg.fs is %g, not %d', label, cfg.fs, FS);
        assert(strcmp(cfg.ode_solver, SOLVER), '%s: cfg.ode_solver is %s, not %s', label, cfg.ode_solver, SOLVER);
        [L, run_dir] = dry_run(cfg, stub_dir);
        % Every selected stage must have reached its stub: a stage that
        % errored is caught by run_all_paper_analyses and would otherwise
        % pass here by being absent.
        plan = paper_stage_plan(cfg);
        for s = intersect([PLUMBED, OWN_RATE, OPT_IN_PLUMBED, OPT_IN_PLAN2, OPT_IN_NO_SIM, ...
                {'paper_illustrations', 'mc_gain_sweep', 'noise_switch'}], plan.stages)
            assert(any(strcmp({L.stage}, s{1})), '%s: selected stage %s never ran (see FAILED in the log)', label, s{1});
        end
        if ismember('sweeps', plan.stages)
            assert(any(ismember({L.stage}, SWEEPS)), '%s: no sweep ran', label);
        end
        if ~strcmp(bundle{1}, 'mu7chain_gainbias_config')
            assert(all(ismember([OPT_IN_PLUMBED, OPT_IN_NO_SIM], plan.stages)), ...
                '%s: the spectral stages and the 80-neuron spectrum are this bundle''s default', label);
        end
        for k = find(ismember({L.stage}, [SWEEPS, PLUMBED, OPT_IN_PLUMBED, OPT_IN_PLAN2, {'paper_illustrations', 'noise_switch'}]))
            assert(isequal(L(k).fs, FS) && strcmp(L(k).solver, SOLVER), ...
                '%s: %s was handed fs %s and ode_solver ''%s'', not %d and %s', ...
                label, L(k).stage, mat2str(L(k).fs), L(k).solver, FS, SOLVER);
        end
        for k = find(strcmp({L.stage}, 'mc_gain_sweep_run'))
            assert(isequal(L(k).fs, FS) && strcmp(L(k).solver, SOLVER), ...
                '%s: an mc_gain_sweep run is planned at protocol.fs %s, ode_solver ''%s''', ...
                label, mat2str(L(k).fs), L(k).solver);
        end
        if ismember('mc_gain_sweep', plan.stages)
            assert(any(strcmp({L.stage}, 'mc_gain_sweep_run')), '%s: mc_gain_sweep planned no run', label);
        end
        k = find(strcmp({L.stage}, 'noise_switch'));
        if ~isempty(k)
            assert(mod(FS, L(k).store_fs) == 0 && mod(L(k).store_fs, L(k).example_fs) == 0, ...
                '%s: noise_switch store_fs %g / example_fs %g do not divide %d', ...
                label, L(k).store_fs, L(k).example_fs, FS);
        end
        for s = OWN_RATE
            k = find(strcmp({L.stage}, s{1}));
            if isempty(k); continue; end
            assert(isempty(L(k).fs) && isempty(L(k).solver), ...
                '%s: stage %s was handed fs or ode_solver, but keeps its own', label, s{1});
        end
        M = load(fullfile(run_dir, 'run_manifest.mat'));
        assert(isfield(M.run_manifest, 'fs') && M.run_manifest.fs == FS, '%s: run_manifest lacks fs = %d', label, FS);
        assert(isfield(M.run_manifest, 'ode_solver') && strcmp(M.run_manifest.ode_solver, SOLVER), ...
            '%s: run_manifest lacks ode_solver = %s', label, SOLVER);
        rmdir(run_dir, 's');
        % The registry entries that simulate while drawing.
        for f = cfg.figures
            if ismember(f{1}.name, {'fig_introductory_concepts', 'fig_stim_engages_adaptation'})
                a = f{1}.args;
                i_fs = find(strcmp(a(1:2:end), 'fs'), 1, 'last');
                i_s  = find(strcmp(a(1:2:end), 'ode_solver'), 1, 'last');
                assert(~isempty(i_fs) && isequal(a{2*i_fs}, FS) && ~isempty(i_s) && strcmp(a{2*i_s}, SOLVER), ...
                    '%s: figure entry %s does not carry fs %d and ode_solver %s', label, f{1}.name, FS, SOLVER);
            end
        end
        ran = unique({L.stage});
        fprintf('  %-44s fs %d, %s everywhere taken; called: %s\n', label, FS, SOLVER, strjoin(ran, ', '));
    end
end

% Without cfg.fs and cfg.ode_solver: nothing is handed on, the sweeps keep
% the mode table's rate and solver. Memory capacity still takes the bundle's
% cfg.mc_ode_solver, which is its own field.
cfg = rmfield(mu7chain_config('fast', '3cond'), {'fs', 'ode_solver'});
[L, run_dir] = dry_run(cfg, stub_dir);
assert(any(strcmp({L.stage}, 'eig_heatmap')) && any(ismember({L.stage}, SWEEPS)), 'no cfg.fs: stages missing');
for s = PLUMBED
    k = find(strcmp({L.stage}, s{1}));
    if isempty(k); continue; end
    assert(isempty(L(k).fs), 'no cfg.fs: stage %s was still handed fs', s{1});
    if ~ismember(s{1}, {'memory_capacity', 'mc_example'})
        assert(isempty(L(k).solver), 'no cfg.ode_solver: stage %s was still handed ode_solver', s{1});
    end
end
% The opt-in stages keep their own 500 Hz and name no integrator (the
% model's builder then picks sra1 for their noisy presets).
for s = OPT_IN_PLUMBED
    k = find(strcmp({L.stage}, s{1}));
    assert(isscalar(k) && isequal(L(k).fs, 500) && isempty(L(k).solver), ...
        'no cfg.fs: opt-in stage %s is planned at %s with ''%s'', not its own 500 Hz and no named integrator', ...
        s{1}, mat2str(L(k).fs), L(k).solver);
end
for s = SWEEPS
    want = analysis_run_config(s{1}, 'fast', srnn_param_preset(cfg.preset_name));
    for k = find(strcmp({L.stage}, s{1}))
        assert(isequal(L(k).fs, want.model.fs) && strcmp(L(k).solver, want.model.ode_solver), ...
            'no cfg.fs: sweep %s ran at %g with %s, not the table''s %g with %s', ...
            s{1}, L(k).fs, L(k).solver, want.model.fs, want.model.ode_solver);
    end
end
M = load(fullfile(run_dir, 'run_manifest.mat'));
assert(~isfield(M.run_manifest, 'fs') && ~isfield(M.run_manifest, 'ode_solver'), ...
    'no cfg.fs: run_manifest has an fs or ode_solver field');
rmdir(run_dir, 's');
fprintf('  %-44s nothing handed on; sweeps at the mode table''s rate and solver\n', ...
    'mu7chain_config(''fast''), both removed');

% The opt-in plans follow the bundle's rate and integrator, whatever they
% are, and a stage's own fs knob overrides the bundle's (plans only).
c = mu7chain_config('fast', '3cond');
c.fs = 400; c.ode_solver = 'heun';
p80 = lyapunov_spectrum_n80_stage_plan(c); po = outlier_spectra_stage_plan(c); ps = outlier_structures_stage_plan(c);
[f80, s80] = fs_solver_of(p80.args);
assert(isequal([f80, p80.fs, po.fs, ps.fs], [400 400 400 400]) && ...
    isequal({s80, p80.ode_solver, po.ode_solver, ps.ode_solver}, {'heun', 'heun', 'heun', 'heun'}), ...
    'the opt-in plans do not follow cfg.fs = 400 and cfg.ode_solver = ''heun''');
c.outlier_spectra_fs = 250; c.outlier_structures_fs = 1000;
po = outlier_spectra_stage_plan(c); ps = outlier_structures_stage_plan(c);
assert(po.fs == 250 && ps.fs == 1000, 'a stage''s own fs knob must override the bundle''s cfg.fs');
% ospec_run integrates in pieces: an adaptive integrator is refused in the plan.
c = mu7chain_config('fast', '3cond'); c.ode_solver = 'ode45';
for pf = {@outlier_spectra_stage_plan, @outlier_structures_stage_plan}
    try
        pf{1}(c);
        error('test_bundle_fs_plumbing:NoError', '%s accepted ode45', func2str(pf{1}));
    catch ME
        assert(endsWith(ME.identifier, ':BadConfig'), '%s: unexpected error %s', func2str(pf{1}), ME.identifier);
    end
end
fprintf('  %-44s follow cfg.fs and cfg.ode_solver; own fs knob overrides; ode45 refused\n', 'the opt-in plans');

%% PART 2: the stages use what they are handed
% From here the stage stubs hand over to the real stages (stub_active).
FS_PLUMB_STAGES_REAL = true;
cfg = mu7chain_config('fast', '3cond');
tmp = tempname; mkdir(tmp);
deep = { ...
    'run_eig_heatmap',   @() run_eig_heatmap('preset_name', cfg.preset_name, 'run_mode', 'fast', ...
        'verbose', 'near-none', 'out_dir', tmp, 'fs', FS, 'ode_solver', SOLVER); ...
    'run_local_lyapunov', @() run_local_lyapunov('preset_name', cfg.local_lyapunov_preset, 'run_mode', 'fast', ...
        'verbose', 'near-none', 'out_dir', tmp, 'fs', FS, 'ode_solver', SOLVER, 'T', 40, ...
        'step_duration_s', 5, 'accumulation_start_s', 1, 'warmup_s', 15, 'simulation_start_s', -15); ...
    'run_paper_illustrations', @() run_paper_illustrations(setfield(cfg, 'verbose', 'near-none'), tmp); ... %#ok<SFLD>
    'fig_introductory_concepts', @() fig_introductory_concepts('preset_name', cfg.sompolinsky_preset, ...
        'save', false, 'visible', false, 'fs', FS, 'ode_solver', SOLVER); ...
    'fig_stim_engages_adaptation', @() fig_stim_engages_adaptation('preset_name', cfg.bursting_preset, ...
        'save', false, 'visible', false, 'fs', FS, 'ode_solver', SOLVER)};
for k = 1:size(deep, 1)
    FS_PLUMB_LOG = struct('stage', {}, 'fs', {}, 'solver', {}, 'store_fs', {}, 'example_fs', {});
    assert(strcmp(fileparts(which('build_from_preset')), stub_dir), 'build_from_preset is not stubbed; refusing');
    try
        evalc('deep{k, 2}();');
        error('test_bundle_fs_plumbing:NoBuild', '%s returned without building a model', deep{k, 1});
    catch ME
        if ~strcmp(ME.identifier, 'STUB:Stop'); rethrow(ME); end
    end
    L = FS_PLUMB_LOG;
    assert(~isempty(L), '%s: no model was built', deep{k, 1});
    assert(isequal(L(1).fs, FS) && strcmp(L(1).solver, SOLVER), ...
        '%s built its model with fs %s and ode_solver ''%s'', not %d and %s', ...
        deep{k, 1}, mat2str(L(1).fs), L(1).solver, FS, SOLVER);
    fprintf('  %-44s builds its model at fs %d with %s\n', deep{k, 1}, FS, SOLVER);
end
rmdir(tmp, 's');

%% stage_fs and stage_solver
assert(stage_fs(0, 400, 't') == 400 && stage_fs([], 200, 't') == 200 && stage_fs(500, 400, 't') == 500);
for bad = {-1, NaN, Inf, [400 500], 'x'}
    try
        stage_fs(bad{1}, 400, 't');
        error('test_bundle_fs_plumbing:NoError', 'stage_fs accepted %s', mat2str(bad{1}));
    catch ME
        assert(strcmp(ME.identifier, 't:BadFs'), 'stage_fs: unexpected error %s', ME.identifier);
    end
end
assert(strcmp(stage_solver('', 'rk4', 't'), 'rk4') && isempty(stage_solver('', '', 't')) && ...
    strcmp(stage_solver('sra1', 'rk4', 't'), 'sra1') && strcmp(stage_solver("heun", '', 't'), 'heun'));
for bad = {'nosuch', 3}
    try
        stage_solver(bad{1}, 'rk4', 't');
        error('test_bundle_fs_plumbing:NoError', 'stage_solver accepted %s', mat2str(bad{1}));
    catch ME
        assert(~strcmp(ME.identifier, 'test_bundle_fs_plumbing:NoError'), 'stage_solver accepted a bad name');
    end
end
fprintf('test_bundle_fs_plumbing: passed\n');
end

%% ------------------------------------------------------------------------
function [L, run_dir] = dry_run(cfg, stub_dir)
global FS_PLUMB_LOG %#ok<GVMIS>
FS_PLUMB_LOG = struct('stage', {}, 'fs', {}, 'solver', {}, 'store_fs', {}, 'example_fs', {});
cfg.run_dir = tempname;                 % absolute, absent: a fresh run directory
cfg.verbose = 'near-none';
run_dir = cfg.run_dir;
% SAFETY: never start the real stages. After the entry points' own
% setup_paths call, every simulating function must still resolve to a stub.
setup_paths();
for name = {'run_sensitivity_analysis', 'run_tau_sensitivity_analysis', 'run_param_space_analysis', ...
        'run_eig_heatmap', 'run_lyapunov_spectrum', 'run_transient_gain', 'run_local_lyapunov', ...
        'run_dc_lle_analysis', 'run_numerics_verification', 'run_memory_capacity', ...
        'run_memory_capacity_example', 'run_paper_illustrations', 'run_mc_gain_sweep_stage', ...
        'run_noise_switch_stage', 'run_lyapunov_spectrum_n80_stage', 'run_outlier_spectra_stage', ...
        'run_frequency_response_stage', 'run_eigenvalue_density_stage', 'run_outlier_structures_stage', ...
        'run_concept_examples_stage', 'run_uncoupled_reference_stage', ...
        'pool_restart', 'build_from_preset'}
    w = which(name{1});
    assert(strcmp(fileparts(w), stub_dir), ...
        'test_bundle_fs_plumbing:NotStubbed', '%s resolves to %s, not to its stub; refusing to run', name{1}, w);
end
evalc('run_all_paper_analyses(cfg);');
L = FS_PLUMB_LOG;
assert(~isempty(L), 'dry run of %s recorded no stage', cfg.run_mode);
end

function [f, s] = fs_solver_of(a)
% The 'fs' and 'ode_solver' of a name-value list ([] and '' = none).
f = []; s = '';
i = find(strcmp(a(1:2:end), 'fs'), 1, 'last'); if ~isempty(i); f = a{2*i}; end
i = find(strcmp(a(1:2:end), 'ode_solver'), 1, 'last'); if ~isempty(i); s = char(a{2*i}); end
end

function cleanup(stub_dir, old_dir)
global FS_PLUMB_LOG %#ok<GVMIS>
cd(old_dir);
rmpath(stub_dir);
rmdir(stub_dir, 's');
clear global FS_PLUMB_LOG FS_PLUMB_STAGES_REAL
rehash;
end

function write_stubs(d)
% Name-value stages: record the 'fs' and 'ode_solver' they were given
% ([] and '' = none). In part 2 (FS_PLUMB_STAGES_REAL true) each hands
% its call to the real stage instead (stub_call_real).
nv = {'run_eig_heatmap', 'eig_heatmap'; 'run_lyapunov_spectrum', 'lyapunov_spectrum'; ...
      'run_transient_gain', 'transient_gain'; 'run_local_lyapunov', 'local_lyapunov'; ...
      'run_dc_lle_analysis', 'dc_lle'; 'run_numerics_verification', 'numerics'; ...
      'run_memory_capacity', 'memory_capacity'; 'run_memory_capacity_example', 'mc_example'};
for k = 1:size(nv, 1)
    put(d, nv{k, 1}, { ...
        sprintf('function out = %s(varargin)', nv{k, 1}), ...
        'global FS_PLUMB_LOG', ...
        'out = ''stub'';', ...
        'if ~stub_active(); out = stub_call_real(mfilename, varargin{:}); return; end', ...
        '[f, s] = stub_fs_solver(varargin);', ...
        sprintf('stub_log(''%s'', f, s);', nv{k, 2}), ...
        'end'});
end
% A stage stub hands over to the real function once part 1 is done
% (stub_active false): the real file is found by name in the project.
root = fileparts(which('setup_paths'));
put(d, 'stub_active', {'function tf = stub_active()', 'global FS_PLUMB_STAGES_REAL %#ok<GVMIS>', ...
    'tf = isempty(FS_PLUMB_STAGES_REAL) || ~FS_PLUMB_STAGES_REAL;', 'end'});
put(d, 'stub_call_real', {'function out = stub_call_real(name, varargin)', ...
    sprintf('hits = dir(fullfile(''%s'', ''src'', ''**'', [name ''.m'']));', root), ...
    'assert(isscalar(hits), ''stub_call_real: %s not found once'', name);', ...
    'here = cd(hits.folder); back = onCleanup(@() cd(here));', ...
    'fh = str2func(name); out = fh(varargin{:});', 'end'});
put(d, 'stub_fs_solver', {'function [f, s] = stub_fs_solver(a)', ...
    'f = []; s = '''';', ...
    'i = find(strcmp(a(1:2:end), ''fs''), 1, ''last''); if ~isempty(i); f = a{2*i}; end', ...
    'i = find(strcmp(a(1:2:end), ''ode_solver''), 1, ''last''); if ~isempty(i); s = char(a{2*i}); end', 'end'});
put(d, 'stub_log', {'function stub_log(stage, f, s, store_fs, example_fs)', 'global FS_PLUMB_LOG %#ok<GVMIS>', ...
    'if nargin < 4; store_fs = []; end', 'if nargin < 5; example_fs = []; end', ...
    'FS_PLUMB_LOG(end+1) = struct(''stage'', stage, ''fs'', f, ''solver'', s, ''store_fs'', store_fs, ''example_fs'', example_fs);', 'end'});
% The sweeps: record the rate and solver the resolved context would build models with.
sw = {'run_sensitivity_analysis', 'sensitivity'; 'run_tau_sensitivity_analysis', 'tau_sensitivity'; ...
      'run_param_space_analysis', 'param_space'};
for k = 1:size(sw, 1)
    put(d, sw{k, 1}, { ...
        sprintf('function %s(ctx)', sw{k, 1}), ...
        'assert(ctx.model_defaults.fs == ctx.cfg.model.fs);', ...
        'assert(strcmp(ctx.model_defaults.ode_solver, ctx.cfg.model.ode_solver));', ...
        sprintf('stub_log(''%s'', ctx.model_defaults.fs, ctx.model_defaults.ode_solver);', sw{k, 2}), ...
        'end'});
end
put(d, 'run_paper_illustrations', { ...
    'function out = run_paper_illustrations(cfg, run_dir)', ...
    'if ~stub_active(); out = stub_call_real(mfilename, cfg, run_dir); return; end', ...
    'f = []; if isfield(cfg, ''fs''); f = cfg.fs; end', ...
    's = ''''; if isfield(cfg, ''ode_solver''); s = cfg.ode_solver; end', ...
    'stub_log(''paper_illustrations'', f, s);', 'out = run_dir;', 'end'});
% Block C: one record per planned run, its solver and its protocol's fs.
put(d, 'run_mc_gain_sweep_stage', {'function out = run_mc_gain_sweep_stage(cfg, run_dir, plan)', ...
    'stub_log(''mc_gain_sweep'', [], '''');', ...
    'for r = plan.runs', ...
    '    [~, s] = stub_fs_solver(r.args);', ...
    '    i = find(strcmp(r.args(1:2:end), ''protocol''), 1, ''last''); f = [];', ...
    '    if ~isempty(i) && isfield(r.args{2*i}, ''fs''); f = r.args{2*i}.fs; end', ...
    '    stub_log(''mc_gain_sweep_run'', f, s);', 'end', 'out = ''stub'';', 'end'});
put(d, 'run_noise_switch_stage', {'function out = run_noise_switch_stage(cfg, run_dir, plan)', ...
    '[f, s] = stub_fs_solver(plan.args); a = plan.args;', ...
    'i = find(strcmp(a(1:2:end), ''store_fs''), 1, ''last''); st = 200; if ~isempty(i); st = a{2*i}; end', ...
    'i = find(strcmp(a(1:2:end), ''example_fs''), 1, ''last''); ex = 100; if ~isempty(i); ex = a{2*i}; end', ...
    'stub_log(''noise_switch'', f, s, st, ex);', 'out = ''stub'';', 'end'});
% The opt-in stages built after the plumbing (FD-032, FD-040 to FD-043): the
% rate and integrator of their plans. Never run for real here: at the
% bundles' sizes they simulate 500 neurons.
put(d, 'run_lyapunov_spectrum_n80_stage', {'function out = run_lyapunov_spectrum_n80_stage(cfg, run_dir, plan)', ...
    '[f, s] = stub_fs_solver(plan.args);', ...
    'assert(isequal(f, plan.fs) && strcmp(s, plan.ode_solver), ''n80 plan: args and fields disagree'');', ...
    'stub_log(''lyapunov_spectrum_n80'', f, s);', 'out = ''stub'';', 'end'});
for st = {'outlier_spectra', 'outlier_structures'}
    put(d, ['run_' st{1} '_stage'], {sprintf('function out = run_%s_stage(cfg, run_dir, plan)', st{1}), ...
        sprintf('stub_log(''%s'', plan.fs, plan.ode_solver);', st{1}), 'out = ''stub'';', 'end'});
end
% Plan 2 (2026-11-21): concept_examples (its plan's fs and ode_solver) and
% uncoupled_reference (the joint sample's resolved model settings). Without
% these stubs the dry runs ran both for real at 500 neurons, six times (29 min).
put(d, 'run_concept_examples_stage', {'function out = run_concept_examples_stage(cfg, run_dir, plan)', ...
    'stub_log(''concept_examples'', plan.fs, plan.ode_solver);', 'out = ''stub'';', 'end'});
put(d, 'run_uncoupled_reference_stage', {'function out = run_uncoupled_reference_stage(cfg, run_dir, plan)', ...
    'stub_log(''uncoupled_reference'', plan.ctx.model_defaults.fs, plan.ctx.model_defaults.ode_solver);', 'out = ''stub'';', 'end'});
for st = {'frequency_response', 'eigenvalue_density'}
    put(d, ['run_' st{1} '_stage'], {sprintf('function out = run_%s_stage(cfg, run_dir, plan)', st{1}), ...
        sprintf('stub_log(''%s'', [], '''');', st{1}), 'out = ''stub'';', 'end'});
end
% Part 2: the model builder records what it was handed and stops.
put(d, 'build_from_preset', {'function [model, info] = build_from_preset(preset_name, condition_name, varargin)', ...
    'global FS_PLUMB_STAGES_REAL %#ok<GVMIS>', ...
    'assert(~isempty(FS_PLUMB_STAGES_REAL) && FS_PLUMB_STAGES_REAL, ''build_from_preset stub reached in part 1'');', ...
    '[f, s] = stub_fs_solver(varargin);', 'stub_log(''build'', f, s);', ...
    'model = []; info = [];', 'error(''STUB:Stop'', ''stub build: stop'');', 'end'});
put(d, 'replot_sensitivity', {'function out = replot_sensitivity(run_dir)', 'out = run_dir;', 'end'});
put(d, 'assemble_sensitivity_figure', {'function assemble_sensitivity_figure(varargin)', 'end'});
put(d, 'write_run_parameters_md', {'function write_run_parameters_md(varargin)', 'end'});
% setup_paths: the entry points call it, and the real one puts the project
% folders back in FRONT of the stubs (addpath prepends), which would run the
% real stages. This one adds the same folders, then the stubs in front again.
put(d, 'setup_paths', { ...
    'function projectRoot = setup_paths()', ...
    sprintf('projectRoot = ''%s'';', root), ...
    'addpath(projectRoot); addpath(genpath(fullfile(projectRoot, ''src'')));', ...
    'addpath(genpath(fullfile(projectRoot, ''scripts'')));', ...
    sprintf('addpath(''%s'', ''-begin'');', d), 'end'});
put(d, 'capture_git_provenance', {'function info = capture_git_provenance(varargin)', ...
    'info = struct(''commit'', '''', ''commit_short'', '''', ''branch'', '''', ''is_dirty'', false, ''hostname'', '''', ''platform'', computer());', 'end'});
put(d, 'pool_restart', {'function pool_restart(varargin)', 'end'});
put(d, 'pool_config', {'function pc = pool_config(varargin)', 'pc = struct(''is_set'', true, ''how_to_set'', '''');', 'end'});
end

function put(d, name, lines)
fid = fopen(fullfile(d, [name '.m']), 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
