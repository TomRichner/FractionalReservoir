function test_bundle_fs_plumbing()
% TEST_BUNDLE_FS_PLUMBING The bundle config's cfg.fs reaches every stage that takes it (FD-031 item 14).
%
% Tom (2026-11-17): 500 Hz for runs from now on, stated once in the bundle
% config (mu7chain_config's cfg.fs, inherited by mu7chain_AC_config). This is
% a DRY RUN of run_all_paper_analyses: the real entry point, the real
% run_all_analyses and resolve_run_context, with every stage that simulates
% replaced by a stub (a folder put at the front of the path) that records
% what it was called with. Nothing is simulated.
%
% Checked, at 'fast' and 'production', for mu7chain_config and
% mu7chain_AC_config:
%   * the three sweeps' resolved model_defaults.fs is cfg.fs (500), i.e. the
%     analysis_run_config mode table's rate is replaced;
%   * eig_heatmap, lyapunov_spectrum, transient_gain, local_lyapunov and
%     dc_lle are called with 'fs', 500;
%   * paper_illustrations gets the bundle cfg carrying fs = 500 (it reads
%     cfg.fs itself);
%   * the stages that keep their own rate by design (numerics,
%     memory_capacity, mc_example; the opt-in stages read cfg) are NOT
%     called with 'fs';
%   * run_manifest records fs = 500.
% And with cfg.fs removed, no stage gets 'fs' and the sweeps keep the mode
% table's rate (each stage's own default), as before the field existed.
% Finally stage_fs itself: 0 and [] give the default, a bad value errors.
%
% Quick: a few seconds.

global FS_PLUMB_LOG %#ok<GVMIS>   the stubs append here

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

PLUMBED = {'eig_heatmap', 'lyapunov_spectrum', 'transient_gain', 'local_lyapunov', 'dc_lle'};
OWN_RATE = {'numerics', 'memory_capacity', 'mc_example'};
SWEEPS = {'sensitivity', 'tau_sensitivity', 'param_space'};

for bundle = {'mu7chain_config', 'mu7chain_AC_config'}
    for mode = {'fast', 'production'}
        cfg = feval(bundle{1}, mode{1}, '3cond');
        assert(isequal(cfg.fs, 500), '%s(%s): cfg.fs is %g, not 500', bundle{1}, mode{1}, cfg.fs);
        [L, run_dir] = dry_run(cfg, stub_dir);
        label = sprintf('%s(''%s'')', bundle{1}, mode{1});
        % Every selected stage must have reached its stub: a stage that
        % errored is caught by run_all_paper_analyses and would otherwise
        % pass here by being absent.
        plan = paper_stage_plan(cfg);
        for s = intersect([PLUMBED, OWN_RATE, {'paper_illustrations'}], plan.stages)
            assert(any(strcmp({L.stage}, s{1})), '%s: selected stage %s never ran (see FAILED in the log)', label, s{1});
        end
        if ismember('sweeps', plan.stages)
            assert(any(ismember({L.stage}, SWEEPS)), '%s: no sweep ran', label);
        end
        for s = SWEEPS
            for k = find(strcmp({L.stage}, s{1}))
                assert(isequal(L(k).fs, 500), '%s: sweep %s ran at fs %g, not 500', label, s{1}, L(k).fs);
            end
        end
        for s = PLUMBED
            k = find(strcmp({L.stage}, s{1}));
            if isempty(k); continue; end          % a stage this bundle does not select
            assert(isequal(L(k).fs, 500), '%s: stage %s was handed fs %s, not 500', label, s{1}, mat2str(L(k).fs));
        end
        k = find(strcmp({L.stage}, 'paper_illustrations'));
        if ~isempty(k)
            assert(isequal(L(k).fs, 500), '%s: paper_illustrations cfg.fs is %s', label, mat2str(L(k).fs));
        end
        for s = OWN_RATE
            k = find(strcmp({L.stage}, s{1}));
            if isempty(k); continue; end
            assert(isempty(L(k).fs), '%s: stage %s was handed fs, but keeps its own rate', label, s{1});
        end
        M = load(fullfile(run_dir, 'run_manifest.mat'));
        assert(isfield(M.run_manifest, 'fs') && M.run_manifest.fs == 500, '%s: run_manifest lacks fs = 500', label);
        rmdir(run_dir, 's');
        ran = unique({L.stage});
        fprintf('  %-38s fs = 500 where taken; stages called: %s\n', label, strjoin(ran, ', '));
    end
end

% Without cfg.fs: nothing is handed on, the sweeps keep the mode table's rate.
cfg = rmfield(mu7chain_config('fast', '3cond'), 'fs');
[L, run_dir] = dry_run(cfg, stub_dir);
assert(any(strcmp({L.stage}, 'eig_heatmap')) && any(ismember({L.stage}, SWEEPS)), 'no cfg.fs: stages missing');
for s = PLUMBED
    k = find(strcmp({L.stage}, s{1}));
    if ~isempty(k); assert(isempty(L(k).fs), 'no cfg.fs: stage %s was still handed fs', s{1}); end
end
for s = SWEEPS
    want = analysis_run_config(s{1}, 'fast', srnn_param_preset(cfg.preset_name));
    for k = find(strcmp({L.stage}, s{1}))
        assert(isequal(L(k).fs, want.model.fs), 'no cfg.fs: sweep %s ran at %g, not the table''s %g', ...
            s{1}, L(k).fs, want.model.fs);
    end
end
M = load(fullfile(run_dir, 'run_manifest.mat'));
assert(~isfield(M.run_manifest, 'fs'), 'no cfg.fs: run_manifest has an fs field');
rmdir(run_dir, 's');
fprintf('  %-38s no fs handed on; sweeps at the mode table''s rate\n', 'mu7chain_config(''fast''), cfg.fs removed');

% stage_fs
assert(stage_fs(0, 400, 't') == 400 && stage_fs([], 200, 't') == 200 && stage_fs(500, 400, 't') == 500);
for bad = {-1, NaN, Inf, [400 500], 'x'}
    try
        stage_fs(bad{1}, 400, 't');
        error('test_bundle_fs_plumbing:NoError', 'stage_fs accepted %s', mat2str(bad{1}));
    catch ME
        assert(strcmp(ME.identifier, 't:BadFs'), 'stage_fs: unexpected error %s', ME.identifier);
    end
end
fprintf('test_bundle_fs_plumbing: passed\n');
end

%% ------------------------------------------------------------------------
function [L, run_dir] = dry_run(cfg, stub_dir)
global FS_PLUMB_LOG %#ok<GVMIS>
FS_PLUMB_LOG = struct('stage', {}, 'fs', {});
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
        'run_noise_switch_stage', 'pool_restart'}
    w = which(name{1});
    assert(strcmp(fileparts(w), stub_dir), ...
        'test_bundle_fs_plumbing:NotStubbed', '%s resolves to %s, not to its stub; refusing to run', name{1}, w);
end
evalc('run_all_paper_analyses(cfg);');
L = FS_PLUMB_LOG;
assert(~isempty(L), 'dry run of %s recorded no stage', cfg.run_mode);
end

function cleanup(stub_dir, old_dir)
global FS_PLUMB_LOG %#ok<GVMIS>
cd(old_dir);
rmpath(stub_dir);
rmdir(stub_dir, 's');
clear global FS_PLUMB_LOG
rehash;
end

function write_stubs(d)
% Name-value stages: record the 'fs' they were given ([] = none).
nv = {'run_eig_heatmap', 'eig_heatmap'; 'run_lyapunov_spectrum', 'lyapunov_spectrum'; ...
      'run_transient_gain', 'transient_gain'; 'run_local_lyapunov', 'local_lyapunov'; ...
      'run_dc_lle_analysis', 'dc_lle'; 'run_numerics_verification', 'numerics'; ...
      'run_memory_capacity', 'memory_capacity'; 'run_memory_capacity_example', 'mc_example'};
for k = 1:size(nv, 1)
    put(d, nv{k, 1}, { ...
        sprintf('function out = %s(varargin)', nv{k, 1}), ...
        'global FS_PLUMB_LOG', ...
        'i = find(strcmp(varargin(1:2:end), ''fs''), 1); f = [];', ...
        'if ~isempty(i); f = varargin{2*i}; end', ...
        sprintf('FS_PLUMB_LOG(end+1) = struct(''stage'', ''%s'', ''fs'', f);', nv{k, 2}), ...
        'out = ''stub'';', 'end'});
end
% The sweeps: record the rate the resolved context would build models at.
sw = {'run_sensitivity_analysis', 'sensitivity'; 'run_tau_sensitivity_analysis', 'tau_sensitivity'; ...
      'run_param_space_analysis', 'param_space'};
for k = 1:size(sw, 1)
    put(d, sw{k, 1}, { ...
        sprintf('function %s(ctx)', sw{k, 1}), ...
        'global FS_PLUMB_LOG', ...
        'assert(ctx.model_defaults.fs == ctx.cfg.model.fs);', ...
        sprintf('FS_PLUMB_LOG(end+1) = struct(''stage'', ''%s'', ''fs'', ctx.model_defaults.fs);', sw{k, 2}), ...
        'end'});
end
put(d, 'run_paper_illustrations', { ...
    'function out = run_paper_illustrations(cfg, run_dir)', 'global FS_PLUMB_LOG', ...
    'f = []; if isfield(cfg, ''fs''); f = cfg.fs; end', ...
    'FS_PLUMB_LOG(end+1) = struct(''stage'', ''paper_illustrations'', ''fs'', f);', 'out = run_dir;', 'end'});
put(d, 'run_mc_gain_sweep_stage', {'function out = run_mc_gain_sweep_stage(varargin)', 'out = ''stub'';', 'end'});
put(d, 'run_noise_switch_stage', {'function out = run_noise_switch_stage(varargin)', 'out = ''stub'';', 'end'});
put(d, 'replot_sensitivity', {'function out = replot_sensitivity(run_dir)', 'out = run_dir;', 'end'});
put(d, 'assemble_sensitivity_figure', {'function assemble_sensitivity_figure(varargin)', 'end'});
put(d, 'write_run_parameters_md', {'function write_run_parameters_md(varargin)', 'end'});
% setup_paths: the entry points call it, and the real one puts the project
% folders back in FRONT of the stubs (addpath prepends), which would run the
% real stages. This one adds the same folders, then the stubs in front again.
root = fileparts(which('setup_paths'));
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
