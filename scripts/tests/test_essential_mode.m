% TEST_ESSENTIAL_MODE The 'essential' run mode and the stage selection (FD-023).
%
% Quick and compute-free: no stage, no sweep and no figure runs. Checks that
%   * run_mode_names holds 'essential', first;
%   * analysis_run_config('...', 'essential') is 4 levels x 2 reps on the 1-D
%     sweeps and 32 joint points, with fast's solver, fs, window and
%     Lyapunov settings in all three analyses, noise on and off;
%   * every OTHER mode's analysis_run_config output is what it was before the
%     mode was added (the literals below were read off the function at
%     0bc35ca, 2026-11-12);
%   * every stage's cost table names 'essential' beside 'fast' (from the
%     stage files' text: the tables are local and reaching one means running
%     the stage; see test_run_modes for the same limit);
%   * paper_stage_plan returns every stage and sweep for an absent or empty
%     selection, the selected ones in the entry point's order for a
%     selection, and errors on an unknown stage or sweep name;
%   * run_all_paper_analyses' stage table is paper_stage_names' list (from
%     its text), and its pre-flight refuses an unknown stage or sweep before
%     the run directory exists; run_all_analyses refuses an unknown sweep,
%     and with an empty sweep list makes the directory and a run_manifest
%     that records the stages and sweeps skipped, simulating nothing;
%   * mu7chain_config('essential') has the run directory, figure root,
%     stages, sweeps, eig-heatmap example, figure-set default and registry
%     specified, and every registry entry is the fast registry's own;
%   * mu7chain_config('fast') and ('medium') select nothing.
%
% Asserts; silent until it throws, then prints one PASS line.
% Assumes setup_paths has run.
% See also: run_mode_names, analysis_run_config, paper_stage_plan,
%           paper_stage_names, mu7chain_config, test_run_modes

%% The mode exists, first in the list
modes = run_mode_names();
assert(isequal(modes, {'essential', 'fast', 'medium', 'medium2', 'production'}), ...
    'run_mode_names: %s', strjoin(modes, ', '));

%% analysis_run_config: essential's sizes, fast's everything else
[noisy, ~, ~] = srnn_param_preset(mu7chain_config('essential').preset_name);
assert(noisy.sigma_u_noise > 0);
for p = {struct(), noisy}
    es = analysis_run_config('sensitivity', 'essential', p{1});
    fs = analysis_run_config('sensitivity', 'fast', p{1});
    assert(es.n_levels == 4 && es.n_reps == 2, 'sensitivity: %d levels x %d reps', es.n_levels, es.n_reps);
    assert(isequal(es.model, fs.model), 'essential sensitivity must run fast''s solver, fs and windows');
    assert(same_but(es, fs, {'n_levels', 'n_reps', 'run_mode'}));

    ep = analysis_run_config('param_space', 'essential', p{1});
    fp = analysis_run_config('param_space', 'fast', p{1});
    assert(ep.n_grid_points == 32 && fp.n_grid_points == 128);
    assert(ep.n_levels == fp.n_levels && ~isfield(ep, 'n_reps'));
    assert(isequal(ep.model, fp.model), 'essential param_space must run fast''s solver, fs and windows');
    assert(same_but(ep, fp, {'n_grid_points', 'run_mode'}));

    et = analysis_run_config('tau_sensitivity', 'essential', p{1});
    ft = analysis_run_config('tau_sensitivity', 'fast', p{1});
    assert(same_but(et, ft, {'run_mode'}), 'essential tau_sensitivity must be fast''s row');
end
% Stated outright as well, so a change to fast's row cannot carry essential's with it unseen.
es = analysis_run_config('sensitivity', 'essential');
assert(strcmp(es.model.ode_solver, 'rk4') && es.model.fs == 200 && ...
    isequal(es.model.T_range, [0 10]) && isequal(es.model.lya_T_interval, [5 10]));
ep = analysis_run_config('param_space', 'essential');
assert(strcmp(ep.model.ode_solver, 'rk4') && ep.model.fs == 200 && ...
    isequal(ep.model.T_range, [0 20]) && isequal(ep.model.lya_T_interval, [10 20]));
es = analysis_run_config('sensitivity', 'essential', noisy);
assert(strcmp(es.model.ode_solver, 'sra1') && es.is_stochastic);

%% Every existing mode's numbers are unchanged
% {analysis, mode, n_levels, n_reps ([] = no reps axis), det solver, fs,
%  T_range, lya_T_interval ([] = class default), n_grid_points ([] = none)}
was = { ...
    'sensitivity',     'fast',       5,  5,  'rk4',   200, [0 10], [5 10],    []; ...
    'sensitivity',     'medium',     11, 15, 'rk4',   400, [0 20], [10 20],   []; ...
    'sensitivity',     'medium2',    13, 12, 'rk4',   800, [0 25], [12.5 25], []; ...
    'sensitivity',     'production', 25, 25, 'ode45', 400, [0 50], [20 50],   []; ...
    'tau_sensitivity', 'fast',       7,  7,  'rk4',   200, [0 20], [10 20],   []; ...
    'tau_sensitivity', 'medium',     11, 25, 'rk4',   400, [0 50], [],        []; ...
    'tau_sensitivity', 'medium2',    13, 15, 'rk4',   800, [0 50], [],        []; ...
    'tau_sensitivity', 'production', 25, 25, 'ode45', 400, [0 50], [],        []; ...
    'param_space',     'fast',       3,  [], 'rk4',   200, [0 20], [10 20],   128; ...
    'param_space',     'medium',     4,  [], 'rk4',   400, [0 20], [10 20],   128; ...
    'param_space',     'medium2',    5,  [], 'rk4',   800, [0 25], [12.5 25], 128; ...
    'param_space',     'production', 5,  [], 'ode45', 400, [0 50], [],        256};
assert(size(was, 1) == 3 * (numel(modes) - 1), 'the table must cover every mode but essential');
for k = 1:size(was, 1)
    [an, md, nl, nr, solver, fs_hz, T, lya, ngp] = was{k, :};
    for stochastic = [false true]
        want = struct('n_levels', nl);
        if ~isempty(nr); want.n_reps = nr; end
        want.model = struct('ode_solver', solver, 'fs', fs_hz, 'T_range', T);
        if ~isempty(lya); want.model.lya_T_interval = lya; end
        want.model.lya_method = 'topk'; want.model.lya_K = 15; want.model.lya_K_auto = true;
        want.model.lya_K_max = 60;      want.model.lya_dt = 0.05;
        want.sde_solver = 'sra1';
        if ~isempty(ngp); want.n_grid_points = ngp; end
        want.analysis = an; want.run_mode = md; want.is_stochastic = stochastic;
        if stochastic
            want.model.ode_solver = 'sra1';
            got = analysis_run_config(an, md, noisy);
        else
            got = analysis_run_config(an, md);
        end
        assert(isequal(got, want), 'analysis_run_config(''%s'', ''%s'') changed (stochastic = %d)', ...
            an, md, stochastic);
    end
end

%% Every stage's cost table names 'essential' beside 'fast'
% Text check. A stage that lacked the case would throw :badMode at run time,
% after the sweeps had finished (run_mode_names tells that story).
root = fileparts(which('setup_paths'));
stage_files = {'run_memory_capacity', 'run_memory_capacity_example', 'run_eig_heatmap', ...
    'run_numerics_verification', 'run_lyapunov_spectrum', 'run_transient_gain', ...
    'run_local_lyapunov', 'run_dc_lle_analysis'};
% A stage either shares 'fast's case (case {'essential', 'fast'}) or has a
% case of its own for 'essential'; in the second form a case for 'fast'
% alone beside it is right. Memory capacity is the one stage with its own:
% three trials where 'fast' has five (Tom, 2026-11-12).
own_case = {'run_memory_capacity'};
for k = 1:numel(stage_files)
    txt = fileread(which(stage_files{k}));
    code = regexprep(txt, '%[^\n]*', '');                       % comments off
    n_shared = numel(regexp(code, 'case\s*\{''essential'',\s*''fast''\}'));
    n_own    = numel(regexp(code, 'case\s*''essential'''));
    n_fast   = numel(regexp(code, 'case\s*''fast'''));
    assert(n_shared + n_own >= 1, '%s has no case for ''essential''', stage_files{k});
    if ismember(stage_files{k}, own_case)
        assert(n_own == 1 && n_shared == 0, '%s is expected to have its own case for ''essential''', stage_files{k});
        assert(n_own == n_fast, '%s: a case for ''fast'' alone needs a case for ''essential'' beside it', stage_files{k});
    else
        assert(n_fast == 0, '%s still has a case for ''fast'' alone, which ''essential'' would miss', stage_files{k});
    end
    assert(isempty(regexp(code, 'strcmp\([^)]*run_mode,\s*''fast''\)', 'once')), ...
        '%s tests run_mode against ''fast'' alone', stage_files{k});
end
% Memory capacity at 'essential': three trials, and the durations it shares
% with every mode untouched (read from the cost table's text; the table is
% a local function).
mc_code = regexprep(fileread(which('run_memory_capacity')), '%[^\n]*', '');
tok = regexp(mc_code, 'case\s*''essential''\s*cfg\.n_trials\s*=\s*(\d+)\s*;', 'tokens', 'once');
assert(~isempty(tok) && str2double(tok{1}) == 3, 'run_memory_capacity: ''essential'' must be 3 trials');
% The training window is a protocol constant, the same in every mode. Since
% the gain-by-bias stage (FD-010) it is set in mc_protocol, which both
% run_memory_capacity and run_mc_gain_sweep call; it must be set there once
% and nowhere in run_memory_capacity.
proto_code = regexprep(fileread(which('mc_protocol')), '%[^\n]*', '');
assert(numel(regexp(proto_code, 'T_train_sec\s*=')) == 1, ...
    'mc_protocol: the training window must be set once, for every mode');
assert(isempty(regexp(mc_code, 'T_train_sec\s*=\s*\d', 'once')), ...
    'run_memory_capacity must take the training window from mc_protocol, not set its own');
% The illustrations stage has no cost table to teach.
assert(~contains(fileread(which('run_paper_illustrations')), 'run_mode'));

%% paper_stage_plan
[all_stages, all_sweeps] = paper_stage_names();
assert(isequal(all_stages, {'sweeps', 'memory_capacity', 'mc_example', 'eig_heatmap', 'numerics', ...
    'lyapunov_spectrum', 'transient_gain', 'local_lyapunov', 'dc_lle', 'paper_illustrations'}));
assert(isequal(all_sweeps, {'sensitivity', 'tau', 'param_space'}));

base = struct('illustrations', true);
% Absent, and empty in each spelling: everything.
for c = {base, setfield(base, 'stages', {}), setfield(base, 'stages', []), ...
         setfield(setfield(base, 'stages', ''), 'sweeps', {})} %#ok<SFLD>
    pl = paper_stage_plan(c{1});
    assert(isequal(pl.stages, all_stages) && isempty(pl.stages_skipped));
    assert(isequal(pl.sweeps, all_sweeps) && isempty(pl.sweeps_skipped));
    assert(~pl.is_selection && all(struct2array(pl.runs)));
end
% A bundle without illustrations has nine stages, none of them skipped.
pl = paper_stage_plan(struct());
assert(isequal(pl.stages, all_stages(1:9)) && isempty(pl.stages_skipped) && ~pl.is_selection);
assert(~pl.runs.paper_illustrations);
pl = paper_stage_plan(struct('illustrations', false));
assert(isequal(pl.stages, all_stages(1:9)));

% A selection: the entry point's order, whatever order the config gives.
c = base; c.stages = {'paper_illustrations', 'local_lyapunov', 'sweeps', 'eig_heatmap'};
c.sweeps = {'param_space', 'sensitivity'};
pl = paper_stage_plan(c);
assert(isequal(pl.stages, {'sweeps', 'eig_heatmap', 'local_lyapunov', 'paper_illustrations'}));
assert(isequal(pl.stages_skipped, {'memory_capacity', 'mc_example', 'numerics', ...
    'lyapunov_spectrum', 'transient_gain', 'dc_lle'}));
assert(isequal(pl.sweeps, {'sensitivity', 'param_space'}) && isequal(pl.sweeps_skipped, {'tau'}));
assert(pl.is_selection && pl.runs.sweeps && pl.runs.eig_heatmap && ~pl.runs.transient_gain);
% Sweeps alone selected: every stage still runs.
pl = paper_stage_plan(setfield(base, 'sweeps', {'tau'})); %#ok<SFLD>
assert(isequal(pl.stages, all_stages) && isequal(pl.sweeps, {'tau'}) && pl.is_selection);
% One name as char; a repeated name counts once.
pl = paper_stage_plan(struct('stages', 'dc_lle'));
assert(isequal(pl.stages, {'dc_lle'}) && isempty(pl.sweeps) && isequal(pl.sweeps_skipped, all_sweeps));
pl = paper_stage_plan(struct('stages', {{'sweeps', 'sweeps'}}, 'sweeps', "sensitivity"));
assert(isequal(pl.stages, {'sweeps'}) && isequal(pl.sweeps, {'sensitivity'}));

% Errors, each naming the valid list.
assert_throws(@() paper_stage_plan(struct('stages', {{'sweeps', 'transient_gian'}})), ...
    'paper_stage_plan:UnknownStage', [all_stages, {'transient_gian'}]);
assert_throws(@() paper_stage_plan(struct('sweeps', {{'tau_sensitivity'}})), ...
    'paper_stage_plan:UnknownSweep', [all_sweeps, {'tau_sensitivity'}]);
% An unknown sweep is refused even when the sweeps stage is not selected.
assert_throws(@() paper_stage_plan(struct('stages', {{'dc_lle'}}, 'sweeps', {{'nope'}})), ...
    'paper_stage_plan:UnknownSweep', {'nope'});
assert_throws(@() paper_stage_plan(struct('stages', 3)), 'paper_stage_plan:BadSelection', {'stages'});
assert_throws(@() paper_stage_plan(struct('stages', {{'paper_illustrations'}})), ...
    'paper_stage_plan:IllustrationsOff', {'cfg.illustrations'});

%% The entry point's stage table is paper_stage_names' list
% Read from the file: the table is built after the sweeps (it needs run_dir),
% so a stage added there and not to paper_stage_names would otherwise be
% found by the entry point's own assert only after the sweeps had run.
code = regexprep(fileread(which('run_all_paper_analyses')), '%[^\n]*', '');
tok = regexp(code, '''(\w+)'',\s*@\(\)\s*run_', 'tokens');
table_names = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
assert(isequal(table_names, all_stages(2:end)), 'run_all_paper_analyses'' stage table is: %s', ...
    strjoin(table_names, ', '));
% The opt-in stages (FD-010) are not written in that table: the entry point
% appends the ones a bundle selects from paper_opt_in_stages, which is also
% where paper_stage_names takes their names from. Each has a plan and a run
% function that exist.
[~, ~, opt_in] = paper_stage_names();
registry = paper_opt_in_stages();
assert(isequal(opt_in, {registry.name}) && isequal(opt_in, {'mc_gain_sweep'}));
assert(isempty(intersect(opt_in, all_stages)), 'an opt-in stage must not also be in the first list');
for k = 1:numel(registry)
    assert(exist(func2str(registry(k).plan), 'file') == 2 && exist(func2str(registry(k).run), 'file') == 2, ...
        'opt-in stage %s: its plan or run function is missing', registry(k).name);
end

%% The pre-flight refuses a bad selection before anything exists
probe_rel = 'data/__essential_preflight_probe__';
probe = fullfile(root, probe_rel);
assert(~isfolder(probe), 'the probe directory %s is in the way', probe);
bad = mu7chain_config('essential'); bad.run_dir = probe_rel;
bad1 = bad; bad1.stages = [bad.stages, {'no_such_stage'}];
assert_throws(@() run_all_paper_analyses(bad1), 'paper_stage_plan:UnknownStage', {'no_such_stage'});
bad2 = bad; bad2.sweeps = {'sensitivity', 'no_such_sweep'};
assert_throws(@() run_all_paper_analyses(bad2), 'paper_stage_plan:UnknownSweep', {'no_such_sweep'});
bad3 = bad; bad3.eig_heatmap_examples = struct('label', 'reference');
assert_throws(@() run_all_paper_analyses(bad3), 'run_all_paper_analyses:BadEigExamples', {'overrides'});
assert(~isfolder(probe), 'a refused run must not create its run directory');
% run_all_analyses takes the sweeps as given and refuses an unknown one, too.
assert_throws(@() run_all_analyses(bad.preset_name, 'essential', 'output_dir', probe, ...
    'sweeps', {'tau_sensitivity'}), 'run_all_analyses:BadSweep', all_sweeps);
assert(~isfolder(probe), 'run_all_analyses must refuse a bad sweep before creating its directory');

%% run_all_analyses with no sweep: the directory, and a manifest that says what was skipped
% No compute: with an empty sweep list it creates the directory, the
% provenance record, run_manifest.mat and parameters.md, and runs nothing.
% This is the path run_all_paper_analyses takes for a bundle that selects
% only later stages, and the one place the manifest's selection fields are
% written.
tmp = tempname; cleanup_tmp = onCleanup(@() rmdir_if(tmp));
extra = struct('stages_run', {{'dc_lle'}}, 'stages_skipped', {{'sweeps', 'numerics'}});
out = run_all_analyses(bad.preset_name, 'essential', 'output_dir', tmp, 'sweeps', {}, ...
    'manifest_fields', extra, 'verbose', 'near-none', 'save_figs', false);
assert(strcmp(out, tmp) && isfile(fullfile(tmp, 'run_manifest.mat')));
M = load(fullfile(tmp, 'run_manifest.mat')); M = M.run_manifest;
assert(strcmp(M.run_mode, 'essential') && strcmp(M.preset_name, bad.preset_name));
assert(isequal(M.stages_run, {'dc_lle'}) && isequal(M.stages_skipped, {'sweeps', 'numerics'}));
assert(isempty(M.sweeps_run) && isequal(M.sweeps_skipped, all_sweeps));
listing = dir(tmp);
assert(~any([listing.isdir] & ~ismember({listing.name}, {'.', '..'})), 'no sweep folder may exist');
% A manifest field of the run's own cannot be replaced by the caller.
tmp2 = tempname; cleanup_tmp2 = onCleanup(@() rmdir_if(tmp2));
assert_throws(@() run_all_analyses(bad.preset_name, 'essential', 'output_dir', tmp2, 'sweeps', {}, ...
    'manifest_fields', struct('run_mode', 'fast'), 'verbose', 'near-none'), ...
    'run_all_analyses:ManifestFieldClash', {'run_mode'});

%% mu7chain_config('essential')
cfg = mu7chain_config('essential');
assert(strcmp(cfg.run_mode, 'essential') && strcmp(cfg.figure_set, '4cond'));
assert(strcmp(cfg.run_dir, 'data/mu7chain_essential'));
assert(strcmp(cfg.fig_root, 'figs/mu7chain_essential_4cond'));
assert(isequal(cfg.figure_sets, {'4cond'}), 'essential draws the four-condition set only');
assert(isequal(cfg.conditions, {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}));
assert(isequal(cfg.stages, {'sweeps', 'eig_heatmap', 'local_lyapunov', 'paper_illustrations'}));
assert(isequal(cfg.sweeps, {'sensitivity', 'param_space'}));
ex = cfg.eig_heatmap_examples;
assert(isstruct(ex) && isscalar(ex) && strcmp(ex.label, 'reference') && ...
    isstruct(ex.overrides) && isempty(fieldnames(ex.overrides)), 'the reference example only');
pl = paper_stage_plan(cfg);
assert(isequal(pl.stages, cfg.stages) && isequal(pl.sweeps, cfg.sweeps) && ...
    isequal(pl.sweeps_skipped, {'tau'}) && numel(pl.stages_skipped) == 6);

% The registry, exactly, in order.
want_figs = {'fig_sensitivity_analysis_allStd', 'fig_sensitivity_medians', 'fig_param_space_allStd', ...
    'fig_EI_weights_param_space', 'fig_lle_vs_rate', 'fig_eig_heatmap', ...
    'fig_main2_grouped', 'fig_main3_grouped', 'fig_main4_grouped'};
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
assert(isequal(names, want_figs), 'the essential registry is: %s', strjoin(names, ', '));
% Each entry is the fast registry's own, except that a grouped entry reads
% ITS set's fig_root (the grouped registry's rule).
fast = mu7chain_config('fast', '4cond');
fast_names = cellfun(@(f) f.name, fast.figures, 'UniformOutput', false);
for k = 1:numel(cfg.figures)
    e = cfg.figures{k}; f = fast.figures{strcmp(fast_names, e.name)};
    assert(isequal(e.fn, f.fn) && e.in_paper == f.in_paper, '%s is not the fast entry', e.name);
    if startsWith(e.name, 'fig_main')
        a = e.args; src = a{find(strcmp(a(1:2:end), 'source_fig_root')) * 2};
        assert(strcmp(src, cfg.fig_root), '%s reads %s, not its own %s', e.name, src, cfg.fig_root);
        a{find(strcmp(a(1:2:end), 'source_fig_root')) * 2} = fast.fig_root;
        assert(isequal(a, f.args));
    else
        assert(isequal(e.args, f.args), '%s: args differ from the fast entry''s', e.name);
    end
end
% The 3cond set can still be asked for; nothing makes it by default.
c3 = mu7chain_config('essential', '3cond');
assert(strcmp(c3.fig_root, 'figs/mu7chain_essential_3cond') && strcmp(c3.run_dir, cfg.run_dir));
assert(isequal(c3.stages, cfg.stages) && isequal(c3.figure_sets, {'4cond'}));
assert(isequal(cellfun(@(f) f.name, c3.figures, 'UniformOutput', false), want_figs));
% Everything else is the fast config's, field for field.
differs = {'run_mode', 'run_dir', 'fig_root', 'figures', 'stages', 'sweeps', ...
    'eig_heatmap_examples', 'figure_sets'};
assert(isequal(rmfield(cfg, differs), rmfield(fast, differs)), ...
    'essential differs from fast in more than: %s', strjoin(differs, ', '));

%% The other modes select nothing
for md = {'fast', 'medium', 'medium2', 'production'}
    for st = {'3cond', '4cond'}
        c = mu7chain_config(md{1}, st{1});
        assert(isempty(c.stages) && isempty(c.sweeps) && isempty(c.eig_heatmap_examples), ...
            '%s must select nothing', md{1});
        assert(isequal(c.figure_sets, {'3cond', '4cond'}));
        pl = paper_stage_plan(c);
        assert(~pl.is_selection && isequal(pl.stages, all_stages) && isequal(pl.sweeps, all_sweeps));
        assert(numel(c.figures) == numel(fast.figures));
    end
end

disp('PASS: test_essential_mode -- the mode, its sizes, the unchanged modes, the stage plan, the pre-flight, the essential bundle.');

%% ------------------------------------------------------------------------
function tf = same_but(a, b, fields)
% a and b are equal once the named fields (those present) are removed.
fa = fields(isfield(a, fields)); fb = fields(isfield(b, fields));
tf = isequal(rmfield(a, fa), rmfield(b, fb));
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end

function assert_throws(fn, want_id, must_name)
% fn throws want_id, and its message names every string in must_name.
try
    fn();
catch err
    assert(strcmp(err.identifier, want_id), 'expected %s, got %s: %s', want_id, err.identifier, err.message);
    for k = 1:numel(must_name)
        assert(contains(err.message, must_name{k}), '%s does not name ''%s'': %s', ...
            want_id, must_name{k}, err.message);
    end
    return
end
error('test_essential_mode:NoThrow', 'expected %s, and nothing was thrown', want_id);
end
