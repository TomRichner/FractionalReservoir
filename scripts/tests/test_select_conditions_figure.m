%% TEST_SELECT_CONDITIONS_FIGURE The 'conditions' option end to end (FD-006).
%
% Part 1, fig_lle_vs_rate (the first figure converted), against an existing
% regression run, read-only (three conditions: no_adaptation / sfa1_std1 /
% sfa3_std2, produced by regression_fast_config):
%   - with no 'conditions' it draws what it drew before the option existed:
%     one panel per condition the run holds, in the run's order, titled by
%     name;
%   - with a subset it draws exactly that subset, in the requested order, and
%     its table lists them in that order;
%   - a name the run lacks is an error, not a missing panel.
%
% Part 2, make_all_paper_figures with a small registry: the converted figure
% and a stand-in that does NOT declare 'conditions'. Checks
% results.conditions_applied ('yes' / 'no' / ''), the minimal line naming the
% entry that did not take the selection, the manifest's Conditions row and
% per-entry column, that an empty selection passes nothing, and that a
% selection naming a condition the run lacks FAILS the converted entry rather
% than falling back to all conditions.
%
% The run is found at <project root>/data/regression_fast/reference_post_refactor,
% else at the same path in the main checkout (so it runs from a git worktree,
% whose data/ is empty); the test SKIPs if neither exists. Figures are
% invisible and written to a temporary directory; nothing is written into the
% run. include_sweeps is false throughout: the weight-balance rebuild is the
% cost (about 20 s per call at regression size) and the sweeps add nothing to
% what is checked. Part 2 needs a pool allocation (make_all_paper_figures
% refuses without one); pool_set(0) is enough.

n_fail = 0;

run_dir = find_regression_run();
if isempty(run_dir)
    fprintf('SKIP: no regression run directory (data/regression_fast/reference_post_refactor)\n');
    return
end
fprintf('run: %s\n', run_dir);

st  = manuscript_style();
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));
base = {'run_dir', run_dir, 'visible', false, 'verbose', 'near-none', 'include_sweeps', false};

% The run's own condition order, from its param_space PSA.
ps = dir(fullfile(run_dir, 'param_space_*'));
psa = ParamSpaceAnalysis2.from_dir(fullfile(ps(1).folder, ps(1).name));
run_conds = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
n_fail = n_fail + check(isequal(run_conds, {'no_adaptation', 'sfa1_std1', 'sfa3_std2'}), ...
    'regression run holds the three expected conditions');

%% Part 1: fig_lle_vs_rate
% (a) no 'conditions': every condition, the run's order, titled by name
out = fig_lle_vs_rate(base{:}, 'save', false);
t = panel_titles(out.figs);
want = cellfun(@(c) st.condition_title(c), run_conds, 'UniformOutput', false);
n_fail = n_fail + check(isequal(t, want), sprintf('no selection: %d panels, run order, titled by name', numel(t)));
close(out.figs);

% (b) a subset, out of the run's order
sel = {'sfa3_std2', 'no_adaptation'};
out = fig_lle_vs_rate(base{:}, 'save', true, 'out_dir', fullfile(tmp, 'sub'), 'conditions', sel);
t = panel_titles(out.figs);
want = cellfun(@(c) st.condition_title(c), sel, 'UniformOutput', false);
n_fail = n_fail + check(isequal(t, want), 'subset: exactly the selected panels, in the requested order');
tbl = fileread(fullfile(tmp, 'sub', 'Fig_LLE_vs_Rate_table.md'));
rows = regexp(tbl, '^\| (\w+) \| [+-]', 'tokens', 'lineanchors');
rows = cellfun(@(c) c{1}, rows, 'UniformOutput', false);
n_fail = n_fail + check(isequal(rows, sel) && ...
    contains(tbl, 'Conditions drawn (selected, in this order): sfa3_std2, no_adaptation'), ...
    'subset: table rows in the requested order, selection stated');
n_fail = n_fail + check(~isempty(out.files), 'subset: files written');
close(out.figs);

% (c) a name the run lacks is an error, not a missing panel
n_open = numel(findall(0, 'Type', 'figure'));
id = error_of(@() fig_lle_vs_rate(base{:}, 'save', false, ...
    'conditions', {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:NoSuchCondition'), ...
    'a condition the run lacks -> NoSuchCondition');
n_fail = n_fail + check(numel(findall(0, 'Type', 'figure')) == n_open, 'no figure left open by the error');

%% Part 2: make_all_paper_figures records whether the selection was applied
pc = pool_config();
if ~pc.is_set
    fprintf('SKIP Part 2: no pool allocation in this session (pool_set(0) is enough)\n');
else
    stub_dir = fullfile(tmp, 'stub'); mkdir(stub_dir);
    write_stub(fullfile(stub_dir, 'fig_selcond_stub.m'));
    addpath(stub_dir);
    cleanup_path = onCleanup(@() rmpath(stub_dir));

    e_fig  = struct('name', 'lle_vs_rate', 'fn', @fig_lle_vs_rate, 'in_paper', true, ...
                    'args', {{'include_sweeps', false}});
    e_stub = struct('name', 'stub', 'fn', @fig_selcond_stub, 'in_paper', false, 'args', {{}});

    cfg = struct();
    cfg.preset_name = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
    cfg.run_mode    = 'fast';
    cfg.run_dir     = run_dir;
    cfg.verbose     = 'minimal';
    cfg.visible_figures = false;

    % (a) with a selection: applied to the figure, not to the stub, and said so
    cfg.figures    = {e_fig, e_stub};
    cfg.conditions = {'sfa1_std1', 'no_adaptation'};
    cfg.fig_root   = fullfile(tmp, 'pass_sel');
    txt = evalc('res = make_all_paper_figures(cfg);');
    n_fail = n_fail + check(all([res.ok]), 'selection pass: both entries built');
    n_fail = n_fail + check(isequal({res.conditions_applied}, {'yes', 'no'}), ...
        'selection pass: applied yes for the converted figure, no for the stub');
    n_fail = n_fail + check(contains(txt, 'NOT applied') && contains(txt, ': stub'), ...
        'selection pass: one minimal line names the entry that did not take it');
    man = fileread(fullfile(cfg.fig_root, 'manifest.md'));
    n_fail = n_fail + check(contains(man, '| Conditions | `sfa1_std1`, `no_adaptation` (in this order) |'), ...
        'manifest: requested conditions row');
    n_fail = n_fail + check(~isempty(regexp(man, '\| `lle_vs_rate` \|[^\n]*\| yes \| ok \|', 'once')) && ...
        ~isempty(regexp(man, '\| `stub` \|[^\n]*\| \*\*no\*\* \(all drawn\) \| ok \|', 'once')) && ...
        contains(man, 'NOT applied to 1 entry'), ...
        'manifest: per-entry applied column and the not-applied note');
    tbl = fileread(fullfile(cfg.fig_root, 'lle_vs_rate', 'Fig_LLE_vs_Rate_table.md'));
    n_fail = n_fail + check(contains(tbl, 'Conditions drawn (selected, in this order): sfa1_std1, no_adaptation'), ...
        'selection pass: the figure drew the selection');
    stub_note = fileread(fullfile(cfg.fig_root, 'stub', 'stub.txt'));
    n_fail = n_fail + check(contains(stub_note, ['run_dir=' run_dir]), ...
        'stub rebuilt without the option, run_dir still passed');

    % (b) without a selection: nothing passed, nothing reported (the stub
    % alone is enough to show the option is not passed)
    cfg.figures    = {e_stub};
    cfg.conditions = {};
    cfg.fig_root   = fullfile(tmp, 'pass_all');
    txt = evalc('res = make_all_paper_figures(cfg);');
    n_fail = n_fail + check(res.ok && strcmp(res.conditions_applied, ''), ...
        'no selection: option not passed, nothing recorded as applied');
    n_fail = n_fail + check(~contains(txt, 'NOT applied'), 'no selection: no not-applied line');
    man = fileread(fullfile(cfg.fig_root, 'manifest.md'));
    n_fail = n_fail + check(contains(man, '| Conditions | all the run holds (none selected) |') && ...
        ~contains(man, 'NOT applied'), 'manifest without a selection');
    cfg = rmfield(cfg, 'conditions');   % absent field = all, too
    cfg.fig_root = fullfile(tmp, 'pass_absent');
    evalc('res = make_all_paper_figures(cfg);');
    n_fail = n_fail + check(res.ok && strcmp(res.conditions_applied, ''), 'absent cfg.conditions = all');

    % (c) a selection naming a condition the run lacks FAILS the converted entry
    cfg.figures    = {e_fig, e_stub};
    cfg.conditions = {'sfa3_stdc2'};
    cfg.fig_root   = fullfile(tmp, 'pass_bad');
    evalc('res = make_all_paper_figures(cfg);');
    n_fail = n_fail + check(~res(1).ok && contains(res(1).err, 'select_conditions:NoSuchCondition') && ...
        res(2).ok && strcmp(res(2).conditions_applied, 'no'), ...
        'missing condition: the converted entry FAILS (no fallback); the stub builds');
    clear cleanup_path
end
close(findall(0, 'Type', 'figure'));
clear cleanup_tmp

%% Verdict
if n_fail == 0
    fprintf('\nALL select_conditions FIGURE CHECKS PASSED\n');
else
    error('test_select_conditions_figure:Failed', '%d check(s) FAILED', n_fail);
end

%% ------------------------------------------------------------------------
function d = find_regression_run()
% <project root>/data/..., else the main checkout's (a worktree's data/ is empty).
rel = fullfile('data', 'regression_fast', 'reference_post_refactor');
root = fileparts(which('setup_paths'));
d = fullfile(root, rel);
if isfolder(d); return; end
[s, common] = system(sprintf('git -C "%s" rev-parse --path-format=absolute --git-common-dir', root));
d = '';
if s == 0
    main_root = fileparts(strtrim(common));   % <main checkout>/.git -> <main checkout>
    if isfolder(fullfile(main_root, rel)); d = fullfile(main_root, rel); end
end
end

function f = check(ok, what)
if ok
    fprintf('PASS  %s\n', what);
    f = 0;
else
    fprintf('FAIL  %s\n', what);
    f = 1;
end
end

function [id, msg] = error_of(fn)
id = ''; msg = '';
try
    fn();
catch ME
    id = ME.identifier; msg = ME.message;
end
end

function t = panel_titles(fig)
% Axes titles, left to right by tile.
ax = findobj(fig, 'Type', 'axes');
tiles = arrayfun(@(a) a.Layout.Tile, ax);
[~, o] = sort(tiles);
ax = ax(o);
t = arrayfun(@(a) char(a.Title.String), ax, 'UniformOutput', false)';
end

function write_stub(path)
% A figure that does NOT declare 'conditions': records what it was given.
fid = fopen(path, 'w');
fprintf(fid, '%s\n', ...
    'function out = fig_selcond_stub(cfg)', ...
    'arguments', ...
    '    cfg.run_dir (1,:) char = ''''', ...
    '    cfg.out_dir (1,:) char = ''''', ...
    '    cfg.save (1,1) logical = true', ...
    '    cfg.visible (1,1) logical = false', ...
    '    cfg.verbose (1,:) char = ''minimal''', ...
    'end', ...
    'if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end', ...
    'p = fullfile(cfg.out_dir, ''stub.txt'');', ...
    'fid = fopen(p, ''w''); fprintf(fid, ''run_dir=%s\n'', cfg.run_dir); fclose(fid);', ...
    'out = struct(''figs'', gobjects(0), ''files'', {{p}}, ''source'', cfg.run_dir);', ...
    'end');
fclose(fid);
end
