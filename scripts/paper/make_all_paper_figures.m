function results = make_all_paper_figures(cfg)
% MAKE_ALL_PAPER_FIGURES Rebuild every manuscript figure from one run directory.
%
%   results = MAKE_ALL_PAPER_FIGURES()
%   results = MAKE_ALL_PAPER_FIGURES(paper_config('run_dir', d))
%   results = MAKE_ALL_PAPER_FIGURES(paper_config('fig_root', 'figs/try'))
%
% Takes ONE argument, a config struct -- not name-value pairs. (This header used
% to advertise an 'only' option the arguments block never accepted.)
%
% ENTRY POINT TWO OF TWO. The cheap half: pure replots, analytic figures, and a
% handful that simulate but finish in seconds to minutes. Everything measured in
% hours lives in run_all_paper_analyses.
%
% ONE RUN, RESOLVED ONCE. The run directory is resolved here and passed to every
% figure, so the manuscript cannot end up built from two different runs. Before
% this, five figures hardcoded a data_root and they did not agree: one pointed
% at an old SRNNModel2 run while the others pointed at a SRNNCellTypePairs run.
%
% PER-FIGURE ERROR ISOLATION. One figure failing does not cost the rest. Each
% result records whether files actually landed on disk -- not merely that the
% call returned -- because save_some_figs_to_folder_2 catches export failures,
% warns, and carries on, so a figure can "succeed" having written nothing.
%
% THE MASTER CLOSES EACH FIGURE'S HANDLES, once that figure has been verified.
% No figure function calls `close all force` any more -- in a loop that destroys
% the PREVIOUS entry's output before it can be checked -- but something must
% still close them, and the caller is the only party that knows when an entry is
% finished.
%
% This is not tidiness. replot_sensitivity and replot_param_space_analysis save
% ALL OPEN FIGURES into their prep folder, so leaving earlier entries' figures
% open pollutes those folders: observed directly on the first full run, where
% replot_sensitivity swept up 15 figures including the bursting panels from four
% entries earlier. The final sheets were still correct (they are built from
% explicitly named handles), but the intermediate folder was wrong, and a future
% figure that globs a prep folder would inherit the mess.
%
% See also: paper_config, run_all_paper_analyses, resolve_run_dir

arguments
    cfg = paper_config()
end

setup_paths();

%% Resolve the run once
% Figures that read no run directory (the analytic ones) ignore this; passing it
% to them is harmless and keeps the call uniform.
%
% cfg.run_dir is the SAME field run_all_paper_analyses writes to: set, it is
% read here; empty, resolve_run_dir searches data/param_space for the newest
% run_all_* whose manifest names this preset.
%
% The relative-path resolution matters. resolve_run_dir's search only ever
% looks in data/param_space/run_all_*, so a config naming its own directory --
% fast4_config's 'data/fast_4' -- wrote somewhere the search will never find.
% Passing it through as an explicit run_dir is what makes such a config work
% end to end; without it the figures died with NoMatchingRun.
want_run_dir = cfg.run_dir;
if ~isempty(want_run_dir) && ~is_absolute_path(want_run_dir)
    want_run_dir = fullfile(fileparts(which('setup_paths')), want_run_dir);
end
try
    run_dir = resolve_run_dir('run_dir', want_run_dir, 'preset_name', cfg.preset_name);
    fprintf('Run directory: %s\n', run_dir);
catch ME
    warning('make_all_paper_figures:NoRun', ...
        ['No run directory matched preset ''%s'':\n  %s\n' ...
         'Replot figures will fail; analytic and self-simulating ones will not.'], ...
        cfg.preset_name, ME.message);
    run_dir = '';
end

%% Resolve where the output goes
% Two modes, per cfg.fig_root (see paper_config): a configured path is STABLE
% and overwritten, so the manuscript can cite something that never moves; empty
% auto-names a dated folder, so a one-off cannot clobber the paper's set. Same
% convention as run_all_analyses' output_dir.
%
% Every entry gets <fig_root>/<entry name>. The per-entry subdirectory is not
% cosmetic: save_figure_stable deletes <out_dir>/<fig_tag>* before saving, and
% 'Fig_single_neuron_SFA_STD' is a strict PREFIX of
% 'Fig_single_neuron_SFA_STD_STF'. Sharing one directory, the first tag's save
% destroys the second's outputs -- which is live today in fig_adaptation_methods
% and survives only because the registry happens to run them in that order.
project_root = fileparts(which('setup_paths'));
if isempty(cfg.fig_root)
    dt_str   = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
    fig_root = fullfile(project_root, 'figs', sprintf('figures_%s', dt_str));
    root_mode = 'auto-named';
else
    fig_root = cfg.fig_root;
    if ~is_absolute_path(fig_root)
        fig_root = fullfile(project_root, fig_root);   % relative -> project root
    end
    root_mode = 'configured';
end
if ~isfolder(fig_root); mkdir(fig_root); end
fprintf('Figure root:   %s  (%s)\n', fig_root, root_mode);

figs = cfg.figures;
n = numel(figs);

fprintf('\n========================================================\n');
fprintf('PAPER FIGURES  (%d entries, preset %s)\n', n, cfg.preset_name);
fprintf('  start: %s\n', datetime('now'));
fprintf('========================================================\n');
t_all = tic;

results = struct('name', {}, 'ok', {}, 'in_paper', {}, 'n_files', {}, ...
                 'seconds', {}, 'files', {}, 'err', {});

for k = 1:n
    f = figs{k};
    fprintf('\n---- [%d/%d] %s%s\n', k, n, f.name, paper_tag(f.in_paper));
    t0 = tic;
    try
        args = [{'run_dir', run_dir, 'out_dir', fullfile(fig_root, f.name), ...
                 'save', true, 'visible', false}, f.args];
        out  = call_figure(f.fn, args);

        n_files = numel(out.files);
        ok = n_files > 0;
        if ~ok
            % A figure that returns without writing anything is a FAILURE here,
            % even though the call succeeded. See the header.
            err = 'returned no files (export may have failed silently)';
            fprintf(2, '     no files written\n');
        else
            err = '';
            fprintf('     %d file(s): %s\n', n_files, strjoin(out.files, ', '));
        end
        results(end+1) = struct('name', f.name, 'ok', ok, 'in_paper', f.in_paper, ...
            'n_files', n_files, 'seconds', toc(t0), 'files', {out.files}, ...
            'err', err); %#ok<AGROW>

        % Close only THIS entry's handles, now that it is verified. See header.
        close_figs(out.figs);
    catch ME
        fprintf(2, '     FAILED %s: %s\n', ME.identifier, ME.message);
        for s = 1:min(3, numel(ME.stack))
            fprintf(2, '       at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
        end
        results(end+1) = struct('name', f.name, 'ok', false, 'in_paper', f.in_paper, ...
            'n_files', 0, 'seconds', toc(t0), 'files', {{}}, ...
            'err', sprintf('%s: %s', ME.identifier, ME.message)); %#ok<AGROW>
    end
end

%% Close anything the replot helpers left behind
% Each entry's OWN handles are closed as it finishes (above). This catches the
% prep figures that replot_sensitivity / replot_param_space_analysis create and
% do not return -- they are invisible to the per-entry cleanup, and left open
% they get swept into the NEXT pass's prep folder and produce
% MATLAB:class:InvalidHandle warnings when a stale handle is saved. Safe here
% because every entry has already been verified.
stray = findobj(0, 'Type', 'figure');
if ~isempty(stray)
    fprintf('\n     closing %d stray prep figure(s) left by the replot helpers\n', ...
        numel(stray));
    close(stray);
end

%% Provenance and manifest, once at the root
% Replaces the per-figure README_*.txt and the sixteen near-identical
% git_provenance.txt files. capture_git_provenance also writes
% working_changes.patch when the tree is dirty, which is what lets a run made
% from uncommitted work be reproduced.
try
    capture_git_provenance(fig_root, project_root);
    write_figure_manifest(fig_root, cfg, run_dir, root_mode, results, toc(t_all)/60);
catch ME
    warning('make_all_paper_figures:ManifestFailed', ...
        'Could not write the manifest: %s', ME.message);
end

%% Summary
n_ok    = sum([results.ok]);
n_paper = sum([results.in_paper]);
n_paper_ok = sum([results.ok] & [results.in_paper]);

fprintf('\n========================================================\n');
fprintf('PAPER FIGURES COMPLETE in %.1f min\n', toc(t_all)/60);
fprintf('  %d/%d succeeded   (%d/%d of the in-paper figures)\n', ...
    n_ok, n, n_paper_ok, n_paper);
fprintf('  output: %s\n', fig_root);
fprintf('--------------------------------------------------------\n');
for k = 1:numel(results)
    r = results(k);
    if r.ok; tag = 'OK    '; else; tag = 'FAILED'; end
    fprintf('  %s  %-32s %2d files %6.1f s%s\n', tag, r.name, r.n_files, ...
        r.seconds, paper_tag(r.in_paper));
end
bad = results(~[results.ok]);
if ~isempty(bad)
    fprintf('--------------------------------------------------------\n');
    for k = 1:numel(bad)
        fprintf(2, '  %s: %s\n', bad(k).name, bad(k).err);
    end
end
fprintf('========================================================\n');
end

%% ------------------------------------------------------------------------
function tf = is_absolute_path(p)
% Absolute on either platform: a leading separator, or a Windows drive letter.
%
% Deliberately NOT isfolder(): whether a relative path happens to exist under
% the current directory must not decide where output lands, or the same config
% would write to two different places depending on cwd.
tf = startsWith(p, filesep) || startsWith(p, '/') || ...
     ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end

function out = call_figure(fn, args)
% Call a figure function, tolerating one that does not accept run_dir.
%
% Every converted figure declares run_dir in its arguments block precisely so
% this call can be uniform -- but a newly added one might not, and failing the
% whole sheet over an unused argument would be a poor trade.
try
    out = fn(args{:});
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction') || ...
       contains(ME.message, 'run_dir')
        keep = true(1, numel(args));
        idx = find(strcmp(args, 'run_dir'), 1);
        if ~isempty(idx); keep(idx:idx+1) = false; end
        out = fn(args{keep});
    else
        rethrow(ME);
    end
end
end

function s = paper_tag(in_paper)
if in_paper; s = '   [in paper]'; else; s = ''; end
end

function close_figs(h)
% Close a figure array, tolerating handles a figure function already closed.
if isempty(h); return; end
for k = 1:numel(h)
    if isgraphics(h(k), 'figure')
        close(h(k));
    end
end
end
