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
% WHICH CONDITIONS ARE DRAWN: cfg.conditions (FD-006). A cell of condition
% names, in the order the figures lay them out; absent or empty means every
% condition the run holds. When non-empty it is passed to every entry as
% 'conditions', and each figure resolves it with select_conditions (an error on
% a name the run lacks). An entry whose function does not declare the option is
% built without it, drawing every condition; that is printed as one minimal
% line after the pass and recorded per entry in the manifest
% (results(k).conditions_applied = 'yes' | 'no' | '').
%
% See also: paper_config, run_all_paper_analyses, resolve_run_dir, select_conditions

arguments
    cfg = paper_config()
end

setup_paths();

%% PREFLIGHT: this session must have a pool allocation
% Several figures simulate, some through a parfor. With no allocation set
% those loops would run serially without anyone having chosen that, so the
% pass is refused before anything is resolved or created. pool_set(0) is an
% allocation (serial, on purpose) and is accepted. Same check as
% run_all_paper_analyses.
pc = pool_config();
if ~pc.is_set
    error('make_all_paper_figures:NoPoolAllocation', ...
        ['No parallel-pool allocation is set for this MATLAB session, so the figure pass is ' ...
         'refused.\n%s'], pc.how_to_set);
end

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
% The verbosity level: needed from the first line that prints, so it is resolved
% here, before the run directory is reported.
if ~isfield(cfg, 'verbose'), cfg.verbose = 'minimal'; end
verbose = cfg.verbose;

% Which conditions the figures draw (FD-006). cfg.conditions is a cell of
% condition names, in the left-to-right order the figures lay them out; absent
% or empty means every condition the run holds, which is what every figure
% drew before the option existed. It is passed to each entry as 'conditions'
% ONLY when non-empty, so a bundle without it calls every figure exactly as
% before. Each figure resolves it through select_conditions, which errors on a
% name the run lacks. An entry whose function does not declare the option is
% built without it (all conditions) and reported, never silently.
if ~isfield(cfg, 'conditions'), cfg.conditions = {}; end
want_conditions = cfg.conditions;
if isempty(want_conditions) || ((ischar(want_conditions) || isstring(want_conditions)) && ...
        all(strlength(want_conditions) == 0))
    want_conditions = {};                    % {}, [], '' and "" mean all, as in select_conditions
end
want_conditions = cellstr(want_conditions);  % char / string / cell -> cellstr
cfg.conditions = want_conditions;            % the manifest reads the normalised list

want_run_dir = cfg.run_dir;
if ~isempty(want_run_dir) && ~is_absolute_path(want_run_dir)
    want_run_dir = fullfile(fileparts(which('setup_paths')), want_run_dir);
end
try
    run_dir = resolve_run_dir('run_dir', want_run_dir, 'preset_name', cfg.preset_name);
    vprintf(verbose, 'minimal', 'Run directory: %s\n', run_dir);
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
% The transcript of the figure pass is saved beside the manifest (TR,
% 2026-09-14), through vlog since 2026-09-15 (diary lost everything under the
% MATLAB MCP server; see vlog). Process-global; the onCleanup closes it
% however this function exits.
log_guard = vlog('open', fullfile(fig_root, 'command_window.log')); %#ok<NASGU>
vprintf(verbose, 'minimal', '[log] transcript -> %s\n', fullfile(fig_root, 'command_window.log'));
vprintf(verbose, 'minimal', 'Figure root:   %s  (%s)\n', fig_root, root_mode);

%% Keep the figures off screen unless asked otherwise
% cfg.visible_figures, default false. Each entry is passed `visible` below --
% that is the knob the fig_* contract already has -- but BOTH are needed:
%
%   the per-figure argument   covers the figure each fig_* creates itself
%   the root default          covers everything else drawn during the pass --
%                             the prep figures the replot_* helpers create and
%                             do not return, and any plotting primitive that
%                             calls figure() without asking
%
% with_graphics_defaults returns an onCleanup, so the session's own setting is
% restored when this function returns, including on the error path. It must stay
% in scope for the whole pass; do not clear it.
visible_figures = isfield(cfg, 'visible_figures') && cfg.visible_figures;
% The verbosity level, threaded into every figure as its `verbose` argument;
% older callers' cfg structs may lack the field.
if ~visible_figures
    fig_visibility = with_graphics_defaults('DefaultFigureVisible', 'off'); %#ok<NASGU>
end

figs = cfg.figures;
n = numel(figs);

vprintf(verbose, 'verbose', '\n========================================================\n');
vprintf(verbose, 'minimal', 'PAPER FIGURES  (%d entries, preset %s)  start %s\n', n, cfg.preset_name, datetime('now'));
if isempty(want_conditions)
    vprintf(verbose, 'minimal', 'Conditions:    all the run holds (none selected)\n');
else
    vprintf(verbose, 'minimal', 'Conditions:    %s\n', strjoin(want_conditions, ', '));
end
vprintf(verbose, 'verbose', '========================================================\n');
t_all = tic;

% conditions_applied: 'yes' (the entry took the selection), 'no' (its function
% does not declare 'conditions', so it drew every condition), or '' (no
% selection requested, or the entry failed before that was known).
results = struct('name', {}, 'ok', {}, 'in_paper', {}, 'n_files', {}, ...
                 'seconds', {}, 'files', {}, 'err', {}, 'conditions_applied', {});

for k = 1:n
    f = figs{k};
    vprintf(verbose, 'minimal', '---- [%d/%d] %s%s\n', k, n, f.name, paper_tag(f.in_paper));
    t0 = tic;
    try
        args = {'run_dir', run_dir, 'out_dir', fullfile(fig_root, f.name), ...
                'save', true, 'visible', visible_figures, 'verbose', verbose};
        if ~isempty(want_conditions)
            args = [args, {'conditions', want_conditions}]; %#ok<AGROW>
        end
        args = [args, f.args]; %#ok<AGROW>
        [out, applied] = call_figure(f.fn, args);
        if isempty(applied)
            cond_status = '';
        elseif applied
            cond_status = 'yes';
        else
            cond_status = 'no';
            vprintf(verbose, 'verbose', '     does not declare ''conditions'': built with every condition\n');
        end

        n_files = numel(out.files);
        ok = n_files > 0;
        if ~ok
            % A figure that returns without writing anything is a FAILURE here,
            % even though the call succeeded. See the header.
            err = 'returned no files (export may have failed silently)';
            vfail('     no files written\n');
        else
            err = '';
            vprintf(verbose, 'minimal', '     %d file(s): %s\n', n_files, strjoin(out.files, ', '));
        end
        results(end+1) = struct('name', f.name, 'ok', ok, 'in_paper', f.in_paper, ...
            'n_files', n_files, 'seconds', toc(t0), 'files', {out.files}, ...
            'err', err, 'conditions_applied', cond_status); %#ok<AGROW>

        % Close only THIS entry's handles, now that it is verified. See header.
        close_figs(out.figs);
    catch ME
        vfail('     FAILED %s: %s\n', ME.identifier, ME.message);
        for s = 1:min(3, numel(ME.stack))
            vfail('       at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
        end
        results(end+1) = struct('name', f.name, 'ok', false, 'in_paper', f.in_paper, ...
            'n_files', 0, 'seconds', toc(t0), 'files', {{}}, ...
            'err', sprintf('%s: %s', ME.identifier, ME.message), ...
            'conditions_applied', ''); %#ok<AGROW>
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
    vprintf(verbose, 'verbose', '\n     closing %d stray prep figure(s) left by the replot helpers\n', ...
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

%% The scrollable report, last
% report.md links every PNG under fig_root in explorer order, and report.pdf is
% pandoc's render of it. AFTER the manifest, because the report echoes the
% manifest's provenance rows. Its own try/catch: it reads what is on disk and
% can only describe the run, so a failure here must not change the result
% count above -- and the PDF step already degrades to a warning on its own.
try
    write_figure_report(fig_root, 'pdf', true);
catch ME
    warning('make_all_paper_figures:ReportFailed', ...
        'Could not write the figure report: %s', ME.message);
end

%% Summary
n_ok    = sum([results.ok]);
n_paper = sum([results.in_paper]);
n_paper_ok = sum([results.ok] & [results.in_paper]);

vprintf(verbose, 'verbose', '\n========================================================\n');
vprintf(verbose, 'near-none', 'PAPER FIGURES COMPLETE in %.1f min: %d/%d succeeded (%d/%d in-paper), output %s\n', ...
    toc(t_all)/60, n_ok, n, n_paper_ok, n_paper, fig_root);
vprintf(verbose, 'minimal', '--------------------------------------------------------\n');
for k = 1:numel(results)
    r = results(k);
    if r.ok; tag = 'OK    '; else; tag = 'FAILED'; end
    vprintf(verbose, 'minimal', '  %s  %-32s %2d files %6.1f s%s\n', tag, r.name, r.n_files, ...
        r.seconds, paper_tag(r.in_paper));
end
% The selection is never dropped silently: one line naming every entry that
% was built with all conditions because its function does not declare the
% option (the manifest records the same per entry).
not_applied = {results(strcmp({results.conditions_applied}, 'no')).name};
if ~isempty(not_applied)
    vprintf(verbose, 'minimal', ['  conditions {%s} NOT applied (figure does not declare ' ...
        '''conditions''; drew every condition) in %d entr%s: %s\n'], ...
        strjoin(want_conditions, ', '), numel(not_applied), ...
        plural_y(numel(not_applied)), strjoin(not_applied, ', '));
end
bad = results(~[results.ok]);
if ~isempty(bad)
    vprintf(verbose, 'minimal', '--------------------------------------------------------\n');
    for k = 1:numel(bad)
        vfail('  %s: %s\n', bad(k).name, bad(k).err);
    end
end
vprintf(verbose, 'minimal', '========================================================\n');
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

function [out, applied] = call_figure(fn, args)
% Call a figure function, tolerating one that does not accept run_dir, or one
% that does not (yet) declare 'conditions'.
%
% Every converted figure declares run_dir in its arguments block precisely so
% this call can be uniform -- but a newly added one might not, and failing the
% whole sheet over an unused argument would be a poor trade.
%
% 'conditions' (FD-006) is passed only when a selection was requested. A figure
% whose arguments block does not declare it is called again WITHOUT it and so
% draws every condition the run holds; `applied` comes back false so the
% caller can say so (one minimal line, and the manifest). The retry is taken
% only on MATLAB's own rejection of the argument NAME ("Invalid argument name
% 'conditions'", raised before the function body runs), never on an error the
% figure raises itself -- a figure that declares the option and is asked for a
% condition the run lacks must FAIL (select_conditions:NoSuchCondition), not
% fall back to drawing everything. A wrapper that forwards to an inner figure
% that does not declare it is caught the same way, and is likewise recorded as
% not applied.
%
% `applied` is true when 'conditions' was passed and accepted, false when it
% had to be dropped, and [] when it was not passed (nothing requested).
has_cond = any(strcmp(args(1:2:end), 'conditions'));
if has_cond; applied = true; else; applied = []; end
try
    out = call_tolerating_run_dir(fn, args);
catch ME
    if has_cond && rejects_argument_name(ME, 'conditions')
        keep = true(1, numel(args));
        idx = find(strcmp(args, 'conditions'), 1);
        keep(idx:idx+1) = false;
        out = call_tolerating_run_dir(fn, args(keep));
        applied = false;
    else
        rethrow(ME);
    end
end
end

function out = call_tolerating_run_dir(fn, args)
% The original call_figure: retry once without run_dir. The one addition: a
% rejection of 'conditions' is passed up to call_figure untouched (its message
% lists the declared names, run_dir among them, which would otherwise read as
% a run_dir problem here).
try
    out = fn(args{:});
catch ME
    if rejects_argument_name(ME, 'conditions')
        rethrow(ME);
    elseif strcmp(ME.identifier, 'MATLAB:UndefinedFunction') || ...
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

function tf = rejects_argument_name(ME, name)
% True when ME is an arguments-block (or inputParser) rejection of `name`.
% Measured in R2025: identifier MATLAB:TooManyInputs, message "Invalid
% argument name 'conditions'. Use one of: ...".
tf = (strcmp(ME.identifier, 'MATLAB:TooManyInputs') && ...
      contains(ME.message, sprintf('argument name ''%s''', name))) || ...
     (strcmp(ME.identifier, 'MATLAB:InputParser:UnmatchedParameter') && ...
      contains(ME.message, sprintf('''%s''', name)));
end

function s = plural_y(k)
if k == 1; s = 'y'; else; s = 'ies'; end
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
