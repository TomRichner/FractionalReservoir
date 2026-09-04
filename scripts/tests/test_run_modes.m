% TEST_RUN_MODES One list of run modes, and every stage validates against it.
%
% THE BUG THIS GUARDS. analysis_run_config accepted five run modes; the four
% later stages -- run_memory_capacity, run_memory_capacity_example,
% run_eig_heatmap, run_dc_lle_analysis -- each carried their own switch accepting
% three. A run_mode valid for the sweeps was therefore invalid downstream, and
% because run_all_paper_analyses wraps each stage, that surfaced as stages
% quietly producing nothing AFTER the sweeps had already run.
%
% On 2026-09-03 single_multi_TS ran at 'fast2': sweeps completed, three stages
% threw "Unknown run_mode", two figures failed on missing data, and
% fig_memory_capacity fell back to a twelve-day-old .mat and reported SUCCESS.
%
% WHAT THIS COVERS: that run_mode_names() is the one list, that
% analysis_run_config accepts exactly it, and that every later stage REJECTS an
% unknown mode while naming that same list in the message.
%
% WHAT IT DOES NOT COVER, and why. It cannot cheaply prove that each stage
% IMPLEMENTS every canonical mode, because the cost tables are local functions
% and reaching one means running the analysis -- minutes at 'fast', hours at
% 'production'. A missing case would fall through to `otherwise` and throw
% badMode at run time. Two things stand in for that: the preflight in
% run_all_paper_analyses fails in microseconds rather than after the sweeps, and
% every stage's `otherwise` branch now points here. Making the cost tables
% callable in their own right is the fix if this gap ever bites.
%
% Assumes setup_paths has run.
%
% See also: run_mode_names, analysis_run_config, run_all_paper_analyses

all_passed = true;
modes = run_mode_names();

fprintf('\n-- the canonical list --\n');
all_passed = check(sprintf('run_mode_names = {%s}', strjoin(modes, ', ')), ...
    iscellstr(modes) && ~isempty(modes)) && all_passed; %#ok<ISCLSTR>
all_passed = check('''fast2'' is gone (removed 2026-09-03)', ...
    ~ismember('fast2', modes)) && all_passed;

%% analysis_run_config accepts exactly the canonical list
fprintf('\n-- analysis_run_config accepts exactly it --\n');
analyses = {'sensitivity', 'tau_sensitivity', 'param_space'};
ok = true; bad = {};
for i = 1:numel(analyses)
    for j = 1:numel(modes)
        try
            c = analysis_run_config(analyses{i}, modes{j});
            if ~isscalar(c.n_levels) || c.n_levels <= 0
                ok = false; bad{end+1} = sprintf('%s/%s: bad n_levels', analyses{i}, modes{j}); %#ok<SAGROW>
            end
        catch ME
            ok = false;
            bad{end+1} = sprintf('%s/%s: %s', analyses{i}, modes{j}, ME.identifier); %#ok<SAGROW>
        end
    end
end
all_passed = check(sprintf('%d analyses x %d modes all resolve', ...
    numel(analyses), numel(modes)), ok) && all_passed;
if ~ok; fprintf(2, '     %s\n', strjoin(bad, '\n     ')); end

all_passed = check('a removed mode is rejected', ...
    throws_id(@() analysis_run_config('sensitivity', 'fast2'), ...
              'analysis_run_config:badMode')) && all_passed;

%% Every later stage validates against the SAME list
% Each is called with an unknown mode; the error must name every canonical mode,
% which is what proves the stage reads run_mode_names rather than a private copy.
fprintf('\n-- later stages name the canonical list --\n');
stages = { ...
    'run_memory_capacity',         @() run_memory_capacity('run_mode', 'no_such_mode'); ...
    'run_memory_capacity_example', @() run_memory_capacity_example('run_mode', 'no_such_mode'); ...
    'run_eig_heatmap',             @() run_eig_heatmap('run_mode', 'no_such_mode'); ...
    'run_dc_lle_analysis',         @() run_dc_lle_analysis('run_mode', 'no_such_mode')};
for k = 1:size(stages, 1)
    [threw, msg, id] = capture(stages{k, 2});
    names_all = threw && all(cellfun(@(m) contains(msg, m), modes));
    is_bad_mode = threw && endsWith(id, ':badMode');
    all_passed = check(sprintf('%-30s rejects and lists all %d modes', ...
        stages{k, 1}, numel(modes)), names_all && is_bad_mode) && all_passed;
    if ~(names_all && is_bad_mode)
        fprintf(2, '     id=%s msg=%s\n', id, msg);
    end
end

%% The preflight fires before any compute
fprintf('\n-- run_all_paper_analyses preflight --\n');
bad_cfg = paper_config();
bad_cfg.run_mode = 'fast2';
all_passed = check('a removed mode is refused up front', ...
    throws_id(@() run_all_paper_analyses(bad_cfg), ...
              'run_all_paper_analyses:BadRunMode')) && all_passed;
% It must fire BEFORE the run directory is touched, or the check is worthless:
% a refused run that has already created output is not a refused run.
bad_cfg2 = paper_config('run_dir', 'data/__preflight_probe__');
bad_cfg2.run_mode = 'fast2';
[~] = capture(@() run_all_paper_analyses(bad_cfg2));
probe = fullfile(fileparts(which('setup_paths')), 'data', '__preflight_probe__');
all_passed = check('...and creates no output directory', ~isfolder(probe)) && all_passed;
if isfolder(probe); rmdir(probe, 's'); end

%% ------------------------------------------------------------------------
if all_passed
    fprintf('\ntest_run_modes: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_run_modes: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function tf = throws_id(fn, want)
try
    fn();
    tf = false;
catch ME
    tf = strcmp(ME.identifier, want);
end
end

function [threw, msg, id] = capture(fn)
try
    fn();
    threw = false; msg = ''; id = '';
catch ME
    threw = true; msg = ME.message; id = ME.identifier;
end
end
