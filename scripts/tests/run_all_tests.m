function results = run_all_tests(selection, varargin)
% RUN_ALL_TESTS  Run the test_*.m files in scripts/tests and summarise them.
%
%   results = run_all_tests()            % every test (slow ones included)
%   results = run_all_tests('quick')     % every test not on the slow list
%   results = run_all_tests('slow')      % only the slow list
%   results = run_all_tests('test_psa*') % wildcard on the test name
%   results = run_all_tests('jacobian')  % no wildcard: regexp on the name
%   results = run_all_tests({'test_c_over_K', 'test_route_scale'})
%
%   Name-value options:
%     'show_output'  'failures' (default) | 'always' | 'never'
%                    whether a test's captured command-window text is
%                    printed after its status line
%     'keep_figures' false (default) | true
%                    figures a test opened are closed after it finishes;
%                    figures that existed before the run are never touched
%                    by this function (a test that itself calls
%                    `close all` still closes them)
%
% The tests come in two styles and this copes with both. Some `assert`
% and are silent until they throw; others print PASS / FAIL per check with
% a final banner and never throw. So a test's status is
%
%   ERROR  it threw
%   FAIL   it returned, but its output contains a failure marker
%          (the word FAIL, FAILED or FAILURES; see FAIL_PATTERN below)
%   pass   neither
%
% Each test runs in its own function workspace through evalc, so a test
% that starts with `clear` or `close all; clear; clc` cannot wipe the
% runner's state, and an error in one test does not stop the rest. The
% working directory, the warning state and the global random stream are
% restored after every test.
%
% RESULTS is a table with one row per test: name, status, seconds, slow,
% message (first line of the error, or the first failing line) and output
% (everything the test printed).
%
% Assumes setup_paths has run. Contains no path code, like the tests.
%
% SLOW TESTS. The list in slow_test_names() below is by measurement
% (2026-11-10): anything that took 30 s or more. 'quick' is the rest.
%
% Two things this cannot protect against, because they are what a test
% script itself does: a test that calls `close all` closes the user's
% figures too (test_pairs_single_celltype and test_sensitivity_refactor
% still do, after plotting), and a test that calls `clc` clears the command
% window, including this function's progress lines so far (the returned
% table and the final summary are unaffected).

    if nargin < 1 || isempty(selection); selection = 'all'; end
    p = inputParser;
    p.addParameter('show_output', 'failures', ...
        @(s) any(strcmp(s, {'failures', 'always', 'never'})));
    p.addParameter('keep_figures', false, @(x) islogical(x) && isscalar(x));
    p.parse(varargin{:});
    opts = p.Results;

    test_dir = fileparts(mfilename('fullpath'));
    listing = dir(fullfile(test_dir, 'test_*.m'));
    all_names = sort(erase({listing.name}, '.m'));
    slow_names = slow_test_names();
    names = select_tests(all_names, slow_names, selection);
    if isempty(names)
        error('run_all_tests:NoTests', 'No test matches the selection.');
    end

    n = numel(names);
    status = strings(n, 1);
    seconds = zeros(n, 1);
    message = strings(n, 1);
    output = strings(n, 1);
    is_slow = ismember(names(:), slow_names);

    fprintf('Running %d test(s) from %s\n\n', n, test_dir);
    t_all = tic;
    for k = 1:n
        fprintf('[%2d/%2d] %-45s ', k, n, names{k});
        [status(k), seconds(k), message(k), output(k)] = ...
            run_one_test(names{k}, opts.keep_figures);
        fprintf('%-5s %7.1f s\n', status(k), seconds(k));
        show = strcmp(opts.show_output, 'always') || ...
            (strcmp(opts.show_output, 'failures') && status(k) ~= "pass");
        if show
            fprintf('%s\n', indent_text(output(k)));
        end
    end
    total_seconds = toc(t_all);

    results = table(string(names(:)), status, seconds, is_slow, message, output, ...
        'VariableNames', {'name', 'status', 'seconds', 'slow', 'message', 'output'});

    fprintf('\n%-45s %-6s %9s\n', 'test', 'status', 'seconds');
    fprintf('%s\n', repmat('-', 1, 62));
    for k = 1:n
        tag = '';
        if is_slow(k); tag = '  (slow)'; end
        fprintf('%-45s %-6s %9.1f%s\n', names{k}, status(k), seconds(k), tag);
        if status(k) ~= "pass" && strlength(message(k)) > 0
            fprintf('      %s\n', message(k));
        end
    end
    fprintf('%s\n', repmat('-', 1, 62));
    n_pass = sum(status == "pass");
    n_fail = sum(status == "FAIL");
    n_err = sum(status == "ERROR");
    fprintf('%d test(s): %d pass, %d FAIL, %d ERROR, %.1f s total\n', ...
        n, n_pass, n_fail, n_err, total_seconds);
    if n_fail + n_err == 0
        fprintf('ALL TESTS PASSED\n');
    else
        fprintf('NOT ALL TESTS PASSED\n');
    end
end

% ------------------------------------------------------------------------
function names = slow_test_names()
% Tests left out of run_all_tests('quick'): the ones that took 30 s or more
% when every test was timed on 2026-11-10 (seconds in the comments; 14-core
% desktop; the stage and sweep tests start a parallel pool if none is open,
% which adds its start-up time once). Everything else took under 30 s each.
% Re-time and edit this list when a test is added or changes weight.
    names = {
        'test_regression_bundle'            % ~50 min  the FULL regression: every stage and figure (32 + 18 min)
        'test_numerics_verification_stage'  % 737 s  the numerics stage at 'fast', five checks x 3 conditions
        'test_regression_stages'            % ~3 min  the MIDDLE regression: three stages and their figures
        'test_lyapunov_topk'                % 165 s  top-K against the ode45 full-spectrum QR reference
        'test_sensitivity_refactor'         %  61 s  two sensitivity sweeps on a parallel pool
        'test_benettin_vs_qr'               %  50 s  SRNNModel2, n = 300, both methods
        'test_RMTBlocks_predictions'        %  37 s  eigenvalue statistics over many draws
        'test_eig_heatmap_stage'            %  33 s  analysis stage + three figures
        };
end

% ------------------------------------------------------------------------
function names = select_tests(all_names, slow_names, selection)
    if iscell(selection) || (isstring(selection) && numel(selection) > 1)
        wanted = erase(cellstr(selection), '.m');
        missing = setdiff(wanted, all_names);
        if ~isempty(missing)
            error('run_all_tests:UnknownTest', 'No such test: %s', ...
                strjoin(missing, ', '));
        end
        names = wanted(:)';
        return
    end
    selection = erase(char(selection), '.m');
    switch selection
        case 'all'
            names = all_names;
        case 'quick'
            names = setdiff(all_names, slow_names, 'stable');
        case 'slow'
            names = intersect(all_names, slow_names, 'stable');
        otherwise
            if contains(selection, '*') || contains(selection, '?')
                pattern = regexptranslate('wildcard', selection);
                hit = ~cellfun(@isempty, regexp(all_names, ['^' pattern '$'], 'once'));
            else
                hit = ~cellfun(@isempty, regexp(all_names, selection, 'once'));
            end
            names = all_names(hit);
    end
end

% ------------------------------------------------------------------------
function [status, seconds, message, output] = run_one_test(name, keep_figures)
    ERROR_MARKER = '##run_all_tests:error##';
    % A printed failure: the bare word FAIL ("x: FAIL", "[FAIL] x",
    % "FAIL  x"), or FAILED / FAILURES in a banner. Case-sensitive, so
    % prose such as "failed to converge" does not count. "0 FAILED" and
    % "FAILED: 0" style counts are removed first.
    FAIL_PATTERN = '(?<![A-Za-z_])(FAIL|FAILED|FAILURES)(?![A-Za-z_])';
    ZERO_COUNT = '(\<0 (checks? )?(FAIL|FAILED|FAILURES)\>)|((FAIL|FAILED|FAILURES)\s*[:=]\s*0\>)';

    figs_before = findall(groot, 'Type', 'figure');
    cwd_before = pwd;
    warn_before = warning;
    rng_before = rng;

    t0 = tic;
    try
        output = isolated_evalc(name, ERROR_MARKER);
        runner_error = '';
    catch ME      % only if evalc itself could not run the command
        output = '';
        runner_error = ME.message;
    end
    seconds = toc(t0);

    cd(cwd_before);
    warning(warn_before);
    rng(rng_before);
    if ~keep_figures
        figs_new = setdiff(findall(groot, 'Type', 'figure'), figs_before);
        delete(figs_new(isvalid(figs_new)));
    end

    message = "";
    marker_at = strfind(output, ERROR_MARKER);
    if ~isempty(runner_error)
        status = "ERROR";
        message = string(first_line(runner_error));
    elseif ~isempty(marker_at)
        status = "ERROR";
        report = output(marker_at(end) + numel(ERROR_MARKER):end);
        message = string(first_line(report));
        output = strrep(output, ERROR_MARKER, 'ERROR: ');
    else
        scrubbed = regexprep(output, ZERO_COUNT, '');
        lines = regexp(scrubbed, '\r?\n', 'split');
        hit = ~cellfun(@isempty, regexp(lines, FAIL_PATTERN, 'once'));
        if any(hit)
            status = "FAIL";
            message = string(strtrim(lines{find(hit, 1)}));
        else
            status = "pass";
        end
    end
    output = string(output);
end

% ------------------------------------------------------------------------
function run_all_tests_output__ = isolated_evalc(run_all_tests_name__, run_all_tests_marker__)
% The only workspace a test script can see. A script that calls `clear`
% wipes the two inputs, which is harmless: the command string has already
% been handed to evalc, and the error is reported through the captured
% text rather than through a variable.
    run_all_tests_output__ = evalc(sprintf([ ...
        'try, %s; catch run_all_tests_ME__, ' ...
        'fprintf(''\\n%s%%s\\n'', getReport(run_all_tests_ME__, ''extended'', ''hyperlinks'', ''off'')); end'], ...
        run_all_tests_name__, run_all_tests_marker__));
end

% ------------------------------------------------------------------------
function s = first_line(txt)
    parts = regexp(strtrim(char(txt)), '\r?\n', 'split');
    parts = parts(~cellfun(@(c) isempty(strtrim(c)), parts));
    if isempty(parts)
        s = '';
    elseif numel(parts) > 1 && startsWith(strtrim(parts{1}), 'Error using')
        s = strtrim(parts{2});      % the message, not "Error using f (line n)"
    else
        s = strtrim(parts{1});
    end
    if numel(s) > 200; s = [s(1:197) '...']; end
end

% ------------------------------------------------------------------------
function s = indent_text(txt)
    s = regexprep(char(txt), '(^|\n)', '$1      | ');
end
