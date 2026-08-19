% test_srnn_param_preset_equivalence.m
% Prove the FLATTENED srnn_param_preset still returns exactly what the old
% CHAINED implementation did, for the presets that existed at the flattening.
%
% srnn_param_preset used to build ten of its fourteen presets by calling itself
% recursively, up to seven hops deep. Flattening it into self-contained cases is
% a pure readability change, so the only thing that must be established is that
% nothing moved. scripts/tests/srnn_param_preset_old.m is the frozen chained
% implementation; this script compares the two output-for-output.
%
% SCOPE -- READ THIS BEFORE ADDING A PRESET. Only the K presets the frozen copy
% knows about are compared. A preset added to srnn_param_preset.m AFTER the
% freeze has no old counterpart to compare against, so it is skipped and NAMED
% in the output rather than treated as a failure. That is what keeps this test
% from taxing every future preset. The requirement it does impose is ordering:
% the frozen K must remain the FIRST K, in their original order, so APPEND NEW
% PRESETS AT THE BOTTOM. Reordering, renaming or deleting one of the original K
% fails, by design -- those are exactly the changes worth catching.
%
% The comparison is only meaningful if the frozen copy is genuinely independent,
% so the first check greps it for calls to the live function. Without that guard
% a mis-frozen copy would delegate to the very code under test and this script
% would pass unconditionally.
%
% Run with the matlab MCP (run_matlab_file), not matlab -batch.
%
% See also: srnn_param_preset, srnn_param_preset_old, test_srnn_param_preset,
%           srnn_adaptation_conditions

fprintf('=== srnn_param_preset: flattened vs. frozen chained implementation ===\n\n');
all_passed = true;

%% The frozen copy must not delegate to the live function
old_path = which('srnn_param_preset_old');
all_passed = check('srnn_param_preset_old is on the path', ...
    ~isempty(old_path)) && all_passed;

% `srnn_param_preset\s*\(` matches a call to the LIVE function but not to
% srnn_param_preset_old( or srnn_param_preset_names( -- in both of those an
% underscore follows the name, not a paren.
%
% Comment lines are dropped first. The frozen file's own header explains this
% very guard and so contains the pattern as prose; a mention in a comment is
% not a delegation, and matching it would make the guard cry wolf.
code_lines = strsplit(fileread(old_path), newline);
is_comment = ~cellfun(@isempty, regexp(code_lines, '^\s*%', 'once'));
code_only = strjoin(code_lines(~is_comment), newline);
delegations = regexp(code_only, 'srnn_param_preset\s*\(', 'once');
all_passed = check('frozen copy never calls the live srnn_param_preset', ...
    isempty(delegations)) && all_passed;
if ~isempty(delegations)
    fprintf(['      The frozen copy delegates to the function under test, so ' ...
        'every comparison below is vacuous.\n']);
end

% Positive control for the guard itself. Now that comments are stripped and the
% name has a look-alike sibling, a guard that silently stopped matching would
% look identical to a guard that passed -- so check it still fires on a line
% that IS a delegation, and still ignores the two legitimate call shapes.
bites = ~isempty(regexp('    [d, mc] = srnn_param_preset(''default'');', ...
    'srnn_param_preset\s*\(', 'once'));
spares_old = isempty(regexp('    [d, mc] = srnn_param_preset_old(''default'');', ...
    'srnn_param_preset\s*\(', 'once'));
spares_names = isempty(regexp('    n = srnn_param_preset_names();', ...
    'srnn_param_preset\s*\(', 'once'));
all_passed = check('the delegation guard fires on a real delegation', ...
    bites && spares_old && spares_names) && all_passed;

%% Scope: the K presets the frozen copy knows about, and only those
% K is whatever srnn_param_preset_old enumerates -- not a hardcoded number, so
% this needs no maintenance. The frozen list must still be a PREFIX of the live
% one: presets appended to the bottom of srnn_param_preset.m are new work with
% no old counterpart and are skipped, while REORDERING, RENAMING or DELETING
% one of the original K is a failure, because that is a change to something this
% test is supposed to be guarding.
%
% Both lists are recovered from the UnknownPreset message rather than kept as a
% third list here, so nothing has to be updated in this file when a preset is
% added. (Same idiom as test_srnn_param_preset.m.)
%
% ==> APPEND NEW PRESETS AT THE BOTTOM of the switch in srnn_param_preset.m and
%     at the end of its srnn_param_preset_names list. Inserting one in the
%     middle breaks the prefix and fails here.
names_old = preset_names_from_error(@srnn_param_preset_old);
names_new = preset_names_from_error(@srnn_param_preset);
K = numel(names_old);

enough = numel(names_new) >= K;
all_passed = check('live function still defines all K frozen presets', enough) && all_passed;
if ~enough
    fprintf('      frozen has %d presets, live has %d -- one was deleted.\n', ...
        K, numel(names_new));
    fprintf('      missing: %s\n', strjoin(setdiff(names_old, names_new), ', '));
end

prefix_ok = enough && isequal(names_new(1:K), names_old);
all_passed = check('the K frozen presets are still the first K, in order', ...
    prefix_ok) && all_passed;
if enough && ~prefix_ok
    for i = 1:K
        if ~strcmp(names_old{i}, names_new{i})
            fprintf('      position %d: frozen ''%s'', live ''%s''\n', ...
                i, names_old{i}, names_new{i});
        end
    end
    fprintf(['      Append new presets at the END of the list; do not insert ' ...
        'or reorder.\n']);
end

added = names_new(K+1:end);
if isempty(added)
    fprintf('\n  %d presets under comparison; none added since the freeze.\n\n', K);
else
    % Not a failure -- just say plainly what this test does NOT cover, so a
    % green banner is never read as "every preset is verified".
    fprintf('\n  %d presets under comparison. NOT COVERED (added since the freeze): %s\n', ...
        K, strjoin(added, ', '));
    fprintf(['  Those have no frozen counterpart; test_srnn_param_preset.m is ' ...
        'what validates them.\n\n']);
end

% Only the frozen K are compared below.
names_old = names_old(1:K);

%% Every preset, all three outputs
for i = 1:numel(names_old)
    nm = names_old{i};
    [d_old, mc_old, cond_old] = srnn_param_preset_old(nm);
    [d_new, mc_new, cond_new] = srnn_param_preset(nm);

    % --- model_class ----------------------------------------------------
    ok = strcmp(mc_old, mc_new);
    all_passed = check(sprintf('%s: model_class', nm), ok) && all_passed;
    if ~ok
        fprintf('      old = %s, new = %s\n', mc_old, mc_new);
    end

    % --- the override struct --------------------------------------------
    % isequaln, not isequal: mu_S_c / sigma_S_c are [] in six presets and
    % NaN would compare unequal to itself under isequal.
    ok = isequaln(d_old, d_new);
    all_passed = check(sprintf('%s: model_defaults struct', nm), ok) && all_passed;
    if ~ok
        report_struct_diff(d_old, d_new);
    end

    % Field ORDER is not required by isequaln (MATLAB compares structs
    % order-independently) and nothing downstream depends on it, but a
    % divergence means the flattened case was written in a different order
    % than the chain produced -- worth knowing while reviewing the diff.
    ok = isequal(fieldnames(d_old), fieldnames(d_new));
    all_passed = check(sprintf('%s: field order', nm), ok) && all_passed;
    if ~ok
        fprintf('      old: %s\n', strjoin(fieldnames(d_old)', ', '));
        fprintf('      new: %s\n', strjoin(fieldnames(d_new)', ', '));
    end

    % --- the adaptation conditions --------------------------------------
    % A 1x4 cell of plain structs (no function handles anywhere), so isequaln
    % is a complete comparison and no bespoke recursive comparator is needed.
    % This is what catches a std_routes that was mis-transcribed: the routes
    % only reach the model through these conditions.
    ok = isequaln(cond_old, cond_new);
    all_passed = check(sprintf('%s: adaptation conditions', nm), ok) && all_passed;
    if ~ok
        report_conditions_diff(cond_old, cond_new);
    end
end

%% The 1- and 2-output call forms
% Five live callers use these (run_*_analysis.m take one output,
% run_overnight_queue.m and the two Fig_sensitivity_analysis_allStd scripts take
% two), and nargout is visible inside the function, so they are worth pinning
% separately rather than assuming they fall out of the 3-output check.
one_ok = true; two_ok = true;
for i = 1:numel(names_old)
    nm = names_old{i};
    d1_old = srnn_param_preset_old(nm);
    d1_new = srnn_param_preset(nm);
    one_ok = one_ok && isequaln(d1_old, d1_new);

    [d2_old, mc2_old] = srnn_param_preset_old(nm);
    [d2_new, mc2_new] = srnn_param_preset(nm);
    two_ok = two_ok && isequaln(d2_old, d2_new) && strcmp(mc2_old, mc2_new);
end
all_passed = check('1-output call form agrees for every preset', one_ok) && all_passed;
all_passed = check('2-output call form agrees for every preset', two_ok) && all_passed;

%% The unknown-preset error
[threw_old, err_old] = capture_error(@() srnn_param_preset_old('no_such_preset'));
[threw_new, err_new] = capture_error(@() srnn_param_preset('no_such_preset'));
all_passed = check('both reject an unknown preset', ...
    threw_old && threw_new) && all_passed;
all_passed = check('unknown-preset identifier unchanged', ...
    threw_old && threw_new && strcmp(err_old.identifier, err_new.identifier)) && all_passed;

% Compare the message STRUCTURE, not the whole string. The message ends with
% the full preset list, which legitimately grows every time a preset is
% appended -- comparing it verbatim would make adding a preset a failure, which
% is exactly what this test is meant not to do. So: the part before the list
% must match exactly, and the frozen names must still all be in the live list.
marker = 'Valid presets:';
ok = threw_old && threw_new && ...
    strcmp(extractBefore(err_old.message, marker), ...
           extractBefore(err_new.message, marker));
all_passed = check('unknown-preset message preamble unchanged', ok) && all_passed;
if ~ok && threw_old && threw_new
    fprintf('      old: %s\n', err_old.message);
    fprintf('      new: %s\n', err_new.message);
end

missing = setdiff(names_old, names_new, 'stable');
all_passed = check('unknown-preset message still lists every frozen preset', ...
    isempty(missing)) && all_passed;
if ~isempty(missing)
    fprintf('      absent from the live message: %s\n', strjoin(missing, ', '));
end

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function names = preset_names_from_error(fcn)
% Every preset FCN knows about, read out of its UnknownPreset message.
%
% srnn_param_preset_names is local to each file and cannot be called from here,
% but the error enumerates it, and that message format is itself asserted below.
try
    fcn('__definitely_not_a_preset__');
    error('test_srnn_param_preset_equivalence:NoError', ...
        'A bogus preset name was accepted; the list cannot be recovered.');
catch err
    if ~strcmp(err.identifier, 'srnn_param_preset:UnknownPreset')
        rethrow(err);
    end
end
marker = 'Valid presets:';
tail = extractAfter(err.message, marker);
names = strtrim(split(strtrim(erase(tail, '.')), ','))';
end

function report_struct_diff(a, b)
% Name the fields responsible, so a failure says what to fix.
fa = fieldnames(a);
fb = fieldnames(b);
only_a = setdiff(fa, fb);
only_b = setdiff(fb, fa);
if ~isempty(only_a); fprintf('      only in old: %s\n', strjoin(only_a', ', ')); end
if ~isempty(only_b); fprintf('      only in new: %s\n', strjoin(only_b', ', ')); end
for k = 1:numel(fa)
    f = fa{k};
    if ~isfield(b, f); continue; end
    if ~isequaln(a.(f), b.(f))
        fprintf('      %s differs:\n', f);
        fprintf('        old: %s\n', compact_value(a.(f)));
        fprintf('        new: %s\n', compact_value(b.(f)));
    end
end
end

function report_conditions_diff(a, b)
% The conditions are a 1x4 cell keyed by .name; report per condition.
if numel(a) ~= numel(b)
    fprintf('      condition count: old %d, new %d\n', numel(a), numel(b));
    return;
end
for k = 1:numel(a)
    if ~isequaln(a{k}, b{k})
        fprintf('      condition ''%s'' differs\n', a{k}.name);
        if isfield(a{k}, 'synapse_config') && isfield(b{k}, 'synapse_config')
            fprintf('        (compare synapse_config / std_routes for this preset)\n');
        end
    end
end
end

function s = compact_value(v)
% One-line rendering good enough to spot a transcription slip.
if isempty(v)
    s = '[]';
elseif ischar(v)
    s = ['''' v ''''];
elseif isnumeric(v) || islogical(v)
    s = ['[' strjoin(compose('%g', double(v(:)))', ' ') ']'];
elseif iscell(v)
    s = ['{' strjoin(cellfun(@compact_value, v(:)', 'UniformOutput', false), ' ') '}'];
elseif isstruct(v)
    s = ['struct(' strjoin(fieldnames(v)', ', ') ')'];
else
    s = class(v);
end
end

function [threw, err] = capture_error(fcn)
threw = false;
err = [];
try
    fcn();
catch err
    threw = true;
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
