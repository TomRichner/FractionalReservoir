function test_preset_golden()
% TEST_PRESET_GOLDEN srnn_param_preset's outputs are unchanged.
%
% Compares every output of srnn_param_preset against a fixture captured from
% pre-refactor code (see make_preset_golden), for all 10 presets and all 42
% conditions:
%
%   d             exact isequaln
%   model_class   exact
%   conditions    exact isequaln -- NO exclusions
%   errors        the identifiers for 13 retired names and one unknown name
%
% WHY EXACT, WITH NOTHING EXCLUDED. The refactor that this guards -- moving
% condition definitions out of srnn_adaptation_conditions and into each preset --
% is supposed to change NOTHING about what the function returns. An exclusion is
% a hole in that claim, so there are none. In particular conditions carry no
% `title` field: titles stay in srnn_condition_titles, keyed by name, which is
% what lets a saved run be retitled without recomputing it.
%
% Two properties make this sound:
%   * srnn_param_preset is PURE -- one char argument, no file reads, no RNG. The
%     generator re-checks this before capturing.
%   * isequaln on structs is field-order-INSENSITIVE, so the refactor is free to
%     write condition fields in whatever order reads best.
%
% IF THIS FAILS, the refactor changed behaviour. Fix the code, not the fixture.
% Regenerating a golden file to make a test pass leaves something that still
% reads like evidence but no longer is. Genuine intended changes are made by
% deleting the .mat, re-running make_preset_golden, and saying in the commit
% message why the expected values moved.
%
% See also: make_preset_golden, srnn_param_preset, test_c_over_K

setup_paths();
all_passed = true;

fixture = fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
    'golden_preset_outputs.mat');
if ~isfile(fixture)
    error('test_preset_golden:NoFixture', ...
        ['Fixture not found:\n  %s\n\n' ...
         'It is gitignored as a *.mat and must be force-added. Regenerate with ' ...
         'make_preset_golden -- but NOTE that a fixture captured from current ' ...
         'code proves nothing about a refactor already applied.'], fixture);
end
L = load(fixture);
G = L.G;

fprintf('\n=== preset outputs vs fixture captured %s ===\n', G.captured.when);

%% Values
fprintf('\n-- outputs unchanged --\n');
n_cond_checked = 0;
for k = 1:numel(G.presets)
    want = G.presets(k);
    try
        [d, model_class, conditions] = srnn_param_preset(want.name);
    catch ME
        all_passed = check(sprintf('%-46s', want.name), false) && all_passed;
        fprintf(2, '      threw %s: %s\n', ME.identifier, ME.message);
        continue
    end

    ok_d   = isequaln(d, want.d);
    ok_cls = isequaln(model_class, want.model_class);
    ok_c   = isequaln(conditions, want.conditions);
    ok     = ok_d && ok_cls && ok_c;
    n_cond_checked = n_cond_checked + numel(want.conditions);

    all_passed = check(sprintf('%-46s %d cond', want.name, ...
        numel(want.conditions)), ok) && all_passed;

    if ~ok
        if ~ok_cls
            fprintf(2, '      model_class: got ''%s'' want ''%s''\n', ...
                model_class, want.model_class);
        end
        if ~ok_d;  report_struct_diff('d', d, want.d);  end
        if ~ok_c;  report_condition_diff(conditions, want.conditions);  end
    end
end

%% Error behaviour
fprintf('\n-- retired presets still error --\n');
bad_retired = {};
for k = 1:numel(G.retired)
    got = error_id(@() srnn_param_preset(G.retired(k).name));
    if ~strcmp(got, G.retired(k).identifier)
        bad_retired{end+1} = sprintf('%s (got ''%s'')', ...
            G.retired(k).name, got); %#ok<AGROW>
    end
end
all_passed = check(sprintf('%d retired names -> %s', numel(G.retired), ...
    G.retired(1).identifier), isempty(bad_retired)) && all_passed;
if ~isempty(bad_retired)
    fprintf(2, '      %s\n', strjoin(bad_retired, '\n      '));
end

got_unknown = error_id(@() srnn_param_preset(G.unknown.name));
all_passed = check(sprintf('unknown name -> %s', G.unknown.identifier), ...
    strcmp(got_unknown, G.unknown.identifier)) && all_passed;

% The valid-preset list is embedded in the unknown-name message, so freezing the
% message freezes the list -- srnn_param_preset_names() is a subfunction and
% cannot be called from here.
got_msg = error_msg(@() srnn_param_preset(G.unknown.name));
all_passed = check('valid-preset list unchanged', ...
    strcmp(got_msg, G.unknown.message)) && all_passed;
if ~strcmp(got_msg, G.unknown.message)
    fprintf(2, '      got  %s\n      want %s\n', got_msg, G.unknown.message);
end

%% Summary
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL PASSED (%d presets, %d conditions)\n', ...
        numel(G.presets), n_cond_checked);
else
    fprintf(2, 'SOME TESTS FAILED -- fix the code, not the fixture\n');
end
fprintf('========================================\n');
end

%% ------------------------------------------------------------------------
function ok = check(name, tf)
if tf
    fprintf('PASS  %s\n', name);
else
    fprintf(2, 'FAIL  %s\n', name);
end
ok = tf;
end

function report_struct_diff(label, got, want)
% Name the fields that actually differ, rather than printing two structs.
fg = fieldnames(got); fw = fieldnames(want);
only_got  = setdiff(fg, fw);
only_want = setdiff(fw, fg);
if ~isempty(only_got)
    fprintf(2, '      %s: extra field(s) %s\n', label, strjoin(only_got', ', '));
end
if ~isempty(only_want)
    fprintf(2, '      %s: missing field(s) %s\n', label, strjoin(only_want', ', '));
end
for i = 1:numel(fw)
    f = fw{i};
    if isfield(got, f) && ~isequaln(got.(f), want.(f))
        fprintf(2, '      %s.%s differs\n', label, f);
    end
end
end

function report_condition_diff(got, want)
if numel(got) ~= numel(want)
    fprintf(2, '      conditions: got %d, want %d\n', numel(got), numel(want));
    return
end
for i = 1:numel(want)
    if isequaln(got{i}, want{i}); continue; end
    gname = 'unnamed'; if isfield(got{i}, 'name'); gname = got{i}.name; end
    if ~strcmp(gname, want{i}.name)
        fprintf(2, '      condition %d: name got ''%s'' want ''%s''\n', ...
            i, gname, want{i}.name);
    end
    if ~isscalar(got{i})
        % The struct('tau_a', taus) trap: a cell value builds a 1xN struct ARRAY.
        fprintf(2, '      condition %d (%s): got a %dx%d struct ARRAY, not a scalar\n', ...
            i, want{i}.name, size(got{i}, 1), size(got{i}, 2));
        continue
    end
    report_struct_diff(sprintf('condition %d (%s)', i, want{i}.name), ...
        got{i}, want{i});
end
end

function id = error_id(fn)
try
    fn();
    id = '';
catch ME
    id = ME.identifier;
end
end

function msg = error_msg(fn)
try
    fn();
    msg = '';
catch ME
    msg = ME.message;
end
end
