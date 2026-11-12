%% TEST_SELECT_CONDITIONS Every rule in select_conditions' header (FD-006).
%
% Empty request = all, unchanged; requested order wins and idx indexes
% `available`; a name the run lacks raises NoSuchCondition listing what is
% available; a duplicate raises DuplicateCondition; a char or a string names
% one condition; malformed input raises BadInput. Pure, no data, instant.
%
% The end-to-end check of the option (a figure and make_all_paper_figures) is
% test_select_conditions_figure.

n_fail = 0;
avail = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};

%% 1. empty request = everything, unchanged
empties = {{}, '', "", [], strings(0), cell(0, 1)};
for k = 1:numel(empties)
    [nm, ix] = select_conditions(avail, empties{k});
    n_fail = n_fail + check(isequal(nm, avail) && isequal(ix, 1:4), ...
        sprintf('empty request #%d (%s) returns available unchanged', k, class(empties{k})));
end
[nm, ix] = select_conditions({}, {});
n_fail = n_fail + check(isempty(nm) && isempty(ix), 'nothing available, nothing requested');

%% 2. requested order wins
req = {'sfa3_stdc2', 'no_adaptation', 'sfa1_std1'};
[nm, ix] = select_conditions(avail, req);
n_fail = n_fail + check(isequal(nm, req) && isequal(ix, [3 1 2]) && isequal(avail(ix), nm), ...
    'result follows the requested order; idx indexes available');
[nm, ix] = select_conditions(avail, {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'});
n_fail = n_fail + check(isequal(nm, {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}) && isequal(ix, 1:3), ...
    'the paper''s three-condition set');
[nm, ix] = select_conditions(avail, avail);
n_fail = n_fail + check(isequal(nm, avail) && isequal(ix, 1:4), 'the four-condition set');
[nm, ix] = select_conditions(avail', req');   % column inputs
n_fail = n_fail + check(isequal(nm, req) && isequal(ix, [3 1 2]) && isrow(nm) && isrow(ix), ...
    'column inputs give row outputs');

%% 3. a name the run lacks is an error that lists what is available
[id, msg] = error_of(@() select_conditions({'no_adaptation', 'sfa1_std1', 'sfa3_std2'}, ...
    {'no_adaptation', 'sfa3_stdc2'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:NoSuchCondition') && ...
    contains(msg, 'sfa3_stdc2') && contains(msg, 'no_adaptation, sfa1_std1, sfa3_std2'), ...
    'missing name -> NoSuchCondition, naming it and listing the available ones');
id = error_of(@() select_conditions(avail, {'SFA1_STD1'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:NoSuchCondition'), 'matching is case-sensitive');
id = error_of(@() select_conditions(avail, {''}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:NoSuchCondition'), ...
    'an empty NAME inside a list is not "all": NoSuchCondition');

%% 4. a duplicate in requested is an error
[id, msg] = error_of(@() select_conditions(avail, {'sfa1_std1', 'no_adaptation', 'sfa1_std1'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:DuplicateCondition') && contains(msg, 'sfa1_std1'), ...
    'duplicate request -> DuplicateCondition');
id = error_of(@() select_conditions(avail, {'nope', 'nope'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:DuplicateCondition'), ...
    'duplicate is reported before a missing name');

%% 5. a char or a string names one condition
[nm, ix] = select_conditions(avail, 'sfa3_std2');
n_fail = n_fail + check(isequal(nm, {'sfa3_std2'}) && ix == 4, 'char names one condition');
[nm, ix] = select_conditions(avail, "sfa1_std1");
n_fail = n_fail + check(isequal(nm, {'sfa1_std1'}) && ix == 2, 'string names one condition');
[nm, ix] = select_conditions(avail, ["sfa3_std2", "no_adaptation"]);
n_fail = n_fail + check(isequal(nm, {'sfa3_std2', 'no_adaptation'}) && isequal(ix, [4 1]), ...
    'string array is a list');
[nm, ix] = select_conditions(avail, {"sfa3_std2", 'no_adaptation'});
n_fail = n_fail + check(isequal(nm, {'sfa3_std2', 'no_adaptation'}) && isequal(ix, [4 1]), ...
    'cell mixing string and char');
[nm, ix] = select_conditions('no_adaptation', {});
n_fail = n_fail + check(isequal(nm, {'no_adaptation'}) && ix == 1, 'char available is one condition');
[nm, ix] = select_conditions(string(avail), {'sfa1_std1'});
n_fail = n_fail + check(isequal(nm, {'sfa1_std1'}) && ix == 2 && iscell(nm) && all(cellfun(@ischar, nm)), ...
    'string-array available; names come back as cellstr');

%% Malformed input
id = error_of(@() select_conditions(avail, {'sfa1_std1', 3}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:BadInput'), 'non-text request entry -> BadInput');
id = error_of(@() select_conditions({'a', 'a'}, {}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:BadInput'), 'duplicate in available -> BadInput');
id = error_of(@() select_conditions(avail, 7));
n_fail = n_fail + check(strcmp(id, 'select_conditions:BadInput'), 'numeric request -> BadInput');

%% Verdict
if n_fail == 0
    fprintf('\nALL select_conditions CHECKS PASSED\n');
else
    error('test_select_conditions:Failed', '%d check(s) FAILED', n_fail);
end

%% ------------------------------------------------------------------------
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
