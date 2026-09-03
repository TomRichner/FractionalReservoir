% TEST_LOG_LADDER The log-spaced ladder helper, and the logspace trap it exists for.
%
% log_ladder(lo, hi, n) is generic numerics living in src/util/, so it gets its
% own test rather than riding along in test_adaptation_conditions -- where these
% checks used to live, back when the helper was srnn_sfa_timescales and knew one
% hardcoded ladder. A util function should fail independently of the conditions
% logic that happens to be its biggest caller.
%
% THE TRAP. MATLAB's logspace(a, b, 1) returns 10^b -- the LAST endpoint. Asking
% for one value out of a range silently gives the top of it, and nothing reports
% the substitution. That produced a real defect: single-timescale SFA conditions
% got the SLOW 10 s constant when they wanted the fast 0.25 s one, and the model
% built and ran without complaint.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: log_ladder, srnn_param_preset, srnn_adaptation_conditions

all_passed = true;
fprintf('\n=== test_log_ladder ===\n');

%% ------------------------------------------------------------------------
fprintf('\n-- the trap this exists for --\n');
% Asserted rather than assumed. If MATLAB ever changed this, the n = 1 branch
% would be dead code and the reader should be told rather than left guessing.
all_passed = check('logspace(a,b,1) really does return the LAST endpoint', ...
    abs(logspace(log10(0.25), log10(10), 1) - 10) < 1e-12) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- the three cases --\n');
all_passed = check('n >= 2 matches logspace exactly', ...
    isequal(log_ladder(0.25, 10, 3), logspace(log10(0.25), log10(10), 3)) && ...
    isequal(log_ladder(0.25, 10, 5), logspace(log10(0.25), log10(10), 5))) && all_passed;
all_passed = check('n = 1 returns lo, not hi', ...
    isequal(log_ladder(0.25, 10, 1), 0.25)) && all_passed;
all_passed = check('n = 0 is a 1x0 row, not []', ...
    isequal(size(log_ladder(0.25, 10, 0)), [1 0]) && ...
    isnumeric(log_ladder(0.25, 10, 0))) && all_passed;
all_passed = check('n = 2 is exactly the two endpoints', ...
    isequal(log_ladder(0.25, 10, 2), [0.25 10])) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- n = 1 is POSITIONAL, not by magnitude --\n');
% A rule about argument ORDER stays predictable where a rule about MAGNITUDE
% would silently flip when someone passes a descending range.
all_passed = check('log_ladder(10, 0.25, 1) returns 10, the FIRST argument', ...
    isequal(log_ladder(10, 0.25, 1), 10)) && all_passed;
all_passed = check('a descending range still spans its endpoints', ...
    isequal(log_ladder(10, 0.25, 2), [10 0.25])) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- bad input is refused, not silently coerced --\n');
all_passed = check('lo = 0 errors (log of zero)', ...
    throws(@() log_ladder(0, 10, 3))) && all_passed;
all_passed = check('a negative endpoint errors', ...
    throws(@() log_ladder(-1, 10, 3))) && all_passed;
all_passed = check('Inf endpoint errors', ...
    throws(@() log_ladder(0.25, Inf, 3))) && all_passed;
all_passed = check('negative n errors', ...
    throws(@() log_ladder(0.25, 10, -1))) && all_passed;
all_passed = check('non-integer n errors', ...
    throws(@() log_ladder(0.25, 10, 2.5))) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- the fast end is DERIVED, not restated --\n');
% srnn_sfa_timescales wrote 0.25 twice: once as a literal for K = 1 and once
% inside log10() for the rest, so retuning the fast end would have missed one.
% Here n = 1 must track lo for ANY endpoints.
all_passed = check('n = 1 tracks lo across different ladders', ...
    isequal(log_ladder(0.5, 15, 1), 0.5) && ...
    isequal(log_ladder(0.1, 10, 1), 0.1) && ...
    isequal(log_ladder(2, 200, 1), 2)) && all_passed;
all_passed = check('n = 1 equals the first element of the n = 3 ladder', ...
    isequal(log_ladder(0.1, 10, 1), subsref(log_ladder(0.1, 10, 3), ...
        substruct('()', {1})))) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n');
if all_passed
    fprintf('test_log_ladder: ALL TESTS PASSED\n');
else
    fprintf(2, 'test_log_ladder: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function passed = check(name, condition)
if condition
    fprintf('  PASS  %s\n', name);
    passed = true;
else
    fprintf(2, '  FAIL  %s\n', name);
    passed = false;
end
end

function tf = throws(fn)
try
    fn();
    tf = false;
catch
    tf = true;
end
end
