function [file, committed] = regression_tol_reference_file(name)
% REGRESSION_TOL_REFERENCE_FILE Where a named reference of the tolerance tier lives.
%
%   file = regression_tol_reference_file('anchor')    % scripts/tests/fixtures/regression_tol_anchor.mat
%   file = regression_tol_reference_file('chain')     % scripts/tests/fixtures/regression_tol_chain.mat
%   file = regression_tol_reference_file('my_name')   % data/regression_tol/my_name.mat (gitignored)
%   [file, committed] = regression_tol_reference_file(name)
%
% One definition, shared by make_regression_tol_reference (which writes)
% and regression_tol_check (which reads), so the two cannot disagree. The
% file is not required to exist. COMMITTED is true for the two names whose
% reference is a fixture kept in git: 'anchor' (the base suite, the model
% class before Plan A touched it) and 'chain' (the chain suite; see
% regression_tol_run). Every other name is a pre-edit reference.
%
% See also: make_regression_tol_reference, regression_tol_check
arguments
    name (1,:) char {mustBeNonempty}
end
project_root = fileparts(which('setup_paths'));
committed = any(strcmp(name, {'anchor', 'chain'}));
if committed
    file = fullfile(project_root, 'scripts', 'tests', 'fixtures', ['regression_tol_' name '.mat']);
else
    file = fullfile(project_root, 'data', 'regression_tol', [name '.mat']);
end
end
