function out_file = make_regression_tol_reference(name, opts)
% MAKE_REGRESSION_TOL_REFERENCE Capture a reference for the tolerance regression tier.
%
%   out_file = make_regression_tol_reference('anchor')
%   out_file = make_regression_tol_reference('chain')
%   out_file = make_regression_tol_reference('before_jacobian_edit')
%   out_file = make_regression_tol_reference('before_jacobian_edit', 'overwrite', true)
%   out_file = make_regression_tol_reference('before_facilitation', 'suite', 'chain')
%
% Runs regression_tol_run in capture mode (about half a minute, in the
% client, no pool needed) and saves what it returned with a description of
% the code and machine it came from.
%
% THE TIER HAS TWO SUITES of cases (regression_tol_run): 'base', the
% networks the model class had before Plan A, and 'chain', networks with
% chained vesicle-pool depression. A reference holds one suite and records
% which (meta.suite); regression_tol_check reads it from there.
%
% THREE KINDS OF REFERENCE (Plan A, sections 2.2 and 7.3 step 5):
%
%   'anchor'   scripts/tests/fixtures/regression_tol_anchor.mat, kept in git
%              (force-added). The BASE suite. It describes the model class
%              BEFORE Plan A touched it and is compared after every merge,
%              so that a chain of re-captures cannot walk away from the
%              original. It is written ONCE: this refuses to overwrite it
%              and refuses a tree with uncommitted changes to tracked
%              files. Unlike the byte-exact references it is not to be
%              deleted and regenerated when results move; once the code it
%              describes is gone it cannot be made again. (The one
%              legitimate re-capture is from a checkout of the commit
%              recorded in it.)
%
%   'chain'    scripts/tests/fixtures/regression_tol_chain.mat, kept in git
%              (force-added). The CHAIN suite, captured when the chain was
%              first merged, before chained facilitation was built into the
%              same route model: with facilitation off the chain must go on
%              agreeing with it. Like the byte-exact references it is
%              replaced only deliberately: this refuses to overwrite it, so
%              delete the .mat, run this again from a clean committed tree,
%              and say in the commit message WHY the expected values moved
%              (a new case in the suite is a reason; a check that fails is
%              not).
%
%   any other name   data/regression_tol/<name>.mat, gitignored, never
%              committed. A PRE-EDIT reference of either suite ('suite',
%              default 'base'): capture one before a critical edit (the
%              Jacobian, say), make the edit, then
%              regression_tol_check('<name>'). Uncommitted changes are
%              allowed and recorded. Refuses to overwrite unless asked.
%
% WHAT IS STORED: data (the struct from regression_tol_run) and meta:
%   name, suite, made_on, seconds   which reference, of which suite, when,
%                            how long the capture took
%   commit, dirty            HEAD, and `git status --porcelain` of tracked files
%   blobs                    the git blob hashes of the working copies of
%                            src/model/SRNNCellTypePairs.m,
%                            src/model/SRNN_ESN_reservoir.m,
%                            src/analysis/ParamSpaceAnalysis2.m and the
%                            three route-model files (git hash-object, so
%                            they equal `git ls-tree` on a clean tree)
%   matlab_release, machine, cores, threads   where it ran; THREADS is
%                            maxNumCompThreads, which moves results at
%                            round-off
%   settings                 regression_tol_settings() as in force
% (The anchor predates meta.suite and the route files' blobs; a reference
% without meta.suite is the base suite.)
%
% The two fixtures need `git add -f` (*.mat is gitignored repo-wide).
%
% See also: regression_tol_run, regression_tol_check, regression_tol_settings
arguments
    name (1,:) char {mustBeNonempty}
    opts.overwrite (1,1) logical = false
    opts.suite (1,:) char {mustBeMember(opts.suite, {'', 'base', 'chain'})} = ''
end
if isempty(regexp(name, '^[A-Za-z0-9_\-]+$', 'once'))
    error('make_regression_tol_reference:BadName', ...
        'A reference name may hold letters, digits, _ and - only: ''%s''.', name);
end
project_root = fileparts(which('setup_paths'));
[out_file, committed] = regression_tol_reference_file(name);
% The two committed references each belong to one suite.
suite = opts.suite;
fixed_suite = struct('anchor', 'base', 'chain', 'chain');
if committed
    if ~isempty(suite) && ~strcmp(suite, fixed_suite.(name))
        error('make_regression_tol_reference:BadSuite', ...
            'The reference ''%s'' holds the ''%s'' suite; it cannot be captured from ''%s''.', ...
            name, fixed_suite.(name), suite);
    end
    suite = fixed_suite.(name);
elseif isempty(suite)
    suite = 'base';
end
if isfile(out_file) && (committed || ~opts.overwrite)
    if strcmp(name, 'anchor')
        error('make_regression_tol_reference:Exists', ...
            ['%s already exists.\nThe anchor describes the model class before Plan A touched it and is ' ...
             'written once. It is not regenerated to make a check pass; see this function''s header.'], out_file);
    elseif committed
        error('make_regression_tol_reference:Exists', ...
            ['%s already exists.\nReplacing it is deliberate: delete the file, run this again from a ' ...
             'clean committed tree, and record in the commit message why the expected values moved.'], out_file);
    end
    error('make_regression_tol_reference:Exists', ...
        '%s already exists. Pass ''overwrite'', true to replace a pre-edit reference.', out_file);
end

[status, commit] = system(sprintf('git -C "%s" rev-parse HEAD', project_root));
if status ~= 0; error('make_regression_tol_reference:NoGit', 'git rev-parse failed: %s', commit); end
commit = strtrim(commit);
[~, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', project_root));
dirty = strtrim(dirty);
if committed
    % Replacing the chain reference means deleting the tracked fixture
    % first, which git reports as a change; that one line is this
    % function's own business. (Not so for the anchor: a deleted anchor
    % leaves the tree dirty and is refused, as it always was.)
    if ~strcmp(name, 'anchor')
        [~, fixture, ext] = fileparts(out_file);
        lines = strtrim(splitlines(string(dirty)));
        lines = lines(strlength(lines) > 0 & ~endsWith(lines, [fixture ext]));
        dirty = char(strjoin(lines, newline));
    end
    if ~isempty(dirty)
        error('make_regression_tol_reference:DirtyTree', ...
            'Tracked files have uncommitted changes; a committed reference must describe a commit.\n%s', dirty);
    end
end
tracked = {'src/model/SRNNCellTypePairs.m', 'src/model/SRNN_ESN_reservoir.m', ...
    'src/analysis/ParamSpaceAnalysis2.m', 'src/model/routes/route_model.m', ...
    'src/model/routes/route_product.m', 'src/model/routes/route_pool.m'};
blobs = struct('file', tracked, 'blob', '');
for k = 1:numel(tracked)
    [status, h] = system(sprintf('git -C "%s" hash-object "%s"', project_root, tracked{k}));
    if status ~= 0; error('make_regression_tol_reference:NoGit', 'git hash-object failed: %s', h); end
    blobs(k).blob = strtrim(h);
end

[data, seconds] = regression_tol_run([], suite);

meta = struct();
meta.name = name;
meta.suite = suite;
meta.made_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
meta.seconds = seconds;
meta.commit = commit;
meta.dirty = dirty;
meta.blobs = blobs;
meta.matlab_release = version;
meta.machine = computer;
meta.cores = feature('numcores');
meta.threads = data.info.threads;
meta.settings = regression_tol_settings();

if ~isfolder(fileparts(out_file)); mkdir(fileparts(out_file)); end
save(out_file, 'data', 'meta', '-v7');
d = dir(out_file);
fprintf('Wrote %s\n  suite ''%s'', %.1f MB, captured in %.0f s on %d thread(s)\n', out_file, suite, ...
    d.bytes / 2^20, seconds, meta.threads);
fprintf('  commit %s%s\n', commit, repmat(' (with uncommitted changes)', 1, ~isempty(dirty)));
for k = 1:numel(blobs)
    fprintf('  %s  %s\n', blobs(k).blob, blobs(k).file);
end
if committed
    fprintf('  *.mat is gitignored -- commit with:  git add -f %s\n', strrep(out_file, '\', '/'));
end
end
