function rows = regression_tol_mutate(mutations, opts)
% REGRESSION_TOL_MUTATE Apply deliberate errors to the source, one at a time, and record what the tolerance tier catches.
%
%   rows = regression_tol_mutate(M)                          % against the anchor
%   rows = regression_tol_mutate(M, 'reference', 'before_edit')
%   rows = regression_tol_mutate(M(3), 'save_dir', folder)
%
% WHY. The tolerance regression tier (regression_tol_check) gates rewrites
% of the model's arithmetic, so it has to be shown to catch errors before
% it is trusted (FD-020; Plan A, step 0c). This function takes a list of
% small deliberate errors AS DATA, and for each one: edits the source file,
% runs the tier against a named reference, restores the file, and returns
% one row saying whether the tier failed, where, and by how much. EVERY
% MUTATION SHOULD BE CAUGHT; one that is not means a tolerance is too
% loose or a layer is missing.
%
% THE MUTATED CODE IS NEVER COMMITTED AND NEVER LEFT ON DISK. The file is
% restored with `git checkout -- <file>` under an onCleanup, so it is
% restored however the run ends (an error, Ctrl-C), and the function then
% checks that `git hash-object` of the file equals the blob it had in HEAD
% and that `git status --porcelain -- <file>` is empty. It refuses to
% start unless `git status --porcelain -- src/` is empty, and each
% mutation's file must be tracked and unmodified.
%
% M is a struct array (or one struct) with the fields
%   name     a valid short label, unique in the list
%   file     the file to edit, relative to the repository root
%   find     the text to replace. IT MUST OCCUR EXACTLY ONCE in the file;
%            the function errors otherwise and touches nothing. It may span
%            lines: write the line break as newline (char(10)); it is
%            matched against the file's own line ending (CRLF or LF).
%   replace  the replacement text
%   what     (optional) one sentence saying what the error is, for the note
% The lists live in functions that return such an array, one per state of
% the code (regression_tol_mutations_pre_refactor is the first), because
% the text to find changes when the code is rewritten and this harness
% does not.
%
% THE CLASS CACHE. MATLAB keeps a loaded classdef and the persistent
% variables of functions, so after the file is written (and again after it
% is restored) this function runs `clear classes` IN THE BASE WORKSPACE.
% That wipes the base workspace of the session it runs in. It is run in
% the base workspace because `clear classes` inside a function clears that
% function's own variables, the onCleanup among them.
%
% EACH ROW (a struct; struct2table(rows) gives the table):
%   name, file, what
%   caught            the tier did not pass
%   caught_by         'values'       at least one compared leaf is FAIL
%                     'throw only'   no leaf failed by value: the tier failed
%                                    only because cases threw (their leaves
%                                    are missing). WEAKER EVIDENCE.
%                     'check threw'  regression_tol_check itself threw
%                     'NOT CAUGHT'
%   exact_changed     an exact leaf (discrete, integer, logical, text) changed
%   n_moved           float leaves that moved at all (pass or fail)
%   n_errors, errors  the cases that threw, with their messages
%   layer             1 x L struct, one per layer (pointwise, stable):
%                       n_fail, n_missing, n_fail_default (failing leaves
%                       under the layer's own tolerance, not an override),
%                       fail_min_factor / fail_min_key: the failing float
%                       leaf that exceeded its tolerance by the SMALLEST
%                       factor (max_rel / rtol), fail_max_factor /
%                       fail_max_key: by the largest
%   near_factor, near_key   among float leaves that MOVED but PASSED, the one
%                     closest to its tolerance (factor < 1)
%   group             1 x G struct, one per leaf group (opts.groups): n,
%                     n_fail, n_missing, max_factor over its float leaves
%                     (failing or not), its key. The groups split each
%                     layer by what a leaf is (the parameter view, dS, J,
%                     JY, trajectories, top-K, each overridden group, ...),
%                     so a row says which KIND of leaf saw the error.
%   only_override     every failing or missing leaf is under an override
%   only_stable       nothing in the pointwise layer failed or went missing
%   only_pointwise    nothing in the stable layer failed or went missing
%   blob_restored     the file's blob after the restore equals HEAD's
%   seconds, log      the tier's duration and its captured printout
%
% Options:
%   reference  the reference to compare with (default 'anchor'); see
%              regression_tol_check
%   settings   tolerance settings (default regression_tol_settings())
%   groups     n x 2 cell {name, regular expression}; a leaf belongs to the
%              first group whose expression matches its key. Default:
%              default_groups() below, written for regression_tol_run.
%   save_dir   when given, each row is also saved there as <name>.mat (with
%              the table of leaves that failed, went missing or moved), so a
%              long list survives a dropped connection
%   verbose    print one line per mutation (default true)
%   dry_run    only check the list (every file tracked and clean, every text
%              found exactly once) and return; nothing is edited or run.
%              The same check is made on the WHOLE list before the first
%              mutation runs, so a bad entry costs no compute.
%
% See also: regression_tol_mutations_pre_refactor, regression_tol_check,
%           regression_tol_compare, regression_tol_settings
arguments
    mutations struct
    opts.reference (1,:) char = 'anchor'
    opts.settings (1,1) struct = regression_tol_settings()
    opts.groups cell = default_groups()
    opts.save_dir (1,:) char = ''
    opts.verbose (1,1) logical = true
    opts.dry_run (1,1) logical = false
end
root = fileparts(which('setup_paths'));
required = {'name', 'file', 'find', 'replace'};
missing = setdiff(required, fieldnames(mutations));
if ~isempty(missing)
    error('regression_tol_mutate:BadMutation', 'A mutation needs the field ''%s''.', missing{1});
end
names = {mutations.name};
if numel(unique(names)) ~= numel(names)
    error('regression_tol_mutate:BadMutation', 'Mutation names must be unique.');
end
dirty = git(root, 'status --porcelain -- src/');
if ~isempty(dirty)
    error('regression_tol_mutate:DirtyTree', ...
        'src/ has uncommitted changes; a mutation is restored with git checkout and would destroy them:\n%s', dirty);
end
if ~isempty(opts.save_dir) && ~isfolder(opts.save_dir)
    mkdir(opts.save_dir);
end

for k = 1:numel(mutations)
    prepare(mutations(k), root);         % errors on a bad entry, edits nothing
end
if opts.dry_run
    rows = struct('name', names, 'file', {mutations.file});
    if opts.verbose
        fprintf('regression_tol_mutate: %d mutation(s) checked, each text found exactly once\n', numel(mutations));
    end
    return;
end

rows = cell(1, numel(mutations));
for k = 1:numel(mutations)
    rows{k} = run_one(mutations(k), root, opts);
    if opts.verbose
        r = rows{k};
        fprintf('%-34s %-11s fail/missing: pointwise %d/%d, stable %d/%d; threw %d; restored %d\n', ...
            r.name, r.caught_by, r.layer(1).n_fail, r.layer(1).n_missing, ...
            r.layer(2).n_fail, r.layer(2).n_missing, r.n_errors, r.blob_restored);
    end
end
rows = [rows{:}];
end

%% ------------------------------------------------------------------------
function row = run_one(mut, root, opts)
[file, rel, blob_head, mutated] = prepare(mut, root);

row = struct('name', mut.name, 'file', mut.file, 'what', '');
if isfield(mut, 'what') && ~isempty(mut.what); row.what = mut.what; end

% From here on the file is restored however this function ends.
restore = onCleanup(@() restore_file(root, rel));
fid = fopen(file, 'w');
fwrite(fid, mutated);
fclose(fid);
drop_cached_code();

check_error = '';
ok = false; leaves = table(); info = struct();
try
    [log, ok, leaves, info] = evalc( ...
        'regression_tol_check(opts.reference, ''settings'', opts.settings, ''verbose'', true)');
catch err
    check_error = sprintf('%s (%s)', err.message, err.identifier);
    log = check_error;
end

delete(restore);                         % git checkout and a second clear, now
blob_now = strtrim(git(root, sprintf('hash-object -- "%s"', rel)));
row.blob_restored = strcmp(blob_now, blob_head) && ...
    isempty(git(root, sprintf('status --porcelain -- "%s"', rel)));
if ~row.blob_restored
    error('regression_tol_mutate:NotRestored', ...
        '%s: %s was NOT restored (blob %s, HEAD has %s). Restore it by hand: git checkout -- %s', ...
        mut.name, mut.file, blob_now, blob_head, rel);
end

row = summarise(row, ok, leaves, info, check_error, opts);
row.log = log;
if ~isempty(opts.save_dir)
    changed = table();
    if ~isempty(leaves)
        changed = leaves(~strcmp(leaves.status, 'pass') | leaves.max_rel > 0, :);
    end
    save(fullfile(opts.save_dir, [mut.name '.mat']), 'row', 'changed', 'mut');
end
end

function [file, rel, blob_head, mutated] = prepare(mut, root)
% Check one mutation and return the mutated file as bytes. Edits nothing.
file = fullfile(root, mut.file);
if ~isfile(file)
    error('regression_tol_mutate:BadMutation', '%s: no file %s.', mut.name, mut.file);
end
rel = strrep(mut.file, '\', '/');
tree = git(root, sprintf('ls-tree HEAD -- "%s"', rel));
tok = regexp(tree, 'blob ([0-9a-f]{40})', 'tokens', 'once');
if isempty(tok)
    error('regression_tol_mutate:BadMutation', '%s: %s is not tracked by git.', mut.name, mut.file);
end
blob_head = tok{1};
if ~isempty(git(root, sprintf('status --porcelain -- "%s"', rel)))
    error('regression_tol_mutate:DirtyTree', '%s: %s is modified; refusing to mutate it.', mut.name, mut.file);
end

% The file as bytes, so its encoding and line endings come back unchanged.
fid = fopen(file, 'r');
text = fread(fid, [1 Inf], 'uint8=>char');   % one char per byte, no decoding
fclose(fid);
needle = char(mut.find);
patch = char(mut.replace);
if contains(text, sprintf('\r\n'))
    needle = strrep(strrep(needle, sprintf('\r\n'), newline), newline, sprintf('\r\n'));
    patch  = strrep(strrep(patch,  sprintf('\r\n'), newline), newline, sprintf('\r\n'));
end
hits = strfind(text, needle);
if numel(hits) ~= 1
    error('regression_tol_mutate:BadMutation', ...
        '%s: the text to find occurs %d times in %s; it must occur exactly once.', ...
        mut.name, numel(hits), mut.file);
end
if strcmp(needle, patch)
    error('regression_tol_mutate:BadMutation', '%s: the replacement equals the text it replaces.', mut.name);
end
mutated = uint8([text(1:hits - 1), patch, text(hits + numel(needle):end)]);
end

function restore_file(root, rel)
git(root, sprintf('checkout -- "%s"', rel));
drop_cached_code();
end

function drop_cached_code()
% See THE CLASS CACHE in the header. rehash makes MATLAB look at the file
% again even where its change notification is late. Warnings are off for
% the one statement: MATLAB warns that classes with live objects (this
% function's onCleanup, the datetime of a pool record) cannot be cleared;
% the model classes have no live object here and are cleared.
state = warning('off', 'all');
evalc('evalin(''base'', ''clear classes'')');
warning(state);
rehash;
end

function out = git(root, args)
[status, out] = system(sprintf('git -C "%s" %s', root, args));
if status ~= 0
    error('regression_tol_mutate:Git', 'git %s failed:\n%s', args, out);
end
out = strtrim(out);
end

%% ------------------------------------------------------------------------
function row = summarise(row, ok, leaves, info, check_error, opts)
layer_names = {'pointwise', 'stable'};
blank_layer = struct('name', '', 'n_fail', 0, 'n_missing', 0, 'n_fail_default', 0, ...
    'fail_min_factor', NaN, 'fail_min_key', '', 'fail_max_factor', NaN, 'fail_max_key', '');
blank_group = struct('name', '', 'n', 0, 'n_fail', 0, 'n_missing', 0, 'max_factor', 0, 'max_key', '');
row.caught = ~ok;
row.caught_by = 'NOT CAUGHT';
row.exact_changed = false;
row.n_moved = 0;
row.n_errors = 0;
row.errors = {};
row.check_error = check_error;
row.layer = repmat(blank_layer, 1, numel(layer_names));
for q = 1:numel(layer_names); row.layer(q).name = layer_names{q}; end
row.near_factor = 0;
row.near_key = '';
row.group = repmat(blank_group, 1, size(opts.groups, 1));
for g = 1:size(opts.groups, 1); row.group(g).name = opts.groups{g, 1}; end
row.only_override = false;
row.only_stable = false;
row.only_pointwise = false;
row.seconds = NaN;
if ~isempty(check_error)
    row.caught_by = 'check threw';
    return;
end
row.seconds = info.seconds;
row.errors = info.data.info.errors(:)';
row.n_errors = numel(row.errors);

is_fail = strcmp(leaves.status, 'FAIL');
is_missing = strcmp(leaves.status, 'missing');
is_float = strcmp(leaves.kind, 'float');
factor = leaves.max_rel ./ leaves.rtol;          % how many tolerances a leaf moved
factor(~is_float | ~(leaves.max_rel > 0)) = 0;   % also the NaN of a non-finite failure
% A float leaf that failed without a number (its size, class or non-finite
% pattern changed) has no factor; it counts as failed, not in the factors.
has_factor = is_float & isfinite(leaves.max_rel) & leaves.max_rel > 0;

if any(is_fail)
    row.caught_by = 'values';
elseif ~ok
    row.caught_by = 'throw only';
end
row.exact_changed = any(is_fail & ~is_float);
row.n_moved = nnz(has_factor);

% Under an override: the leaf's tolerance is not its layer's own.
overridden = false(height(leaves), 1);
for k = 1:size(opts.settings.override, 1)
    overridden = overridden | ~cellfun(@isempty, regexp(leaves.key, opts.settings.override{k, 1}, 'once'));
end

for q = 1:numel(layer_names)
    in = strcmp(leaves.layer, layer_names{q});
    L = row.layer(q);
    L.n_fail = nnz(in & is_fail);
    L.n_missing = nnz(in & is_missing);
    L.n_fail_default = nnz(in & is_fail & ~overridden);
    pick = find(in & is_fail & has_factor);
    if ~isempty(pick)
        [L.fail_min_factor, w] = min(factor(pick));
        L.fail_min_key = leaves.key{pick(w)};
        [L.fail_max_factor, w] = max(factor(pick));
        L.fail_max_key = leaves.key{pick(w)};
    end
    row.layer(q) = L;
end

pick = find(strcmp(leaves.status, 'pass') & has_factor);
if ~isempty(pick)
    [row.near_factor, w] = max(factor(pick));
    row.near_key = leaves.key{pick(w)};
end

group_of = zeros(height(leaves), 1);
for g = size(opts.groups, 1):-1:1                 % so that the first match wins
    group_of(~cellfun(@isempty, regexp(leaves.key, opts.groups{g, 2}, 'once'))) = g;
end
for g = 1:size(opts.groups, 1)
    in = group_of == g;
    G = row.group(g);
    G.n = nnz(in);
    G.n_fail = nnz(in & is_fail);
    G.n_missing = nnz(in & is_missing);
    pick = find(in & has_factor);
    if ~isempty(pick)
        [G.max_factor, w] = max(factor(pick));
        G.max_key = leaves.key{pick(w)};
    end
    row.group(g) = G;
end

bad = is_fail | is_missing;
row.only_override  = any(bad) && all(overridden(bad));
row.only_stable    = any(bad) && ~any(bad & strcmp(leaves.layer, 'pointwise'));
row.only_pointwise = any(bad) && ~any(bad & strcmp(leaves.layer, 'stable'));
end

function groups = default_groups()
% What kind of leaf saw the error. The first matching row wins, so the four
% overridden groups of regression_tol_settings come first.
groups = { ...
    'ov_orth_defect',   '(^|\.)orth_defect_max$'; ...
    'ov_fd_jacobian',   '(^|\.)finite_difference_J$'; ...
    'ov_qr',            '^stable\.estimators\.\w+\.qr\.'; ...
    'ov_mc',            '^stable\.mc\.\w+\.(MC|R2)$'; ...
    'pw_view',          '^pointwise\..*\.view\.'; ...
    'pw_dS',            '^pointwise\..*\.dS$'; ...
    'pw_J',             '^pointwise\..*\.J(\{\d+\}|#size)$'; ...
    'pw_JY',            '^pointwise\..*\.JY$'; ...
    'pw_rate_theta',    '^pointwise\..*\.(rate|synaptic_output)\.'; ...
    'pw_other',         '^pointwise\.'; ...
    'st_trajectory',    '^stable\.(model\.\w+\.\w+|options)\.(S0|S_kept|t_kept|W|plot_data\.|tau_a_clipped)'; ...
    'st_topk',          '^stable\.(model\.\w+\.\w+|options)\.(lya_results|lya_summary)\.'; ...
    'st_gain',          '^stable\.extra\.transient_gain\.'; ...
    'st_benettin',      '^stable\.estimators\.\w+\.benettin\.'; ...
    'st_sweep',         '^stable\.sweep\.'; ...
    'st_other',         '^stable\.'};
end
