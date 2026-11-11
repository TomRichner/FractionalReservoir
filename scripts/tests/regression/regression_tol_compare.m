function [ok, leaves, summary] = regression_tol_compare(ref, new, opts)
% REGRESSION_TOL_COMPARE Compare two outputs of regression_tol_run element by element, within a tolerance.
%
%   [ok, leaves, summary] = regression_tol_compare(ref, new)
%   ... = regression_tol_compare(ref, new, 'rtol', 0)            % every tolerance zero: exact
%   ... = regression_tol_compare(ref, new, 'settings', T, 'verbose', false)
%
% REF and NEW are structs from regression_tol_run (or any nested structs
% and cells of arrays). Both are flattened to LEAVES, one per array, keyed
% by where the array sits ('pointwise.pairs_steps.sfa3_std2.sra1_noise_on.J{3}').
% Only the branches 'pointwise' and 'stable' are compared when REF has
% them; 'info' (the thread count and the like) is not a result.
%
% WHY NOT regression_compare. That function compares a digest: seven
% statistics and a hash per array. A sum over a Jacobian can hide a wrong
% entry behind a large one and cannot say where a difference is. This one
% compares every element.
%
% PER LEAF it reports, on a pass as well as a failure, so that a pass says
% how close it was:
%   max_abs       the largest |NEW - REF|
%   max_rel       the largest relative difference under the tier's rule
%                 (regression_tol_settings): |NEW_i - REF_i| divided by
%                 max(|REF_i|, floor_frac * the leaf's largest |REF|).
%                 The leaf passes when max_rel <= rtol.
%   max_rel_elem  the largest purely element-wise |NEW_i - REF_i| / |REF_i|
%                 over the entries where REF is not zero (information only)
%   where         the subscripts of the element that gave max_rel
%   rtol          the tolerance that applied
%
% WHAT IS EXACT whatever the tolerance: sizes; classes; the positions (and
% signs) of non-finite values; integer, logical and text leaves; and the
% float leaves named DISCRETE in the settings (counts and indices).
%
% SPARSE MATRICES compare by value over the union of the two patterns. A
% changed sparsity pattern (an entry stored in one and absent in the
% other) is noted in DETAIL and does not fail by itself: if the values
% agree within tolerance the entry that appeared or vanished is within
% tolerance of zero.
%
% STATUS of a leaf:
%   pass      within tolerance (or exactly equal, for the exact kinds)
%   FAIL      outside it; DETAIL says how
%   missing   in REF, not in NEW: fails
%   new       in NEW, not in REF: REPORTED, does not fail. The refactor and
%             later models add fields (to route_fields, to lya_summary);
%             a reference cannot hold what did not exist when it was made.
%
% OK is true when no leaf is FAIL or missing. LEAVES is a table of every
% leaf (key, layer, kind, status, n, max_abs, max_rel, max_rel_elem, rtol,
% where, detail). SUMMARY has one entry per layer: the number of leaves,
% the failures, the worst max_rel and its key, and the tolerance.
%
% Options: settings (default regression_tol_settings()), rtol (a scalar
% that replaces EVERY tolerance, overrides included; 0 is the exact
% comparison), verbose (default true), max_print (default 40).
%
% See also: regression_tol_run, regression_tol_settings, regression_tol_check
arguments
    ref (1,1) struct
    new (1,1) struct
    opts.settings (1,1) struct = regression_tol_settings()
    opts.rtol double {mustBeScalarOrEmpty, mustBeNonnegative} = []
    opts.verbose (1,1) logical = true
    opts.max_print (1,1) double = 40
end
T = opts.settings;

layers = {'pointwise', 'stable'};
if all(isfield(ref, layers))
    ref = struct('pointwise', ref.pointwise, 'stable', ref.stable);
    if all(isfield(new, layers))
        new = struct('pointwise', new.pointwise, 'stable', new.stable);
    end
end
[ref_keys, ref_vals] = flatten(ref, '');
[new_keys, new_vals] = flatten(new, '');
[~, ia, ib] = intersect(ref_keys, new_keys, 'stable');
only_ref = setdiff(ref_keys, new_keys, 'stable');
only_new = setdiff(new_keys, ref_keys, 'stable');

n_rows = numel(ia) + numel(only_ref) + numel(only_new);
rows = repmat(blank_row(), n_rows, 1);
for k = 1:numel(ia)
    key = ref_keys{ia(k)};
    [rtol, layer] = tolerance_for(key, T, opts.rtol);
    discrete = any(~cellfun(@isempty, regexp(key, T.discrete, 'once')));
    row = compare_leaf(ref_vals{ia(k)}, new_vals{ib(k)}, rtol, T.floor_frac, discrete);
    row.key = key;
    row.layer = layer;
    rows(k) = row;
end
at = numel(ia);
for k = 1:numel(only_ref)
    row = blank_row();
    row.key = only_ref{k};
    [~, row.layer] = tolerance_for(row.key, T, opts.rtol);
    row.status = 'missing';
    row.detail = 'in the reference, not in this run';
    rows(at + k) = row;
end
at = at + numel(only_ref);
for k = 1:numel(only_new)
    row = blank_row();
    row.key = only_new{k};
    [~, row.layer] = tolerance_for(row.key, T, opts.rtol);
    row.status = 'new';
    row.detail = 'in this run, not in the reference (reported, not failed)';
    rows(at + k) = row;
end
leaves = struct2table(rows, 'AsArray', true);
leaves = leaves(:, {'key', 'layer', 'kind', 'status', 'n', 'max_abs', 'max_rel', ...
    'max_rel_elem', 'rtol', 'where', 'detail'});

bad = ismember(leaves.status, {'FAIL', 'missing'});
ok = ~any(bad);

layer_names = unique(leaves.layer, 'stable');
summary = struct('layer', {}, 'n_leaves', {}, 'n_fail', {}, 'n_missing', {}, 'n_new', {}, ...
    'worst_max_rel', {}, 'worst_key', {}, 'rtol', {});
for q = 1:numel(layer_names)
    in = strcmp(leaves.layer, layer_names{q});
    L = leaves(in, :);
    s = struct('layer', layer_names{q}, 'n_leaves', height(L), ...
        'n_fail', nnz(strcmp(L.status, 'FAIL')), 'n_missing', nnz(strcmp(L.status, 'missing')), ...
        'n_new', nnz(strcmp(L.status, 'new')), 'worst_max_rel', 0, 'worst_key', '', 'rtol', NaN);
    is_float = strcmp(L.kind, 'float') & isfinite(L.max_rel);
    if any(is_float)
        F = L(is_float, :);
        [s.worst_max_rel, w] = max(F.max_rel);
        s.worst_key = F.key{w};
        s.rtol = min(F.rtol);
    end
    summary(end + 1) = s; %#ok<AGROW>
end

if opts.verbose
    fprintf('regression_tol_compare: %d leaves in the reference, %d in this run, %d in common\n', ...
        numel(ref_keys), numel(new_keys), numel(ia));
    for q = 1:numel(summary)
        s = summary(q);
        fprintf('  %-10s %4d leaves, rtol %-7.1e  worst relative difference %.2e', ...
            s.layer, s.n_leaves, s.rtol, s.worst_max_rel);
        if s.worst_max_rel > 0
            fprintf('\n             at %s', s.worst_key);
            if s.rtol > 0
                fprintf('\n             (tolerance / worst = %.1e)', s.rtol / s.worst_max_rel);
            end
        end
        fprintf('\n             %d FAIL, %d missing, %d new\n', s.n_fail, s.n_missing, s.n_new);
    end
    noted = leaves(strcmp(leaves.status, 'pass') & ~cellfun(@isempty, leaves.detail), :);
    for k = 1:min(height(noted), opts.max_print)
        fprintf('    note    %s\n            %s\n', noted.key{k}, noted.detail{k});
    end
    fresh = leaves(strcmp(leaves.status, 'new'), :);
    for k = 1:min(height(fresh), opts.max_print)
        fprintf('    new     %s\n', fresh.key{k});
    end
    if height(fresh) > opts.max_print
        fprintf('    ... and %d more new leaves (see the returned table)\n', height(fresh) - opts.max_print);
    end
    failed = leaves(bad, :);
    [~, order] = sort(fillmissing(failed.max_rel, 'constant', Inf), 'descend');
    failed = failed(order, :);
    for k = 1:min(height(failed), opts.max_print)
        fprintf('    %-7s %s\n            %s\n', failed.status{k}, failed.key{k}, failed.detail{k});
    end
    if height(failed) > opts.max_print
        fprintf('    ... and %d more (see the returned table)\n', height(failed) - opts.max_print);
    end
    if ok
        fprintf('  PASS: every leaf of the reference is within tolerance\n');
    else
        fprintf('  FAIL: %d leaf(s) outside tolerance or missing\n', nnz(bad));
    end
end
end

%% ------------------------------------------------------------------------
function row = blank_row()
row = struct('key', '', 'layer', '', 'kind', '', 'status', '', 'n', NaN, 'max_abs', NaN, ...
    'max_rel', NaN, 'max_rel_elem', NaN, 'rtol', NaN, 'where', '', 'detail', '');
end

function [rtol, layer] = tolerance_for(key, T, forced)
layer = strtok(key, '.{(');
if isfield(T.rtol, layer)
    rtol = T.rtol.(layer);
else
    % A struct that is not a regression_tol_run output: the tightest layer.
    rtol = T.rtol.pointwise;
    layer = 'other';
end
for k = 1:size(T.override, 1)
    if ~isempty(regexp(key, T.override{k, 1}, 'once'))
        rtol = T.override{k, 2};
        break;
    end
end
if ~isempty(forced)
    rtol = forced;
end
end

function [keys, vals] = flatten(v, prefix)
% One leaf per array. Scalar structs recurse by field, struct arrays by
% element then field, cells by element; an empty struct or cell is a leaf.
keys = {}; vals = {};
if isstruct(v) && isscalar(v)
    f = fieldnames(v);
    if isempty(f)
        keys = {prefix}; vals = {v};
        return;
    end
    for k = 1:numel(f)
        if isempty(prefix); sub = f{k}; else; sub = [prefix '.' f{k}]; end
        [kk, vv] = flatten(v.(f{k}), sub);
        keys = [keys; kk]; vals = [vals; vv]; %#ok<AGROW>
    end
elseif (isstruct(v) || iscell(v)) && ~isempty(v)
    for k = 1:numel(v)
        if iscell(v)
            [kk, vv] = flatten(v{k}, sprintf('%s{%d}', prefix, k));
        else
            [kk, vv] = flatten(v(k), sprintf('%s(%d)', prefix, k));
        end
        keys = [keys; kk]; vals = [vals; vv]; %#ok<AGROW>
    end
    % The container's own shape is a leaf too, so 2 x 3 against 3 x 2 is seen.
    keys = [keys; {[prefix '#size']}]; vals = [vals; {int32(size(v))}];
else
    keys = {prefix}; vals = {v};
end
end

function row = compare_leaf(a, b, rtol, floor_frac, discrete)
% A is the reference's array, B this run's.
row = blank_row();
row.status = 'pass';
row.n = numel(a);
if ~strcmp(class(a), class(b))
    row.kind = class(a);
    row.status = 'FAIL';
    row.detail = sprintf('class %s -> %s', class(a), class(b));
    return;
end
if ~isequal(size(a), size(b))
    row.kind = class(a);
    row.status = 'FAIL';
    row.detail = sprintf('size %s -> %s', size_text(a), size_text(b));
    return;
end

if ~isfloat(a) || discrete
    % Exact kinds: integers, logicals, text, discrete floats, and anything
    % else (compared with isequaln; an object that cannot be compared is
    % noted and passed, since only its class and size are known to match).
    if discrete; row.kind = 'discrete'; else; row.kind = class(a); end
    try
        same = isequaln(a, b);
    catch
        same = true;
        row.detail = sprintf('%s not comparable; class and size match', class(a));
    end
    if ~same
        row.status = 'FAIL';
        if isnumeric(a) || islogical(a)
            d = abs(double(full(b(:))) - double(full(a(:))));
            [row.max_abs, w] = max(d);
            row.where = subscript_text(size(a), w);
            row.detail = sprintf('exact leaf differs: %d of %d element(s), largest |difference| %g at %s', ...
                nnz(d > 0 | isnan(d)), numel(d), row.max_abs, row.where);
        elseif ischar(a) || isstring(a)
            row.detail = sprintf('text ''%s'' -> ''%s''', clip(char(join(string(a(:)'), ''))), ...
                clip(char(join(string(b(:)'), ''))));
        else
            row.detail = sprintf('%s leaf differs', class(a));
        end
    else
        row.max_abs = 0; row.max_rel = 0;
    end
    return;
end

row.kind = 'float';
row.rtol = rtol;
note = '';
if issparse(a) || issparse(b)
    if issparse(a) ~= issparse(b)
        note = 'storage changed between sparse and full';
    end
    pa = a ~= 0; pb = b ~= 0;
    if ~isequal(pa, pb)
        note = strtrim(sprintf('%s sparsity pattern changed (nnz %d -> %d; %d entries in one only)', ...
            note, nnz(pa), nnz(pb), nnz(xor(pa, pb))));
    end
    idx = find(pa | pb);
    av = full(a(idx)); bv = full(b(idx));
else
    idx = (1:numel(a))';
    av = a(:); bv = b(:);
end
odd = ~isfinite(av) | ~isfinite(bv);
if any(odd)
    if ~isequaln(av(odd), bv(odd))
        row.status = 'FAIL';
        w = find(odd & ~arrayfun(@isequaln, av, bv), 1);
        row.where = subscript_text(size(a), idx(w));
        row.detail = sprintf('non-finite values differ (%d in the reference, %d in this run), first at %s', ...
            nnz(~isfinite(av)), nnz(~isfinite(bv)), row.where);
        return;
    end
    idx = idx(~odd); av = av(~odd); bv = bv(~odd);
end
if isempty(av)
    row.max_abs = 0; row.max_rel = 0; row.max_rel_elem = 0;
    row.detail = note;
    return;
end
scale = max(abs(av));
d = abs(bv - av);
unit = max(abs(av), floor_frac * scale);
rel = d ./ unit;
rel(d == 0) = 0;                     % 0/0 where the leaf is identically zero
[row.max_rel, w] = max(rel);
row.max_abs = max(d);
row.where = subscript_text(size(a), idx(w));
nz = av ~= 0;
if any(nz)
    row.max_rel_elem = max(d(nz) ./ abs(av(nz)));
else
    row.max_rel_elem = 0;
end
if any(d > rtol .* unit)
    row.status = 'FAIL';
    row.detail = strtrim(sprintf('relative difference %.3e > rtol %.1e (|difference| %.3e, leaf scale %.3e) at %s; %d of %d element(s) outside. %s', ...
        row.max_rel, rtol, d(w), scale, row.where, nnz(d > rtol .* unit), numel(d), note));
else
    row.detail = note;
end
end

function t = size_text(v)
t = strjoin(arrayfun(@(n) sprintf('%d', n), size(v), 'UniformOutput', false), 'x');
end

function t = subscript_text(sz, linear_index)
if isempty(linear_index)
    t = '';
    return;
end
subs = cell(1, numel(sz));
[subs{:}] = ind2sub(sz, linear_index);
t = ['(' strjoin(cellfun(@(n) sprintf('%d', n), subs, 'UniformOutput', false), ',') ')'];
end

function t = clip(t)
if numel(t) > 60; t = [t(1:57) '...']; end
end
