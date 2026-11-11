function [ok, diffs] = regression_compare(ref, new, opts)
% REGRESSION_COMPARE Compare two regression digests and list what differs.
%
%   [ok, diffs] = regression_compare(ref, new)
%   [ok, diffs] = regression_compare(ref, new, 'rtol', 1e-9, 'ignore', patterns)
%
% REF and NEW are digests from regression_digest. OK is true when nothing
% differs outside the ignored keys. DIFFS is a table, one row per
% difference, most serious kind first:
%
%   missing      a key in REF that NEW does not have
%   extra        a key in NEW that REF does not have
%   shape        the array changed size, or its count of non-finite values
%   value        a statistic moved by more than the tolerance
%   text         a text value changed
%   file         a file is in one run directory and not the other
%   figure       a figure file is in one figure folder and not the other
%
% TOLERANCE. A statistic s (sum, weighted sum, min, max) differs when
%   |s_new - s_ref| > atol + rtol * scale
% where scale is the array's sum of absolute values (for the two sums) or the
% larger absolute end (for min and max). Scaling the sums by sum_abs, not by
% the sum itself, keeps a quantity that sums to nearly zero from producing
% huge relative differences. The defaults, rtol = 0 and atol = 0, are EXACT
% mode: besides the statistics, the MD5 of each array's bytes is compared,
% so any changed bit is reported. On one machine the pipeline was found
% reproducible to the last bit (see the reference's notes), so in exact mode
% any difference is a change. With a tolerance the hashes are not used.
% DIFFS.rel reports |difference| / scale for value rows, so the size of a
% change can be read off.
%
% IGNORE is a cell of regular expressions matched against the key (or the
% file name for file and figure rows). The reference carries its own list,
% ref.ignore, of keys known to differ between identical runs (timings and
% the like); opts.ignore adds to it.
%
% See also: regression_digest, regression_check
arguments
    ref (1,1) struct
    new (1,1) struct
    opts.rtol (1,1) double {mustBeNonnegative} = 0
    opts.atol (1,1) double {mustBeNonnegative} = 0
    opts.ignore cell = {}
    opts.verbose (1,1) logical = true
    opts.max_print (1,1) double = 40
end
ignore = opts.ignore(:)';
if isfield(ref, 'ignore'); ignore = [ref.ignore(:)', ignore]; end

rows = cell(0, 5);      % kind, key, detail, rel, ignored
% Exact mode (no tolerance) also compares the MD5 of each array's bytes.
exact = opts.rtol == 0 && opts.atol == 0 && isfield(ref, 'hash') && isfield(new, 'hash');

%% keys
[~, ia, ib] = intersect(ref.keys, new.keys, 'stable');
only_ref = setdiff(ref.keys, new.keys, 'stable');
only_new = setdiff(new.keys, ref.keys, 'stable');
for k = 1:numel(only_ref); rows(end+1, :) = {'missing', only_ref{k}, 'in the reference, not in this run', NaN, false}; end %#ok<AGROW>
for k = 1:numel(only_new); rows(end+1, :) = {'extra',   only_new{k}, 'in this run, not in the reference', NaN, false}; end %#ok<AGROW>

%% leaves present in both
names = {'n', 'n_nonfinite', 'sum', 'sum_abs', 'weighted_sum', 'min', 'max'};
for k = 1:numel(ia)
    i = ia(k); j = ib(k);
    key = ref.keys{i};
    if ~strcmp(ref.kind{i}, new.kind{j})
        rows(end+1, :) = {'shape', key, sprintf('kind %s -> %s', ref.kind{i}, new.kind{j}), NaN, false}; %#ok<AGROW>
        continue;
    end
    if strcmp(ref.kind{i}, 'txt')
        if ~strcmp(ref.text{i}, new.text{j})
            rows(end+1, :) = {'text', key, sprintf('''%s'' -> ''%s''', clip(ref.text{i}), clip(new.text{j})), NaN, false}; %#ok<AGROW>
        end
        continue;
    end
    a = ref.stats(i, :); b = new.stats(j, :);
    if ~strcmp(ref.shape{i}, new.shape{j}) || a(1) ~= b(1) || a(2) ~= b(2)
        rows(end+1, :) = {'shape', key, sprintf('%s (n %d, non-finite %d) -> %s (n %d, non-finite %d)', ...
            ref.shape{i}, a(1), a(2), new.shape{j}, b(1), b(2)), NaN, false}; %#ok<AGROW>
        continue;
    end
    scale_sum = max(a(4), b(4));
    scale_ext = max(abs([a(6:7), b(6:7)]), [], 'omitnan');
    if isempty(scale_ext) || isnan(scale_ext); scale_ext = 0; end
    scale = [scale_sum, scale_sum, scale_sum, scale_ext, scale_ext];   % for stats 3:7
    d = abs(b(3:7) - a(3:7));
    d(isnan(a(3:7)) & isnan(b(3:7))) = 0;
    d(isnan(d)) = Inf;                                  % NaN on one side only
    bad = d > opts.atol + opts.rtol .* scale;
    rel = d ./ max(scale, realmin);
    rel(d == 0) = 0;
    if any(bad)
        [worst, w] = max(rel .* bad);
        rows(end+1, :) = {'value', key, sprintf('%s %.10g -> %.10g', names{2 + w}, a(2 + w), b(2 + w)), worst, false}; %#ok<AGROW>
    elseif exact && ~strcmp(ref.hash{i}, new.hash{j})
        % every statistic agrees to the last digit, yet the bytes differ
        rows(end+1, :) = {'value', key, 'bits differ; all statistics agree to double precision', max(rel), false}; %#ok<AGROW>
    end
end

%% files and figures
rows = [rows; list_diff('file', ref.files, new.files)];
if isfield(ref, 'fig_files') && isfield(new, 'fig_files') && ~isempty(ref.fig_files) && ~isempty(new.fig_files)
    rows = [rows; list_diff('figure', ref.fig_files, new.fig_files)];
end

%% ignored keys
for k = 1:size(rows, 1)
    for p = 1:numel(ignore)
        if ~isempty(regexp(rows{k, 2}, ignore{p}, 'once'))
            rows{k, 5} = true;
            break;
        end
    end
end

% Built column by column so that an empty result is still a typed table
% (cell2table on a 0 x 5 cell gives double columns).
diffs = table(rows(:, 1), rows(:, 2), rows(:, 3), ...
    reshape(cell2mat(rows(:, 4)), [], 1), reshape(logical(cell2mat(rows(:, 5))), [], 1), ...
    'VariableNames', {'kind', 'key', 'detail', 'rel', 'ignored'});
order = {'missing', 'extra', 'shape', 'value', 'text', 'file', 'figure'};
if height(diffs) > 0
    [~, rank] = ismember(diffs.kind, order);
    [~, idx] = sortrows([double(diffs.ignored), rank, -fillmissing(diffs.rel, 'constant', 0)]);
    diffs = diffs(idx, :);
end
live = diffs(~diffs.ignored, :);
ok = isempty(live);

if opts.verbose
    fprintf('regression_compare: %d leaves in the reference, %d in this run, %d in common\n', ...
        numel(ref.keys), numel(new.keys), numel(ia));
    fprintf('  tolerance: rtol %g, atol %g;  %d ignore pattern(s)\n', opts.rtol, opts.atol, numel(ignore));
    if ok
        fprintf('  NO DIFFERENCES (%d ignored)\n', nnz(diffs.ignored));
    else
        for q = 1:numel(order)
            c = nnz(strcmp(live.kind, order{q}));
            if c > 0; fprintf('  %-8s %d\n', order{q}, c); end
        end
        fprintf('  (%d ignored)\n', nnz(diffs.ignored));
        for k = 1:min(height(live), opts.max_print)
            if isnan(live.rel(k))
                fprintf('    %-7s %s\n            %s\n', live.kind{k}, live.key{k}, live.detail{k});
            else
                fprintf('    %-7s %s\n            %s  (rel %.2e)\n', live.kind{k}, live.key{k}, live.detail{k}, live.rel(k));
            end
        end
        if height(live) > opts.max_print
            fprintf('    ... and %d more (see the returned table)\n', height(live) - opts.max_print);
        end
    end
end
end

%% ------------------------------------------------------------------------
function rows = list_diff(kind, a, b)
rows = cell(0, 5);
gone = setdiff(a, b, 'stable');
came = setdiff(b, a, 'stable');
for k = 1:numel(gone); rows(end+1, :) = {kind, gone{k}, 'in the reference, not in this run', NaN, false}; end %#ok<AGROW>
for k = 1:numel(came); rows(end+1, :) = {kind, came{k}, 'in this run, not in the reference', NaN, false}; end %#ok<AGROW>
end

function t = clip(t)
if numel(t) > 60; t = [t(1:57) '...']; end
end
