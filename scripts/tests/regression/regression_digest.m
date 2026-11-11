function D = regression_digest(run_dir, opts)
% REGRESSION_DIGEST Reduce a run directory to a small, comparable summary of every number in it.
%
%   D = regression_digest(run_dir)
%   D = regression_digest(run_dir, 'fig_root', fig_root)
%
% WHAT IT IS FOR. A run directory of the regression bundle is tens of
% megabytes of .mat files in dated subfolders. To tell whether a code change
% moved any result, two runs have to be compared number for number, and the
% reference has to be small enough to keep in git. This walks every .mat file
% under RUN_DIR and reduces each numeric array to seven statistics; two
% digests are then compared by regression_compare.
%
% It is deliberately GENERIC. It knows nothing about growth rates or memory
% capacity: it records whatever the stages saved. So a stage that starts
% saving a new quantity is covered without touching this file, and a change
% anywhere in the pipeline shows up, not only in the quantities someone
% thought to list.
%
% WHAT IS RECORDED, one row ("leaf") per array:
%   key     where it is: the file's path relative to RUN_DIR, then the
%           variable and field path, e.g.
%           'memory_capacity/mc_results.mat::results.sfa3_std2.MC_total'
%   kind    'num' or 'txt'
%   shape   the array's size, as text ('3x40'); for a stacked leaf, the
%           number of elements stacked and the total count
%   hash    for 'num' leaves, the MD5 of the array's bytes. An exact
%           comparison uses it, so any changed bit is seen.
%   stats   [n, n_nonfinite, sum, sum_abs, weighted_sum, min, max] over the
%           finite values. Used when a tolerance is allowed, and to say how
%           large a change is. The weighted sum uses position-dependent weights
%           in [1, 2), so two arrays holding the same numbers in a different
%           order do not compare equal.
%   text    for 'txt' leaves, the text (paths under RUN_DIR and under the
%           project root are replaced by <run_dir> and <root>)
%
% STACKING. A sweep saves one result struct per job, hundreds per file. A
% cell or struct array of more than STACK_MIN elements that all have the same
% fields is recorded as ONE leaf per field, over the values of all elements
% in order, not one leaf per element per field. That is what keeps the
% digest to a few thousand rows. Nothing is lost for comparison purposes: a
% change in any one job moves the sum, and a reordering moves the weighted
% sum.
%
% DATED NAMES. The sweeps name their subfolders with the launch time
% ('..._nLevs_5_nov_10_26_19_05'). Those suffixes are stripped from the keys
% so that two runs of the same code have the same keys.
%
% NOT RECORDED: objects other than tables (psa_object.mat holds the same
% results as the summary files), function handles, datetimes and durations
% (their class is noted, not their value), and files that are not .mat. The
% list of every file under RUN_DIR (and under FIG_ROOT when given) IS
% recorded, with dated suffixes stripped, so a stage or a figure that
% silently writes nothing is caught.
%
% See also: regression_compare, regression_check, make_regression_reference
arguments
    run_dir (1,:) char
    opts.fig_root (1,:) char = ''
    opts.verbose (1,1) logical = true
end
STACK_MIN = 8;          % collections larger than this are stacked
TEXT_MAX  = 2000;       % longer text is recorded by length and checksum only

project_root = fileparts(which('setup_paths'));
if ~isfolder(run_dir)
    if isfolder(fullfile(project_root, run_dir))
        run_dir = fullfile(project_root, run_dir);
    else
        error('regression_digest:NoRunDir', 'Run directory not found: %s', run_dir);
    end
end
run_dir = char(java.io.File(run_dir).getCanonicalPath());

% Accumulators (grown geometrically; a digest is a few thousand rows)
cap = 4096; n_leaf = 0;
keys = cell(cap, 1); kind = cell(cap, 1); shape = cell(cap, 1);
stats = nan(cap, 7); txts = cell(cap, 1); hashes = cell(cap, 1);

listing = dir(fullfile(run_dir, '**', '*'));
listing = listing(~[listing.isdir]);
rel = cell(numel(listing), 1);
raw = cell(numel(listing), 1);
for k = 1:numel(listing)
    full_path = fullfile(listing(k).folder, listing(k).name);
    raw{k} = strrep(full_path(numel(run_dir) + 2:end), '\', '/');
    rel{k} = strip_dates(raw{k});
end
% Sort by the date-stripped name, then by the original name, so two folders
% that differ only in their launch time (a stage run twice into one run
% directory, e.g. the sensitivity replot made by the analyses and again by
% the figures) come in time order. The later ones are then told apart with
% '#2', '#3', ...: without that they would share every key.
[~, order] = sortrows([rel, raw]);
rel = rel(order); listing = listing(order);
for k = 2:numel(rel)
    base = regexprep(rel{k - 1}, '#\d+$', '');
    if strcmp(rel{k}, base)
        prev = regexp(rel{k - 1}, '#(\d+)$', 'tokens', 'once');
        if isempty(prev); n_dup = 2; else; n_dup = str2double(prev{1}) + 1; end
        rel{k} = sprintf('%s#%d', rel{k}, n_dup);
    end
end

is_mat = endsWith({listing.name}', '.mat');
warn_state = warning('off', 'all');            % loading old objects is noisy
restore_warn = onCleanup(@() warning(warn_state));
basic = {'double', 'single', 'logical', 'char', 'string', 'cell', 'struct', 'table', ...
         'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'};
for k = find(is_mat)'
    full_path = fullfile(listing(k).folder, listing(k).name);
    try
        info = whos('-file', full_path);
    catch err
        add_text([rel{k} '::<unreadable>'], err.identifier);
        continue;
    end
    for v = 1:numel(info)
        key = [rel{k} '::' info(v).name];
        if ~ismember(info(v).class, basic)
            add_text([key '<class>'], info(v).class);
            continue;
        end
        S = load(full_path, info(v).name);
        walk(key, S.(info(v).name));
    end
    if opts.verbose && mod(find(find(is_mat) == k), 25) == 0
        fprintf('  [digest] %d of %d .mat files, %d leaves\n', ...
            find(find(is_mat) == k), nnz(is_mat), n_leaf);
    end
end

D = struct();
D.version   = 2;
D.keys      = keys(1:n_leaf);
D.kind      = kind(1:n_leaf);
D.shape     = shape(1:n_leaf);
D.stats     = stats(1:n_leaf, :);
D.text      = txts(1:n_leaf);
D.hash      = hashes(1:n_leaf);
D.stat_names = {'n', 'n_nonfinite', 'sum', 'sum_abs', 'weighted_sum', 'min', 'max'};
D.files     = rel;
D.fig_files = {};
if ~isempty(opts.fig_root)
    fr = opts.fig_root;
    if ~isfolder(fr) && isfolder(fullfile(project_root, fr)); fr = fullfile(project_root, fr); end
    if isfolder(fr)
        fr = char(java.io.File(fr).getCanonicalPath());
        fl = dir(fullfile(fr, '**', '*'));
        fl = fl(~[fl.isdir]);
        ff = cell(numel(fl), 1);
        for k = 1:numel(fl)
            full_path = fullfile(fl(k).folder, fl(k).name);
            ff{k} = strip_dates(strrep(full_path(numel(fr) + 2:end), '\', '/'));
        end
        D.fig_files = sort(ff);
    end
end
if opts.verbose
    fprintf('  [digest] %s: %d files (%d .mat), %d leaves, %d figure files\n', ...
        run_dir, numel(rel), nnz(is_mat), n_leaf, numel(D.fig_files));
end

%% ---- nested: the walk -------------------------------------------------
    function walk(key, v)
        if isnumeric(v) || islogical(v)
            add_num(key, v, size_text(v));
        elseif ischar(v)
            add_text(key, v);
        elseif isstring(v)
            add_text(key, char(strjoin(v(:)', '|')));
        elseif iscell(v)
            if numel(v) > STACK_MIN
                walk_many([key '{*}'], v(:)');
            else
                for i = 1:numel(v)
                    walk(sprintf('%s{%d}', key, i), v{i});
                end
            end
        elseif isstruct(v)
            f = fieldnames(v);
            if isscalar(v)
                for i = 1:numel(f)
                    walk([key '.' f{i}], v.(f{i}));
                end
            elseif numel(v) > STACK_MIN
                walk_many([key '(*)'], num2cell(v(:)'));
            else
                for j = 1:numel(v)
                    for i = 1:numel(f)
                        walk(sprintf('%s(%d).%s', key, j, f{i}), v(j).(f{i}));
                    end
                end
            end
        elseif istable(v)
            names = v.Properties.VariableNames;
            for i = 1:numel(names)
                walk([key '.' names{i}], v.(names{i}));
            end
        else
            add_text([key '<class>'], class(v));
        end
    end

    function walk_many(key, vals)
        % VALS: a 1 x m cell of the values one slot takes across a collection.
        m = numel(vals);
        if m == 0
            add_text([key '<empty>'], '');
            return;
        end
        is_num = cellfun(@(x) isnumeric(x) || islogical(x), vals);
        if all(is_num)
            cat_vals = cellfun(@(x) double(x(:)), vals, 'UniformOutput', false);
            add_num(key, vertcat(cat_vals{:}), sprintf('stack %d', m));
            return;
        end
        is_txt = cellfun(@(x) ischar(x) || isstring(x), vals);
        if all(is_txt)
            % a list of names: record the character codes, in order
            codes = cellfun(@(x) [double(char(strjoin(string(x(:)'), '|'))), 0]', ...
                vals, 'UniformOutput', false);
            add_num([key '<chars>'], vertcat(codes{:}), sprintf('stack %d text', m));
            return;
        end
        is_struct = cellfun(@(x) isstruct(x) && isscalar(x), vals);
        if all(is_struct)
            f = fieldnames(vals{1});
            same = all(cellfun(@(x) isequal(fieldnames(x), f), vals));
            if same
                for i = 1:numel(f)
                    walk_many([key '.' f{i}], cellfun(@(x) x.(f{i}), vals, 'UniformOutput', false));
                end
                return;
            end
        end
        % Empty placeholders ([]) among structs or cells: note how many, then
        % stack the rest.
        is_empty_num = cellfun(@(x) (isnumeric(x) || islogical(x)) && isempty(x), vals);
        if any(is_empty_num) && ~all(is_empty_num)
            add_num([key '<n_empty>'], find(is_empty_num(:)), sprintf('stack %d', m));
            walk_many(key, vals(~is_empty_num));
            return;
        end
        is_cell = cellfun(@iscell, vals);
        if all(is_cell)
            flat = cellfun(@(x) x(:)', vals, 'UniformOutput', false);
            add_num([key '<cell_counts>'], cellfun(@numel, flat(:)), sprintf('stack %d', m));
            walk_many([key '{*}'], [flat{:}]);
            return;
        end
        % Mixed content that does not stack: walk the first few, note the rest.
        classes = strjoin(unique(cellfun(@class, vals, 'UniformOutput', false)), ',');
        add_text([key '<mixed>'], sprintf('%d elements: %s', m, classes));
        for i = 1:min(m, STACK_MIN)
            walk(sprintf('%s[%d]', key, i), vals{i});
        end
    end

%% ---- nested: leaves ---------------------------------------------------
    function add_num(key, v, shp)
        x = full(double(v(:)));
        if ~isreal(x)
            x = [real(x); imag(x)];
            shp = [shp ' complex'];
        end
        n = numel(x);
        fin = isfinite(x);
        xs = x(fin);
        w = 1 + mod((1:n)', 1009) / 1009;       % position weights in [1, 2)
        ws = w(fin);
        if isempty(xs)
            row = [n, n, 0, 0, 0, NaN, NaN];
        else
            row = [n, n - numel(xs), sum(xs), sum(abs(xs)), sum(xs .* ws), min(xs), max(xs)];
        end
        add_row(key, 'num', shp, row, '', md5_hex(x));
    end

    function add_text(key, t)
        t = char(t);
        t = strrep(t, '\', '/');
        t = strrep(t, strrep(run_dir, '\', '/'), '<run_dir>');
        t = strrep(t, strrep(project_root, '\', '/'), '<root>');
        t = strip_dates(t);
        if numel(t) > TEXT_MAX
            c = double(t(:));
            t = sprintf('<%d chars, checksum %.0f>', numel(c), sum(c .* (1 + mod((1:numel(c))', 1009))));
        end
        add_row(key, 'txt', sprintf('1x%d', numel(t)), nan(1, 7), t, '');
    end

    function add_row(key, k, shp, row, t, h)
        n_leaf = n_leaf + 1;
        if n_leaf > cap
            cap = 2 * cap;
            keys{cap} = []; kind{cap} = []; shape{cap} = []; txts{cap} = []; hashes{cap} = [];
            stats(cap, 7) = NaN; stats(n_leaf:cap, :) = NaN;
        end
        keys{n_leaf} = key; kind{n_leaf} = k; shape{n_leaf} = shp;
        stats(n_leaf, :) = row; txts{n_leaf} = t; hashes{n_leaf} = h;
    end
end

%% ------------------------------------------------------------------------
function s = size_text(v)
s = strjoin(arrayfun(@(d) sprintf('%d', d), size(v), 'UniformOutput', false), 'x');
if issparse(v); s = [s ' sparse']; end
end

function h = md5_hex(x)
% MD5 of the array's bytes, so an exact comparison sees ANY changed bit. The
% statistics cannot: one element of a million changing in its last few digits
% moves the sum by less than the sum's own rounding.
if isempty(x); h = 'empty'; return; end
x(isnan(x)) = NaN;                           % one NaN bit pattern
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast(x, 'uint8'));
h = sprintf('%02x', typecast(md.digest(), 'uint8'));
end

function s = strip_dates(s)
% Launch-time suffixes the stages put in folder and file names:
%   _nov_10_26_19_05     (mon_dd_yy_HH_MM, the sweeps)
%   _20261110_190512     (yyyymmdd_HHMMSS)
%   2026-11-10 19:05:12, 10-Nov-2026 19:05:12   (inside text)
s = regexprep(s, '_(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)_\d{1,2}_\d{2}_\d{1,2}_\d{2}', '', 'ignorecase');
s = regexprep(s, '_?\d{8}[_T]\d{6}', '');
s = regexprep(s, '\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}(:\d{2}(\.\d+)?)?', '<time>');
s = regexprep(s, '\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d{4}( \d{2}:\d{2}:\d{2})?', '<time>');
end
