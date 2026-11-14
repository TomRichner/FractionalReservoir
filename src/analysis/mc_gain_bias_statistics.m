function [S, files, SB, SG] = mc_gain_bias_statistics(cfg)
% MC_GAIN_BIAS_STATISTICS Robustness statistics of memory capacity over the map of gain by the sweep's second axis (block C).
%
%   S = MC_GAIN_BIAS_STATISTICS('data_file', f)
%   S = MC_GAIN_BIAS_STATISTICS('data', D_or_table, 'out_dir', d)
%   [S, files, SB, SG] = ...
%
% One row per condition (cfg.conditions through select_conditions; empty =
% every condition the result holds). Every per-cell value is the median
% over seeds of that grid point (mc_gain_grid), except the map
% distribution, which pools every job. With out_dir given, writes
% <tag>.csv (S), <tag>_by_bias.csv (SB), <tag>_by_gain.csv (SG) and one
% <tag>.md holding all three.
%
% THE SECOND AXIS IS THE RESULT'S (D.second_axis): the bias S_c, or the
% property the sweep ran, for example the self-feedback k. Everything below
% is written for the bias, and holds for any axis with its column name in
% place of "bias": the columns are then peak_self_feedback,
% gain_oct_fitself_feedback, self_feedback_range_best, ..., SB's key
% column is self_feedback and its file <tag>_by_self_feedback.csv, and the
% markdown names the axis and its unit in mc_sweep_quantity's wording. The
% code computes with axis-neutral names and renames the columns once, at
% the end (name_columns).
%
% THE FIRST AXIS IS THE RESULT'S TOO (D.first_axis): the gain, or another
% property, for example the self-coupling s. Below it is written for the
% gain, and "gain" in a column name stands for the first axis's column
% (peak_self_coupling, self_coupling_oct_fit..., SG's key column and the
% file <tag>_by_self_coupling.csv). Widths and fits along the first axis
% are taken logarithmically when every level is positive (the gain on any
% grid): in octaves and log10 units, as below. On an axis with a level at
% or below zero (a self-coupling that changes sign, a noise level from 0)
% they are in the axis's own units, and the columns say so: "oct" becomes "range"
% (self_coupling_range_best, width_used_range), SB's octaves column
% becomes range, and the _log10 suffix is dropped (fit_width, fit_fwhm,
% fwhm_meas); the Gaussian is then fitted against the axis itself.
%
% S, per condition (Plan C, section 5):
%   n_jobs, n_failed         rows of the condition; failed rows are skipped
%   mc_median, mc_iqr,       memory capacity over the WHOLE map: every grid
%   mc_p05, mc_p95           point and every seed (each job one value)
%   share_half_own           share of the grid points whose median is >= half
%                            of THIS condition's best cell
%   share_half_any           ... >= half of the best cell of ANY condition
%                            selected (best_any, the same for every row;
%                            the denominator is the grid points with at
%                            least one successful seed)
%   peak_mc, peak_gain,      the best cell and where it is
%   peak_bias
%
%   TUNING WIDTH, the answer to "does no adaptation need precise tuning":
%   the range over which memory capacity is >= half_any = best_any / 2 (half
%   the best ANY condition reaches, so one bar for every condition). Along
%   a row or column of the grid it is the contiguous stretch around that
%   row's best cell, its two ends by linear interpolation of the medians
%   (in log2 gain along gain, in S_c along bias); a stretch that reaches the
%   grid's edge is open (the width is then a lower bound); a row that never
%   reaches half_any has width 0. SB and SG hold every row and column.
%   gain_oct_fitbias, _open  along gain, in OCTAVES (log2 of the gain
%                            ratio), at cfg.second_fit (default the preset's
%                            bias; the nearest level)
%   gain_oct_best, _bias,    the widest along gain over all biases, and at
%   _open                    which bias
%   bias_range_presetgain,   along bias, in S_c units, at the preset's gain
%   _open                    (nearest level)
%   bias_range_best, _gain,  the widest along bias over all gains, and at
%   _open                    which gain
%
%   GAUSSIAN FIT along gain at second_fit (Plan C decision 5):
%   fit_*                    a Gaussian in log10(gain) fitted to the cell
%                            medians, mc = h exp(-(x - c)^2 / (2 w^2)), no
%                            offset, least squares from several starts (the
%                            best of each of the three highest cells as
%                            centre x three widths). fit_height h,
%                            fit_center_gain 10^c, fit_width_log10 w,
%                            fit_fwhm_log10 = 2.355 w, fit_r2 = 1 -
%                            SS_res / SS_tot on the medians.
%   fit_status,              'accepted' | 'rejected' | 'none' (fewer than three
%   fit_reject_reason        finite cells: nothing fitted). A fit is accepted
%                            only when it describes a peak the data hold:
%                              R^2 >= cfg.fit_r2_min (0.8);
%                              its center lies inside the gains of the row;
%                              its height is positive and at most
%                              cfg.fit_height_max (2) times the row's largest
%                              cell;
%                              its FWHM is at most the span of the row's gains
%                              (in log10 gain).
%                            R^2 alone does not do this. A row that only
%                            falls (its best cell at the lowest gain, then a
%                            collapse) is matched closely by the flank of a
%                            Gaussian centered far outside the grid: in the
%                            fast run of 2026-11-14 the no-adaptation row at
%                            self-coupling 0 gave height 248150, center
%                            5e-08, width 10.5 octaves at R^2 0.997, and two
%                            rows with a peak between two gain levels gave
%                            heights 179 and 100 where the largest cells
%                            were 15.4 and 13.6. The fitted values are kept
%                            in the fit_* columns as fitted; the reason (every
%                            rule that failed) says why they are not used.
%   fwhm_meas_log10, _open   the MEASURED full width at half of the row's own
%                            maximum, same interpolation (open = lower bound)
%   width_used_oct,          the width to quote, in octaves: the fit's FWHM
%   width_source             when the fit is accepted, else the measured one
%                            ('fit' | 'measured'). A curve that is flat and
%                            then collapses (chained MTS) is not a Gaussian;
%                            the fit is rejected and the measured width is
%                            used.
%
%   PAIRED with the reference (default sfa3_stdc2, the chained MTS):
%   share_ref_higher         share of grid points where the median over seeds
%                            of (ref - this, same gain, bias, seed) is > 0;
%                            NaN on the reference's own row or when the
%                            reference is not in the result
%   median_paired_diff       median over all paired jobs of (ref - this)
%   n_pairs                  the number of paired jobs behind it
%
% SB, one row per (condition, bias): row_peak_mc, row_peak_gain, gain_lo,
% gain_hi, octaves, open (along gain, at half_any).
% SG, one row per (condition, gain): col_peak_mc, col_peak_bias, bias_lo,
% bias_hi, range, open (along bias, at half_any).
%
% See also: read_mc_gain_sweep, mc_gain_grid, fig_mc_robustness

arguments
    cfg.data_file  (1,:) char = ''
    cfg.data                  = []      % read_mc_gain_sweep output, or anything it reads
    cfg.conditions            = {}
    cfg.second_fit            = []      % level of the second axis for the fit; empty: D.preset_second
    cfg.reference  (1,:) char = 'sfa3_stdc2'
    cfg.fit_r2_min (1,1) double = 0.8
    cfg.fit_height_max (1,1) double {mustBePositive} = 2   % the fit's height, at most this many times the row's largest cell
    cfg.out_dir   (1,:) char = ''
    cfg.tag        (1,:) char = 'MC_gain_bias_statistics'
    cfg.verbose    (1,:) char = 'minimal'
end

if ~isempty(cfg.data)
    D = read_mc_gain_sweep(cfg.data);
elseif ~isempty(cfg.data_file)
    D = read_mc_gain_sweep(cfg.data_file);
else
    error('mc_gain_bias_statistics:NoData', 'Give data_file or data.');
end
conds = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
levels1 = D.first_levels;                       % the first axis: its levels,
column1 = D.first_axis.column;                  % its column in the table,
Q1 = mc_sweep_quantity(column1, D);             % and its wording
levels2 = D.second_levels;                      % the same for the second axis
column2 = D.second_axis.column;
Q2 = mc_sweep_quantity(column2, D);
second_fit = cfg.second_fit;
if isempty(second_fit); second_fit = D.preset_second; end
[~, ib_fit] = min(abs(levels2 - second_fit));
% Along the first axis, widths and the fit are taken as the axis is drawn:
% for a log axis widths in octaves (log2) and the fit in log10, for a
% linear one both in the axis's own units. to_w / from_w map a level to
% the width coordinate and back, to_f / from_f to the fit coordinate.
X = first_axis_coordinates(levels1);
[~, ig_pre] = min(abs(X.to_w(levels1) - X.to_w(D.preset_first)));
x_fit = X.to_f(levels1);
x_w = X.to_w(levels1);

% Best cell of any selected condition
M = cell(1, nc);
for c = 1:nc; M{c} = mc_gain_grid(D, conds{c}, 'mc_total'); end
best_any = max(cellfun(@(m) max([m(:); -Inf], [], 'omitnan'), M));
half_any = best_any / 2;

has_ref = ismember(cfg.reference, D.conditions);
Tref = D.T(strcmp(D.T.condition, cfg.reference), :);

rows = cell(nc, 1); sb = {}; sg = {};
for c = 1:nc
    T = D.T(strcmp(D.T.condition, conds{c}), :);
    ok = ~T.failed;
    v = T.mc_total(ok);
    m = M{c};
    valid = isfinite(m);
    r = struct();
    r.condition = conds{c};
    r.n_jobs    = height(T);
    r.n_failed  = nnz(T.failed);
    r.mc_median = med(v);
    r.mc_iqr    = pct(v, 75) - pct(v, 25);
    r.mc_p05    = pct(v, 5);
    r.mc_p95    = pct(v, 95);
    best_own = max(m(:), [], 'omitnan');
    r.share_half_own = nnz(m(valid) >= 0.5 * best_own) / max(1, nnz(valid));
    r.share_half_any = nnz(m(valid) >= half_any) / max(1, nnz(valid));
    r.best_any = best_any;
    [r.peak_mc, k] = max(m(:), [], 'omitnan');
    if isempty(k) || isnan(r.peak_mc)
        r.peak_first = NaN; r.peak_second = NaN;
    else
        [bk, gk] = ind2sub(size(m), k);
        r.peak_first = levels1(gk); r.peak_second = levels2(bk);
    end

    % Tuning width along the first axis (octaves, or its own units), every
    % level of the second axis
    oct = zeros(1, numel(levels2)); oct_open = false(1, numel(levels2));
    for b = 1:numel(levels2)
        [w, lo, hi, op] = above_range(x_w, m(b, :), half_any);
        oct(b) = w; oct_open(b) = op;
        [pk, kk] = max(m(b, :), [], 'omitnan');
        pg = NaN; if ~isempty(kk) && isfinite(pk); pg = levels1(kk); end
        sb{end+1, 1} = struct('condition', conds{c}, 'second', levels2(b), 'row_peak_mc', pk, ...
            'row_peak_first', pg, 'first_lo', X.from_w(lo), 'first_hi', X.from_w(hi), ...
            'octaves', w, 'open', op); %#ok<AGROW>
    end
    r.first_oct_fitsecond = oct(ib_fit); r.first_oct_fitsecond_open = oct_open(ib_fit);
    [r.first_oct_best, kb] = max(oct);
    r.first_oct_best_second = levels2(kb); r.first_oct_best_open = oct_open(kb);

    % Tuning width along the second axis (in its own units), every level of
    % the first
    rng_b = zeros(1, numel(levels1)); rng_open = false(1, numel(levels1));
    for g = 1:numel(levels1)
        [w, lo, hi, op] = above_range(levels2, m(:, g)', half_any);
        rng_b(g) = w; rng_open(g) = op;
        [pk, kk] = max(m(:, g), [], 'omitnan');
        pb = NaN; if ~isempty(kk) && isfinite(pk); pb = levels2(kk); end
        sg{end+1, 1} = struct('condition', conds{c}, 'first', levels1(g), 'col_peak_mc', pk, ...
            'col_peak_second', pb, 'second_lo', lo, 'second_hi', hi, 'range', w, 'open', op); %#ok<AGROW>
    end
    r.second_range_presetfirst = rng_b(ig_pre); r.second_range_presetfirst_open = rng_open(ig_pre);
    [r.second_range_best, kg] = max(rng_b);
    r.second_range_best_first = levels1(kg); r.second_range_best_open = rng_open(kg);

    % Along the first axis at the fit level: Gaussian fit and measured width
    y = m(ib_fit, :);
    r.fit_second = levels2(ib_fit);
    [r.fit_height, r.fit_center_first, r.fit_width_log10, r.fit_r2] = gauss_fit(x_fit, y);
    r.fit_center_first = X.from_f(r.fit_center_first);   % fitted in the fit coordinate
    r.fit_fwhm_log10 = 2 * sqrt(2 * log(2)) * r.fit_width_log10;
    [r.fwhm_meas_log10, ~, ~, r.fwhm_open] = above_range(x_fit, y, max(y, [], 'omitnan') / 2);
    [r.fit_status, r.fit_reject_reason] = judge_fit(r, x_fit, y, cfg, X, Q1);
    if strcmp(r.fit_status, 'accepted')
        r.width_used_oct = r.fit_fwhm_log10 * X.f_to_w; r.width_source = 'fit';
    else
        r.width_used_oct = r.fwhm_meas_log10 * X.f_to_w; r.width_source = 'measured';
    end

    % Paired with the reference, same grid point and seed
    r.paired_ref = cfg.reference;
    r.share_ref_higher = NaN; r.median_paired_diff = NaN; r.n_pairs = 0;
    if has_ref && ~strcmp(conds{c}, cfg.reference)
        A = T(ok, {column1, column2, 'seed', 'mc_total'});
        B = Tref(~Tref.failed, {column1, column2, 'seed', 'mc_total'});
        B = renamevars(B, 'mc_total', 'mc_ref');
        P = innerjoin(A, B, 'Keys', {column1, column2, 'seed'});
        if height(P) > 0
            P.diff = P.mc_ref - P.mc_total;
            G = groupsummary(P, {column1, column2}, 'median', 'diff');
            r.share_ref_higher   = mean(G.median_diff > 0);
            r.median_paired_diff = median(P.diff);
            r.n_pairs            = height(P);
        end
    end
    rows{c} = r;
end
S = struct2table([rows{:}], 'AsArray', true);
S.condition    = cellstr(S.condition);
S.paired_ref   = cellstr(S.paired_ref);
S.width_source = cellstr(S.width_source);
S.fit_status   = cellstr(S.fit_status);
S.fit_reject_reason = cellstr(S.fit_reject_reason);
SB = struct2table([sb{:}], 'AsArray', true); SB.condition = cellstr(SB.condition);
SG = struct2table([sg{:}], 'AsArray', true); SG.condition = cellstr(SG.condition);

vprintf(cfg.verbose, 'minimal', '[mc_gain_bias_statistics] %d conditions, %d %s x %d %s x %d seeds; fit at %s; half the best of any = %.3g\n', ...
    nc, numel(levels1), Q1.plural, numel(levels2), Q2.plural, numel(D.seeds), Q2.at(levels2(ib_fit)), half_any);

% The markdown is written from the tables as computed; the tables then
% leave, and are saved, under the axes' own column names.
files = {};
write_files = ~isempty(cfg.out_dir);
if write_files
    if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end
    f_md = fullfile(cfg.out_dir, [cfg.tag '.md']);
    write_md(f_md, S, SB, SG, D, Q1, Q2, X, second_fit, cfg, half_any, levels1(ig_pre));
end
S  = name_columns(S, column1, column2, X.is_log);
SB = name_columns(SB, column1, column2, X.is_log);
SG = name_columns(SG, column1, column2, X.is_log);
if write_files
    f_csv = fullfile(cfg.out_dir, [cfg.tag '.csv']);
    writetable(S, f_csv);
    f_sb = fullfile(cfg.out_dir, [cfg.tag '_by_' column2 '.csv']);
    writetable(SB, f_sb);
    f_sg = fullfile(cfg.out_dir, [cfg.tag '_by_' column1 '.csv']);
    writetable(SG, f_sg);
    files = {f_csv, f_md, f_sb, f_sg};
end
end

%% ------------------------------------------------------------------------
function X = first_axis_coordinates(levels)
% The coordinates along the first axis in which widths and the Gaussian
% fit are taken: logarithmic when every level is positive (the gain,
% whether its grid is log- or evenly spaced, as these statistics always
% took it), the axis's own units when a level is zero or negative (a
% self-coupling that changes sign, a noise axis from 0). to_w / from_w for
% widths (log2, octaves, when logarithmic), to_f / from_f for the fit
% (log10), f_to_w the factor from a fit-coordinate width to a
% width-coordinate one (log2(10)).
X.is_log = all(levels > 0);
if X.is_log
    X.to_w = @log2;  X.from_w = @(u) 2.^u;
    X.to_f = @log10; X.from_f = @(u) 10.^u;
    X.f_to_w = log2(10);
else
    X.to_w = @(v) v; X.from_w = @(u) u;
    X.to_f = @(v) v; X.from_f = @(u) u;
    X.f_to_w = 1;
end
end

function T = name_columns(T, column1, column2, is_log)
% The statistics are computed with the axes called "first" and "second"
% (peak_first, peak_second, second_range_best, ...). A table leaves with
% the axes' own column names in their place: peak_gain and peak_bias for
% the default axes, as these columns were always called, peak_self_coupling
% for the self-coupling. On a linear first axis the units in the names
% follow the axis: "oct" becomes "range", octaves becomes range, and the
% _log10 suffix goes.
names = T.Properties.VariableNames;
names = regexprep(names, '(^|_|fit)second(_|$)', ['$1' column2 '$2']);
names = regexprep(names, '(^|_|preset)first(_|$)', ['$1' column1 '$2']);
if ~is_log
    names = regexprep(names, '(^|_)oct(_|$)', '$1range$2');
    names = regexprep(names, '^octaves$', 'range');
    names = regexprep(names, '_log10$', '');
end
T.Properties.VariableNames = names;
end

%% ------------------------------------------------------------------------
function v = med(x)
if isempty(x); v = NaN; else; v = median(x); end
end

function v = pct(x, p)
if isempty(x); v = NaN; else; v = prctile(x, p); end
end

function [h, c, w, r2] = gauss_fit(x, y)
% Least squares on the finite medians; parameters (h, c, log w), the best
% of several starts (centre at each of the three highest cells, three
% widths), so a curve with a plateau does not trap the fit in a poor local
% minimum. x is the fit coordinate (log10 gain); c and w are in it.
ok = isfinite(y);
h = NaN; c = NaN; w = NaN; r2 = NaN;
if nnz(ok) < 3; return; end
x = x(ok); y = y(ok);
model = @(p, x) p(1) * exp(-(x - p(2)).^2 / (2 * exp(2 * p(3))));
cost  = @(p) sum((y - model(p, x)).^2);
[~, order] = sort(y, 'descend');
span = max(x) - min(x);
best = Inf; p = [];
for k = order(1:min(3, numel(order)))
    for w0 = [0.1 0.25 0.5] * max(span, eps)
        q = fminsearch(cost, [y(k), x(k), log(w0)], ...
            optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000));
        if cost(q) < best; best = cost(q); p = q; end
    end
end
h = p(1); c = p(2); w = exp(p(3));
ss_tot = sum((y - mean(y)).^2);
r2 = 1 - best / max(ss_tot, eps);
end

function [status, reason] = judge_fit(r, x, y, cfg, X, Q1)
% Whether the Gaussian fitted to one row describes a peak the row holds.
% r: the row's fit_* values; x, y: the fit coordinate of the first axis
% (log10 gain) and the cell medians it was fitted to; X, Q1: the axis's
% coordinates and wording. status 'accepted' | 'rejected' | 'none'
% (nothing was fitted); reason: every rule that failed, joined by '; ',
% empty when accepted.
% The rules and why R^2 alone is not enough: the header, under fit_status.
if ~isfinite(r.fit_r2)
    status = 'none'; reason = 'fewer than three finite cells'; return;
end
ok = isfinite(y);
first_lo = X.from_f(min(x(ok))); first_hi = X.from_f(max(x(ok)));
span = max(x(ok)) - min(x(ok));
y_max = max(y(ok));
plural = strrep(Q1.plural, '\_', '_');
span_unit = '';
if X.is_log; span_unit = ' (log10)'; end
why = {};
if r.fit_r2 < cfg.fit_r2_min
    why{end+1} = sprintf('R^2 %.2f < %.2f', r.fit_r2, cfg.fit_r2_min);
end
if ~(r.fit_center_first >= first_lo && r.fit_center_first <= first_hi)
    why{end+1} = sprintf('center %.3g outside the %s %.3g to %.3g', r.fit_center_first, plural, first_lo, first_hi);
end
if ~(r.fit_height > 0 && r.fit_height <= cfg.fit_height_max * y_max)
    why{end+1} = sprintf('height %.3g not within 0 to %.3g x the largest cell %.3g', r.fit_height, cfg.fit_height_max, y_max);
end
if ~(r.fit_fwhm_log10 <= span)
    why{end+1} = sprintf('FWHM %.3g > the span of the %s %.3g%s', r.fit_fwhm_log10, plural, span, span_unit);
end
reason = strjoin(why, '; ');
if isempty(why); status = 'accepted'; else; status = 'rejected'; end
end

function [w, lo, hi, open] = above_range(x, y, thr)
% The contiguous stretch around the maximum of y with y >= thr, its ends by
% linear interpolation in x. open: it reaches an edge (a lower bound).
% No value reaches thr: w = 0, lo = hi = NaN.
w = 0; lo = NaN; hi = NaN; open = false;
ok = isfinite(y);
if ~any(ok) || ~isfinite(thr); return; end
x = x(ok); y = y(ok);
[ym, k] = max(y);
if ym < thr; return; end
i = k; while i > 1 && y(i - 1) >= thr; i = i - 1; end
j = k; while j < numel(y) && y(j + 1) >= thr; j = j + 1; end
if i == 1
    lo = x(1); open = true;
else
    lo = x(i - 1) + (thr - y(i - 1)) * (x(i) - x(i - 1)) / (y(i) - y(i - 1));
end
if j == numel(y)
    hi = x(end); open = true;
else
    hi = x(j) + (y(j) - thr) * (x(j + 1) - x(j)) / (y(j) - y(j + 1));
end
w = hi - lo;
end

function write_md(file, S, SB, SG, D, Q1, Q2, X, second_fit, cfg, half_any, first_at_preset)
% The three tables as one markdown file, from the tables under their
% axis-neutral column names. Q1, Q2: the axes' wording (mc_sweep_quantity);
% X: the first axis's coordinates (first_axis_coordinates). second_fit:
% the level asked for the fit. first_at_preset: the level of the first
% axis nearest the preset's, where the width along the second axis is
% quoted.
fid = fopen(file, 'w');
if fid < 0; warning('mc_gain_bias_statistics:Write', 'Cannot write %s', file); return; end
c = onCleanup(@() fclose(fid));
titles = titles_by_name(S.condition');
levels1 = D.first_levels;
levels2 = D.second_levels;
name1   = plain(Q1.symbol);                      % 'gain', 'self-coupling'
short1  = plain(Q1.variable);                    % its symbol in a column head: 'g', 's'
if isempty(short1); short1 = D.first_axis.file_token; end
plural1 = plain(Q1.plural);                      % 'gains'
if X.is_log                                      % how widths and the fit are measured along it
    width1 = 'octaves';
    fit1   = ['log10 ' name1];
else
    width1 = plain(Q1.unit);
    if isempty(width1); width1 = [name1 ' units']; end
    fit1   = name1;
end
axis_name = plain(Q2.symbol);                    % 'bias', 'self-feedback'
axis_var  = plain(Q2.variable);                  % 'S_c', 'k'
axis_unit = plain(Q2.unit);                      % '', 'R_W'
if isempty(axis_var); axis_var = D.second_axis.file_token; end
if isempty(axis_unit); axis_unit = axis_var; end % the bias is quoted in S_c
fprintf(fid, '# Memory capacity across %s and %s: robustness statistics\n\n', name1, axis_name);
fprintf(fid, 'Source: `%s`\n\n', D.source);
fprintf(fid, ['%d %s (%.3g to %.3g, %s) x %d %s (%s %.3g to %.3g) x %d seed(s). ' ...
    'Map statistics pool every grid point and seed; shares, peaks and widths use the median over seeds per grid point. ' ...
    'Half the best cell of any condition: %.3f. ' ...
    'Gaussian fit against %s at %s %.3g (the nearest level to %.3g). ' ...
    'Paired differences are %s minus the row''s condition, same %s, %s and seed.\n\n'], ...
    numel(levels1), plural1, min(levels1), max(levels1), D.first_axis.property, ...
    numel(levels2), plain(Q2.plural), axis_var, min(levels2), max(levels2), ...
    numel(D.seeds), half_any, fit1, axis_name, S.fit_second(1), second_fit, cfg.reference, name1, axis_name);
fprintf(fid, '## The map\n\n');
fprintf(fid, '| Condition | jobs (failed) | MC median | IQR | 5th | 95th | share >= half own best | share >= half best of any | peak MC | at %s | at %s |\n', name1, axis_name);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    fprintf(fid, '| %s | %d (%d) | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.3f | %.3g | %.3g |\n', ...
        titles{i}, S.n_jobs(i), S.n_failed(i), S.mc_median(i), S.mc_iqr(i), ...
        S.mc_p05(i), S.mc_p95(i), S.share_half_own(i), S.share_half_any(i), S.peak_mc(i), ...
        S.peak_first(i), S.peak_second(i));
end
fprintf(fid, ['\n## Tuning width: where MC >= %.3f (half the best of any condition)\n\n' ...
    'Contiguous stretch around the row''s (column''s) best cell; "+" marks a stretch that reaches the grid''s edge (a lower bound).\n\n'], half_any);
fprintf(fid, '| Condition | along %s at %s %.3g (%s) | widest along %s (%s) | at %s | along %s at %s %.3g (%s) | widest along %s (%s) | at %s |\n', ...
    name1, axis_name, S.fit_second(1), width1, name1, width1, axis_name, axis_name, name1, first_at_preset, ...
    axis_unit, axis_name, axis_unit, name1);
fprintf(fid, '|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    fprintf(fid, '| %s | %.2f%s | %.2f%s | %.3g | %.3f%s | %.3f%s | %.3g |\n', titles{i}, ...
        S.first_oct_fitsecond(i), open_mark(S.first_oct_fitsecond_open(i)), S.first_oct_best(i), open_mark(S.first_oct_best_open(i)), ...
        S.first_oct_best_second(i), S.second_range_presetfirst(i), open_mark(S.second_range_presetfirst_open(i)), ...
        S.second_range_best(i), open_mark(S.second_range_best_open(i)), S.second_range_best_first(i));
end
fprintf(fid, '\n### Along %s, every %s (%s)\n\n| Condition |%s\n|---|%s\n', name1, axis_name, width1, ...
    sprintf([' ' axis_var ' %.3g |'], levels2), repmat('---|', 1, numel(levels2)));
for i = 1:height(S)
    R = SB(strcmp(SB.condition, S.condition{i}), :);
    fprintf(fid, '| %s |%s\n', titles{i}, ...
        strjoin(arrayfun(@(k) sprintf(' %.2f%s |', R.octaves(k), open_mark(R.open(k))), 1:height(R), 'UniformOutput', false), ''));
end
fprintf(fid, '\n### Along %s, every %s (%s units)\n\n| Condition |%s\n|---|%s\n', axis_name, name1, axis_unit, ...
    sprintf([' ' short1 ' %.3g |'], levels1), repmat('---|', 1, numel(levels1)));
for i = 1:height(S)
    R = SG(strcmp(SG.condition, S.condition{i}), :);
    fprintf(fid, '| %s |%s\n', titles{i}, ...
        strjoin(arrayfun(@(k) sprintf(' %.3f%s |', R.range(k), open_mark(R.open(k))), 1:height(R), 'UniformOutput', false), ''));
end
fprintf(fid, ['\n## Width at half of its own maximum, along %s at %s %.3g\n\n' ...
    'The width quoted is the Gaussian fit''s when the fit is accepted, else the measured one. ' ...
    'A fit is accepted when its R^2 >= %.2f, its center lies inside the %s of the row, ' ...
    'its height is positive and at most %.3g times the row''s largest cell, and its FWHM is at most the span of the %s; ' ...
    'a rejected fit''s values are listed as fitted, with every rule it failed.\n\n'], ...
    name1, axis_name, S.fit_second(1), cfg.fit_r2_min, plural1, cfg.fit_height_max, plural1);
fprintf(fid, '| Condition | fit height | fit center (%s) | fit FWHM (%s) | fit R^2 | fit | measured FWHM (%s) | width quoted (%s) | from | share of grid points %s higher | median paired diff | pairs |\n', ...
    name1, fit1, fit1, width1, cfg.reference);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    open_txt = ''; if S.fwhm_open(i); open_txt = ' (open, lower bound)'; end
    fit_txt = S.fit_status{i};
    if ~isempty(S.fit_reject_reason{i}); fit_txt = sprintf('%s: %s', fit_txt, S.fit_reject_reason{i}); end
    fprintf(fid, '| %s | %.3g | %.3g | %.3f | %.3f | %s | %.3f%s | %.2f | %s | %.2f | %+.3f | %d |\n', ...
        titles{i}, S.fit_height(i), S.fit_center_first(i), S.fit_fwhm_log10(i), ...
        S.fit_r2(i), fit_txt, S.fwhm_meas_log10(i), open_txt, S.width_used_oct(i), S.width_source{i}, ...
        S.share_ref_higher(i), S.median_paired_diff(i), S.n_pairs(i));
end
end

function s = open_mark(tf)
if tf; s = '+'; else; s = ''; end
end

function s = plain(s)
% A label written for the tex interpreter, as plain text for markdown.
s = strrep(s, '\_', '_');
end
