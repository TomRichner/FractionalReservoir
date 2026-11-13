function [S, files, SB, SG] = mc_gain_bias_statistics(cfg)
% MC_GAIN_BIAS_STATISTICS Robustness statistics of memory capacity over the gain-by-bias map (block C).
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
%                            ratio), at bias_fit (the preset's bias, nearest
%                            level)
%   gain_oct_best, _bias,    the widest along gain over all biases, and at
%   _open                    which bias
%   bias_range_presetgain,   along bias, in S_c units, at the preset's gain
%   _open                    (nearest level)
%   bias_range_best, _gain,  the widest along bias over all gains, and at
%   _open                    which gain
%
%   GAUSSIAN FIT along gain at bias_fit (Plan C decision 5):
%   fit_*                    a Gaussian in log10(gain) fitted to the cell
%                            medians, mc = h exp(-(x - c)^2 / (2 w^2)), no
%                            offset, least squares from several starts (the
%                            best of each of the three highest cells as
%                            centre x three widths). fit_height h,
%                            fit_center_gain 10^c, fit_width_log10 w,
%                            fit_fwhm_log10 = 2.355 w, fit_r2 = 1 -
%                            SS_res / SS_tot on the medians.
%   fwhm_meas_log10, _open   the MEASURED full width at half of the row's own
%                            maximum, same interpolation (open = lower bound)
%   width_used_oct,          the width to quote, in octaves: the fit's FWHM
%   width_source             when fit_r2 >= cfg.fit_r2_min (0.8), else the
%                            measured one ('fit' | 'measured'). A curve that
%                            is flat and then collapses (chained MTS) is
%                            not a Gaussian; its R^2 shows it and the
%                            measured width is used.
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
    cfg.bias_fit              = []      % empty: D.preset_bias
    cfg.reference  (1,:) char = 'sfa3_stdc2'
    cfg.fit_r2_min (1,1) double = 0.8
    cfg.out_dir    (1,:) char = ''
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
bias_fit = cfg.bias_fit;
if isempty(bias_fit); bias_fit = D.preset_bias; end
[~, ib_fit] = min(abs(D.biases - bias_fit));
[~, ig_pre] = min(abs(log(D.gains) - log(D.preset_gain)));
x_log = log10(D.gains);
x_oct = log2(D.gains);

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
        r.peak_gain = NaN; r.peak_bias = NaN;
    else
        [bk, gk] = ind2sub(size(m), k);
        r.peak_gain = D.gains(gk); r.peak_bias = D.biases(bk);
    end

    % Tuning width along gain (octaves), every bias
    oct = zeros(1, numel(D.biases)); oct_open = false(1, numel(D.biases));
    for b = 1:numel(D.biases)
        [w, lo, hi, op] = above_range(x_oct, m(b, :), half_any);
        oct(b) = w; oct_open(b) = op;
        [pk, kk] = max(m(b, :), [], 'omitnan');
        pg = NaN; if ~isempty(kk) && isfinite(pk); pg = D.gains(kk); end
        sb{end+1, 1} = struct('condition', conds{c}, 'bias', D.biases(b), 'row_peak_mc', pk, ...
            'row_peak_gain', pg, 'gain_lo', 2^lo, 'gain_hi', 2^hi, 'octaves', w, 'open', op); %#ok<AGROW>
    end
    r.gain_oct_fitbias = oct(ib_fit); r.gain_oct_fitbias_open = oct_open(ib_fit);
    [r.gain_oct_best, kb] = max(oct);
    r.gain_oct_best_bias = D.biases(kb); r.gain_oct_best_open = oct_open(kb);

    % Tuning width along bias (S_c), every gain
    rng_b = zeros(1, numel(D.gains)); rng_open = false(1, numel(D.gains));
    for g = 1:numel(D.gains)
        [w, lo, hi, op] = above_range(D.biases, m(:, g)', half_any);
        rng_b(g) = w; rng_open(g) = op;
        [pk, kk] = max(m(:, g), [], 'omitnan');
        pb = NaN; if ~isempty(kk) && isfinite(pk); pb = D.biases(kk); end
        sg{end+1, 1} = struct('condition', conds{c}, 'gain', D.gains(g), 'col_peak_mc', pk, ...
            'col_peak_bias', pb, 'bias_lo', lo, 'bias_hi', hi, 'range', w, 'open', op); %#ok<AGROW>
    end
    r.bias_range_presetgain = rng_b(ig_pre); r.bias_range_presetgain_open = rng_open(ig_pre);
    [r.bias_range_best, kg] = max(rng_b);
    r.bias_range_best_gain = D.gains(kg); r.bias_range_best_open = rng_open(kg);

    % Along gain at the fit bias: Gaussian fit and measured width
    y = m(ib_fit, :);
    r.fit_bias = D.biases(ib_fit);
    [r.fit_height, r.fit_center_gain, r.fit_width_log10, r.fit_r2] = gauss_fit(x_log, y);
    r.fit_fwhm_log10 = 2 * sqrt(2 * log(2)) * r.fit_width_log10;
    [r.fwhm_meas_log10, ~, ~, r.fwhm_open] = above_range(x_log, y, max(y, [], 'omitnan') / 2);
    if isfinite(r.fit_r2) && r.fit_r2 >= cfg.fit_r2_min
        r.width_used_oct = r.fit_fwhm_log10 * log2(10); r.width_source = 'fit';
    else
        r.width_used_oct = r.fwhm_meas_log10 * log2(10); r.width_source = 'measured';
    end

    % Paired with the reference, same gain, bias and seed
    r.paired_ref = cfg.reference;
    r.share_ref_higher = NaN; r.median_paired_diff = NaN; r.n_pairs = 0;
    if has_ref && ~strcmp(conds{c}, cfg.reference)
        A = T(ok, {'gain', 'bias', 'seed', 'mc_total'});
        B = Tref(~Tref.failed, {'gain', 'bias', 'seed', 'mc_total'});
        B = renamevars(B, 'mc_total', 'mc_ref');
        P = innerjoin(A, B, 'Keys', {'gain', 'bias', 'seed'});
        if height(P) > 0
            P.diff = P.mc_ref - P.mc_total;
            G = groupsummary(P, {'gain', 'bias'}, 'median', 'diff');
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
SB = struct2table([sb{:}], 'AsArray', true); SB.condition = cellstr(SB.condition);
SG = struct2table([sg{:}], 'AsArray', true); SG.condition = cellstr(SG.condition);

vprintf(cfg.verbose, 'minimal', '[mc_gain_bias_statistics] %d conditions, %d gains x %d biases x %d seeds; fit at bias %.3g; half the best of any = %.3g\n', ...
    nc, numel(D.gains), numel(D.biases), numel(D.seeds), D.biases(ib_fit), half_any);

files = {};
if ~isempty(cfg.out_dir)
    if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end
    f_csv = fullfile(cfg.out_dir, [cfg.tag '.csv']);
    writetable(S, f_csv);
    f_sb = fullfile(cfg.out_dir, [cfg.tag '_by_bias.csv']);
    writetable(SB, f_sb);
    f_sg = fullfile(cfg.out_dir, [cfg.tag '_by_gain.csv']);
    writetable(SG, f_sg);
    f_md = fullfile(cfg.out_dir, [cfg.tag '.md']);
    write_md(f_md, S, SB, SG, D, cfg, half_any, D.gains(ig_pre));
    files = {f_csv, f_md, f_sb, f_sg};
end
end

%% ------------------------------------------------------------------------
function v = med(x)
if isempty(x); v = NaN; else; v = median(x); end
end

function v = pct(x, p)
if isempty(x); v = NaN; else; v = prctile(x, p); end
end

function [h, c_gain, w, r2] = gauss_fit(x, y)
% Least squares on the finite medians; parameters (h, c, log w), the best
% of several starts (centre at each of the three highest cells, three
% widths), so a curve with a plateau does not trap the fit in a poor local
% minimum.
ok = isfinite(y);
h = NaN; c_gain = NaN; w = NaN; r2 = NaN;
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
h = p(1); c_gain = 10^p(2); w = exp(p(3));
ss_tot = sum((y - mean(y)).^2);
r2 = 1 - best / max(ss_tot, eps);
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

function write_md(file, S, SB, SG, D, cfg, half_any, gain_at_preset)
% The three tables as one markdown file. gain_at_preset: the gain level
% nearest the preset's, where the width along bias is quoted.
fid = fopen(file, 'w');
if fid < 0; warning('mc_gain_bias_statistics:Write', 'Cannot write %s', file); return; end
c = onCleanup(@() fclose(fid));
titles = titles_by_name(S.condition');
fprintf(fid, '# Memory capacity across gain and bias: robustness statistics\n\n');
fprintf(fid, 'Source: `%s`\n\n', D.source);
fprintf(fid, ['%d gains (%.3g to %.3g, level_of_chaos) x %d biases (S_c %.3g to %.3g) x %d seed(s). ' ...
    'Map statistics pool every grid point and seed; shares, peaks and widths use the median over seeds per grid point. ' ...
    'Half the best cell of any condition: %.3f. ' ...
    'Gaussian fit against log10 gain at bias %.3g (the nearest level to %.3g). ' ...
    'Paired differences are %s minus the row''s condition, same gain, bias and seed.\n\n'], ...
    numel(D.gains), min(D.gains), max(D.gains), numel(D.biases), min(D.biases), max(D.biases), ...
    numel(D.seeds), half_any, S.fit_bias(1), ifempty(cfg.bias_fit, D.preset_bias), cfg.reference);
fprintf(fid, '## The map\n\n');
fprintf(fid, '| Condition | jobs (failed) | MC median | IQR | 5th | 95th | share >= half own best | share >= half best of any | peak MC | at gain | at bias |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    fprintf(fid, '| %s | %d (%d) | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.3f | %.3g | %.3g |\n', ...
        titles{i}, S.n_jobs(i), S.n_failed(i), S.mc_median(i), S.mc_iqr(i), ...
        S.mc_p05(i), S.mc_p95(i), S.share_half_own(i), S.share_half_any(i), S.peak_mc(i), ...
        S.peak_gain(i), S.peak_bias(i));
end
fprintf(fid, ['\n## Tuning width: where MC >= %.3f (half the best of any condition)\n\n' ...
    'Contiguous stretch around the row''s (column''s) best cell; "+" marks a stretch that reaches the grid''s edge (a lower bound).\n\n'], half_any);
fprintf(fid, '| Condition | along gain at bias %.3g (octaves) | widest along gain (octaves) | at bias | along bias at gain %.3g (S_c) | widest along bias (S_c) | at gain |\n', ...
    S.fit_bias(1), gain_at_preset);
fprintf(fid, '|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    fprintf(fid, '| %s | %.2f%s | %.2f%s | %.3g | %.3f%s | %.3f%s | %.3g |\n', titles{i}, ...
        S.gain_oct_fitbias(i), open_mark(S.gain_oct_fitbias_open(i)), S.gain_oct_best(i), open_mark(S.gain_oct_best_open(i)), ...
        S.gain_oct_best_bias(i), S.bias_range_presetgain(i), open_mark(S.bias_range_presetgain_open(i)), ...
        S.bias_range_best(i), open_mark(S.bias_range_best_open(i)), S.bias_range_best_gain(i));
end
fprintf(fid, '\n### Along gain, every bias (octaves)\n\n| Condition |%s\n|---|%s\n', ...
    sprintf(' S_c %.3g |', D.biases), repmat('---|', 1, numel(D.biases)));
for i = 1:height(S)
    R = SB(strcmp(SB.condition, S.condition{i}), :);
    fprintf(fid, '| %s |%s\n', titles{i}, ...
        strjoin(arrayfun(@(k) sprintf(' %.2f%s |', R.octaves(k), open_mark(R.open(k))), 1:height(R), 'UniformOutput', false), ''));
end
fprintf(fid, '\n### Along bias, every gain (S_c units)\n\n| Condition |%s\n|---|%s\n', ...
    sprintf(' g %.3g |', D.gains), repmat('---|', 1, numel(D.gains)));
for i = 1:height(S)
    R = SG(strcmp(SG.condition, S.condition{i}), :);
    fprintf(fid, '| %s |%s\n', titles{i}, ...
        strjoin(arrayfun(@(k) sprintf(' %.3f%s |', R.range(k), open_mark(R.open(k))), 1:height(R), 'UniformOutput', false), ''));
end
fprintf(fid, ['\n## Width at half of its own maximum, along gain at bias %.3g\n\n' ...
    'The width quoted is the Gaussian fit''s when its R^2 >= %.2f, else the measured one.\n\n'], S.fit_bias(1), cfg.fit_r2_min);
fprintf(fid, '| Condition | fit height | fit centre (gain) | fit FWHM (log10 gain) | fit R^2 | measured FWHM (log10 gain) | width quoted (octaves) | from | share of grid points %s higher | median paired diff | pairs |\n', cfg.reference);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    open_txt = ''; if S.fwhm_open(i); open_txt = ' (open, lower bound)'; end
    fprintf(fid, '| %s | %.3f | %.3g | %.3f | %.3f | %.3f%s | %.2f | %s | %.2f | %+.3f | %d |\n', ...
        titles{i}, S.fit_height(i), S.fit_center_gain(i), S.fit_fwhm_log10(i), ...
        S.fit_r2(i), S.fwhm_meas_log10(i), open_txt, S.width_used_oct(i), S.width_source{i}, ...
        S.share_ref_higher(i), S.median_paired_diff(i), S.n_pairs(i));
end
end

function s = open_mark(tf)
if tf; s = '+'; else; s = ''; end
end

function v = ifempty(a, b)
if isempty(a); v = b; else; v = a; end
end
