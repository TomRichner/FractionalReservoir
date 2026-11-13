function [S, files] = mc_gain_bias_statistics(cfg)
% MC_GAIN_BIAS_STATISTICS Robustness statistics of memory capacity over the gain-by-bias map (block C).
%
%   S = MC_GAIN_BIAS_STATISTICS('data_file', f)
%   S = MC_GAIN_BIAS_STATISTICS('data', D_or_table, 'out_dir', d)
%   [S, files] = ...
%
% One row per condition (cfg.conditions through select_conditions; empty =
% every condition the result holds). The table is written as
% <out_dir>/<tag>.csv and <tag>.md when out_dir is given. The quantities
% (Plan C, section 5):
%
%   n_jobs, n_failed         rows of the condition; failed rows are skipped
%   mc_median, mc_iqr,       memory capacity over the WHOLE map: every grid
%   mc_p05, mc_p95           point and every seed (each job one value)
%   share_half_own           share of the map's grid points whose median over
%                            seeds is >= half of THIS condition's best cell
%   share_half_any           ... >= half of the best cell of ANY condition
%                            selected (the denominator is the grid points
%                            with at least one successful seed)
%   peak_mc, peak_gain,      the best cell (median over seeds) and where it is
%   peak_bias
%   fit_*                    along gain at bias_fit (default the preset's
%                            bias; the nearest level, reported as fit_bias):
%                            a Gaussian in log10(gain) fitted to the cell
%                            medians, mc = h exp(-(x - c)^2 / (2 w^2)), no
%                            offset. fit_height h, fit_center_gain 10^c,
%                            fit_width_log10 w, fit_fwhm_log10 = 2.355 w,
%                            fit_r2 (1 - SS_res / SS_tot on the medians: a
%                            poor fit shows here, not hidden)
%   fwhm_meas_log10          the MEASURED full width at half maximum of the
%                            same medians, linear interpolation in log10 gain,
%                            the contiguous stretch around the peak;
%   fwhm_open                true when that stretch reaches the grid's edge,
%                            so the measured width is a lower bound
%   paired_ref               the reference condition (default sfa3_stdc2, the
%                            chained MTS) of the paired comparison
%   share_ref_higher         share of grid points where the median over seeds
%                            of the paired difference (ref - this, same seed)
%                            is > 0; NaN on the reference's own row or when
%                            the reference is not in the result
%   median_paired_diff       median over all paired jobs of (ref - this)
%   n_pairs                  the number of paired jobs behind it
%
% See also: read_mc_gain_sweep, mc_gain_grid, fig_mc_robustness

arguments
    cfg.data_file  (1,:) char = ''
    cfg.data                  = []      % read_mc_gain_sweep output, or anything it reads
    cfg.conditions            = {}
    cfg.bias_fit              = []      % empty: D.preset_bias
    cfg.reference  (1,:) char = 'sfa3_stdc2'
    cfg.out_dir    (1,:) char = ''
    cfg.tag        (1,:) char = 'MC_gain_bias_statistics'
    cfg.verbose    (1,:) char = 'minimal'
end

if ~isempty(cfg.data)
    if isstruct(cfg.data) && isfield(cfg.data, 'gains'); D = cfg.data; else; D = read_mc_gain_sweep(cfg.data); end
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
x_log = log10(D.gains);

% Best cell of any selected condition
M = cell(1, nc);
for c = 1:nc; M{c} = mc_gain_grid(D, conds{c}, 'mc_total'); end
best_any = max(cellfun(@(m) max(m(:), [], 'omitnan'), M));

has_ref = ismember(cfg.reference, D.conditions);
Tref = D.T(strcmp(D.T.condition, cfg.reference), :);

rows = cell(nc, 1);
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
    r.share_half_any = nnz(m(valid) >= 0.5 * best_any) / max(1, nnz(valid));
    [r.peak_mc, k] = max(m(:), [], 'omitnan');
    if isempty(k) || isnan(r.peak_mc)
        r.peak_gain = NaN; r.peak_bias = NaN;
    else
        [bk, gk] = ind2sub(size(m), k);
        r.peak_gain = D.gains(gk); r.peak_bias = D.biases(bk);
    end

    % Along gain at the fit bias
    y = m(ib_fit, :);
    r.fit_bias = D.biases(ib_fit);
    [r.fit_height, r.fit_center_gain, r.fit_width_log10, r.fit_r2] = gauss_fit(x_log, y);
    r.fit_fwhm_log10 = 2 * sqrt(2 * log(2)) * r.fit_width_log10;
    [r.fwhm_meas_log10, r.fwhm_open] = measured_fwhm(x_log, y);

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
S.condition  = cellstr(S.condition);
S.paired_ref = cellstr(S.paired_ref);

vprintf(cfg.verbose, 'minimal', '[mc_gain_bias_statistics] %d conditions, %d gains x %d biases x %d seeds; fit at bias %.3g\n', ...
    nc, numel(D.gains), numel(D.biases), numel(D.seeds), D.biases(ib_fit));

files = {};
if ~isempty(cfg.out_dir)
    if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end
    f_csv = fullfile(cfg.out_dir, [cfg.tag '.csv']);
    writetable(S, f_csv);
    f_md = fullfile(cfg.out_dir, [cfg.tag '.md']);
    write_md(f_md, S, D, cfg);
    files = {f_csv, f_md};
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
% Least squares on the finite medians; parameters (h, c, log w).
ok = isfinite(y);
h = NaN; c_gain = NaN; w = NaN; r2 = NaN;
if nnz(ok) < 3; return; end
x = x(ok); y = y(ok);
[h0, k] = max(y);
p0 = [h0, x(k), log(max(0.1, 0.5 * (max(x) - min(x)) / 2))];
model = @(p, x) p(1) * exp(-(x - p(2)).^2 / (2 * exp(2 * p(3))));
cost  = @(p) sum((y - model(p, x)).^2);
p = fminsearch(cost, p0, optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000));
h = p(1); c_gain = 10^p(2); w = exp(p(3));
ss_tot = sum((y - mean(y)).^2);
r2 = 1 - cost(p) / max(ss_tot, eps);
end

function [fw, open] = measured_fwhm(x, y)
% Contiguous stretch around the peak with y >= max/2, crossings by linear
% interpolation in x. Open (a lower bound) when it reaches either edge.
fw = NaN; open = false;
ok = isfinite(y);
if nnz(ok) < 2; return; end
x = x(ok); y = y(ok);
[ym, k] = max(y); half = ym / 2;
i = k; while i > 1 && y(i - 1) >= half; i = i - 1; end
j = k; while j < numel(y) && y(j + 1) >= half; j = j + 1; end
if i == 1
    xl = x(1); open = true;
else
    xl = x(i - 1) + (half - y(i - 1)) * (x(i) - x(i - 1)) / (y(i) - y(i - 1));
end
if j == numel(y)
    xr = x(end); open = true;
else
    xr = x(j) + (y(j) - half) * (x(j + 1) - x(j)) / (y(j) - y(j + 1));
end
fw = xr - xl;
end

function write_md(file, S, D, cfg)
fid = fopen(file, 'w');
if fid < 0; warning('mc_gain_bias_statistics:Write', 'Cannot write %s', file); return; end
c = onCleanup(@() fclose(fid));
titles = srnn_condition_titles();
fprintf(fid, '# Memory capacity across gain and bias: robustness statistics\n\n');
fprintf(fid, 'Source: `%s`\n\n', D.source);
fprintf(fid, ['%d gains (%.3g to %.3g, level_of_chaos) x %d biases (S_c %.3g to %.3g) x %d seeds. ' ...
    'Map statistics pool every grid point and seed; shares and peaks use the median over seeds per grid point. ' ...
    'Gaussian fit against log10 gain at bias %.3g (the nearest level to %.3g). ' ...
    'Paired differences are %s minus the row''s condition, same gain, bias and seed.\n\n'], ...
    numel(D.gains), min(D.gains), max(D.gains), numel(D.biases), min(D.biases), max(D.biases), ...
    numel(D.seeds), S.fit_bias(1), ifempty(cfg.bias_fit, D.preset_bias), cfg.reference);
fprintf(fid, '| Condition | jobs (failed) | MC median | IQR | 5th | 95th | share >= half own best | share >= half best of any | peak MC | at gain | at bias |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    fprintf(fid, '| %s | %d (%d) | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.3f | %.3g | %.3g |\n', ...
        ttl(titles, S.condition{i}), S.n_jobs(i), S.n_failed(i), S.mc_median(i), S.mc_iqr(i), ...
        S.mc_p05(i), S.mc_p95(i), S.share_half_own(i), S.share_half_any(i), S.peak_mc(i), ...
        S.peak_gain(i), S.peak_bias(i));
end
fprintf(fid, '\n| Condition | fit height | fit centre (gain) | fit FWHM (log10 gain) | measured FWHM (log10 gain) | fit R^2 | share of grid points %s higher | median paired diff | pairs |\n', cfg.reference);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
for i = 1:height(S)
    open_txt = ''; if S.fwhm_open(i); open_txt = ' (open, lower bound)'; end
    fprintf(fid, '| %s | %.3f | %.3g | %.3f | %.3f%s | %.3f | %.2f | %+.3f | %d |\n', ...
        ttl(titles, S.condition{i}), S.fit_height(i), S.fit_center_gain(i), S.fit_fwhm_log10(i), ...
        S.fwhm_meas_log10(i), open_txt, S.fit_r2(i), S.share_ref_higher(i), ...
        S.median_paired_diff(i), S.n_pairs(i));
end
end

function t = ttl(titles, nm)
if isKey(titles, nm); t = titles(nm); else; t = nm; end
end

function v = ifempty(a, b)
if isempty(a); v = b; else; v = a; end
end
