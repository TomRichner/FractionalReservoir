% OW_EXTRACT Keep, from every run of the copied fastmedium outlier_spectra stage, only what the window study needs.
%
% Reads data/explorations/outlier_window/fastmedium_copy/outlier_spectra
% one run file at a time (ospec_read_index, ospec_load_run) and writes
% data/explorations/outlier_window/ow_series.mat with a struct array X, one
% element per index row:
%   network, condition, seed, window_s, dt, N, seconds (R.seconds: sim,
%   fast, null, full, dense, all), seconds_index (the index's)
%   t, rate_top                 the frozen-gain rightmost series
%   full_t, full_re, full_dt    real part of the full Jacobian's rightmost
%                               network mode (NaN: none found, below zero)
%   full_sec                    wall time of each full-pass state (s)
%   pre                         window-wide quantities of the paper figures
%                               on the first L s, L = lengths (20 50 100 all)
%
% Nothing under the main checkout is read here; nothing is simulated.

setup_paths();
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % the worktree
src = fullfile(root, 'data', 'explorations', 'outlier_window', 'fastmedium_copy', 'outlier_spectra');
out = fullfile(root, 'data', 'explorations', 'outlier_window', 'ow_series.mat');
lengths = [20 50 100 200];
[runs, settings] = ospec_read_index(src);
n_runs = height(runs);
X = struct([]);
t0 = tic;
for k = 1:n_runs
    if ~runs.ok(k); continue; end
    R = ospec_load_run(src, runs(k, :));
    x = struct();
    x.network = R.network; x.condition = R.condition; x.seed = R.seed;
    x.window_s = R.window_s; x.dt = R.dt; x.N = R.N; x.n = R.n;
    x.seconds = R.seconds; x.seconds_index = runs.seconds(k);
    x.t = R.t; x.rate_top = R.rate_top;
    x.full_t = []; x.full_re = []; x.full_dt = NaN; x.full_sec = [];
    if R.has_full_series && ~isempty(R.full)
        x.full_t = R.full.t; x.full_re = real(R.full.net); x.full_dt = R.full.dt;
        x.full_sec = R.full.seconds;
    end
    % Window-wide quantities drawn in the paper (ospec_measures' row), on prefixes
    E = 1; I = min(2, size(R.dbar, 1));
    rt = R.ratio_theory;
    gr = squeeze(R.dbar(E, E, :))' ./ squeeze(R.dbar(E, I, :))';
    P = struct('L', {}, 'ratio_med', {}, 'ratio_p05', {}, 'ratio_p95', {}, 'frac_on_E', {}, ...
        'frac_on_I', {}, 'gain_ratio_EI', {}, 'rate_top_med', {}, 'frac_unstable', {});
    for j = 1:numel(lengths)
        L = min(lengths(j), R.window_s);
        m = min(numel(R.t), round(L / R.dt));
        P(j).L = L;
        P(j).ratio_med = median(rt(1:m), 'omitnan');
        P(j).ratio_p05 = pctl(rt(1:m), 0.05);
        P(j).ratio_p95 = pctl(rt(1:m), 0.95);
        P(j).frac_on_E = median(R.frac_on(E, 1:m));
        P(j).frac_on_I = median(R.frac_on(I, 1:m));
        P(j).gain_ratio_EI = median(gr(1:m), 'omitnan');
        P(j).rate_top_med = median(R.rate_top(1:m));
        P(j).frac_unstable = mean(R.rate_top(1:m) > 0);
    end
    x.pre = P;
    x.cal_t = R.cal.t;
    if isempty(X); X = x; else; X(end + 1) = x; end %#ok<SAGROW>
    clear R
    if mod(k, 10) == 0; fprintf('%d/%d runs read (%.0f s)\n', k, n_runs, toc(t0)); end
end
save(out, 'X', 'settings', 'lengths', '-v7');
fprintf('wrote %s: %d runs (%.0f s)\n', out, numel(X), toc(t0));

function v = pctl(x, p)
% p-quantile ignoring NaN, linear between order statistics at (i - 0.5)/n.
s = sort(x(~isnan(x))); n = numel(s);
if n == 0; v = NaN; return; end
pos = min(max(p * n + 0.5, 1), n); lo = floor(pos); hi = ceil(pos);
v = s(lo) + (pos - lo) * (s(hi) - s(lo));
end
