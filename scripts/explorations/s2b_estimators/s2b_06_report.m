function s2b_06_report(out)
%S2B_06_REPORT Print the comparison between seeds for each case of s2b_06.
if nargin < 1; L = load(fullfile(s2b_out_dirs(), 's2b_06_paper_runs.mat')); out = L.out; end
pr = nchoosek(1:numel(out(1).C), 2);
sm = @(x) movmean(x, 10);                       % 0.5 s, about the 2-Hz low-pass Figure 2 G draws
for k = 1:numel(out)
    C = out(k).C;  X = [C.loc];  ns = size(X, 2);
    fprintf('\n=== %s, %s: window [%g, %g] s, %g s of alignment before it ===\n', out(k).cond, out(k).case, out(k).interval, out(k).warmup);
    fprintf('accumulated lambda_1 per seed: %s  (range %.4f)\n', mat2str([C.LLE], 4), range([C.LLE]));
    fprintf('five leading exponents, seed 1: %s ; gap 1-2 per seed: %s\n', mat2str(C(1).spec(1:5), 3), mat2str(arrayfun(@(c) c.spec(1) - c.spec(2), C), 2));
    fprintf('raw column behind the leading rate, per seed: %s\n', mat2str([C.raw_col]));
    cs = arrayfun(@(i) abs(C(pr(i, 1)).q_start' * C(pr(i, 2)).q_start), 1:size(pr, 1));
    ce = arrayfun(@(i) abs(C(pr(i, 1)).q_end' * C(pr(i, 2)).q_end), 1:size(pr, 1));
    fprintf('pairwise |cos| of the leading vectors: start of window %.2f..%.2f (median %.2f), end %.2f..%.2f (median %.2f)\n', ...
        min(cs), max(cs), median(cs), min(ce), max(ce), median(ce));
    d = arrayfun(@(i) rms(X(:, pr(i, 1)) - X(:, pr(i, 2))), 1:size(pr, 1));
    ds = arrayfun(@(i) rms(sm(X(:, pr(i, 1))) - sm(X(:, pr(i, 2)))), 1:size(pr, 1));
    fprintf('local rate: std of a trace %.2f..%.2f s^-1; pairwise rms difference %.2f..%.2f (median %.2f) = %.0f%% of the std; 0.5-s smoothed: %.2f (median) against std %.2f = %.0f%%\n', ...
        min(std(X)), max(std(X)), min(d), max(d), median(d), 100 * median(d) / median(std(X)), median(ds), median(std(sm(X))), 100 * median(ds) / median(std(sm(X))));
    h = floor(size(X, 1) / 2);
    fprintf('share of samples > 0 per seed: %s (range %.3f); seed 1 first half %.3f, second half %.3f\n', ...
        mat2str(round(mean(X > 0), 3)), range(mean(X > 0)), mean(X(1:h, 1) > 0), mean(X(h + 1:end, 1) > 0));
    qs = quantile(X, [0.25 0.5 0.75]);
    q1 = quantile(X(1:h, 1), [0.25 0.5 0.75]);  q2 = quantile(X(h + 1:end, 1), [0.25 0.5 0.75]);
    fprintf('quartiles (25, 50, 75%%) per seed:\n');  disp(round(qs', 3));
    fprintf('   largest difference between seeds %s; between halves of seed 1 %s\n', mat2str(round(range(qs, 2)', 3)), mat2str(round(abs(q1 - q2), 3)));
    fprintf('   (%d seeds, %d samples each)\n', ns, size(X, 1));
end
end
