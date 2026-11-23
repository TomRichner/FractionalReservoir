% OW_COST_FD057 Estimated wall time of outlier_spectra on 12 workers after FD-057: the reference at 200 s, E-to-E x 1.5 at 20 s.
%
% From the per-run timings of the fastmedium runs (ow_series.mat, ow_extract).
% Reference runs: as measured. E-to-E x 1.5 runs: their measured parts
% rescaled to a 20 s window at the short window's steps (frozen 0.01 s:
% 2000 states; full 0.1 s: 200 states): simulation by 80/260, the frozen
% pass by its per-state time times 2000, the full pass by its per-state time
% times 200, null and dense unchanged. The other networks: as measured.
% Scheduling: run_outlier_spectra_stage's order (cost key window/step +
% 50 window for full runs, descending, stable), one run at a time on each of
% 12 workers. EXTRAPOLATED for the E-to-E x 1.5 runs; the scheduling model
% gives 393 min for the measured times against 387.8 min logged.

setup_paths();
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
S0 = load(fullfile(root, 'data', 'explorations', 'outlier_window', 'ow_series.mat'));
X = S0.X;
n = numel(X);
t_run = [X.seconds_index];
key = zeros(1, n);
for r = 1:n
    x = X(r); s = x.seconds;
    switch x.network
        case 'reference'
            key(r) = 200 / 0.05 + 50 * 200;
        case 'EEx1p5'
            other = x.seconds_index - (s.sim + s.fast + s.null + s.full + s.dense);
            t_run(r) = s.sim * 80 / 260 + s.fast / numel(x.t) * 2000 + s.null + ...
                s.full / numel(x.full_t) * 200 + s.dense + other;
            key(r) = 20 / 0.01 + 50 * 20;
        otherwise
            key(r) = 20 / 0.01;
    end
end
[~, order] = sort(key, 'descend');
free = zeros(1, 12); finish = zeros(1, n);
for r = order
    [t0, w] = min(free);
    free(w) = t0 + t_run(r);
    finish(r) = free(w);
end
[t_long, i_long] = max(t_run);
nets = {X.network};
fprintf('worker-min: reference %.0f, EEx1p5 %.0f, others %.0f, all %.0f\n', sum(t_run(strcmp(nets, 'reference'))) / 60, ...
    sum(t_run(strcmp(nets, 'EEx1p5'))) / 60, sum(t_run(~ismember(nets, {'reference', 'EEx1p5'}))) / 60, sum(t_run) / 60);
fprintf('estimated wall time on 12 workers: %.0f min (lower bound sum/12 = %.0f min)\n', max(free) / 60, sum(t_run) / 12 / 60);
fprintf('longest run: %s %s [%d %d], %.0f min, ends at %.0f min\n', X(i_long).network, X(i_long).condition, X(i_long).seed, ...
    t_long / 60, finish(i_long) / 60);
for c = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
    sel = strcmp(nets, 'EEx1p5') & strcmp({X.condition}, c{1});
    fprintf('EEx1p5 %s at 20 s: median %.1f min per run (was %.1f at 200 s)\n', c{1}, median(t_run(sel)) / 60, ...
        median([X(sel).seconds_index]) / 60);
end
