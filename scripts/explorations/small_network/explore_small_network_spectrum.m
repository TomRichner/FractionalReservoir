% EXPLORE_SMALL_NETWORK_SPECTRUM The distribution of the top-K Lyapunov exponents of one condition of the 80-neuron network.
%
% WHY. At 80 neurons the top 20 exponents of the adapting conditions are
% still one tight band (explore_small_network). This script takes many
% more of them, K = 200 as set below, for one condition (the chained
% multiple-timescale condition sfa3_stdc2 by default), so the bands and the
% gaps between them can be seen, and prints the per-neuron time constants
% of the slow variables the bands may belong to.
%
% THE NETWORK AND THE RUN are those of explore_small_network: the network
% from small_network_model (n = 80, alpha = 0.2, F = 0.18634), seed [7 8],
% 15 s at 500 Hz with sra1, topk exponents accumulated over 5 to 15 s
% after a 5 s warm-up, lya_dt 0.05 s. Only K differs.
%
% THE FIGURE (one, left open, saved as PNG to
% figs/explorations/small_network/, gitignored):
%   (a) a histogram of the K exponents, in 1/s, bins bin_width wide
%       (bin_width = [] picks it from the data: about 8 exponents per bin
%       on average over their range, rounded to a multiple of 0.005 1/s)
%   (b) the same exponents sorted, exponent against rank 1..K, as markers:
%       a band shows as a plateau, a gap as a step
% Both panels carry a zero line. No predicted rates are drawn.
%
% WHAT IT PRINTS
%   - all K exponents, lambda_1, the number of positive exponents, and the
%     wall time split into build, the Lyapunov computation (lyapunov_topk's
%     own timer, lya_results.seconds) and the rest of run() (mainly the
%     trajectory)
%   - the largest steps between consecutive sorted exponents and the ranks
%     at which they fall (rank i means between exponent i and i+1)
%   - from the built model, as rates 1/tau (min / median / max over the
%     neurons): the slowest SFA rung of every neuron (tau_a_matrix, last
%     column), and, on the depression routes, the releasable-pool refill
%     tau_rec1, the release time tau_rel and the reserve refill tau_rec2,
%     plus the reserve size n2, per presynaptic neuron (the compiled route
%     parameters of get_params; a neuron's routes share them)
%   - a CONVERGENCE CHECK (print only, no figure): the same condition and K
%     run for 30 s with the window 5 to 30 s, and how far the exponents
%     moved from the 15 s run (max and median absolute change over the K,
%     and the change in lambda_1). Off by default (run_check); it was run
%     at K = 100 (2026-11-19): every exponent rose by about 0.015 to 0.022
%     1/s, lambda_1 from -0.1185 to -0.1035.
%
% POOL: one condition, so no pool. pool_set(0) records that this session
% holds none, which keeps library code serial and quiet; nothing is opened.
%
% PATH: run setup_paths once in the session first, with the current folder
% at the repository root (setup_paths adds scripts/ with genpath, which
% includes this folder).
%
% See also: explore_small_network, small_network_model, lyapunov_topk

%% Settings
condition      = 'sfa3_stdc2';  % any condition small_network_model accepts
lya_K          = 200;           % how many exponents
rng_seeds      = [7 8];         % network seed pair
T_range        = [0 15];        % s
fs             = 500;           % Hz
lya_dt         = 0.05;          % s, re-orthonormalization interval
lya_warmup     = 5;             % s
lya_T_interval = [5 15];        % s, accumulation window
bin_width      = [];            % 1/s, histogram bin width; [] = from the data (see below)
n_gaps         = 8;             % how many of the largest steps to list
run_check      = false;         % also run the 30 s convergence check (done at K = 100)
T_range_check  = [0 30];        % s
lya_T_interval_check = [5 30];  % s
save_fig       = true;
out_dir        = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', 'small_network');

%% No pool
pool_set(0, 'label', 'small_network_spectrum');

%% Build and run
lya_args = {'lya_method', 'topk', 'lya_K', lya_K, 'lya_dt', lya_dt, ...
    'lya_warmup', lya_warmup};
t_start = tic;
m = small_network_model(condition, 'seed', rng_seeds, 'T_range', T_range, ...
    'fs', fs, 'lya_T_interval', lya_T_interval, lya_args{:});
build_s = toc(t_start);
t_run = tic;
m.run();
run_s = toc(t_run);
wall_s = toc(t_start);
lya_s = m.lya_results.seconds;    % timed inside lyapunov_topk
lambda = sort(m.lya_results.LE_spectrum(:), 'descend');

%% Diagnostics: the exponents
fprintf('\n=== n = %d, %s, seed [%d %d], K = %d, %g-%g s ===\n', ...
    m.n, condition, rng_seeds, lya_K, lya_T_interval);
fprintf('  state dimension %d\n', m.N_sys_eqs);
fprintf(['  wall time %.2f s: build %.2f s, run %.2f s ' ...
    '(Lyapunov %.2f s, trajectory and the rest of run() %.2f s)\n'], ...
    wall_s, build_s, run_s, lya_s, run_s - lya_s);
fprintf('  lambda_1 = %+.6f, positive exponents %d of %d\n', ...
    lambda(1), nnz(lambda > 0), numel(lambda));
fprintf('  all exponents (1/s), by rank:\n');
for i0 = 1:10:numel(lambda)
    i1 = min(i0 + 9, numel(lambda));
    fprintf('   %3d-%3d:', i0, i1);
    fprintf(' %+.4f', lambda(i0:i1));
    fprintf('\n');
end

% The gaps: steps between consecutive sorted exponents, largest first.
steps = lambda(1:end-1) - lambda(2:end);
[step_sorted, step_rank] = sort(steps, 'descend');
fprintf('  largest steps between consecutive exponents (rank i: between i and i+1):\n');
for g = 1:min(n_gaps, numel(steps))
    i = step_rank(g);
    fprintf('    rank %3d: %+.4f -> %+.4f, step %.4f\n', ...
        i, lambda(i), lambda(i + 1), step_sorted(g));
end
fprintf('  median step %.4f\n', median(steps));

%% Diagnostics: the slow time constants of the built model, as rates 1/tau
p = m.get_params();
types = m.cell_type_names;
fprintf('\n  time constants as rates 1/tau (1/s): min / median / max over the neurons\n');
tau_slow = [];
for q = 1:numel(types)
    tam = p.tau_a_matrix{q};
    if ~isempty(tam)
        tau_slow = [tau_slow; max(tam, [], 2)]; %#ok<AGROW>
    end
end
if isempty(tau_slow)
    fprintf('    no SFA in this condition\n');
else
    print_rates('slowest SFA rung', tau_slow);
end
% Depression: one value per presynaptic neuron. A neuron's routes share
% their draws, so take each presynaptic type's first route that carries
% depression.
dep = struct('tau_rec1', [], 'tau_rel', [], 'tau_rec2', [], 'n2', []);
dep_names = fieldnames(dep);
for pre = 1:numel(types)
    for post = 1:numel(types)
        P = p.route{pre, post}.P;
        if isfield(P, 'tau_rec1')
            for f = 1:numel(dep_names)
                if isfield(P, dep_names{f})
                    dep.(dep_names{f}) = [dep.(dep_names{f}); P.(dep_names{f})(:)];
                end
            end
            break;
        end
    end
end
for f = 1:3
    if ~isempty(dep.(dep_names{f}))
        print_rates(['depression ' dep_names{f}], dep.(dep_names{f}));
    end
end
if ~isempty(dep.n2)
    fprintf('    %-26s %.3f / %.3f / %.3f  (a size, not a rate)\n', ...
        'chain reserve size n2', min(dep.n2), median(dep.n2), max(dep.n2));
end
if all(cellfun(@isempty, struct2cell(dep)))
    fprintf('    no per-neuron depression parameters exposed for this condition\n');
end

%% Figure
fig = figure('Name', sprintf('Small network (n = 80): %s, top-%d Lyapunov exponents', ...
    condition, lya_K), 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 900 650]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('n = %d, %s, seed [%d %d], K = %d, \\lambda_1 = %+.4f 1/s  (window %g-%g s)', ...
    m.n, strrep(condition, '_', '\_'), rng_seeds, lya_K, lambda(1), lya_T_interval));

ax1 = nexttile(tl);
if isempty(bin_width)
    % From the data: about per_bin exponents per bin on average over the
    % range, rounded to a multiple of 0.005 1/s (at least 0.005). The
    % exponents are spread nearly evenly over their range, so this keeps
    % the bins narrow enough to show structure inside the top band as K
    % grows. (A Freedman-Diaconis width was tried first: at K = 200 it
    % gave 0.045 1/s, which put the whole top band in one or two bins.)
    per_bin = 8;
    range_lambda = lambda(1) - lambda(end);
    bin_width = max(0.005, round(range_lambda * per_bin / numel(lambda) / 0.005) * 0.005);
    fprintf('  histogram bins: range %.4f 1/s, %d exponents, %d per bin on average -> %.3f 1/s\n', ...
        range_lambda, numel(lambda), per_bin, bin_width);
end
edges = (floor(min(lambda) / bin_width) * bin_width):bin_width:(ceil(max(lambda) / bin_width) * bin_width + bin_width);
histogram(ax1, lambda, edges, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'w');
xline(ax1, 0, 'k--');
xlim(ax1, [edges(1), max(edges(end), 0) + bin_width]);
xlabel(ax1, 'Lyapunov exponent (1/s)');
ylabel(ax1, 'Count');
title(ax1, sprintf('(a) Histogram, bins %.3g 1/s', bin_width));
box(ax1, 'off');

ax2 = nexttile(tl);
plot(ax2, 1:numel(lambda), lambda, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.2 0.4 0.7], 'MarkerEdgeColor', 'none');
yline(ax2, 0, 'k--');
xlim(ax2, [0, numel(lambda) + 1]);
xlabel(ax2, 'Rank');
ylabel(ax2, 'Lyapunov exponent (1/s)');
title(ax2, '(b) Sorted exponents');
box(ax2, 'off');
grid(ax2, 'on');

if save_fig
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    png = fullfile(out_dir, sprintf('small_network_spectrum_%s_K%d.png', condition, lya_K));
    exportgraphics(fig, png, 'Resolution', 150);
    fprintf('\n  saved %s\n', png);
end

%% Convergence check: the same run for 30 s (print only)
if run_check
    t_start = tic;
    m_check = small_network_model(condition, 'seed', rng_seeds, ...
        'T_range', T_range_check, 'fs', fs, ...
        'lya_T_interval', lya_T_interval_check, lya_args{:});
    m_check.run();
    wall_check = toc(t_start);
    lambda_check = sort(m_check.lya_results.LE_spectrum(:), 'descend');
    change = abs(lambda_check - lambda);
    fprintf('\n=== convergence check: %g-%g s window against %g-%g s ===\n', ...
        lya_T_interval_check, lya_T_interval);
    fprintf('  wall time %.1f s\n', wall_check);
    fprintf('  lambda_1: %+.6f (15 s) -> %+.6f (30 s), change %+.6f\n', ...
        lambda(1), lambda_check(1), lambda_check(1) - lambda(1));
    fprintf('  |change| over the %d exponents (rank by rank): max %.4f (rank %d), median %.4f\n', ...
        numel(change), max(change), find(change == max(change), 1), median(change));
    fprintf('  positive exponents in the 30 s run: %d\n', nnz(lambda_check > 0));
    fprintf('  30 s exponents (1/s), by rank:\n');
    for i0 = 1:10:numel(lambda_check)
        i1 = min(i0 + 9, numel(lambda_check));
        fprintf('   %3d-%3d:', i0, i1);
        fprintf(' %+.4f', lambda_check(i0:i1));
        fprintf('\n');
    end
end

%% Local function
function print_rates(label, tau)
% PRINT_RATES One line: min / median / max of 1/tau, and the tau range.
rate = 1 ./ tau(:);
fprintf('    %-26s %.4f / %.4f / %.4f   (tau %.3g to %.3g s, %d neurons)\n', ...
    label, min(rate), median(rate), max(rate), min(tau), max(tau), numel(tau));
end
