function T = p6_table()
% P6_TABLE Block B pilot 6, questions 1, 2 and 4: count the hyperexcitable
% networks of the saved mu7chain essential run.
%
% Reads the joint parameter-space sample (32 connectivity points x 4
% conditions) through ParamSpaceAnalysis2, builds one row per network with
% its parameters (effective_param), its realised E:I weight balance (the
% network is rebuilt, not simulated) and lambda_1 / mean rate / n_positive /
% h_KS per condition, classifies, prints the counts with Wilson intervals and
% the cross-tabulations, and writes the table and two figures.
%
% Simulates nothing. Assumes setup_paths has been run. About 20 s.

P = p6_config();
if ~exist(P.out_dir, 'dir'), mkdir(P.out_dir); end

d = dir(fullfile(P.run_dir, 'param_space_*'));
assert(numel(d) == 1, 'expected one param_space directory');
psa = ParamSpaceAnalysis2.from_dir(fullfile(d.folder, d.name));
psa.load_condition_results();

axes_names = {'n', 'f_E', 'level_of_chaos', 'mu_EE_relative', ...
    'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'};
r0 = psa.results.(P.conds{1});
idx = find(~cellfun(@isempty, r0));
nNet = numel(idx);

T = table(idx(:), 'VariableNames', {'config_idx'});
for a = 1:numel(axes_names)
    T.(axes_names{a}) = arrayfun(@(i) psa.effective_param(r0{i}, axes_names{a}), idx(:));
end
T.network_seed = arrayfun(@(i) r0{i}.network_seed, idx(:));

% Realised E:I weight balance, |sum W_E| / (|sum W_E| + |sum W_I|), the
% quantity fig_EI_weights_param_space colours by (0.5 = balanced).
T.B_E = nan(nNet, 1);
for k = 1:nNet
    m = psa.rebuild_model(r0{idx(k)});
    m.verbose = 'near-none';
    m.build();
    ti = m.type_indices;
    SE = full(sum(sum(m.W(:, ti{1}))));
    SI = full(sum(sum(m.W(:, ti{2}))));
    T.B_E(k) = abs(SE) / (abs(SE) + abs(SI));
end

for c = 1:numel(P.conds)
    rc = psa.results.(P.conds{c});
    assert(isequal(find(~cellfun(@isempty, rc)), idx), 'conditions differ in sample');
    for f = {'LLE', 'mean_rate', 'n_positive', 'h_KS_bits', 'D_KY', 'lambda_1_drift', 'frac_local_positive'}
        T.([f{1} '_' P.conds{c}]) = arrayfun(@(i) rc{i}.(f{1}), idx(:));
    end
    assert(all(arrayfun(@(i) rc{i}.success, idx)), 'failed job in %s', P.conds{c});
end

%% Classification
L = zeros(nNet, 4); R = L; FLP = L;
for c = 1:4
    L(:, c) = T.(['LLE_' P.conds{c}]);
    R(:, c) = T.(['mean_rate_' P.conds{c}]);
    FLP(:, c) = T.(['frac_local_positive_' P.conds{c}]);
end
A  = L > P.lle_thr;                 % (a) chaotic
Bh = R > P.rate_hi;                 % (b) high rate, primary threshold
Bl = R > P.rate_lo;                 % (b') high rate, lower threshold
Ch = A & Bh;  Cl = A & Bl;          % (c) both
sat = ~A & R > 0.9;                 % saturated, stable
for c = 1:4
    T.(['chaotic_' P.conds{c}])  = A(:, c);
    T.(['highrate_' P.conds{c}]) = Bh(:, c);
end

fid = fopen(fullfile(P.out_dir, 'p6_counts.txt'), 'w');
assert(fid > 0, 'cannot open the counts file');
out = @(varargin) both(fid, varargin{:});
out('Pilot 6 counts, %d networks, thresholds: lambda_1 > %g, rate > %.2f (second: %.2f)\n', ...
    nNet, P.lle_thr, P.rate_hi, P.rate_lo);
out('Intervals are 95 %% Wilson score intervals.\n\n');
defs = {'(a) lambda_1 > 0', A; ...
    sprintf('(a'') lambda_1 > %.1f', P.lle_margin), L > P.lle_margin; ...
    sprintf('(a") lambda_1 > -%.1f', P.lle_margin), L > -P.lle_margin; ...
    sprintf('(b) rate > %.2f', P.rate_hi), Bh; ...
    sprintf('(b'') rate > %.2f', P.rate_lo), Bl; ...
    sprintf('(c) (a) and rate > %.2f', P.rate_hi), Ch; ...
    sprintf('(c'') (a) and rate > %.2f', P.rate_lo), Cl; ...
    '(a) or (b): chaotic or high rate', A | Bh; ...
    'saturated: stable and rate > 0.9', sat; ...
    'not chaotic, but leading local rate positive part of the time', ~A & FLP > 0.05};
% The last row is a proxy for a stable network that discharges in bursts: the
% stored frac_local_positive (share of the Lyapunov window in which the
% leading direction's local rate is positive). Network 761 under chained MTS,
% the one periodic burster among the six rebuilt, has 0.115; the stable
% non-bursting ones have 0. Validated on that one network only.
counts = struct();
for q = 1:size(defs, 1)
    out('%s\n', defs{q, 1});
    for c = 1:4
        k = sum(defs{q, 2}(:, c));
        [lo, hi] = wilson(k, nNet);
        out('   %-16s %2d / %d = %.2f  [%.2f, %.2f]\n', P.cond_lbl{c}, k, nNet, k / nNet, lo, hi);
    end
end

out('\nCross-tabulation against no adaptation (rows: without; columns: with)\n');
for q = [1 4 6 8]
    for c = 2:4
        x = defs{q, 2};
        tab = [sum(x(:, 1) & x(:, c)), sum(x(:, 1) & ~x(:, c)); ...
               sum(~x(:, 1) & x(:, c)), sum(~x(:, 1) & ~x(:, c))];
        out('%s | %s\n', defs{q, 1}, P.cond_lbl{c});
        out('                 with: yes   no\n');
        out('   without: yes        %3d  %3d\n', tab(1, :));
        out('   without: no         %3d  %3d\n', tab(2, :));
    end
end

out('\nWhere the classes lie (chained MTS), mean of each axis [count]\n');
classes = {'chaotic, rate > 0.40', Ch(:, 3); 'chaotic, rate <= 0.40', A(:, 3) & ~Bh(:, 3); ...
    'stable, rate > 0.40', ~A(:, 3) & Bh(:, 3); 'stable, rate <= 0.40', ~A(:, 3) & ~Bh(:, 3)};
out('%-24s %3s %6s %5s %5s %6s %6s %6s %6s %5s\n', 'class', 'N', 'n', 'f_E', 'gain', 'muEE', 'muEI', 'muIE', 'muII', 'B_E');
for q = 1:size(classes, 1)
    s = classes{q, 2};
    out('%-24s %3d %6.0f %5.2f %5.2f %6.2f %6.2f %6.2f %6.2f %5.2f\n', classes{q, 1}, sum(s), ...
        mean(T.n(s)), mean(T.f_E(s)), mean(T.level_of_chaos(s)), mean(T.mu_EE_relative(s)), ...
        mean(T.mu_EI_relative(s)), mean(T.mu_IE_relative(s)), mean(T.mu_II_relative(s)), mean(T.B_E(s)));
end
out('\nlambda_1 > 0 by level_of_chaos (gain) level, per condition\n');
for g = unique(T.level_of_chaos)'
    s = T.level_of_chaos == g;
    out('   gain %.2f (N = %2d): %s\n', g, sum(s), mat2str(sum(A(s, :), 1)));
end
out('rate > %.2f by f_E level, per condition\n', P.rate_hi);
for g = unique(T.f_E)'
    s = T.f_E == g;
    out('   f_E %.2f (N = %2d): %s\n', g, sum(s), mat2str(sum(Bh(s, :), 1)));
end

%% Question 4: sample size for 50 hyperexcitable networks
out('\nNetworks to sample for 50 hyperexcitable ones (chained MTS), 50 / p at the estimate and at the interval ends\n');
for q = [1 4 6 7 8]
    k = sum(defs{q, 2}(:, 3));
    [lo, hi] = wilson(k, nNet);
    out('   %-34s p = %2d/%d: %5.0f  [%5.0f to %5.0f]\n', defs{q, 1}, k, nNet, ...
        50 / (k / nNet), 50 / hi, 50 / lo);
    counts.(sprintf('def%d', q)) = [k lo hi];
end

%% 1-D sweeps: the same definitions (7 sweeps x 4 levels x 2 reps, others at the preset)
out('\n1-D sensitivity sweeps (everything else at the preset; 16 jobs per sweep and condition)\n');
out('%-18s %-16s %9s %9s %9s\n', 'axis', 'condition', 'lam1>0', sprintf('r>%.2f', P.rate_hi), 'both');
sd = dir(fullfile(P.run_dir, '1D_sensitivity_*'));
tot = zeros(4, 4);
for s = 1:numel(sd)
    ps = ParamSpaceAnalysis2.from_dir(fullfile(sd(s).folder, sd(s).name));
    ps.load_condition_results();
    ax_name = regexprep(sd(s).name, '^1D_sensitivity_sensitivity_(.*)_nLevs.*$', '$1');
    for c = 1:4
        rc = ps.results.(P.conds{c});
        rc = rc(~cellfun(@isempty, rc));
        l = cellfun(@(x) x.LLE, rc); r = cellfun(@(x) x.mean_rate, rc);
        out('%-18s %-16s %4d / %2d %4d / %2d %4d / %2d\n', ax_name, P.cond_lbl{c}, ...
            sum(l > 0), numel(l), sum(r > P.rate_hi), numel(l), sum(l > 0 & r > P.rate_hi), numel(l));
        tot(c, :) = tot(c, :) + [sum(l > 0), sum(r > P.rate_hi), sum(l > 0 & r > P.rate_hi), numel(l)];
    end
end
for c = 1:4
    out('%-18s %-16s %4d / %2d %4d / %2d %4d / %2d\n', 'ALL SWEEPS', P.cond_lbl{c}, ...
        tot(c, 1), tot(c, 4), tot(c, 2), tot(c, 4), tot(c, 3), tot(c, 4));
end
fclose(fid);

writetable(T, fullfile(P.out_dir, 'p6_network_table.csv'));
save(fullfile(P.out_dir, 'p6_network_table.mat'), 'T', 'P', 'counts');

%% Figure 1: lambda_1 against mean rate, one panel per condition
f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 430]);
tl = tiledlayout(f1, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    yl = [-11 8];
    plot(ax, [P.rate_hi P.rate_hi], yl, '--', 'Color', [0.5 0.5 0.5]);
    plot(ax, [P.rate_lo P.rate_lo], yl, ':', 'Color', [0.5 0.5 0.5]);
    plot(ax, [0.005 1.2], [0 0], '-', 'Color', [0 0.6 0]);
    s = Ch(:, c);
    plot(ax, R(~s, c), L(~s, c), 'o', 'Color', P.cond_col(c, :), 'MarkerSize', 7, 'LineWidth', 1.2);
    plot(ax, R(s, c), L(s, c), 'o', 'Color', P.cond_col(c, :), ...
        'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerSize', 9, 'LineWidth', 1.2);
    if c == 3
        for k = find(ismember(T.config_idx, P.rebuild_idx))'
            text(ax, R(k, c) * 1.12, L(k, c) + 0.45, sprintf('%d', T.config_idx(k)), 'FontSize', 8);
        end
    end
    set(ax, 'XScale', 'log', 'XLim', [0.005 1.2], 'YLim', yl, 'FontSize', 11, 'Box', 'off', ...
        'XTick', [0.01 0.1 1], 'XTickLabel', {'0.01', '0.1', '1'});
    title(ax, sprintf('%s\n\\lambda_1>0: %d   rate>%.2f: %d   both: %d', P.cond_lbl{c}, ...
        sum(A(:, c)), P.rate_hi, sum(Bh(:, c)), sum(Ch(:, c))), 'FontSize', 11, 'FontWeight', 'normal');
    xlabel(ax, 'mean rate (0 to 1)');
    if c == 1, ylabel(ax, '\lambda_1 (1/s)'); end
end
title(tl, sprintf(['The 32 sampled networks under each condition. Red fill: chaotic and high rate. ' ...
    'Dashed / dotted: rate thresholds %.2f / %.2f. Numbers: networks rebuilt in question 3.'], ...
    P.rate_hi, P.rate_lo), 'FontSize', 10);
save_fig(f1, P.out_dir, 'p6_lle_vs_rate');

%% Figure 2: where the classes lie in connectivity space (chained MTS and no adaptation)
f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 620]);
tl = tiledlayout(f2, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
rng_state = rng; rng(1); jit = @(v, w) v + w * (rand(size(v)) - 0.5);
xs = {T.B_E, T.f_E, T.level_of_chaos, T.mu_EE_relative + T.mu_EI_relative};
xl = {'E:I weight balance B_E (0.5 = balanced)', 'E fraction f_E', ...
    'gain on all weights (level\_of\_chaos)', '\mu_{EE} + \mu_{EI} (net mean drive onto E)'};
jw = [0 0.08 0.2 1.2];
rows = [1 3];
for rr = 1:2
    c = rows(rr);
    for p = 1:4
        ax = nexttile(tl); hold(ax, 'on');
        x = jit(xs{p}, jw(p));
        cls = {~A(:, c) & ~Bh(:, c), ~A(:, c) & Bh(:, c), A(:, c) & ~Bh(:, c), A(:, c) & Bh(:, c)};
        col = {[0.6 0.6 0.6], [0.9 0.6 0], [0 0.45 0.7], [0.85 0.1 0.1]};
        for q = 1:4
            plot(ax, x(cls{q}), R(cls{q}, c), 'o', 'Color', col{q}, 'MarkerFaceColor', col{q}, 'MarkerSize', 6);
        end
        plot(ax, xlim(ax), [P.rate_hi P.rate_hi], '--', 'Color', [0.5 0.5 0.5]);
        set(ax, 'YScale', 'log', 'YLim', [0.005 1.2], 'FontSize', 10, 'Box', 'off', ...
            'YTick', [0.01 0.1 1], 'YTickLabel', {'0.01', '0.1', '1'});
        if p == 1, ylabel(ax, sprintf('%s\nmean rate', P.cond_lbl{c})); end
        if rr == 2, xlabel(ax, xl{p}); end
        if rr == 1 && p == 4
            legend(ax, {'stable, low rate', 'stable, high rate', 'chaotic, low rate', 'chaotic, high rate'}, ...
                'Location', 'southeast', 'FontSize', 8, 'Box', 'off');
        end
    end
end
rng(rng_state);
title(tl, 'Where the classes lie in connectivity space (points jittered horizontally on the discrete axes)', 'FontSize', 10);
save_fig(f2, P.out_dir, 'p6_where_in_parameter_space');
end

function both(fid, varargin)
fprintf(1, varargin{:});
fprintf(fid, varargin{:});
end

function [lo, hi] = wilson(k, n)
z = 1.959964; p = k / n;
den = 1 + z^2 / n;
ctr = (p + z^2 / (2 * n)) / den;
hw = z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den;
lo = max(0, ctr - hw); hi = min(1, ctr + hw);
end

function save_fig(f, out_dir, tag)
set(f, 'Visible', 'on');   % so the saved .fig opens visible
savefig(f, fullfile(out_dir, [tag '.fig']));
set(f, 'Visible', 'off');
exportgraphics(f, fullfile(out_dir, [tag '.png']), 'Resolution', 110);
close(f);
end
