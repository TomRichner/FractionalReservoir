function out = fig_memory_capacity(cfg)
% FIG_MEMORY_CAPACITY Paper-ready memory-capacity strip, from a saved MC run.
%
%   out = FIG_MEMORY_CAPACITY()
%   out = FIG_MEMORY_CAPACITY('mat_file', f)
%   out = FIG_MEMORY_CAPACITY('run_dir', d)   % looks in <run_dir>/memory_capacity
%
% A 1 x 3 strip assembled from a finished run_memory_capacity ensemble:
%   (a) cumulative memory capacity against delay, bootstrap band, y axis
%       spanning the data (the old fixed [0 10.9] left the paper's curves in
%       the bottom tenth of the panel)
%   (b) per-delay R^2 with the saved bootstrap band (summary.R2_ci95)
%   (c) memory horizon: every paired trial as a connected grey line across the
%       conditions, condition-coloured markers, the median as a thick bar, and
%       the paired sign-flip p and Cohen's d_z for the adjacent pairs printed
%       in the panel -- READ from summary.stats, never recomputed here.
%
% Beside the figure, Fig_Memory_Capacity_table.md: MC and horizon mean
% [bootstrap 95% CI] per condition, every paired test, n trials, readout
% signal, and the source .mat, so a number quoted in the manuscript can be
% traced to the saved summary rather than read off a plot.
%
% No simulation is re-run. Everything is read from `results_all` (see
% run_memory_capacity.m): conditions, H_trials, summary.R2_mean, summary.R2_ci95,
% summary.MC_mean/MC_ci95, summary.H_mean/H_ci95, summary.stats(p).pair /
% p_perm / cohens_dz / exact / n_patterns, settings.delay_s / n_trials /
% readout_signal / R2_threshold_for_horizon.
%
% RESOLUTION: an explicit mat_file, else the newest *_results.mat under
% <run_dir>/memory_capacity -- and NOTHING outside the run directory
% (resolve_data_file). It ERRORS rather than silently plotting nothing.
%
% Until 2026-09-14 this delegated to plot_memory_capacity_combined and opened
% the two working figures of replot_memory_capacity as well; the strip is now
% drawn here so the panel changes above stay local to the paper figure, and the
% working figures are no longer opened (the master closes what a figure
% returns, and those were never returned).
%
% See also: run_memory_capacity, plot_memory_capacity, plot_memory_capacity_combined

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.mat_file    (1,:) char    = ''
    cfg.run_dir     (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));

% INSIDE THE RUN DIRECTORY ONLY. This used to fall back to
% data/memory_capacity/ and then to a paper_ready/ subfolder that has never
% existed. On 2026-09-03 that fallback made this figure plot a .mat from
% 2026-08-22 -- a different network -- and report success, because the run's
% memory_capacity stage had failed and left nothing. To plot a standalone
% analysis, pass its location AS run_dir (fullfile(project_root, 'data')), or
% pass mat_file directly.
mat_file = resolve_data_file(cfg.mat_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'memory_capacity')}, ...
    '*_results.mat', ...
    'Run run_memory_capacity first');
vprintf(cfg.verbose, 'verbose', '[fig_memory_capacity] source: %s\n', mat_file);

S = load(mat_file, 'results_all');
R = S.results_all;

%% Unpack (no recompute)
cond_keys   = R.conditions;
cond_titles = mc_display_names(cond_keys);
n_cond   = numel(cond_keys);
H_trials = R.H_trials;
MC_trials = R.MC_trials;
n_trials = size(H_trials, 1);
sm = R.summary;
st_ = R.settings;
delay_s = st_.delay_s;
R2_mean = sm.R2_mean;  R2_ci = sm.R2_ci95;
R2_cum_mean  = cumsum(R2_mean, 2);
R2_cum_ci_lo = cumsum(R2_ci.lo, 2);
R2_cum_ci_hi = cumsum(R2_ci.hi, 2);
stats = sm.stats;

st = manuscript_style();
style_cleanup = with_graphics_defaults( ...
    'DefaultAxesFontSize',      st.tick_fs, ...
    'DefaultAxesLineWidth',     st.axis_lw, ...
    'DefaultTextInterpreter',   'none', ...
    'DefaultLegendInterpreter', 'none'); %#ok<NASGU>
colors = mc_condition_colors(cond_keys);
xpos = 1:n_cond;
label_fs = 15.4;

%% The strip
fig = figure('Color', 'w', 'Position', [100 82 1089 398]);
tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'loose');
tl.OuterPosition = [0 0.06 1 0.86];

% (a) cumulative MC vs delay, y spanning the data
ax_a = nexttile(tl); hold(ax_a, 'on'); box(ax_a, 'off');
for i = 1:n_cond
    shaded_ci(ax_a, delay_s, R2_cum_ci_lo(i, :), R2_cum_ci_hi(i, :), colors(i, :), 0.12);
    plot(ax_a, delay_s, R2_cum_mean(i, :), '-', 'Color', colors(i, :), 'LineWidth', 2);
end
y_top = 1.1 * max([R2_cum_ci_hi(:); R2_cum_mean(:); 1e-3]);
xlim(ax_a, [0, delay_s(end)]); ylim(ax_a, [0, y_top]);
xlabel(ax_a, 'Delay (s)', 'FontSize', label_fs);
ylabel(ax_a, {'Cumulative', 'Memory Capacity'}, 'FontSize', label_fs);

% (b) per-delay R^2 with the saved bootstrap band; carries the legend
ax_b = nexttile(tl); hold(ax_b, 'on'); box(ax_b, 'off');
h = gobjects(1, n_cond);
for i = 1:n_cond
    shaded_ci(ax_b, delay_s, R2_ci.lo(i, :), R2_ci.hi(i, :), colors(i, :), 0.12);
    h(i) = plot(ax_b, delay_s, R2_mean(i, :), '-', 'Color', colors(i, :), 'LineWidth', 2);
end
yline(ax_b, st_.R2_threshold_for_horizon, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xlim(ax_b, [0, delay_s(end)]); ylim(ax_b, [0, 1]);
xlabel(ax_b, 'Delay (s)', 'FontSize', label_fs);
ylabel(ax_b, '$R^2$', 'Interpreter', 'latex', 'FontSize', label_fs);
set(ax_b, 'YTick', [0 0.5 1]);
legend(ax_b, h, cond_titles, 'Location', 'northeast', 'Box', 'off');

% (c) horizon, paired across trials, median bar, paired tests from the summary
ax_c = nexttile(tl); hold(ax_c, 'on'); box(ax_c, 'off');
jit = 0.08 * (rand(n_trials, 1) - 0.5);   % cosmetic only; same jitter on every condition keeps lines readable
for k = 1:n_trials
    plot(ax_c, xpos + jit(k), H_trials(k, :), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
end
for i = 1:n_cond
    scatter(ax_c, i + jit, H_trials(:, i), 22, 'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
    plot(ax_c, i + [-0.25 0.25], median(H_trials(:, i), 'omitnan') * [1 1], '-', ...
        'Color', colors(i, :), 'LineWidth', 3);
end
xlim(ax_c, [0.5, n_cond + 0.5]);
set(ax_c, 'XTick', xpos, 'XTickLabel', cond_titles);
ylabel(ax_c, 'Memory Horizon (s)', 'FontSize', label_fs);
y_h = max([H_trials(:); 0.1]);
ylim(ax_c, [0, 1.35 * y_h]);
% Adjacent pairs, matched by condition name in stats(p).pair (never by position)
adj = [1:n_cond - 1; 2:n_cond]';
if n_cond >= 3; adj = [adj; 1 n_cond]; end
lines_txt = {};
for a = 1:size(adj, 1)
    p = find_pair(stats, cond_keys{adj(a, 1)}, cond_keys{adj(a, 2)});
    if isempty(p); continue; end
    lines_txt{end + 1} = sprintf('%s vs %s: p = %s, d_z = %.2f', ...
        short_of(st, cond_keys{adj(a, 1)}), short_of(st, cond_keys{adj(a, 2)}), ...
        p_txt(stats(p)), stats(p).cohens_dz); %#ok<AGROW>
end
text(ax_c, 0.03, 0.98, strjoin(lines_txt, newline), 'Units', 'normalized', ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'Interpreter', 'none');
text(ax_c, 0.97, 0.02, sprintf('n = %d paired trials; tests on total MC', n_trials), ...
    'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'Color', [0.3 0.3 0.3]);

AddLetters2Plots(fig, {'(a)', '(b)', '(c)'}, 'FontSize', 18, 'FontWeight', 'normal', ...
    'HShift', -0.04, 'VShift', -0.09);
if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% Table
hdr1 = '| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] | median MC | median horizon |';
rows1 = cell(1, n_cond);
for i = 1:n_cond
    rows1{i} = sprintf('| %s | %.3f [%.3f, %.3f] | %.3f [%.3f, %.3f] | %.3f | %.3f |', cond_titles{i}, ...
        sm.MC_mean(i), sm.MC_ci95.lo(i), sm.MC_ci95.hi(i), sm.H_mean(i), sm.H_ci95.lo(i), sm.H_ci95.hi(i), ...
        median(MC_trials(:, i), 'omitnan'), median(H_trials(:, i), 'omitnan'));
end
hdr2 = '| Pair | mean diff (total MC) | p (sign-flip) | patterns | Cohen''s d_z |';
rows2 = arrayfun(@(p) sprintf('| %s | %+.3f | %s | %d%s | %.2f |', p.pair, p.mean_diff, ...
    p_txt(p), p.n_patterns, tern(p.exact, ' (exact)', ' (Monte Carlo)'), p.cohens_dz), stats, 'UniformOutput', false);
vprintf(cfg.verbose, 'verbose', '%s\n|---|---|---|---|---|\n%s\n\n%s\n|---|---|---|---|---|\n%s\n', ...
    hdr1, strjoin(rows1, newline), hdr2, strjoin(rows2, newline));

%% Save
fig_tag = 'Fig_Memory_Capacity';
out = struct('figs', fig, 'files', {{}}, 'source', mat_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '# Memory capacity\n\nSource: `%s`\n\nPreset `%s`, run mode `%s`, %d paired trials, readout `%s`, horizon threshold R^2 > %.2f, %d bootstrap samples.\n\n', ...
            mat_file, st_.preset_name, st_.run_mode, n_trials, st_.readout_signal, ...
            st_.R2_threshold_for_horizon, st_.n_boot);
        fprintf(fid, '%s\n|---|---|---|---|---|\n%s\n\n', hdr1, strjoin(rows1, newline));
        fprintf(fid, '%s\n|---|---|---|---|---|\n%s\n', hdr2, strjoin(rows2, newline));
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function shaded_ci(ax, x, lo, hi, rgb, alpha_fill)
x = x(:)'; lo = lo(:)'; hi = hi(:)';
fill(ax, [x, fliplr(x)], [hi, fliplr(lo)], rgb, 'FaceAlpha', alpha_fill, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function p = find_pair(stats, a, b)
% Index of the saved pair "a vs b" or "b vs a"; empty if absent.
p = find(strcmp({stats.pair}, sprintf('%s vs %s', a, b)) | ...
         strcmp({stats.pair}, sprintf('%s vs %s', b, a)), 1);
end

function s = p_txt(stat)
if stat.p_perm < 1e-3
    s = sprintf('%.1e', stat.p_perm);
else
    s = sprintf('%.3f', stat.p_perm);
end
end

function v = tern(c, a, b)
if c; v = a; else; v = b; end
end

function s = short_of(st, key)
% Short label when manuscript_style has one, else the full title, else the key.
if isKey(st.condition_short, key)
    s = st.condition_short(key);
elseif isKey(st.condition_title, key)
    s = st.condition_title(key);
else
    s = key;
end
end
