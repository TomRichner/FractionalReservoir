function out = fig_memory_capacity_example(cfg)
% FIG_MEMORY_CAPACITY_EXAMPLE Example memory capacity, one network, 4 conditions.
%
%   out = FIG_MEMORY_CAPACITY_EXAMPLE()
%   out = FIG_MEMORY_CAPACITY_EXAMPLE('data_file', f)
%
% (a) cumulative memory capacity against delay, all four conditions.
% (b) per-delay R^2, all four conditions.
% Below: input reconstruction (target against trained readout) at a few delays
% for the SFA+STD condition, each panel titled with its delay and R^2.
%
% TWO-STEP, no re-simulation at plot time: run_memory_capacity_example runs the
% protocol and saves mc_example_data.mat (gitignored); this renders it, so the
% look can be iterated without re-running the reservoir.
%
% See also: run_memory_capacity_example, fig_memory_capacity, manuscript_style

arguments
    cfg.data_file   (1,:) char    = ''    % '' -> search run_dir, then data/mc_example
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''    % the run whose mc_example data to plot
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));

%% Load precomputed results
% Search the RUN first, then the standalone location. This used to load a
% hardcoded .mat beside this file, with run_dir marked "unused" -- so on
% 2026-08-26 it plotted Aug 22 data while the sweep it was handed was Aug 25.
data_file = resolve_data_file(cfg.data_file, ...
    {fullfile(cfg.run_dir, 'mc_example'), ...
     fullfile(project_root, 'data', 'mc_example')}, ...
    'mc_example_data.mat', ...
    'Run run_memory_capacity_example first');
S = load(data_file);
results         = S.results;
R2              = S.R2; %#ok<NASGU>
delay_s         = S.delay_s;
condition_names = S.condition_names;
n_cond          = numel(condition_names);

%% Style + palette (matched to Fig_memory_capacity.m / plot_memory_capacity_combined.m)
% Root defaults are RESTORED on return. These were bare set(0,...) calls that
% leaked 'none' interpreters into the rest of the session and broke every tex
% label drawn afterwards. style_cleanup must stay in scope for the whole build.
st = manuscript_style();
style_cleanup = with_manuscript_defaults(); %#ok<NASGU>

% Okabe-Ito colorblind-safe palette (same as the combined MC figure):
%   Baseline black, SFA orange, STD sky blue, SFA+STD reddish purple.
% Condition colours come from manuscript_style, keyed BY NAME, so this figure
% and the ensemble MC figure cannot drift apart on them.
colors = [st.condition_color('Baseline');
          st.condition_color('SFA');
          st.condition_color('STD');
          st.condition_color('SFA+STD')];
if size(colors,1) < n_cond
    colors = lines(n_cond);
end

% Short condition labels for the legend (match the combined figure).
if n_cond == 4
    display_names = {'Baseline', 'SFA', 'STD', 'SFA+STD'};
else
    display_names = condition_names;
end

xmax_s       = 7.5;             % (a)/(b) delay-axis limit (s)
recon_cond   = 4;               % reconstruct the SFA+STD condition
recon_delays = [1, 5, 10, 15];  % hold-delay indices to show (labeled in seconds): 0.3, 1.5, 3.0, 4.5 s at T_hold=0.3
recon_ylim   = [-0.6, 0.6];     % shared y-limits across all reconstruction panels
line_w       = 2;               % (a)/(b) curve width (matches the reference figure)

%% ==================== Combined figure ====================
% 24x2 tiled grid (fine rows so the gap below (a)/(b) is tunable): (a)/(b) span
% the top 8 rows; row 9 is left empty as a spacer; the reconstruction block
% (nested layout) fills rows 10-24. Explicit tile indices keep the spans clear.
fig = figure('Color', 'w', 'Position', [100 237 603 652]);   % compact size (hand-tuned)
tl = tiledlayout(24, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% Leave headroom above the tiles so the (a)/(b) panel letters can sit clearly
% above the axes (and their y-labels) without clipping at the figure top. The
% reduced width leaves room at the right for the (c) amplitude scale-bar label.
tl.OuterPosition = [0 0 0.90 0.95];

% (a) Cumulative MC vs delay
ax_a = nexttile(1, [8 1]); hold on; grid off; box off;
for i = 1:n_cond
    plot(delay_s, cumsum(R2{i}), '-', 'Color', colors(i, :), 'LineWidth', line_w);
end
xlim([0 xmax_s]); ylim([0 11]);
xlabel('Delay (s)', 'FontSize', 15.4); ylabel({'Cumulative', 'Memory Capacity'}, 'FontSize', 15.4);
set(ax_a, 'XTick', [0 2.5 5], 'YTick', [0 5 10]);   % match plot_memory_capacity_combined.m
hold off;

% (b) Per-delay R^2 vs delay -- carries the one legend for the figure
ax_b = nexttile(2, [8 1]); hold on; grid off; box off;
for i = 1:n_cond
    plot(delay_s, R2{i}, '-', 'Color', colors(i, :), 'LineWidth', line_w);
end
xlim([0 xmax_s]); ylim([0 1.005]);
xlabel('Delay (s)', 'FontSize', 15.4); ylabel('$R^2$', 'Interpreter', 'latex', 'FontSize', 15.4);
set(ax_b, 'XTick', [0 2.5 5], 'YTick', [0 0.5 1]);   % match plot_memory_capacity_combined.m
lgd = legend(display_names, 'Location', 'northeast', 'Box', 'off');
% Nudge the legend further up and to the right (into the corner headroom).
lgd.Units = 'normalized';
lp = lgd.Position;
lgd.Position = [lp(1) + 0.02, lp(2) + 0.03, lp(3), lp(4)];
hold off;

% Reconstruction rows (SFA+STD): target u(t-d) vs trained readout, delay in s.
% The x-axis (time) is hidden on every row -- a 15 s scale bar on the bottom row
% conveys the timescale instead.
mcr        = results{recon_cond};
nR         = numel(recon_delays);
scalebar_s  = 15;             % length of the horizontal time scale bar (s)
scalebar_au = 1;              % length of the vertical amplitude scale bar (A.U.)
recon_win   = 60;             % time span shown in each reconstruction panel (s)
recon_ax   = gobjects(1, nR);

% Nested tiledlayout for the (c) rows so their vertical spacing can be tightened
% independently of (a)/(b). It occupies grid rows 10-24 of the outer 24x2 layout
% (row 9 above it is the empty spacer that opens the gap below (a)/(b)); 'tight'
% spacing packs the four panels closer than the outer 'compact' while still
% leaving room for each panel's title.
tl_c = tiledlayout(tl, nR, 1, 'TileSpacing', 'tight', 'Padding', 'none');
tl_c.Layout.Tile     = 19;       % top-left tile of the (c) block (row 10, col 1)
tl_c.Layout.TileSpan = [15 2];   % span rows 10-24 x 2 cols
for k = 1:nR
    d = recon_delays(k);
    p = mcr.predictions(d);
    t_delay = mcr.t_pred(p.t_indices);

    recon_ax(k) = nexttile(tl_c); hold on; grid off; box off;
    plot(t_delay, p.y_true, 'k-', 'LineWidth', 0.9);
    plot(t_delay, p.y_pred, '-', 'Color', colors(recon_cond, :), 'LineWidth', 1.4);
    ylim(recon_ylim);
    % Show only the first recon_win seconds of the test window.
    xlim([t_delay(1), t_delay(1) + recon_win]);
    title(sprintf('delay = %.1f seconds,  R^2 = %.2f', delay_s(d), mcr.R2_d(d)), ...
        'FontWeight', 'normal', 'Interpreter', 'tex');
    % Hide both axes entirely (no ticks, no black axis lines) -- the scale bars
    % convey time (x) and amplitude (y) instead.
    set(gca, 'XColor', 'none', 'YColor', 'none');
    if k == nR
        % L-shaped scale bar in the lower-right corner: a horizontal 15 s time bar
        % and a vertical 1 A.U. amplitude bar meeting at a right angle. Clipping is
        % off so the bars/labels can sit beneath and right of the axis limits.
        xl = get(gca, 'XLim');
        xc = xl(2) + 0.02 * diff(xl);                    % corner x (just right of the data)
        x1 = xc - scalebar_s;                            % left end of the time bar
        yb = recon_ylim(1) - 0.18 * diff(recon_ylim);    % corner y (below the trace)
        ya = yb + scalebar_au;                           % top of the amplitude bar

        % Horizontal time bar + vertical amplitude bar, sharing the corner (xc, yb).
        plot([x1 xc], [yb yb], 'k-', 'LineWidth', 2.5, 'Clipping', 'off');
        plot([xc xc], [yb ya], 'k-', 'LineWidth', 2.5, 'Clipping', 'off');

        % Time label centered below the horizontal bar.
        text((x1 + xc) / 2, yb - 0.06 * diff(recon_ylim), sprintf('%d s', scalebar_s), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 12, 'Clipping', 'off');
        % Amplitude label just right of the vertical bar, rotated -90 deg.
        text(xc + 0.03 * diff(xl), (yb + ya) / 2, sprintf('%g AU', scalebar_au), ...
            'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Clipping', 'off');
    end
    hold off;
end

% Panel letters: (a)/(b) on the analysis panels via AddLetters2Plots. The
% remaining reconstruction rows are self-labeled by their delay/R^2 titles.
% Per-axis VShift: a modest upward nudge for (a)/(b) -- the two-line y-label on
% (a) frees up space, so the letters can sit lower than before.
AddLetters2Plots({ax_a, ax_b}, {'(a)', '(b)'}, ...
    'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.06, 'VShift', [-0.02 -0.02]);

% (c) is lettered manually: recon_ax(1) lives in the nested tiledlayout, whose
% reported position confuses AddLetters2Plots (it lands up by (a)). Place it in
% axis-normalized coordinates at the panel's top-left corner instead.
text(recon_ax(1), -0.025, 1.16, '(c)', 'Units', 'normalized', ...
    'FontSize', 18, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Clipping', 'off');

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_MC_Example';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);


end
end




