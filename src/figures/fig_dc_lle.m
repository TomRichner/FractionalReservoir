function out = fig_dc_lle(cfg)
% FIG_DC_LLE Local Lyapunov exponent vs tonic DC drive, one band per regime.
%
%   out = FIG_DC_LLE()
%   out = FIG_DC_LLE('data_file', f)
%
% Mean local largest Lyapunov exponent as a function of the DC level held on
% every neuron, with a +/- 1 std band across network seeds, drawn once per
% adaptation regime. The dashed line at zero is the edge of chaos: a band
% crossing it means tonic drive moved that regime between stable and chaotic.
%
% Reading it: the QUESTION is whether tonic drive changes stability, and whether
% adaptation changes that answer. A flat band means DC does not matter for that
% regime; separation between bands means the regimes respond differently to the
% same drive.
%
% PLOTTING HALF ONLY. run_dc_lle_analysis runs the staircase across seeds and
% conditions and writes dc_lle_results.mat, so the look can be iterated without
% re-simulating -- which matters here, since a production run is 30 seeds x 7
% conditions.
%
% NOT IN THE PAPER (in_paper = false in paper_config). It is the quantitative
% companion to fig_stim_engages_adaptation, which shows the same effect
% qualitatively on one bursting network.
%
% See also: run_dc_lle_analysis, replot_dc_lle, srnn_condition_titles, confplot

arguments
    cfg.data_file   (1,:) char    = ''    % '' -> search run_dir, then data/dc_lle
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''    % the run whose dc_lle data to plot
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

% The RUN first, then the standalone location -- the convention every figure
% follows since two of them plotted four-day-old data from a sibling .mat.
% run_dc_lle_analysis writes into <run_dir>/dc_lle/dc_lle_nSeeds_*/.
data_file = resolve_data_file(cfg.data_file, ...
    [run_dir_candidates(cfg.run_dir), ...
     dir_candidates(fullfile(project_root, 'data', 'dc_lle'))], ...
    'dc_lle_results.mat', ...
    'Run run_dc_lle_analysis first');
D = load(data_file);
r = D.dc_lle_results;

conditions = r.conditions;
n_cond     = numel(conditions);

%% ---- Display names, from the one shared source ----------------------------
titles = srnn_condition_titles();
labels = conditions;
for c = 1:n_cond
    if ischar(conditions{c}) && isKey(titles, conditions{c})
        labels{c} = titles(conditions{c});
    end
end

%% ---- Colours: the manuscript's condition palette where it knows the name ---
colors = lines(n_cond);
if isfield(st, 'condition_color')
    for c = 1:n_cond
        if isKey(st.condition_color, conditions{c})
            colors(c, :) = st.condition_color(conditions{c});
        end
    end
end

%% ---- Draw ------------------------------------------------------------------
fig = figure('Name', 'DC drive vs local Lyapunov exponent', ...
    'Position', [200, 200, 720, 500]);
ax = axes(fig);
hold(ax, 'on');

% confplot draws the band and the line together but returns nothing useful for
% a legend, so a proxy line per condition carries the legend entry.
h = gobjects(1, n_cond);
for c = 1:n_cond
    band = 1 - 0.35 * (1 - colors(c, :));    % lighter, same hue
    confplot(r.dc_levels, r.level_mean(:, c), r.level_std(:, c), ...
        r.level_std(:, c), [colors(c, :); band]);
    h(c) = plot(ax, nan, nan, '-', 'Color', colors(c, :), 'LineWidth', 2);
end

yline(ax, 0, '--k', 'edge of chaos', 'HandleVisibility', 'off', ...
    'LabelHorizontalAlignment', 'left');
hold(ax, 'off');

legend(ax, h, labels, 'Location', 'best', 'Box', 'off');
xlabel(ax, 'DC level (input units)');
ylabel(ax, 'mean local \lambda_1   (\pm 1 std across seeds)');
title(ax, sprintf('Local Lyapunov exponent vs tonic DC drive  (%d seeds, %s)', ...
    r.n_seeds, strrep(r.config.preset_name, '_', '\_')));
grid(ax, 'off');
box(ax, 'off');

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'fig_dc_lle';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

%% ------------------------------------------------------------------------
function c = run_dir_candidates(run_dir)
% Every dc_lle_nSeeds_* folder under <run_dir>/dc_lle, newest first.
% run_dc_lle_analysis timestamps its own subfolder, so the run directory holds a
% level of nesting that resolve_data_file does not search on its own.
c = {};
if isempty(run_dir); return; end
c = dir_candidates(fullfile(run_dir, 'dc_lle'));
end

function c = dir_candidates(parent)
% parent itself, then its dc_lle_nSeeds_* children, newest child first.
c = {};
if isempty(parent) || ~isfolder(parent); return; end
c = {parent};
d = dir(fullfile(parent, 'dc_lle_nSeeds_*'));
d = d([d.isdir]);
if isempty(d); return; end
[~, ord] = sort([d.datenum], 'descend');
d = d(ord);
c = [c, arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false)'];
end
