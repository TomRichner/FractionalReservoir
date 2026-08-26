function out = fig_EI_weights_param_space(cfg)
% FIG_EI_WEIGHTS_PARAM_SPACE Param-space distributions coloured by E:I WEIGHT balance.
%
%   out = FIG_EI_WEIGHTS_PARAM_SPACE()
%   out = FIG_EI_WEIGHTS_PARAM_SPACE('run_dir', d)
%
% The same sheet as fig_EI_param_space, coloured by a different quantity.
%
% WHY IT EXISTS. fig_EI_param_space colours each network by f_E, the fraction of
% neurons that are excitatory. That was the whole story when the only thing
% varying was how many E neurons there were. It is not the whole story now: the
% joint grid also sweeps the four connectivity blocks mu_EE / mu_EI / mu_IE /
% mu_II over -75% to +200%, so two networks with identical f_E can sit at
% opposite ends of the excitation-inhibition balance. A network can be 80%
% excitatory and still inhibition-dominated if its inhibitory synapses are three
% times as strong.
%
% So this figure colours by the balance of the WEIGHTS actually drawn:
%
%   w_frac_E = sum(W(:, E)) / ( sum(W(:, E)) + |sum(W(:, I))| )
%
% summing over presynaptic columns, i.e. total excitatory drive against total
% inhibitory drive. Like f_E it lives in (0, 1) and reads as an E:I ratio, so the
% two sheets use the same colormap and the same tick vocabulary and can be read
% against each other directly.
%
% MEASURED on the medium run of 2026-08-26 (64 grid points): w_frac_E spans
% [0.074, 0.980] where f_E spans [0.2, 0.8], and the two correlate at only 0.675.
% The weight balance reaches roughly 1:12 to 48:1 where the neuron-count balance
% is capped at 1:4 to 4:1 by construction. That gap is the figure's entire point.
%
% WHERE THE WEIGHTS COME FROM. A stored result records scalars (LLE, mean rate)
% and never W, so each grid point's network is REBUILT: psa.rebuild_model(res)
% reconstructs the constructor call the sweep made, including
% rng_seeds = [network_seed, network_seed + 1], and build() then redraws exactly
% the same W. It costs about 9 s for 64 points and is cached per config_idx,
% because the network is a property of the grid point and is shared by all the
% adaptation conditions run at it.
%
% Note the ratio is invariant to both level_of_chaos and rescale_by_abscissa:
% each is a scalar multiplying all of W, so it cancels between numerator and
% denominator. What moves this axis is f and the mu blocks.
%
% See also: fig_EI_param_space, ParamSpaceAnalysis2/rebuild_model,
%           load_and_make_unit_histograms, resolve_run_dir

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

[~, model_class] = srnn_param_preset(cfg.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('fig_EI_weights_param_space:WrongModelClass', ...
        ['This figure needs SRNNCellTypePairs (it indexes W by cell type ' ...
         'through type_indices); preset ''%s'' is a %s preset.'], ...
        cfg.preset_name, model_class);
end

ps_dirs = dir(fullfile(run_dir, 'param_space_*'));
ps_dirs = ps_dirs([ps_dirs.isdir]);
assert(~isempty(ps_dirs), 'No param_space_* subdir found in %s', run_dir);
ps_dir  = fullfile(ps_dirs(1).folder, ps_dirs(1).name);

% The colour axis. FIXED at [0.05, 0.95] = 1:19 to 19:1, chosen from the measured
% spread rather than from whatever this particular run sampled -- the same reason
% fig_EI_param_space pins its own bar. It is wider than that figure's [0.2, 0.8]
% because the weight balance genuinely reaches further than the neuron-count
% balance can; the shared 1:4, 1:2, 1:1, 2:1, 4:1 ticks are what let the two
% sheets still be compared. A handful of networks at the very top of the range
% clamp to the end colour, which load_and_make_unit_histograms warns about.
W_CLIM = [0.05, 0.95];

% Rebuilding a network is the expensive part, so cache it. The network depends
% only on the grid POSITION (network_seed = config_idx*100 + offset), so all
% seven adaptation conditions at one grid point share a single rebuild.
cache = containers.Map('KeyType', 'double', 'ValueType', 'double');
color_fcn = @(psa, res) ei_weight_fraction(psa, res, cache);

fprintf('[fig_EI_weights_param_space] rebuilding networks to weigh E against I...\n');
t_rebuild = tic;
[~, ~] = load_and_make_unit_histograms(ps_dir, ...
    'Metrics', {'lle', 'r'}, 'NormalizeMode', 'probability', 'LLERange', [-1.5, 1.5], ...
    'ColorBy', 'E:I weight balance', 'ColorFcn', color_fcn, 'CLim', W_CLIM, ...
    'ColorLabel', 'excitatory weight fraction');
fprintf('[fig_EI_weights_param_space] %d networks rebuilt in %.1f s\n', ...
    cache.Count, toc(t_rebuild));

lle_fig = findobj(0, 'Type', 'figure', 'Name', 'LLE Unit Histogram');
mr_fig  = findobj(0, 'Type', 'figure', 'Name', 'mean_rate Unit Histogram');
cb_fig  = findobj(0, 'Type', 'figure', 'Name', 'f Value Colorbar');

lle_ax = sort_axes_left_to_right(lle_fig);
mr_ax  = sort_axes_left_to_right(mr_fig);

%% Combined sheet: conditions across, LLE on top and mean rate below, colorbar last
nCols = numel(lle_ax);
nRows = 2;
nGrid = nCols + 1;
src   = {lle_ax, mr_ax};
combined = figure('Color', 'w', 'Position', [100 100 350*nGrid 300*nRows]);
cax = gobjects(nRows, nCols);
for r = 1:nRows
    for c = 1:nCols
        ph = subplot(nRows, nGrid, (r-1)*nGrid + c, 'Parent', combined);
        target_pos = get(ph, 'Position');
        delete(ph);
        cax(r, c) = copyobj(src{r}(c), combined);
        set(cax(r, c), 'Position', target_pos);
    end
end

cbax = gobjects(0);
if ~isempty(cb_fig)
    ph = subplot(nRows, nGrid, nGrid, 'Parent', combined);
    cb_target_pos = get(ph, 'Position');
    delete(ph);
    src_cb = findobj(cb_fig, 'Type', 'axes');
    cbax = copyobj(src_cb(1), combined);
    set(cbax, 'Position', cb_target_pos);
end
close(lle_fig);
close(mr_fig);
if ~isempty(cb_fig); close(cb_fig); end

%% Styling -- kept in step with fig_EI_param_space so the pair reads as a pair
tick_fs  = st.tick_fs;
label_fs = st.label_fs;
title_fs     = 20;
axes_lw      = 1.0;
letter_fs    = 18;
row_shrink   = 0.85;
top_headroom = 0.06;
title_y      = 1.22;
lle_yticks   = [0, 0.5];
rate_yticks  = [0, 0.3];
for r = 1:nRows
    for c = 1:nCols
        ax = cax(r, c);
        set(ax, 'FontSize', tick_fs, 'LineWidth', axes_lw);
        set(ax.YLabel, 'FontSize', label_fs);
        if r == 1
            xlabel(ax, 'Growth Rate', 'FontSize', label_fs);
            set(ax.Title, 'FontWeight', 'normal', 'FontSize', title_fs);
            set(findobj(ax, 'Type', 'constantline'), 'Color', [0 0.7 0]);
            set(ax, 'YTick', lle_yticks);
        else
            set(ax.XLabel, 'FontSize', label_fs);
            title(ax, '');
            set(ax, 'YTick', rate_yticks);
        end
    end
    linkaxes(cax(r, :), 'y');
end

for r = 1:nRows
    for c = 1:nCols
        ax = cax(r, c);
        p  = get(ax, 'Position');
        new_h = p(4) * row_shrink;
        set(ax, 'Position', [p(1), p(2) + (p(4) - new_h) - top_headroom, p(3), new_h]);
    end
end
for c = 1:nCols
    t = get(cax(1, c), 'Title');
    if ~isempty(get(t, 'String'))
        set(t, 'Units', 'normalized', 'Position', [0.5, title_y, 0], ...
            'VerticalAlignment', 'bottom', 'FontSize', title_fs);
    end
end
cb_x_shift = 0.045;
if isgraphics(cbax)
    p = get(cbax, 'Position');
    new_h = p(4) * row_shrink;
    set(cbax, 'Position', [p(1) - cb_x_shift, p(2) + (p(4) - new_h) - top_headroom, p(3), new_h]);
end

%% Vertical dividers between condition columns
pos = cell2mat(get(cax(:), 'Position'));
[~, ~, col_of] = uniquetol(pos(:,1), 0.01);
ncol      = max(col_of);
col_left  = accumarray(col_of, pos(:,1),          [ncol 1], @mean);
col_right = accumarray(col_of, pos(:,1)+pos(:,3), [ncol 1], @mean);
[col_left, ord] = sort(col_left);
col_right = col_right(ord);
y_bot = min(pos(:,2));
y_top = max(pos(:,2) + pos(:,4));
x_shift = 0.012;
for c = 1:ncol-1
    x_div = (col_right(c) + col_left(c+1)) / 2 - x_shift;
    annotation(combined, 'line', [x_div x_div], [y_bot y_top], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 2);
end

letter_axes = cell(1, numel(cax));
k = 0;
for r = 1:nRows
    for c = 1:nCols
        k = k + 1;
        letter_axes{k} = cax(r, c);
    end
end
AddLetters2Plots(letter_axes, panel_letters(numel(cax)), ...
    'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.06);

%% Colorbar ticks as E:I ratios
% 1:4 through 4:1 are shared with fig_EI_param_space on purpose -- they are the
% anchor that lets the two sheets be compared -- with 1:19 and 19:1 added at the
% ends because the weight balance reaches that far and the neuron-count balance
% cannot.
if isgraphics(cbax)
    ei_f   = [0.05,  0.2,   1/3,   0.5,   2/3,   0.8,   0.95];
    ei_lab = {'1:19', '1:4', '1:2', '1:1', '2:1', '4:1', '19:1'};
    ylim_cb = get(cbax, 'YLim');
    keep = ei_f >= ylim_cb(1) - 1e-6 & ei_f <= ylim_cb(2) + 1e-6;
    assert(all(keep), ['E:I weight colorbar lost %d tick(s): CLim is fixed at ' ...
        '[%g %g] so every label should fit.'], sum(~keep), W_CLIM(1), W_CLIM(2));
    set(cbax, 'YTick', ei_f(keep), 'YTickLabel', ei_lab(keep), 'FontSize', tick_fs);
    ylabel(cbax, 'E:I weight ratio', 'FontSize', label_fs);
end

if ~cfg.visible; set(combined, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_EI_Weights_ParamSpace';
out = struct('figs', combined, 'files', {{}}, 'source', ps_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, combined);
    out.files = existing_outputs(out_dir, fig_tag);
    capture_git_provenance(out_dir, project_root);

    [~, ps_name] = fileparts(ps_dir);
    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_EI_weights_param_space', ...
        'title',  'Stability_Manuscript figure: parameter-space distributions, E:I WEIGHT coloured', ...
        'script', 'fig_EI_weights_param_space.m', ...
        'what',   ['The same sheet as fig_EI_param_space -- one column per ' ...
                   'adaptation condition, row 1 the LLE distribution (green ' ...
                   'dashed zero line), row 2 the mean firing rate, each bar a ' ...
                   'stack of per-network patches -- but coloured by the balance ' ...
                   'of the synaptic WEIGHTS rather than by the fraction of ' ...
                   'neurons that are excitatory. The colorbar in the final ' ...
                   'column is labelled as an E:I weight ratio.'], ...
        'how',    ['Presentation replot; no simulation is re-run, but each grid ' ...
                   'point''s NETWORK is rebuilt. A stored result records only ' ...
                   'scalars, never W, so psa.rebuild_model(res) reconstructs the ' ...
                   'constructor call the sweep made -- including ' ...
                   'rng_seeds = [network_seed, network_seed + 1] -- and build() ' ...
                   'redraws the identical W. The colour value is ' ...
                   'sum(W(:,E)) / (sum(W(:,E)) + |sum(W(:,I))|), summed over ' ...
                   'presynaptic columns. Rebuilds are cached per config_idx, ' ...
                   'since the network belongs to the grid point and is shared by ' ...
                   'every condition run there.'], ...
        'source', struct('run_dir', run_dir, 'param_space_subfolder', ps_name, ...
                         'preset', cfg.preset_name, ...
                         'clim', sprintf('[%g %g]', W_CLIM(1), W_CLIM(2))), ...
        'figures', {out.files}, ...
        'notes',  {{['WHY NOT f_E. Its sibling fig_EI_param_space colours by the ' ...
                    'fraction of neurons that are excitatory, which told the ' ...
                    'whole story while that was the only thing varying. The ' ...
                    'joint grid now also sweeps mu_EE / mu_EI / mu_IE / mu_II ' ...
                    'over -75% to +200%, so two networks with the same f_E can ' ...
                    'sit at opposite ends of the E:I balance: an 80% excitatory ' ...
                    'network with inhibitory synapses three times as strong is ' ...
                    'inhibition-dominated.'], ...
                   ['MEASURED on the medium run of 2026-08-26, 64 grid points: ' ...
                    'the weight fraction spans [0.074, 0.980] where f_E spans ' ...
                    '[0.2, 0.8], and the two correlate at only 0.675. Roughly ' ...
                    '1:12 to 48:1 in weights against a construction-capped 1:4 ' ...
                    'to 4:1 in neuron counts.'], ...
                   ['COLOUR LIMITS are fixed at [0.05, 0.95] rather than taken ' ...
                    'from the data, so the bar means the same thing in every ' ...
                    'run -- with a randomly subsampled grid, derived limits ' ...
                    'depend on which networks happened to be drawn. A few ' ...
                    'networks above 0.95 clamp to the end colour; ' ...
                    'load_and_make_unit_histograms warns when that happens.'], ...
                   ['The ratio is invariant to level_of_chaos and to ' ...
                    'rescale_by_abscissa: both are scalars multiplying all of W ' ...
                    'and cancel between numerator and denominator. What moves ' ...
                    'this axis is f and the four mu blocks.']}}));
end
end

%% ------------------------------------------------------------------------
function v = ei_weight_fraction(psa, res, cache)
% Excitatory share of total synaptic weight for the network behind one result.
%
% Cached on config_idx: the network is pinned by the grid position, so every
% adaptation condition at a given point shares one rebuild.
key = res.config_idx;
if isKey(cache, key)
    v = cache(key);
    return
end

m = psa.rebuild_model(res);
% build() prints its spectral radius and dead-state count per call; silenced
% here because this runs once per grid point and would otherwise bury the
% figure's own output.
evalc('m.build()');

ti = m.type_indices;
assert(numel(ti) >= 2, ...
    'ei_weight_fraction needs at least two cell types; got %d.', numel(ti));

% Sum over PRESYNAPTIC columns: total excitatory drive vs total inhibitory
% drive. Type 1 is E and type 2 is I by the class's own convention, which
% SRNNCellTypePairs enforces through the f_E / mu_*_relative aliases.
S_E = full(sum(sum(m.W(:, ti{1}))));
S_I = full(sum(sum(m.W(:, ti{2}))));

denom = abs(S_E) + abs(S_I);
if denom == 0
    v = 0.5;   % a zero W has no balance; centre it rather than divide by zero
else
    v = abs(S_E) / denom;
end
cache(key) = v; %#ok<NASGU>  handle object: this mutates the caller's map
end
