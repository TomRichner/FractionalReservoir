function out = fig_EI_weights_param_space(cfg)
% FIG_EI_WEIGHTS_PARAM_SPACE Param-space distributions coloured by E:I WEIGHT balance.
%
%   out = FIG_EI_WEIGHTS_PARAM_SPACE()
%   out = FIG_EI_WEIGHTS_PARAM_SPACE('run_dir', d)
%
% The same sheets as fig_EI_param_space (one per registry measure, tags
% Fig_EI_Weights_ParamSpace_<stem>), coloured by a different quantity.
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

% The colour axis. FIXED, not derived from whatever this run sampled -- the same
% reason fig_EI_param_space pins its own bar. It is wider than that figure's
% [0.2, 0.8] because the weight balance genuinely reaches further than the
% neuron-count balance can; the shared 1:4, 1:2, 1:1, 2:1, 4:1 ticks are what
% let the two sheets still be compared.
%
% 1:10 to 10:1, i.e. [1/11, 10/11]. Was 1:19 to 19:1, which covered nearly the
% whole measured spread ([0.074, 0.980], about 1:12 to 48:1) but spent most of
% the bar's dynamic range on a handful of extreme networks, flattening the
% colour contrast across the bulk that sits nearer 1:1. Tightening to 1:10 gives
% the crowded middle more of the colormap; the cost is that the outliers clamp
% to the end colours and read as "at least 10:1" rather than showing how far
% past it they go. load_and_make_unit_histograms warns with the count whenever
% that happens, so the clamping is never silent.
W_CLIM = [1/11, 10/11];

% Rebuilding a network is the expensive part, so cache it. The network depends
% only on the grid POSITION (network_seed = config_idx*100 + offset), so all
% seven adaptation conditions at one grid point share a single rebuild.
cache = containers.Map('KeyType', 'double', 'ValueType', 'double');
color_fcn = @(psa, res) ei_weight_fraction(psa, res, cache);

fprintf('[fig_EI_weights_param_space] rebuilding networks to weigh E against I...\n');
t_rebuild = tic;
specs = sweep_metrics();
specs = specs([specs.in_sheets]);
[~, ~] = load_and_make_unit_histograms(ps_dir, ...
    'Metrics', {specs.key}, 'NormalizeMode', 'probability', 'LLERange', [-1.5, 1.5], ...
    'ColorBy', 'E:I weight balance', 'ColorFcn', color_fcn, 'CLim', W_CLIM, ...
    'ColorLabel', 'excitatory weight fraction');
fprintf('[fig_EI_weights_param_space] %d networks rebuilt in %.1f s\n', ...
    cache.Count, toc(t_rebuild));

% One styled sheet per registry measure (ei_metric_sheet holds the layout
% shared with fig_EI_param_space); the colorbar is copied for each sheet.
opts = struct('tick_fs', st.tick_fs, 'label_fs', st.label_fs, 'title_fs', 20, ...
    'axes_lw', 1.0, 'letter_fs', 18, 'row_shrink', 0.85, 'top_headroom', 0.06, ...
    'title_y', 1.22, 'cb_x_shift', 0.045, 'xlabel', '', 'yticks', [], ...
    'zero_color', [0 0.7 0], 'cb_clim', W_CLIM, ...
    'cb_ticks', [1/11, 0.2, 1/3, 0.5, 2/3, 0.8, 10/11], ...
    'cb_labels', {{'1:10', '1:4', '1:2', '1:1', '2:1', '4:1', '10:1'}}, ...
    'cb_ylabel', 'E:I weight ratio');
cb_fig = findobj(0, 'Type', 'figure', 'Name', 'f Value Colorbar');
figs = gobjects(1, numel(specs));
tags = cell(1, numel(specs));
for mi = 1:numel(specs)
    spec = specs(mi);
    src_fig = findobj(0, 'Type', 'figure', 'Name', sprintf('%s Unit Histogram', spec.field));
    assert(isscalar(src_fig), 'fig_EI_weights_param_space:MissingFigure', ...
        'Expected one "%s Unit Histogram" figure, found %d.', spec.field, numel(src_fig));
    o = opts;
    if strcmp(spec.field, 'LLE')
        o.xlabel = 'Growth Rate';  o.yticks = [0, 0.5];
    elseif strcmp(spec.field, 'mean_rate')
        o.yticks = [0, 0.3];
    end
    cb_copy = gobjects(0);
    if isgraphics(cb_fig)
        cb_copy = copyobj(cb_fig, 0); set(cb_copy, 'Visible', 'off');
    end
    figs(mi) = ei_metric_sheet(src_fig, cb_copy, spec, o);
    tags{mi} = sprintf('Fig_EI_Weights_ParamSpace_%s', spec.stem);
    if ~cfg.visible; set(figs(mi), 'Visible', 'off'); end
end
if isgraphics(cb_fig); close(cb_fig); end

out = struct('figs', figs, 'files', {{}}, 'source', ps_dir);
if cfg.save
    for mi = 1:numel(specs)
        save_figure_stable(out_dir, tags{mi}, figs(mi));
        out.files = [out.files, existing_outputs(out_dir, tags{mi})];
    end
end
end

function v = ei_weight_fraction(psa, res, cache)
% The realized E:I weight balance of the network at this grid point, cached
% per config_idx (the network is shared by every condition run at a point).
key = res.config_idx;
if isKey(cache, key)
    v = cache(key);
    return
end
m = psa.rebuild_model(res);
evalc('m.build()');
ti = m.type_indices;
assert(numel(ti) >= 2, ...
    'ei_weight_fraction needs at least two cell types; got %d.', numel(ti));
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
