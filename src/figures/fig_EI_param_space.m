function out = fig_EI_param_space(cfg)
% FIG_EI_PARAM_SPACE Param-space distributions, each bar coloured by E:I ratio.
%
%   out = FIG_EI_PARAM_SPACE()
%   out = FIG_EI_PARAM_SPACE('run_dir', d)
%
% ONE SHEET PER MEASURE (sweep_metrics entries with in_sheets set: lambda_1,
% mean rate, h_KS, D_KY, ...), each 1 x (N+1): one column per adaptation
% condition; each bar is a STACK of per-network patches coloured by the
% fraction excitatory (blue = inhibition-dominated, red = excitation-dominated).
% The last column embeds the colorbar, encoding the fraction as an E:I ratio.
% Tags: Fig_EI_ParamSpace_<stem>. Layout lives in ei_metric_sheet, shared with
% fig_EI_weights_param_space.
%
% The sibling of fig_param_space_allStd, which shows the same distributions in
% plain grey. Both are paper figures: this one answers "does the E:I balance
% explain where in the distribution a network lands", which the grey sheet
% cannot.
%
% TWO THINGS CHANGED IN THE PORT.
%
% 1. IT WAS BUILT FROM A DIFFERENT RUN. The hardcoded data_root pointed at
%    run_all_jul_06_26_22_00 -- an old SRNNModel2 run -- while every sibling
%    figure pointed at an SRNNCellTypePairs run. The manuscript was therefore
%    showing two param-space figures computed from different models. It now
%    resolves the same run as everything else.
%
% 2. THE COLOUR AXIS IS NAMED, not defaulted. load_and_make_unit_histograms
%    colours by its ColorBy option, default 'f'. On SRNNCellTypePairs `f` is a
%    1 x C ROW, not a scalar, so the default breaks the per-network colour
%    assignment outright. The Pairs name for the same quantity is the scalar
%    alias f_E, which is exactly f(1) (SRNNCellTypePairs.m get.f_E), so the
%    colormap and the E:I tick labels carry over unchanged -- only the name
%    read off each result differs. It is resolved through
%    psa.effective_param, NOT res.config: effective_param('f') on a Pairs run
%    returns the class default [0.5 0.5] rather than the swept value, which
%    would have coloured every network identically.
%
% See also: fig_param_space_allStd, load_and_make_unit_histograms, resolve_run_dir

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.color_by    (1,:) char    = ''      % '' -> per model class (f_E / f)
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% The colour axis, per model class. get.f_E asserts exactly two cell types, so
% a three-type run must name its own axis rather than fall through to f_E.
color_by = cfg.color_by;
if isempty(color_by)
    [~, model_class] = srnn_param_preset(cfg.preset_name);
    if strcmp(model_class, 'SRNNCellTypePairs')
        color_by = 'f_E';
    else
        color_by = 'f';
    end
end

ps_dirs = dir(fullfile(run_dir, 'param_space_*'));
ps_dirs = ps_dirs([ps_dirs.isdir]);
assert(~isempty(ps_dirs), 'No param_space_* subdir found in %s', run_dir);
ps_dir  = fullfile(ps_dirs(1).folder, ps_dirs(1).name);

% Start from a clean slate.
% (no 'close all force' -- destroyed sibling figures in a batch; see header)

% 1) Build the f-coloured unit-histogram figures, ONE PER REGISTRY MEASURE
%    with in_sheets set (sweep_metrics: LLE, mean rate, h_KS, D_KY, ...),
%    matching the grey figure's LLE range [-1.5, 1.5] and probability
%    normalization. This also creates a separate 'f Value Colorbar' figure per
%    call, so the loader is called once and the figures are picked up by Name.
% CLim is FIXED at the swept range [0.2, 0.8] = 1:4 to 4:1, not derived from the
% data. Derived limits came from the min/max of the f_E values actually sampled,
% which was fine when the grid enumerated every level but is luck once the grid
% is randomly subsampled: a run that never draws f_E = 0.2 loses the 1:4 tick,
% the bar's span changes, and two runs' colours stop meaning the same thing. The
% axis is a property of the SWEEP, so it belongs in the script, not in the data.
EI_CLIM = [0.2, 0.8];
specs = sweep_metrics();
specs = specs([specs.in_sheets]);
[~, ~] = load_and_make_unit_histograms(ps_dir, ...
    'Metrics', {specs.key}, 'NormalizeMode', 'probability', 'LLERange', [-1.5, 1.5], ...
    'ColorBy', color_by, 'CLim', EI_CLIM);

% 2) One styled sheet per measure (ei_metric_sheet holds the layout that the
%    two E:I figures share). The colorbar figure is consumed by the first
%    sheet; later sheets get a fresh copy of it.
opts = struct('tick_fs', st.tick_fs, 'label_fs', st.label_fs, 'title_fs', 20, ...
    'axes_lw', 1.0, 'letter_fs', 18, 'row_shrink', 0.85, 'top_headroom', 0.06, ...
    'title_y', 1.22, 'cb_x_shift', 0.045, 'xlabel', '', 'yticks', [], ...
    'zero_color', [0 0.7 0], 'cb_clim', EI_CLIM, ...
    'cb_ticks', [0.2, 0.25, 1/3, 0.4, 0.5, 0.6, 2/3, 0.75, 0.8], ...
    'cb_labels', {{'1:4', '1:3', '1:2', '2:3', '1:1', '3:2', '2:1', '3:1', '4:1'}}, ...
    'cb_ylabel', 'E:I ratio');
cb_fig = findobj(0, 'Type', 'figure', 'Name', 'f Value Colorbar');
figs = gobjects(1, numel(specs));
tags = cell(1, numel(specs));
for mi = 1:numel(specs)
    spec = specs(mi);
    src_fig = findobj(0, 'Type', 'figure', 'Name', sprintf('%s Unit Histogram', spec.field));
    assert(isscalar(src_fig), 'fig_EI_param_space:MissingFigure', ...
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
    tags{mi} = sprintf('Fig_EI_ParamSpace_%s', spec.stem);
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
