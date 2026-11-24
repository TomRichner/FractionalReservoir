function cfg = with_spectral_supplement(cfg)
% WITH_SPECTRAL_SUPPLEMENT The spectral stages of 2026-11-20 and their figures, added to a bundle config (FD-040 to FD-046).
%
%   cfg = WITH_SPECTRAL_SUPPLEMENT(cfg)
%
% Tom, 2026-11-20: "make the changes (e.g. with more seeds) discussed so
% far, and integrate these figures into the full run"
% (docs/Plan_spectral_supplement.md). mu7chain_config calls this at its
% end, so it is the DEFAULT of the mu7chain bundles in every run mode
% (mu7chain_AC_config inherits it through mu7chain_config); no other bundle
% calls it, and the regression bundle does not get the stages. Called on a
% config by hand it adds the same to that config only; mu7chain_run
% rebuilds the config per figure set, so a config changed by hand does not
% reach its figure passes.
%
% What it does:
%   * appends to cfg.stages, in this order, the opt-in stages of PARTS
%     below (paper_opt_in_stages runs them in its own order, after the
%     first list and lyapunov_spectrum_n80, before block C's sweep):
%       outlier_spectra      frozen-gain and full-Jacobian spectra of the
%                            bundle's network and four E-to-E doses, per
%                            condition and seed (FD-040; sizes by run mode
%                            in outlier_spectra_stage_plan,
%                            knobs cfg.outlier_spectra_*). Its own four
%                            figures (FD-045): fig_outlier_plane, _dose,
%                            _slow_variables, _dwell. The plane and the
%                            slow-variables figure draw the networks
%                            'reference' and 'EEx1p5', the dwell figure the
%                            reference only (Tom, 2026-11-23, FD-057, over
%                            the long window the stage's long_networks
%                            gives it); they need the dense pass and the
%                            full-Jacobian series, which the stage's
%                            default full_networks gives
%                            (a bundle that changes that knob loses them)
%       frequency_response   at outlier_spectra's calibration states (FD-041;
%                            cfg.frequency_response)
%       eigenvalue_density   at the same states (FD-042; cfg.eigenvalue_density)
%       outlier_structures   fifteen networks with an intended eigenvalue
%                            outside the bulk, the mode followed (FD-043;
%                            sizes in outlier_structures_stage_plan, knobs
%                            cfg.outlier_structures_*). Reads nothing of
%                            outlier_spectra's output. Six figures:
%                            fig_ostruct_summary, fig_ostruct_plane once
%                            per group (rotation, local, blockmean; the
%                            entry's name carries the group),
%                            fig_ostruct_rotation, fig_ostruct_h1
%     An empty cfg.stages (every stage) is written out first
%     (paper_stage_names), since an empty selection never runs an opt-in
%     stage. The two reading stages' plans refuse a bundle that names them
%     without outlier_spectra, so the stages are added together.
%   * appends the figure entries of PARTS to cfg.figures, all 'appendix'
%     (until the supplementary composites of Plan 2 are wired in, which
%     reuse them as panels; earlier: none in_paper, the plan's section 2:
%     supplementary
%     candidates). The density figures are made for two networks: the
%     reference under the function's own name, and E-to-E x 1.5 (the
%     network outlier_spectra's index calls 'EEx1p5') under that name with
%     _EEx1p5 appended, through the figures' 'network' option. Each entry
%     that draws a run's conditions names its own list, the bundle's
%     cfg.conditions that outlier_spectra simulates (its plan's
%     conditions): at 'essential' the four-condition set draws product MTS,
%     which that stage does not run, and an entry's own args override the
%     bundle-wide list (as the single-neuron entry's do).
%   * then the COMPOSITES of the figure plan that read these stages (Plan 2,
%     integration, 2026-11-21; COMPOSITES below), each with the same
%     conditions rule, and each only when cfg.stages names every stage it
%     reads (so 'essential', which skips transient_gain, has no Figure 6):
%       fig_main5_grouped   main        fig_main5_mechanism (outlier_spectra)
%       fig_main6_grouped   main        fig_main6_amplification
%                                       (frequency_response; transient_gain
%                                       and eig_heatmap of the same run)
%       fig_supp6_grouped   supplement  eigenvalue_density, and since Plan 3
%                                       frequency_response (S8's two panels
%                                       are its D and E)
%       fig_supp7_grouped   supplement  outlier_spectra, outlier_structures
%       fig_supp8_grouped   appendix    frequency_response. Since Plan 3
%                                       (2026-11-22) its two panels are
%                                       S6 D and E, so it is no longer a
%                                       supplementary figure; it still draws
%                                       them on their own
%
% TO ADD A STAGE OR FIGURES: a row of PARTS. A new stage is a row with its
% name, its entries and the plan that says which conditions it simulates
% (the stage must also be a row of paper_opt_in_stages); figures of an
% existing stage are entries in that stage's row.
%
% Errors :AlreadyAdded when cfg.stages already names one of the stages or
% cfg.figures already holds one of the entries, :NoFigures when cfg has
% no registry.
%
% See also: mu7chain_config, paper_opt_in_stages, outlier_spectra_stage_plan,
%           frequency_response_stage_plan, eigenvalue_density_stage_plan,
%           outlier_structures_stage_plan,
%           with_switch_spectra, with_ritz_overlay
arguments
    cfg (1,1) struct
end
if ~isfield(cfg, 'figures') || ~iscell(cfg.figures)
    error('with_spectral_supplement:NoFigures', 'The config has no figure registry (cfg.figures).');
end

%% What it adds: one row per stage, {stage, its figure entries, conditions plan}
% Entries are {name, function, extra args, draws_conditions}. The second
% network of the density figures is named as outlier_spectra's index names
% it. The third column is the plan whose conditions the row's entries can
% draw (the two reading stages draw outlier_spectra's runs).
EE = 'EEx1p5';
dens = {@fig_eigdens_radial, @fig_eigdens_gains, @fig_eigdens_measures, @fig_eigdens_plane};
dens_entries = cell(0, 4);
for k = 1:numel(dens)
    name = func2str(dens{k});
    dens_entries(end + 1, :) = {name, dens{k}, {}, true}; %#ok<AGROW>
    dens_entries(end + 1, :) = {[name '_' EE], dens{k}, {'network', EE}, true}; %#ok<AGROW>
end
groups = {'rotation', 'local', 'blockmean'};
plane_entries = cell(numel(groups), 4);
for k = 1:numel(groups)
    plane_entries(k, :) = {['fig_ostruct_plane_' groups{k}], @fig_ostruct_plane, {'group', groups{k}}, true};
end
PARTS = { ...
    'outlier_spectra',    { ...
        'fig_outlier_plane',          @fig_outlier_plane,          {}, true; ...
        'fig_outlier_dose',           @fig_outlier_dose,           {}, true; ...
        'fig_outlier_slow_variables', @fig_outlier_slow_variables, {}, true; ...
        'fig_outlier_dwell',          @fig_outlier_dwell,          {}, true}, @outlier_spectra_stage_plan; ...
    'frequency_response', { ...
        'fig_freqresp_network',       @fig_freqresp_network,       {}, true; ...
        'fig_freqresp_plane',         @fig_freqresp_plane,         {}, true; ...
        'fig_freqresp_phase',         @fig_freqresp_phase,         {}, true; ...
        'fig_freqresp_single_neuron', @fig_freqresp_single_neuron, {}, true}, @outlier_spectra_stage_plan; ...
    'eigenvalue_density', [dens_entries; ...
        {'fig_eigdens_theory_check',  @fig_eigdens_theory_check,   {}, false}], @outlier_spectra_stage_plan; ...
    'outlier_structures', [ ...
        {'fig_ostruct_summary',       @fig_ostruct_summary,        {}, true}; ...
        plane_entries; ...
        {'fig_ostruct_rotation',      @fig_ostruct_rotation,       {}, true; ...
         'fig_ostruct_h1',            @fig_ostruct_h1,             {}, true}], @outlier_structures_stage_plan};

%% Stages
stages = {};
if isfield(cfg, 'stages') && ~isempty(cfg.stages); stages = reshape(cellstr(cfg.stages), 1, []); end
if isempty(stages); stages = paper_stage_names(); end
added_stages = PARTS(:, 1)';
already = intersect(added_stages, stages);
if ~isempty(already)
    error('with_spectral_supplement:AlreadyAdded', 'cfg.stages already names %s.', strjoin(already, ', '));
end
cfg.stages = [stages, added_stages];

%% Figure entries
entries = vertcat(PARTS{:, 2});
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
already = intersect(entries(:, 1), names);
if ~isempty(already)
    error('with_spectral_supplement:AlreadyAdded', 'cfg.figures already holds %s.', strjoin(already, ', '));
end
% The conditions an entry can draw: the bundle's, as far as its stage
% simulates them (the stage's plan; no compute).
has_conditions = isfield(cfg, 'conditions') && ~isempty(cfg.conditions);
for p = 1:size(PARTS, 1)
    entries = PARTS{p, 2};
    draws = {};
    if has_conditions && any([entries{:, 4}])
        plan = PARTS{p, 3}(cfg);
        draws = intersect(cfg.conditions, plan.conditions, 'stable');
    end
    for k = 1:size(entries, 1)
        args = entries{k, 3};
        if entries{k, 4} && ~isempty(draws); args = [args, {'conditions', draws}]; end %#ok<AGROW>
        cfg.figures{end + 1} = paper_figure_entry(entries{k, 1}, entries{k, 2}, 'appendix', args);
    end
end

%% The composites of the figure plan that read these stages
% {name, function, group, args, the stages it reads}. Conditions as above,
% from outlier_spectra's plan (every composite here draws its runs, or
% runs made at its calibration states).
COMPOSITES = { ...
    'fig_main5_grouped', @fig_main5_mechanism,     'main',       {'fig_tag', 'Fig_Main5_Grouped'}, {'outlier_spectra'}; ...
    'fig_main6_grouped', @fig_main6_amplification, 'main',       {'fig_tag', 'Fig_Main6_Grouped'}, ...
        {'frequency_response', 'transient_gain', 'eig_heatmap'}; ...
    'fig_supp6_grouped', @fig_supp6_grouped,       'supplement', {}, {'eigenvalue_density', 'frequency_response'}; ...
    'fig_supp7_grouped', @fig_supp7_grouped,       'supplement', {}, {'outlier_spectra', 'outlier_structures'}; ...
    'fig_supp8_grouped', @fig_supp8_grouped,       'appendix',   {}, {'frequency_response'}};
already = intersect(COMPOSITES(:, 1), names);
if ~isempty(already)
    error('with_spectral_supplement:AlreadyAdded', 'cfg.figures already holds %s.', strjoin(already, ', '));
end
draws = {};
if has_conditions
    plan = outlier_spectra_stage_plan(cfg);
    draws = intersect(cfg.conditions, plan.conditions, 'stable');
end
for k = 1:size(COMPOSITES, 1)
    if ~all(ismember(COMPOSITES{k, 5}, cfg.stages)); continue; end
    args = COMPOSITES{k, 4};
    if ~isempty(draws); args = [args, {'conditions', draws}]; end %#ok<AGROW>
    cfg.figures{end + 1} = paper_figure_entry(COMPOSITES{k, 1}, COMPOSITES{k, 2}, COMPOSITES{k, 3}, args);
end
end
