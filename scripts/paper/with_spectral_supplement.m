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
%                            knobs cfg.outlier_spectra_*)
%       frequency_response   at outlier_spectra's calibration states (FD-041;
%                            cfg.frequency_response)
%       eigenvalue_density   at the same states (FD-042; cfg.eigenvalue_density)
%     An empty cfg.stages (every stage) is written out first
%     (paper_stage_names), since an empty selection never runs an opt-in
%     stage. The two reading stages' plans refuse a bundle that names them
%     without outlier_spectra, so the three are added together.
%   * appends the figure entries of PARTS to cfg.figures, none in_paper
%     (Tom has not chosen among them; the plan's section 2: supplementary
%     candidates). The density figures are made for two networks: the
%     reference under the function's own name, and E-to-E x 1.5 (the
%     network outlier_spectra's index calls 'EEx1p5') under that name with
%     _EEx1p5 appended, through the figures' 'network' option. Each entry
%     that draws a run's conditions names its own list, the bundle's
%     cfg.conditions that outlier_spectra simulates (its plan's
%     conditions): at 'essential' the four-condition set draws product MTS,
%     which that stage does not run, and an entry's own args override the
%     bundle-wide list (as the single-neuron entry's do).
%
% TO ADD A STAGE OR FIGURES (FD-043 outlier_structures; FD-045's four
% fig_outlier_* on FD-040's output): a row of PARTS. A new stage is a row
% with its name and its entries (the stage must also be a row of
% paper_opt_in_stages); figures of an existing stage are entries in that
% stage's row (FD-045's go in outlier_spectra's).
%
% Errors :AlreadyAdded when cfg.stages already names one of the stages or
% cfg.figures already holds one of the entries, :NoFigures when cfg has
% no registry.
%
% See also: mu7chain_config, paper_opt_in_stages, outlier_spectra_stage_plan,
%           frequency_response_stage_plan, eigenvalue_density_stage_plan,
%           with_switch_spectra, with_ritz_overlay
arguments
    cfg (1,1) struct
end
if ~isfield(cfg, 'figures') || ~iscell(cfg.figures)
    error('with_spectral_supplement:NoFigures', 'The config has no figure registry (cfg.figures).');
end

%% What it adds: one row per stage, {stage, its figure entries}
% Entries are {name, function, extra args, draws_conditions}. The second
% network of the density figures is named as outlier_spectra's index names it.
EE = 'EEx1p5';
dens = {@fig_eigdens_radial, @fig_eigdens_gains, @fig_eigdens_measures, @fig_eigdens_plane};
dens_entries = cell(0, 4);
for k = 1:numel(dens)
    name = func2str(dens{k});
    dens_entries(end + 1, :) = {name, dens{k}, {}, true}; %#ok<AGROW>
    dens_entries(end + 1, :) = {[name '_' EE], dens{k}, {'network', EE}, true}; %#ok<AGROW>
end
PARTS = { ...
    'outlier_spectra',    cell(0, 4); ...
    'frequency_response', { ...
        'fig_freqresp_network',       @fig_freqresp_network,       {}, true; ...
        'fig_freqresp_plane',         @fig_freqresp_plane,         {}, true; ...
        'fig_freqresp_phase',         @fig_freqresp_phase,         {}, true; ...
        'fig_freqresp_single_neuron', @fig_freqresp_single_neuron, {}, true}; ...
    'eigenvalue_density', [dens_entries; ...
        {'fig_eigdens_theory_check',  @fig_eigdens_theory_check,   {}, false}]};

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
% The conditions the figures can draw: the bundle's, as far as
% outlier_spectra simulates them (its plan; no compute).
draws = {};
if isfield(cfg, 'conditions') && ~isempty(cfg.conditions)
    plan = outlier_spectra_stage_plan(cfg);
    draws = intersect(cfg.conditions, plan.conditions, 'stable');
end
for k = 1:size(entries, 1)
    args = entries{k, 3};
    if entries{k, 4} && ~isempty(draws); args = [args, {'conditions', draws}]; end %#ok<AGROW>
    cfg.figures{end + 1} = struct('name', entries{k, 1}, 'fn', entries{k, 2}, 'in_paper', false, 'args', {args});
end
end
