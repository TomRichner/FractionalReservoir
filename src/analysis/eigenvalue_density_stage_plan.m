function plan = eigenvalue_density_stage_plan(cfg)
% EIGENVALUE_DENSITY_STAGE_PLAN What a bundle's eigenvalue_density stage will do, validated before any compute.
%
%   plan = EIGENVALUE_DENSITY_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'eigenvalue_density' (FD-042;
% docs/Plan_spectral_supplement.md). Config in, plan out: nothing is
% computed and no file is touched. run_eigenvalue_density_stage carries
% the plan out.
%
% THE QUESTION. Harris et al. (2023) show that unequal column variances
% concentrate the eigenvalues of a random matrix toward the center of its
% disk. In the frozen-gain effective connectivity M = W diag(d) the
% unevenness comes from the neurons' gains d (W alone has a uniform
% density in the paper's network), and adaptation flattens the density by
% un-pinning neurons. The stage measures the radial distribution of the
% spectrum of M at the calibration states of every run of the stage
% outlier_spectra, sets it beside the random-matrix theory evaluated with
% the saved gains, and takes the counts near and beyond the stability
% boundary (eigdens_run_measures). Write-up:
% ../StochasticPlasticDynamicalSystemPaper/docs/future_directions/
% eigenvalue_density_and_expressivity/Report.md.
%
% THE INPUT is the output of the opt-in stage outlier_spectra (FD-040) in
% the SAME run directory: every network, condition and seed its index
% lists, read through eigdens_load_runs and nothing else. So a bundle that
% names this stage must name outlier_spectra too (checked here, so a
% missing input fails in the pre-flight and not after hours of compute;
% the registry of opt-in stages runs outlier_spectra first). The stage has
% no size of its own: the networks, seeds and states are outlier_spectra's
% (its table by run mode), and the cost scales with them (about 1 s per
% run at 500 neurons and 12 states; see run_eigenvalue_density_stage).
%
% THE CONFIG
%   cfg.run_mode     must be one of run_mode_names (the stage is the same
%                    in every mode)
%   cfg.stages       must name 'outlier_spectra' (see above)
%   cfg.eigenvalue_density   optional scalar struct; no field is accepted
%                    today (a field is refused, so a typo cannot be ignored)
%
% THE PLAN
%   plan.folder        where it writes, relative to the run directory
%                      ('eigenvalue_density')
%   plan.input_folder  where it reads ('outlier_spectra')
%   plan.data_file     the file it writes, eigenvalue_density_data.mat
%   plan.grids         the radial grids (u, u2_edges, u_theory)
%   plan.measures      eigdens_measures('settings')
%
% See also: run_eigenvalue_density_stage, eigdens_load_runs,
%           eigdens_run_measures, fig_eigdens_radial, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

ALLOWED = {};
if ~isfield(cfg, 'run_mode') || ~(ischar(cfg.run_mode) || isstring(cfg.run_mode)) || ...
        ~ismember(char(cfg.run_mode), run_mode_names())
    bad_config('The eigenvalue_density stage needs cfg.run_mode, one of %s.', strjoin(run_mode_names(), ', '));
end
stages = {};
if isfield(cfg, 'stages') && ~isempty(cfg.stages); stages = cellstr(cfg.stages); end
if ~ismember('outlier_spectra', stages)
    bad_config(['The eigenvalue_density stage reads the output of the stage outlier_spectra in the same ' ...
        'run directory, and cfg.stages does not name it.']);
end
opts = struct();
if isfield(cfg, 'eigenvalue_density') && ~isempty(cfg.eigenvalue_density)
    opts = cfg.eigenvalue_density;
end
if ~(isstruct(opts) && isscalar(opts))
    bad_config('cfg.eigenvalue_density must be a scalar struct.');
end
unknown = setdiff(fieldnames(opts), ALLOWED);
if ~isempty(unknown)
    bad_config('cfg.eigenvalue_density has unknown field(s) %s; it takes none.', strjoin(unknown', ', '));
end

plan = struct();
plan.folder       = 'eigenvalue_density';
plan.input_folder = 'outlier_spectra';
plan.data_file    = 'eigenvalue_density_data.mat';
plan.grids        = struct('u', 0:0.01:1.4, 'u2_edges', 0:0.04:1.44, 'u_theory', sqrt(linspace(0, 1, 101)));
plan.measures     = eigdens_measures('settings');
end

%% ------------------------------------------------------------------------
function bad_config(varargin)
error('eigenvalue_density_stage_plan:BadConfig', varargin{:});
end
