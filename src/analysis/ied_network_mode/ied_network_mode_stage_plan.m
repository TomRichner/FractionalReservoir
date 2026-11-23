function plan = ied_network_mode_stage_plan(cfg)
%IED_NETWORK_MODE_STAGE_PLAN The network-mode stage of the ied_stim bundle, validated before any compute.
%
%   plan = IED_NETWORK_MODE_STAGE_PLAN(cfg)      cfg = ied_stim_config(mode)
%
% The pre-flight half of the ied_stim bundle's stage 'network_mode' (Plan 2,
% package E; an appendix figure, fig_ied_network_mode, kept separate at
% Tom's word: "make the I and J ... keep it separate as I'm not sure we will
% keep it"). Config in, plan out; nothing computed, no file touched.
% run_ied_network_mode_stage carries it out.
%
% THE QUESTION (Report 5.1, second pilot): in the modular network of
% Figure 8 G to I, does the rightmost NETWORK mode of the full Jacobian
% separate baseline from stimulation, and does it rise before discharges?
% The frozen-gain boundary is no use there (the bulk of M already reaches
% past it), so only the full Jacobian is computed.
%
% WHAT IS RUN, per network seed of the plan: the bundle's first
% no-stimulation block (block 1) and first stimulation block (block
% n_no + 1), re-simulated exactly as ied_stim_blocks runs them (same
% network, noise seed base.noise_seed + (seed - 42) + 1000 * block, the
% whole block, the bundle's detector), with the full state kept. Where the
% bundle's block file exists the re-simulation is compared with it (mean E
% dendritic potential). Then the full Jacobian (1190 state variables) is
% decomposed:
%   - every full_dt s over the first window_s s of each block's analysed
%     part: the rightmost network mode (>= x_thr of the eigenvector's
%     norm^2 in the dendritic rows), its coherence over the ten subnet
%     populations, its share per subnet, the number of unstable network
%     modes;
%   - around the first n_events discharges of the baseline block and as
%     many control times (drawn with rng(seed), at least ctrl_gap s from
%     every discharge), at the lags lags (s) from the discharge peak.
% A DENSE eig is used (full_jacobian_network_mode): ospec_full_spectrum
% switches to eigs above 1000 state variables, and at 1190 its iterative
% solver converged on 2, 12 and 4 of 12 eigenvalues in the pilot, none of
% them a network mode.
%
% THE SIZES, by the bundle's run mode (Tom: the small modes stay short):
%   mode                   networks   window_s   n_events
%   essential                 1          20          3
%   fast                      3          60          6
%   production, stim_blocks_20261117
%                            12         120          6
% COST, measured on the workstation (2026-11-21, 2 pool workers): see
% run_ied_network_mode_stage; per network it is the two block simulations
% (0.31 s per simulated second) plus 2 window_s / full_dt + 2 n_events
% numel(lags) dense eigendecompositions of the 1190 x 1190 Jacobian.
%
% THE KNOBS, cfg.network_mode.<name> (all optional): seeds (default the
% bundle's first n), window_s, full_dt (1 s), n_events, lags
% (-0.6:0.05:0.2), x_thr (0.1), ctrl_gap (1 s).
%
% THE PLAN: plan.folder 'network_mode' (under cfg.run_dir), plan.runs
% (seed, file), the resolved knobs, plan.sim (cfg.sim), plan.run_dir.
%
% See also: run_ied_network_mode_stage, fig_ied_network_mode,
%           ied_stim_config, ied_stim_blocks, full_jacobian_network_mode

arguments
    cfg (1,1) struct
end
SIZES = {'essential', 1, 20, 3; 'fast', 3, 60, 6; 'production', 12, 120, 6; 'stim_blocks_20261117', 12, 120, 6};
KNOBS = {'seeds', 'window_s', 'full_dt', 'n_events', 'lags', 'x_thr', 'ctrl_gap'};
if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, SIZES(:, 1))
    bad_config('The network_mode stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(SIZES(:, 1)', ', '));
end
if ~isfield(cfg, 'sim') || ~all(isfield(cfg.sim, {'seeds', 'amplitude', 'n_no', 'n_stim', 'block_len', 'settle', 'fs'}))
    bad_config('cfg.sim must be the ied_stim bundle''s experiment (ied_stim_config).');
end
if ~isfield(cfg, 'run_dir') || isempty(cfg.run_dir)
    bad_config('cfg.run_dir is not set.');
end
K = struct();
if isfield(cfg, 'network_mode') && ~isempty(cfg.network_mode); K = cfg.network_mode; end
if ~isstruct(K); bad_config('cfg.network_mode must be a struct of knobs.'); end
unknown = setdiff(fieldnames(K), KNOBS);
if ~isempty(unknown)
    bad_config('Unknown setting(s) cfg.network_mode.%s. The knobs are %s.', strjoin(unknown', ', '), strjoin(KNOBS, ', '));
end
row = SIZES(strcmp(SIZES(:, 1), cfg.run_mode), :);
sim = cfg.sim;
P = struct();
P.seeds    = knob(K, 'seeds', sim.seeds(1:min(row{2}, numel(sim.seeds))));
P.window_s = knob(K, 'window_s', min(row{3}, sim.block_len));
P.full_dt  = knob(K, 'full_dt', 1);
P.n_events = knob(K, 'n_events', row{4});
P.lags     = knob(K, 'lags', -0.6:0.05:0.2);
P.x_thr    = knob(K, 'x_thr', 0.1);
P.ctrl_gap = knob(K, 'ctrl_gap', 1);

if ~(isnumeric(P.seeds) && isvector(P.seeds) && all(ismember(P.seeds, sim.seeds)) && numel(unique(P.seeds)) == numel(P.seeds))
    bad_config('cfg.network_mode.seeds must be distinct network seeds of the bundle (%s).', mat2str(sim.seeds));
end
if ~(isscalar(P.window_s) && P.window_s > 0 && P.window_s <= sim.block_len)
    bad_config('cfg.network_mode.window_s must be positive and at most the block length (%g s).', sim.block_len);
end
if ~(isscalar(P.full_dt) && P.full_dt >= 1 / sim.fs && P.full_dt <= P.window_s)
    bad_config('cfg.network_mode.full_dt must lie between 1/fs and window_s.');
end
if ~(isscalar(P.n_events) && P.n_events >= 1 && P.n_events == round(P.n_events))
    bad_config('cfg.network_mode.n_events must be a whole number of at least 1.');
end
if ~(isnumeric(P.lags) && isvector(P.lags) && all(diff(P.lags) > 0))
    bad_config('cfg.network_mode.lags must be increasing.');
end
if ~(isscalar(P.x_thr) && P.x_thr > 0 && P.x_thr < 1)
    bad_config('cfg.network_mode.x_thr must be a fraction between 0 and 1.');
end
if sim.n_stim < 1 || sim.n_no < 1
    bad_config('The bundle must have at least one no-stimulation and one stimulation block.');
end
runs = struct('seed', num2cell(P.seeds(:)'), 'file', arrayfun(@(s) sprintf('ied_netmode_seed%03d.mat', s), ...
    P.seeds(:)', 'UniformOutput', false));
plan = P;
plan.folder  = 'network_mode';
plan.blocks  = [1, sim.n_no + 1];
plan.sim     = sim;
plan.run_dir = cfg.run_dir;
plan.run_mode = cfg.run_mode;
plan.runs    = runs;
end

%% ------------------------------------------------------------------------
function v = knob(K, name, default)
v = default;
if isfield(K, name) && ~isempty(K.(name)); v = K.(name); end
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('ied_network_mode_stage_plan:BadConfig', varargin{:});
end
