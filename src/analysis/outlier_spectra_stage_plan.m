function plan = outlier_spectra_stage_plan(cfg)
% OUTLIER_SPECTRA_STAGE_PLAN The runs a bundle's outlier_spectra stage makes, validated before any compute.
%
%   plan = OUTLIER_SPECTRA_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'outlier_spectra'
% (paper_opt_in_stages; FD-040, docs/Plan_spectral_supplement.md): the
% spectra of the frozen-gain effective connectivity and of the full
% Jacobian along settled trajectories of the bundle's network and of
% copies with a stronger E-to-E mean, under each condition, on several
% seeds. It is the data behind the outlier figures (FD-036 to FD-039) and
% the input of the frequency_response and eigenvalue_density stages.
% Config in, plan out; nothing is computed and no file is touched.
% run_outlier_spectra_stage carries the plan out; ospec_run is one run;
% the header of ospec_load_run documents what a run saves.
%
% THE CONFIG
%   cfg.preset_name   the network: the bundle's NOISY SRNNCellTypePairs
%                     preset with two cell types. Required.
%   cfg.run_mode      sets the size (the table below)
%   cfg.run_conditions   the conditions the bundle simulates
%                     (bundle_run_conditions); the stage's default
%                     conditions are those of its three that are among them
%   cfg.verbose
%   cfg.outlier_spectra_<knob>   optional, one field per knob below. A
%                     field that starts with outlier_spectra_ and is not a
%                     knob is refused, so a typo cannot be ignored.
%
% THE SIZES, by run mode (the knobs a mode sets):
%   mode                      n_seeds  window_s  full_window_s  full_dt_sample
%   essential, fast              2       20 s       20 s           0.01 s
%   fastmedium, medium(2),
%     production                10       20 s      200 s           0.05 s
% and in every mode: 500 neurons (the preset's), 60 s of settling, 12
% calibration states, the full-Jacobian series and the dense spectra for
% reference and E-to-E x 1.5 only. Tom, 2026-11-20: more seeds, and a much
% longer window for the dwell-time figure. THE NUMBERS ABOVE 'fast' ARE THE
% WORKER'S PROPOSAL FROM ESTIMATED COST (FD-040's report), to be confirmed
% by timing on the workstation: nothing at this size was run on the
% MacBook.
%
% THE KNOBS (cfg.outlier_spectra_<name>; default)
%   networks        {'reference', 'EEx1p1', 'EEx1p25', 'EEx1p5', 'EEx1p75'}:
%                   ospec_network's names (the preset, and its E-to-E block
%                   mean times 1.1, 1.25, 1.5, 1.75: the dose)
%   full_networks   {'reference', 'EEx1p5'}: the networks that get the
%                   long window (full_window_s, full_dt_sample), the
%                   full-Jacobian series and the dense spectra at the
%                   calibration states. The others get window_s and
%                   dt_sample and neither pass.
%   conditions      those of {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
%                   the bundle simulates
%   n_seeds         the table's. Seeds are the first n_seeds rows of
%   seeds           [7 8; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 8 9; ...]: the
%                   exploration's seed first (the figures drawn so far are
%                   on it), then the pairs [k k+1] the other stages use
%   T_settle        60 s before the window
%   window_s        20 s: the measured window of the other networks
%   full_window_s   the table's: the measured window of full_networks
%   dt_sample       0.01 s: sampling step of the frozen-gain series
%   full_dt_sample  the table's: the same for full_networks. 0.05 s at
%                   200 s keeps the frozen-gain pass of a long run at the
%                   cost of two short ones; excursions shorter than 50 ms
%                   are then missed or merged in those runs (the
%                   exploration's frozen-gain excursions beyond zero had
%                   medians of 0.4 to 10 s with multiple timescales, and
%                   single 10 ms samples with one timescale at E-to-E x 1.5)
%   n_cal           12 calibration states
%   n_null, n_shuffle   40, 20: null draws and gain shuffles per state
%   n_top, n_keep   3, 30: rightmost modes kept with eigenvectors; largest
%                   real parts and moduli kept
%   dt_vectors      0.1 s: the eigenvectors of the rightmost modes are
%                   saved this often (0: every sample)
%   dt_mode         0.1 s: the intended mode is followed this often (0: not)
%   dt_full         0.1 s: step of the full-Jacobian series
%   k_full          8 rightmost eigenvalues of the full Jacobian per state
%   x_thr           0.1: the dendritic fraction of an eigenvector that
%                   makes it a network mode (ospec_full_spectrum)
%   lya_K           3 Lyapunov exponents (0: none), over the first
%   lya_window_s    20 s of the window
%   fs              500 Hz
%   chunk_s         10 s integrated per piece (memory only; the trajectory
%                   does not depend on it)
%   n_override      0: the preset's size. Tests use 40 to 80.
%
% EVERY RUN OF ONE NETWORK GROUP USES ONE LENGTH. Runs of different length
% on one seed are different noise realizations, so a full_networks run at
% 260 s is not the first 80 s of anything else.
%
% THE PLAN
%   plan.folder    'outlier_spectra', relative to the run directory
%   plan.runs      1 x m struct array, seeds outermost, then network, then
%                  condition: network, condition, seed (1 x 2), window_s,
%                  dt_sample, full (logical), file (the run's .mat name)
%   plan.<knob>    every knob above, as resolved; plan.preset_name,
%                  plan.run_mode, plan.seeds (the rows used), plan.verbose
%
% See also: run_outlier_spectra_stage, ospec_run, ospec_load_run,
%           ospec_network, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

%            mode          n_seeds  window_s  full_window_s  full_dt_sample
SIZES = {   'essential',     2,       20,        20,            0.01; ...
            'fast',          2,       20,        20,            0.01; ...
            'fastmedium',   10,       20,       200,            0.05; ...
            'medium',       10,       20,       200,            0.05; ...
            'medium2',      10,       20,       200,            0.05; ...
            'production',   10,       20,       200,            0.05};
KNOBS = {'networks', 'full_networks', 'conditions', 'n_seeds', 'seeds', 'T_settle', 'window_s', ...
    'full_window_s', 'dt_sample', 'full_dt_sample', 'n_cal', 'n_null', 'n_shuffle', 'n_top', 'n_keep', ...
    'dt_vectors', 'dt_mode', 'dt_full', 'k_full', 'x_thr', 'lya_K', 'lya_window_s', 'fs', 'chunk_s', 'n_override'};
PREFIX = 'outlier_spectra_';
STAGE_CONDITIONS = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};

assert(isequal(sort(SIZES(:, 1)), sort(run_mode_names())'), 'outlier_spectra_stage_plan:TableIncomplete', ...
    'The size table must name every run mode of run_mode_names.');
if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, SIZES(:, 1))
    bad_config('The outlier_spectra stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(SIZES(:, 1)', ', '));
end
preset_name = field_or_empty(cfg, 'preset_name');
if isempty(preset_name)
    bad_config('The outlier_spectra stage runs on cfg.preset_name, which this bundle does not set.');
end
names = fieldnames(cfg);
names = names(startsWith(names, PREFIX));
unknown = setdiff(names, strcat(PREFIX, KNOBS));
if ~isempty(unknown)
    bad_config('Unknown setting(s) %s. The knobs of the outlier_spectra stage are %s%s.', ...
        strjoin(unknown', ', '), PREFIX, strjoin(KNOBS, [', ' PREFIX]));
end
knob = @(name, default) override(cfg, [PREFIX name], default);

%% The network: a noisy two-type SRNNCellTypePairs preset
[preset, model_class, conditions] = srnn_param_preset(preset_name);   % errors on an unknown name
if ~strcmp(model_class, 'SRNNCellTypePairs')
    bad_config('The outlier_spectra stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        preset_name, model_class);
end
if ~(isfield(preset, 'sigma_u_noise') && any(preset.sigma_u_noise(:) > 0))
    bad_config(['The outlier_spectra stage runs the noisy network (cfg.preset_name), integrated in ' ...
        'pieces by a fixed-step scheme; preset ''%s'' has no noise.'], preset_name);
end
if ~(isfield(preset, 'mu_tilde_relative') && isequal(size(preset.mu_tilde_relative), [2 2]))
    bad_config(['The network names of the outlier_spectra stage scale the block means of a two-type ' ...
        'network; preset ''%s'' does not state a 2 x 2 mu_tilde_relative.'], preset_name);
end

%% The size: the mode's row, then the bundle's knobs
row = SIZES(strcmp(SIZES(:, 1), cfg.run_mode), :);
P = struct();
P.networks       = knob_cellstr(cfg, [PREFIX 'networks'], {'reference', 'EEx1p1', 'EEx1p25', 'EEx1p5', 'EEx1p75'});
P.full_networks  = knob_cellstr(cfg, [PREFIX 'full_networks'], intersect({'reference', 'EEx1p5'}, P.networks, 'stable'));
P.n_seeds        = knob('n_seeds', row{2});
P.T_settle       = knob('T_settle', 60);
P.window_s       = knob('window_s', row{3});
P.full_window_s  = knob('full_window_s', row{4});
P.dt_sample      = knob('dt_sample', 0.01);
P.full_dt_sample = knob('full_dt_sample', row{5});
P.n_cal          = knob('n_cal', 12);
P.n_null         = knob('n_null', 40);
P.n_shuffle      = knob('n_shuffle', 20);
P.n_top          = knob('n_top', 3);
P.n_keep         = knob('n_keep', 30);
P.dt_vectors     = knob('dt_vectors', 0.1);
P.dt_mode        = knob('dt_mode', 0.1);
P.dt_full        = knob('dt_full', 0.1);
P.k_full         = knob('k_full', 8);
P.x_thr          = knob('x_thr', 0.1);
P.lya_K          = knob('lya_K', 3);
P.lya_window_s   = knob('lya_window_s', 20);
P.fs             = knob('fs', 500);
P.chunk_s        = knob('chunk_s', 10);
P.n_override     = knob('n_override', 0);

whole_pos = {'n_seeds', 'n_cal', 'n_null', 'n_shuffle', 'n_top', 'n_keep', 'k_full'};
for k = 1:numel(whole_pos)
    v = P.(whole_pos{k});
    if ~(is_whole(v) && v >= 1)
        bad_config('cfg.%s%s must be a whole number of at least 1.', PREFIX, whole_pos{k});
    end
end
positive = {'T_settle', 'window_s', 'full_window_s', 'dt_sample', 'full_dt_sample', 'fs', 'chunk_s', 'lya_window_s'};
for k = 1:numel(positive)
    v = P.(positive{k});
    if ~(is_real_scalar(v) && v > 0)
        bad_config('cfg.%s%s must be a positive scalar.', PREFIX, positive{k});
    end
end
nonneg = {'dt_vectors', 'dt_mode'};
for k = 1:numel(nonneg)
    v = P.(nonneg{k});
    if ~(is_real_scalar(v) && v >= 0)
        bad_config('cfg.%s%s must be a scalar of at least 0 (0 switches it off).', PREFIX, nonneg{k});
    end
end
if ~(is_real_scalar(P.dt_full) && P.dt_full > 0)
    bad_config(['cfg.%sdt_full must be a positive scalar (seconds); to run without the full-Jacobian ' ...
        'series name no network in cfg.%sfull_networks.'], PREFIX, PREFIX);
end
if ~(is_real_scalar(P.x_thr) && P.x_thr > 0 && P.x_thr < 1)
    bad_config('cfg.%sx_thr must be a fraction between 0 and 1.', PREFIX);
end
if ~(is_whole(P.lya_K) && P.lya_K >= 0)
    bad_config('cfg.%slya_K must be a whole number of at least 0.', PREFIX);
end
if ~(is_whole(P.n_override) && (P.n_override == 0 || P.n_override >= 20))
    bad_config('cfg.%sn_override must be 0 (the preset''s size) or a whole number of at least 20.', PREFIX);
end
if min(P.dt_sample, P.full_dt_sample) < 1 / P.fs - 1e-12
    bad_config('The sampling steps (dt_sample %g s, full_dt_sample %g s) cannot be shorter than the integration step 1/fs = %g s.', ...
        P.dt_sample, P.full_dt_sample, 1 / P.fs);
end
if P.chunk_s < 3 / P.fs
    bad_config('cfg.%schunk_s must cover at least three integration steps.', PREFIX);
end

%% The networks
if isempty(P.networks) || numel(unique(P.networks)) ~= numel(P.networks)
    bad_config('cfg.%snetworks must name at least one network, each once.', PREFIX);
end
for k = 1:numel(P.networks)
    try
        ospec_network(P.networks{k}, preset.mu_tilde_relative);
    catch err
        bad_config('cfg.%snetworks: %s', PREFIX, err.message);
    end
end
missing = setdiff(P.full_networks, P.networks);
if ~isempty(missing)
    bad_config('cfg.%sfull_networks names %s, which cfg.%snetworks does not.', PREFIX, strjoin(missing, ', '), PREFIX);
end

%% The conditions
preset_conditions = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
simulated = bundle_run_conditions(cfg);
if isempty(simulated); simulated = preset_conditions; end
P.conditions = knob_cellstr(cfg, [PREFIX 'conditions'], intersect(STAGE_CONDITIONS, simulated, 'stable'));
if isempty(P.conditions)
    bad_config(['The outlier_spectra stage has no condition to run: the bundle simulates none of %s. ' ...
        'Name the conditions in cfg.%sconditions.'], strjoin(STAGE_CONDITIONS, ', '), PREFIX);
end
missing = setdiff(P.conditions, preset_conditions);
if ~isempty(missing)
    bad_config('cfg.%sconditions names %s, which preset ''%s'' does not state.', ...
        PREFIX, strjoin(missing, ', '), preset_name);
end
if numel(unique(P.conditions)) ~= numel(P.conditions)
    bad_config('cfg.%sconditions names a condition twice.', PREFIX);
end

%% The seeds
seeds = knob('seeds', default_seeds(max(P.n_seeds, 10)));
if ~(isnumeric(seeds) && ismatrix(seeds) && size(seeds, 2) == 2 && all(seeds(:) == round(seeds(:))) ...
        && all(seeds(:) >= 0) && size(unique(seeds, 'rows'), 1) == size(seeds, 1))
    bad_config('cfg.%sseeds must be an m x 2 matrix of distinct pairs of whole numbers (rng_seeds).', PREFIX);
end
if size(seeds, 1) < P.n_seeds
    bad_config('cfg.%sn_seeds is %d but cfg.%sseeds has %d row(s).', PREFIX, P.n_seeds, PREFIX, size(seeds, 1));
end
P.seeds = seeds(1:P.n_seeds, :);

%% The runs: seeds outermost, so the first seed is listed before the second
runs = struct('network', {}, 'condition', {}, 'seed', {}, 'window_s', {}, 'dt_sample', {}, 'full', {}, 'file', {});
for s = 1:P.n_seeds
    for w = 1:numel(P.networks)
        is_full = ismember(P.networks{w}, P.full_networks);
        for c = 1:numel(P.conditions)
            r = struct('network', P.networks{w}, 'condition', P.conditions{c}, 'seed', P.seeds(s, :), ...
                'window_s', P.window_s, 'dt_sample', P.dt_sample, 'full', is_full, ...
                'file', ospec_run_name(P.networks{w}, P.conditions{c}, P.seeds(s, :)));
            if is_full
                r.window_s = P.full_window_s;
                r.dt_sample = P.full_dt_sample;
            end
            runs(end + 1) = r; %#ok<AGROW>
        end
    end
end

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
plan = P;
plan.folder      = 'outlier_spectra';
plan.preset_name = preset_name;
plan.run_mode    = cfg.run_mode;
plan.verbose     = verbose;
plan.runs        = runs;
end

%% ------------------------------------------------------------------------
function seeds = default_seeds(m)
% [7 8] (the exploration's), then [k k+1] for k = 1, 2, ... without 7.
k = setdiff(1:(m + 1), 7, 'stable')';
seeds = [7 8; k, k + 1];
seeds = seeds(1:m, :);
end

function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function v = knob_cellstr(s, name, default)
% A list of names: a cellstr, or one name as a char. Unlike the other
% knobs an EMPTY CELL is taken (no network gets the full treatment).
v = default;
if isfield(s, name) && (iscell(s.(name)) || ~isempty(s.(name))); v = s.(name); end
if ischar(v) || (isstring(v) && isscalar(v)); v = {char(v)}; end
if ~iscellstr(v)
    bad_config('cfg.%s must be a cell of names.', name);
end
v = reshape(v, 1, []);
end

function tf = is_real_scalar(v)
tf = isscalar(v) && isnumeric(v) && isreal(v) && isfinite(v);
end

function tf = is_whole(v)
tf = is_real_scalar(v) && v == round(v);
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('outlier_spectra_stage_plan:BadConfig', varargin{:});
end
