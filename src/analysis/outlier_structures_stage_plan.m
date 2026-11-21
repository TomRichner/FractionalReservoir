function plan = outlier_structures_stage_plan(cfg)
% OUTLIER_STRUCTURES_STAGE_PLAN The runs a bundle's outlier_structures stage makes, validated before any compute.
%
%   plan = OUTLIER_STRUCTURES_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'outlier_structures' (FD-043;
% docs/Plan_spectral_supplement.md). Config in, plan out: nothing is
% computed and no file is touched. run_outlier_structures_stage checks
% each network's structural spectrum and carries the plan out; the runs
% are the outlier_spectra stage's (ospec_run, with the intended mode
% followed), so a run's file is documented by the header of
% ospec_load_run, and ostruct_load reads the stage back.
%
% THE QUESTION (Tom, 2026-11-20). Is adaptation's reining in of outlier
% eigenvalues limited to outliers made by a mismatch of the E and I block
% means? Fifteen networks put an intended eigenvalue outside the bulk in
% other ways (scripts/explorations/outliers/explore_outlier_structures.m
% has the reasoning and the hypotheses H1 to H3), and each run follows that
% mode by its eigenvector, whether or not it stays the rightmost. The
% exploration found that with adaptation the intended mode ends at or near
% the bulk edge in every network.
%
% THE NETWORKS (cfg.outlier_structures_networks), an m x 3 cell, one row
% per network: {name, weight_seed, expect}
%   name         ospec_network's name ('lrA1p4W0p35', 'EEx1p5', ...;
%                no noise token: the stage adds those)
%   weight_seed  [] : a SEED of the stage is a new network: rng_seeds is
%                the seed row, so new weights, new added patterns (drawn
%                from rng_seeds(1)) and new noise, all at once.
%                [a b] : the weights are PART OF THE NETWORK'S DEFINITION
%                (a seed chosen because its draw has a local outlier):
%                every run uses rng_seeds [a b], so weights, setpoints,
%                time constants and initial state are fixed, and a seed of
%                the stage is a NOISE REALIZATION on those weights. The
%                first realization is the network's own noise (the class's
%                default noise seed, a + 224737: the exploration's run);
%                realization s >= 2 draws the noise from noise_seed
%                seeds(s, 1), through the network name's token
%                '-noise<seeds(s,1)>' (ospec_network), so it has a file and
%                a row of its own.
%   expect       where the intended eigenvalue should be, in units of the
%                structural bulk radius, for the stage's structural check:
%                []   from the name: a + i w for 'lrA<a>W<w>' (the
%                     eigenvalue the addition is built to have); for a
%                     block-mean network the rightmost eigenvalue of the
%                     theoretical block-mean matrix (the class's lambda_O
%                     over its R)
%                z    complex: this value (the local outliers, whose place
%                     the name cannot say; the exploration's measurement)
%                NaN  no outlier intended (the reference network for
%                     comparison: its intended mode is its rightmost
%                     eigenvalue, a bulk-edge mode); recorded, never skipped
% The DEFAULT is the exploration's fifteen, in three groups (the group
% follows from ospec_network's target):
%   rotation   added low-rank structure (target 'lowrank'): lrA1p5W0p0,
%              lrA1p4W0p15, lrA1p4W0p35, lrA1p4W0p6, lrA1p4W1p0 (about
%              1.4 R, rotation from 0 to 1 R), lrA2p5W0p0, lrA2p5W1p0
%              (about 2.5 R)
%   local      no structure put in (target 'rightmost'): reference on
%              weight seed [27 28] (a real local outlier at 1.33 R),
%              ind200 on weight seed [25 26] ((1.16 + 0.64i) R), and the
%              reference network with new weights per seed (comparison)
%   blockmean  block means (target 'blockmean'): EEx1p5, EIx0p9, IIx1p25,
%              the E-I loops EEx1p6-EIx1p33-IEx1p33 and EEx1p6-EIx1p31-IEx1p31
% A network's KEY (the name the figures and tables use) is its name, or
% '<name>_w<a>' for a network with a fixed weight seed ('reference_w27').
%
% THE STRUCTURAL CHECK (run_outlier_structures_stage, before any run).
% For each network and each distinct weight draw, the structural W is
% built (no integration), its intended eigenvalue lam0 found as the runs
% will find it (ospec_mode_reference), and the check passes when
%   |lam0 / R - expect| <= check_tol |expect|   and   Re lam0 / R >= check_min_re
% with R the structural bulk radius of the drawn matrix. A draw that fails
% is recorded with the reason and its runs are skipped.
%
% THE CONFIG
%   cfg.preset_name   the bundle's NOISY two-type SRNNCellTypePairs preset
%   cfg.run_mode      sets the number of seeds (the table below)
%   cfg.run_conditions   the conditions the bundle simulates
%   cfg.verbose
%   cfg.outlier_structures_<knob>   one field per knob; an unknown one is refused
%
% THE SIZES, by run mode
%   mode                                       n_seeds
%   essential, fast                               1
%   fastmedium, medium, medium2, production       5
% and in every mode the exploration's sampling: 500 neurons, 60 s of
% settling, a 20 s window sampled every 50 ms, the intended mode followed
% every 0.1 s, 8 calibration states, the dense full-Jacobian pass at the
% calibration states, lambda_1 over the window. (Tom, 2026-11-20: "much of
% the remaining figs below in the supplement too", with more seeds.)
%
% THE KNOBS (cfg.outlier_structures_<name>; default)
%   networks      the m x 3 cell above
%   conditions    those of {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
%                 the bundle simulates
%   n_seeds       the table's
%   seeds         [7 8; 1 2; 2 3; 3 4; 4 5; ...] (outlier_spectra's
%                 rows: the exploration's seed first). The first column of
%                 the rows used must be distinct when a network has a fixed
%                 weight seed (it names the noise realizations).
%   T_settle 60, window_s 20, dt_sample 0.05, n_cal 8, n_null 40,
%   n_shuffle 20, n_top 3, n_keep 30, dt_vectors 0.1, dt_mode 0.1 (must be
%   above 0: the stage is about the followed mode), k_full 8, x_thr 0.1,
%   lya_K 3, lya_window_s 20, chunk_s 10, n_override 0 (tests: 40 to 80)
%   fs            the bundle's cfg.fs (stage_fs; 500 Hz in the mu7chain
%                 bundles), else 500 Hz; cfg.outlier_structures_fs, when
%                 set, overrides it for this stage. The integrator is not a
%                 knob: the bundle's cfg.ode_solver (plan.ode_solver), else
%                 ospec_build_model's rule, exactly as in
%                 outlier_spectra_stage_plan (the runs are that stage's).
%   dt_full       window_s: the full-Jacobian series of ospec_run at the
%                 start and the end of the window only (the exploration's
%                 dt_full of 20 s). The stage wants the dense pass, which
%                 ospec_run's stage machinery gives together with the series.
%   check_tol     0.15, check_min_re 1.05 (the structural check)
%
% THE PLAN
%   plan.folder     'outlier_structures'
%   plan.networks   1 x m struct: key, name, group, target, weight_seed
%                   ([] or 1 x 2), expect ([] or a scalar)
%   plan.runs       1 x k struct, seeds outermost, then network, then
%                   condition: key, group, base (the network's name),
%                   network (the name ospec_run is given, with any noise
%                   token), condition, seed (rng_seeds), realization (the
%                   row of plan.seeds), noise_seed (NaN: the class's
%                   default), window_s, dt_sample, full (true), file
%   plan.checks     1 x j struct, the weight draws to check: key, seed
%   plan.<knob>     every knob, resolved; plan.preset_name, plan.run_mode,
%                   plan.verbose; and what run_outlier_spectra_stage reads
%                   from a plan (full_networks, full_window_s, the list of
%                   network names), so that the runs go through it
%
% See also: run_outlier_structures_stage, ostruct_load, ospec_run,
%           ospec_network, outlier_spectra_stage_plan, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

%            mode          n_seeds
SIZES = {   'essential',     1; ...
            'fast',          1; ...
            'fastmedium',    5; ...
            'medium',        5; ...
            'medium2',       5; ...
            'production',    5};
KNOBS = {'networks', 'conditions', 'n_seeds', 'seeds', 'T_settle', 'window_s', 'dt_sample', 'n_cal', ...
    'n_null', 'n_shuffle', 'n_top', 'n_keep', 'dt_vectors', 'dt_mode', 'dt_full', 'k_full', 'x_thr', ...
    'lya_K', 'lya_window_s', 'fs', 'chunk_s', 'n_override', 'check_tol', 'check_min_re'};
PREFIX = 'outlier_structures_';
STAGE_CONDITIONS = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
GROUP_OF = struct('lowrank', 'rotation', 'rightmost', 'local', 'blockmean', 'blockmean');

assert(isequal(sort(SIZES(:, 1)), sort(run_mode_names())'), 'outlier_structures_stage_plan:TableIncomplete', ...
    'The size table must name every run mode of run_mode_names.');
if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, SIZES(:, 1))
    bad_config('The outlier_structures stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(SIZES(:, 1)', ', '));
end
preset_name = field_or_empty(cfg, 'preset_name');
if isempty(preset_name)
    bad_config('The outlier_structures stage runs on cfg.preset_name, which this bundle does not set.');
end
names = fieldnames(cfg);
names = names(startsWith(names, PREFIX));
unknown = setdiff(names, strcat(PREFIX, KNOBS));
if ~isempty(unknown)
    bad_config('Unknown setting(s) %s. The knobs of the outlier_structures stage are %s%s.', ...
        strjoin(unknown', ', '), PREFIX, strjoin(KNOBS, [', ' PREFIX]));
end
knob = @(name, default) override(cfg, [PREFIX name], default);

%% The network: a noisy two-type SRNNCellTypePairs preset
[preset, model_class, conditions] = srnn_param_preset(preset_name);   % errors on an unknown name
if ~strcmp(model_class, 'SRNNCellTypePairs')
    bad_config('The outlier_structures stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        preset_name, model_class);
end
if ~(isfield(preset, 'sigma_u_noise') && any(preset.sigma_u_noise(:) > 0))
    bad_config('The outlier_structures stage runs the noisy network (cfg.preset_name); preset ''%s'' has no noise.', preset_name);
end
if ~(isfield(preset, 'mu_tilde_relative') && isequal(size(preset.mu_tilde_relative), [2 2]))
    bad_config(['The network names of the outlier_structures stage scale the block means of a two-type ' ...
        'network; preset ''%s'' does not state a 2 x 2 mu_tilde_relative.'], preset_name);
end

%% The knobs
row = SIZES(strcmp(SIZES(:, 1), cfg.run_mode), :);
P = struct();
P.n_seeds      = knob('n_seeds', row{2});
P.T_settle     = knob('T_settle', 60);
P.window_s     = knob('window_s', 20);
P.dt_sample    = knob('dt_sample', 0.05);
P.n_cal        = knob('n_cal', 8);
P.n_null       = knob('n_null', 40);
P.n_shuffle    = knob('n_shuffle', 20);
P.n_top        = knob('n_top', 3);
P.n_keep       = knob('n_keep', 30);
P.dt_vectors   = knob('dt_vectors', 0.1);
P.dt_mode      = knob('dt_mode', 0.1);
P.dt_full      = knob('dt_full', P.window_s);
P.k_full       = knob('k_full', 8);
P.x_thr        = knob('x_thr', 0.1);
P.lya_K        = knob('lya_K', 3);
P.lya_window_s = knob('lya_window_s', 20);
P.fs           = knob('fs', stage_fs(field_or_empty(cfg, 'fs'), 500, 'outlier_structures_stage_plan'));   % the bundle's rate, else 500 Hz
P.chunk_s      = knob('chunk_s', 10);
P.n_override   = knob('n_override', 0);
P.check_tol    = knob('check_tol', 0.15);
P.check_min_re = knob('check_min_re', 1.05);

whole_pos = {'n_seeds', 'n_cal', 'n_null', 'n_shuffle', 'n_top', 'n_keep', 'k_full'};
for k = 1:numel(whole_pos)
    v = P.(whole_pos{k});
    if ~(is_whole(v) && v >= 1)
        bad_config('cfg.%s%s must be a whole number of at least 1.', PREFIX, whole_pos{k});
    end
end
positive = {'T_settle', 'window_s', 'dt_sample', 'fs', 'chunk_s', 'lya_window_s', 'dt_mode', 'dt_full', 'check_tol'};
for k = 1:numel(positive)
    v = P.(positive{k});
    if ~(is_real_scalar(v) && v > 0)
        bad_config('cfg.%s%s must be a positive scalar (the intended mode is what this stage follows: dt_mode cannot be 0).', ...
            PREFIX, positive{k});
    end
end
if ~(is_real_scalar(P.dt_vectors) && P.dt_vectors >= 0)
    bad_config('cfg.%sdt_vectors must be a scalar of at least 0 (0: every sample).', PREFIX);
end
if ~is_real_scalar(P.check_min_re)
    bad_config('cfg.%scheck_min_re must be a real scalar (units of the structural bulk radius).', PREFIX);
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
if P.dt_sample < 1 / P.fs - 1e-12
    bad_config('cfg.%sdt_sample (%g s) cannot be shorter than the integration step 1/fs = %g s.', PREFIX, P.dt_sample, 1 / P.fs);
end
if P.chunk_s < 3 / P.fs
    bad_config('cfg.%schunk_s must cover at least three integration steps.', PREFIX);
end

%% The networks
spec = knob('networks', default_networks());
if ~(iscell(spec) && ismatrix(spec) && size(spec, 2) == 3 && size(spec, 1) >= 1)
    bad_config('cfg.%snetworks must be an m x 3 cell, one row {name, weight_seed, expect} per network.', PREFIX);
end
m = size(spec, 1);
networks = struct('key', {}, 'name', {}, 'group', {}, 'target', {}, 'weight_seed', {}, 'expect', {}, 'label', {});
for k = 1:m
    [name, ws_, ex] = spec{k, :};
    if ~(ischar(name) || (isstring(name) && isscalar(name))) || isempty(name)
        bad_config('cfg.%snetworks, row %d: the name must be text.', PREFIX, k);
    end
    name = char(name);
    if any(~cellfun(@isempty, regexp(strsplit(name, '-'), '^noise\d+$', 'once')))
        bad_config('cfg.%snetworks, row %d (%s): no noise token in a name; the stage adds them for the noise realizations.', ...
            PREFIX, k, name);
    end
    try
        net = ospec_network(name, preset.mu_tilde_relative);
    catch err
        bad_config('cfg.%snetworks, row %d: %s', PREFIX, k, err.message);
    end
    if ~(isempty(ws_) || (isnumeric(ws_) && isequal(size(ws_), [1 2]) && all(ws_ == round(ws_)) && all(ws_ >= 0)))
        bad_config('cfg.%snetworks, row %d (%s): the weight seed is [] or a 1 x 2 pair of whole numbers (rng_seeds).', ...
            PREFIX, k, name);
    end
    if ~(isempty(ex) || (isnumeric(ex) && isscalar(ex)))
        bad_config('cfg.%snetworks, row %d (%s): expect is [], a complex number in bulk radii, or NaN.', PREFIX, k, name);
    end
    if isempty(ex) && strcmp(net.target, 'rightmost')
        bad_config(['cfg.%snetworks, row %d (%s): the name puts no eigenvalue anywhere (its intended mode is the ' ...
            'rightmost of the draw), so expect must be given: a complex number, or NaN for no outlier.'], PREFIX, k, name);
    end
    key = name;
    if ~isempty(ws_); key = sprintf('%s_w%d', name, ws_(1)); end
    networks(end + 1) = struct('key', key, 'name', name, 'group', GROUP_OF.(net.target), 'target', net.target, ...
        'weight_seed', ws_, 'expect', ex, 'label', net.label); %#ok<AGROW>
end
if numel(unique({networks.key})) ~= m
    bad_config('cfg.%snetworks names a network twice (the same name and weight seed).', PREFIX);
end

%% The conditions
preset_conditions = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
simulated = bundle_run_conditions(cfg);
if isempty(simulated); simulated = preset_conditions; end
P.conditions = knob_cellstr(cfg, [PREFIX 'conditions'], intersect(STAGE_CONDITIONS, simulated, 'stable'));
if isempty(P.conditions)
    bad_config(['The outlier_structures stage has no condition to run: the bundle simulates none of %s. ' ...
        'Name the conditions in cfg.%sconditions.'], strjoin(STAGE_CONDITIONS, ', '), PREFIX);
end
missing = setdiff(P.conditions, preset_conditions);
if ~isempty(missing)
    bad_config('cfg.%sconditions names %s, which preset ''%s'' does not state.', PREFIX, strjoin(missing, ', '), preset_name);
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
has_fixed = any(arrayfun(@(w) ~isempty(w.weight_seed), networks));
if has_fixed && numel(unique(P.seeds(:, 1))) ~= P.n_seeds
    bad_config(['cfg.%sseeds: the first column of the rows used must be distinct, because it names the noise ' ...
        'realizations of the networks with a fixed weight seed.'], PREFIX);
end

%% The runs: seeds outermost; and the weight draws to check
runs = struct('key', {}, 'group', {}, 'base', {}, 'network', {}, 'condition', {}, 'seed', {}, ...
    'realization', {}, 'noise_seed', {}, 'window_s', {}, 'dt_sample', {}, 'full', {}, 'file', {});
checks = struct('key', {}, 'seed', {});
for s = 1:P.n_seeds
    for w = 1:m
        nw = networks(w);
        if isempty(nw.weight_seed)
            seed = P.seeds(s, :); run_name = nw.name; noise_seed = NaN;
            checks(end + 1) = struct('key', nw.key, 'seed', seed); %#ok<AGROW>
        else
            seed = nw.weight_seed;
            if s == 1
                run_name = nw.name; noise_seed = NaN;
                checks(end + 1) = struct('key', nw.key, 'seed', seed); %#ok<AGROW>
            else
                noise_seed = P.seeds(s, 1);
                run_name = sprintf('%s-noise%d', nw.name, noise_seed);
            end
        end
        for c = 1:numel(P.conditions)
            runs(end + 1) = struct('key', nw.key, 'group', nw.group, 'base', nw.name, 'network', run_name, ...
                'condition', P.conditions{c}, 'seed', seed, 'realization', s, 'noise_seed', noise_seed, ...
                'window_s', P.window_s, 'dt_sample', P.dt_sample, 'full', true, ...
                'file', ospec_run_name(run_name, P.conditions{c}, seed)); %#ok<AGROW>
        end
    end
end
assert(numel(unique({runs.file})) == numel(runs), 'outlier_structures_stage_plan:Files', 'Two runs would share a file.');

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
[ode_solver, problem] = ospec_bundle_solver(cfg, 'outlier_structures_stage_plan');
if ~isempty(problem); bad_config('%s', problem); end
plan = P;
plan.ode_solver    = ode_solver;
plan.folder        = 'outlier_structures';
plan.preset_name   = preset_name;
plan.run_mode      = cfg.run_mode;
plan.verbose       = verbose;
plan.networks      = networks;
plan.runs          = runs;
plan.checks        = checks;
% What run_outlier_spectra_stage reads from a plan, besides the run knobs:
% every run gets the same window and the dense pass.
plan.full_networks = unique({runs.network}, 'stable');
plan.full_window_s = P.window_s;
end

%% ------------------------------------------------------------------------
function spec = default_networks()
% The exploration's fifteen (explore_outlier_structures).
spec = { ...
    'lrA1p5W0p0',  [], []; ...
    'lrA1p4W0p15', [], []; ...
    'lrA1p4W0p35', [], []; ...
    'lrA1p4W0p6',  [], []; ...
    'lrA1p4W1p0',  [], []; ...
    'lrA2p5W0p0',  [], []; ...
    'lrA2p5W1p0',  [], []; ...
    'reference',   [27 28], 1.33; ...
    'ind200',      [25 26], complex(1.16, 0.64); ...
    'reference',   [], NaN; ...
    'EEx1p5',      [], []; ...
    'EIx0p9',      [], []; ...
    'IIx1p25',     [], []; ...
    'EEx1p6-EIx1p33-IEx1p33', [], []; ...
    'EEx1p6-EIx1p31-IEx1p31', [], []};
end

function seeds = default_seeds(m)
% outlier_spectra's: [7 8] (the exploration's), then [k k+1] for k = 1, 2, ... without 7.
k = setdiff(1:(m + 1), 7, 'stable')';
seeds = [7 8; k, k + 1];
seeds = seeds(1:m, :);
end

function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function v = knob_cellstr(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
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
error('outlier_structures_stage_plan:BadConfig', varargin{:});
end
