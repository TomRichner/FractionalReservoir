function detail = run_uncoupled_reference_stage(cfg, run_dir, plan)
% RUN_UNCOUPLED_REFERENCE_STAGE A bundle's uncoupled_reference stage: lambda_1 of the nearly uncoupled network, per condition, into <run_dir>/uncoupled_reference.
%
%   detail = RUN_UNCOUPLED_REFERENCE_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (unused: the plan holds what it asked for)
%   run_dir  the bundle's run directory
%   plan     from uncoupled_reference_stage_plan(cfg)
%   detail   the file written, and what was run, on one line
%
% The run half of the opt-in stage 'uncoupled_reference'
% (paper_opt_in_stages; Plan 2, package C). One ParamSpaceAnalysis2 with the
% joint sample's model settings, the axis level_of_chaos at plan.gains and
% plan.n_reps repetitions (network seeds), every condition; the sweep driver
% is used as it is. It writes, under <run_dir>/uncoupled_reference/:
%
%   uncoupled_reference_nLevs_<G>_<date>/   the driver's own folder
%                                  (psa_object.mat and per-condition results)
%   uncoupled_reference_data.mat   what the figure reads:
%     ref       struct array, one element per condition:
%                 name
%                 gains           1 x G, as plan.gains
%                 lambda_1        R x G, the leading exponent per repetition
%                                 (NaN where a run failed)
%                 lead_frac_x     R x G, the dendritic share of the leading
%                                 perturbation's squared norm
%                 mean_rate       R x G
%                 lambda_1_median, lambda_1_min, lambda_1_max   1 x G over
%                                 the repetitions
%                 n_ok            1 x G, runs that succeeded
%                 reference       lambda_1_median at the FIRST gain: the
%                                 value fig_leading_direction draws
%     settings  preset_name, run_mode, gains, n_reps, conditions, n,
%               T_range, lya_T_interval, fs, ode_solver, psa_dir, minutes
%
% A condition in which no run at the first gain succeeded is an error: a
% reference of NaN would be a missing line that nobody notices.
%
% The repetitions run in the session's pool when it has one (the driver's
% parfor), serially otherwise.
%
% See also: uncoupled_reference_stage_plan, fig_leading_direction,
%           ParamSpaceAnalysis2, run_param_space_analysis

arguments
    cfg     (1,1) struct %#ok<INUSA>
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = plan.verbose;
ctx = plan.ctx;
out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
data_file = fullfile(out_dir, 'uncoupled_reference_data.mat');
if isfile(data_file); delete(data_file); end   % never leave a file of an earlier run beside a failed one

n_gains = numel(plan.gains);
vprintf(verbose, 'minimal', '  uncoupled_reference: %d condition(s) x %d gain(s) (%s) x %d repetition(s)\n', ...
    numel(plan.conditions), n_gains, mat2str(plan.gains), plan.n_reps);
t_stage = tic;

psa = ParamSpaceAnalysis2('n_levels', n_gains, 'batch_size', 25, 'note', '', ...
    'randomize_order', false, 'verbose', verbose);
psa.folder_prefix  = 'uncoupled_reference';
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = ctx.model_defaults;
psa.output_dir     = out_dir;
psa.set_conditions(ctx.conditions);
psa.add_grid_levels('level_of_chaos', plan.gains);
psa.add_grid_levels('reps', 1:plan.n_reps);
psa.run();

%% The summary the figure reads
ref = struct('name', plan.conditions, 'gains', plan.gains, 'lambda_1', [], 'lead_frac_x', [], 'mean_rate', [], ...
    'lambda_1_median', [], 'lambda_1_min', [], 'lambda_1_max', [], 'n_ok', [], 'reference', []);
for c = 1:numel(plan.conditions)
    name = plan.conditions{c};
    R = psa.results.(name);
    lam = nan(plan.n_reps, n_gains); fx = lam; rate = lam;
    for k = 1:numel(R)
        r = R{k};
        if ~(isstruct(r) && isfield(r, 'success') && r.success && isfield(r, 'LLE') && isfinite(r.LLE)); continue; end
        g = find(plan.gains == r.config.level_of_chaos);
        rep = r.config.reps;
        assert(isscalar(g) && rep >= 1 && rep <= plan.n_reps, 'run_uncoupled_reference_stage:Job', ...
            'A result of %s is not on the stage''s grid (level_of_chaos %g, repetition %g).', ...
            name, r.config.level_of_chaos, rep);
        lam(rep, g) = r.LLE; fx(rep, g) = r.lead_frac_x; rate(rep, g) = r.mean_rate;
    end
    ref(c).lambda_1 = lam; ref(c).lead_frac_x = fx; ref(c).mean_rate = rate;
    ref(c).lambda_1_median = median(lam, 1, 'omitnan');
    ref(c).lambda_1_min = min(lam, [], 1); ref(c).lambda_1_max = max(lam, [], 1);
    ref(c).n_ok = sum(isfinite(lam), 1);
    ref(c).reference = ref(c).lambda_1_median(1);
    if ref(c).n_ok(1) == 0
        error('run_uncoupled_reference_stage:NoRun', ...
            'No run of %s succeeded at level_of_chaos = %g: there is no uncoupled reference for it.', ...
            name, plan.gains(1));
    end
    vprintf(verbose, 'minimal', '    %-14s lambda_1 = %s 1/s at gains %s (%s of %d runs)\n', name, ...
        mat2str(ref(c).lambda_1_median, 4), mat2str(plan.gains), mat2str(ref(c).n_ok), plan.n_reps);
end
md = psa.model_defaults;
settings = struct('preset_name', plan.preset_name, 'run_mode', ctx.run_mode, 'gains', plan.gains, ...
    'n_reps', plan.n_reps, 'conditions', {plan.conditions}, 'n', field_or(md, 'n', NaN), ...
    'T_range', field_or(md, 'T_range', []), 'lya_T_interval', field_or(md, 'lya_T_interval', []), ...
    'fs', field_or(md, 'fs', NaN), 'ode_solver', field_or(md, 'ode_solver', ''), ...
    'psa_dir', psa.output_dir, 'minutes', toc(t_stage) / 60);
save(data_file, 'ref', 'settings');
detail = sprintf('%s (%d conditions x %d gains x %d repetitions)', data_file, ...
    numel(plan.conditions), n_gains, plan.n_reps);
end

function v = field_or(s, name, default)
v = default;
if isfield(s, name); v = s.(name); end
end
