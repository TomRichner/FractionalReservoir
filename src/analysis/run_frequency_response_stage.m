function detail = run_frequency_response_stage(cfg, run_dir, plan)
% RUN_FREQUENCY_RESPONSE_STAGE A bundle's frequency_response stage: the network frequency response at the outlier_spectra calibration states, into <run_dir>/frequency_response.
%
%   detail = RUN_FREQUENCY_RESPONSE_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (its verbosity; the plan already holds
%            what the config asked for)
%   run_dir  the bundle's run directory
%   plan     from frequency_response_stage_plan(cfg)
%   detail   the file written, and what was computed, on one line
%
% The run half of the opt-in stage 'frequency_response'
% (paper_opt_in_stages; FD-041). It reads the index of the stage
% outlier_spectra in THIS run directory (fresp_load_states: the one place
% that layout is read; no index is an error, there is no fallback), takes
% every run of the plan's networks, whatever conditions and seeds are
% there, and computes fresp_run_one for each: no simulation, linear algebra
% at the saved states. Then, per seed and per condition with slow
% variables, the single-neuron responses at the matched operating point on
% the first of the plan's networks (fresp_matched: it uses the network's
% parameters, not its states).
%
% OUTPUT  <run_dir>/frequency_response/frequency_response_data.mat
%   settings   struct: the plan's settings (f, f_single, sigma, w, r0,
%              ...), networks, conditions (in the index's order), seeds
%              (n_seeds x 2, in the index's order), input (the index
%              folder), seconds
%   runs       struct array, one per (network, condition, seed): the fields
%              of fresp_run_one
%   matched    struct array, one per (condition with slow variables, seed):
%              network, condition, seed and the fields of fresp_matched
%
% PARALLEL. The runs are independent and go through the session's pool when
% one is up (pool_ensure), serially otherwise. Per run at 500 neurons and
% 4000 state variables the largest objects are the dense full Jacobian and
% its eigenvalue workspace (about 0.4 GB while eig runs) and the gains on a
% chunk of the complex grid (about 0.3 GB): under 1 GB per worker.
%
% TIME, AN ESTIMATE TO BE MEASURED ON THE WORKSTATION. Not timed at 500
% neurons in this form. The exploration it was ported from
% (scripts/explorations/frequency_response, the same computation per run)
% took, serially on a MacBook Air M3 at 500 neurons, 12 states and 61
% frequencies: about 2 minutes a run with one timescale (95 and 127 s) and
% about 4 minutes with multiple timescales (248 and 262 s; 4000 state
% variables); a condition without slow variables takes seconds. So about
% 6 to 7 minutes per network and seed for the three conditions, and about
% 2 hours of single-process time for 2 networks x 3 conditions x 10 seeds,
% divided by the pool's workers (60 independent runs). It is dominated by
% the dense eigenvalue problems: one n x n per frequency and state, one
% N x N per check state.
%
% See also: frequency_response_stage_plan, fresp_run_one, fresp_load_states,
%           fig_freqresp_network, fig_freqresp_plane, fig_freqresp_phase,
%           fig_freqresp_single_neuron

arguments
    cfg     (1,1) struct
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
t_stage = tic;

T = fresp_load_states(run_dir);
absent = setdiff(plan.networks, T.network);
if ~isempty(absent)
    error('run_frequency_response_stage:NoSuchNetwork', ...
        'The outlier_spectra stage of %s saved no run of the network(s) %s (it has %s).', ...
        run_dir, strjoin(absent, ', '), strjoin(unique(T.network, 'stable')', ', '));
end
T = T(ismember(T.network, plan.networks), :);
rows = table2struct(T);
n_runs = numel(rows);
conditions = unique(T.condition, 'stable')';
seeds = unique([T.seed1, T.seed2], 'rows', 'stable');
vprintf(verbose, 'minimal', ...
    '  frequency_response: %d run(s): %d network(s) x %d condition(s) x %d seed(s); %d frequencies, plane %d x %d\n', ...
    n_runs, numel(plan.networks), numel(conditions), size(seeds, 1), numel(plan.opts.f), ...
    numel(plan.opts.sigma), numel(plan.opts.w));

n_workers = 0;
pool = pool_ensure('verbose', verbose);
if ~isempty(pool); n_workers = pool.NumWorkers; end

opts = plan.opts;
runs = cell(1, n_runs);
parfor (k = 1:n_runs, n_workers)
    st = fresp_load_states(run_dir, rows(k));
    runs{k} = fresp_run_one(st, opts);
end
runs = [runs{:}];
for k = 1:n_runs
    vprintf(verbose, 'verbose', '    %-10s %-14s seed [%d %d]: %.0f s\n', runs(k).network, ...
        runs(k).condition, runs(k).seed(1), runs(k).seed(2), runs(k).seconds.all);
end

% The matched operating point: the first network's parameters, every seed,
% every condition with slow variables.
sel = find(strcmp({runs.network}, plan.networks{1}) & [runs.has_slow]);
matched = cell(1, numel(sel));
parfor (q = 1:numel(sel), n_workers)
    st = fresp_load_states(run_dir, rows(sel(q)));         %#ok<PFBNS>  a small struct array
    mt = fresp_matched(st.params, opts.f_single, opts.r0);  %#ok<PFBNS>
    mt.network = st.network; mt.condition = st.condition; mt.seed = st.seed;
    matched{q} = mt;
end
matched = [matched{:}];

settings = plan.opts;
settings.networks = plan.networks;
settings.conditions = conditions;
settings.seeds = seeds;
settings.input = fullfile(run_dir, plan.input_folder);
settings.seconds = toc(t_stage);

out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
mat_file = fullfile(out_dir, plan.file);
save(mat_file, 'settings', 'runs', 'matched', '-v7.3');
detail = sprintf('%s (%d runs: %d networks x %d conditions x %d seeds)', mat_file, n_runs, ...
    numel(plan.networks), numel(conditions), size(seeds, 1));
end
