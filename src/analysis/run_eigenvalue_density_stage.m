function detail = run_eigenvalue_density_stage(cfg, run_dir, plan)
% RUN_EIGENVALUE_DENSITY_STAGE A bundle's eigenvalue_density stage: the radial eigenvalue density of the frozen-gain connectivity, into <run_dir>/eigenvalue_density.
%
%   detail = RUN_EIGENVALUE_DENSITY_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (its verbosity; the plan holds the rest)
%   run_dir  the bundle's run directory, which must hold the output of the
%            stage outlier_spectra
%   plan     from eigenvalue_density_stage_plan(cfg)
%   detail   the file written and what it holds, on one line
%
% The run half of the opt-in stage 'eigenvalue_density' (FD-042). For every
% run the outlier_spectra index lists (all networks, conditions and seeds),
% eigdens_load_runs reads it and eigdens_run_measures measures it; the runs
% are independent and go through a parfor in the session's pool
% (pool_ensure; serial in the client without one). Nothing is simulated.
%
% OUTPUT <run_dir>/eigenvalue_density/eigenvalue_density_data.mat with
%   runs       struct array, one element per index row, in index order
%              (eigdens_run_measures' fields)
%   networks   struct array, one per network: name, s2, mu, R_struct
%              (realized block statistics of W, from its first run), harris
%              (eigdens_harris_density for W alone, with R_class the class's
%              radius), conditions and seeds present
%   settings   the grids, eigdens_measures('settings'), the input index's
%              path, n_runs, seconds, created
% fig_eigdens_radial, _gains, _measures and _plane read it.
%
% COST. The exploration measured the same quantities on twelve 500-neuron
% runs (12 states each) in about three minutes, most of it spent RECOVERING
% the gains by least squares from stored eigenvectors, which this stage
% does not do (it reads them). What remains per run: loading the run file,
% redrawing W (no integration), and per state the block theory (a
% fixed-point iteration over 101 radii) and the radial curves and measures
% of 1 + q + 1 spectra; see the timing in FD-042's report. Memory: one
% run's spectra at a time per worker (the run file's other fields are
% loaded with it; tens of MB at 500 neurons).
%
% See also: eigenvalue_density_stage_plan, eigdens_load_runs,
%           eigdens_run_measures, fig_eigdens_radial

arguments
    cfg     (1,1) struct
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
t_all = tic;
index = eigdens_load_runs(run_dir);
n_runs = height(index);
vprintf(verbose, 'minimal', '  eigenvalue_density: %d run(s) of %s, %d network(s)\n', ...
    n_runs, plan.input_folder, numel(unique(index.network)));

pool = pool_ensure();
n_workers = 0;
if ~isempty(pool); n_workers = pool.NumWorkers; end
grids = plan.grids;
res = cell(1, n_runs);
net = cell(1, n_runs);
parfor (k = 1:n_runs, n_workers)
    X = eigdens_load_runs(run_dir, index, k);
    res{k} = eigdens_run_measures(X, grids);
    net{k} = struct('name', X.network, 's2', X.s2, 'mu', X.mu, 'R_struct', X.R_struct, 'harris_in', X.harris);
end
runs = [res{:}];
net = [net{:}];

% One entry per network, from its first run (W is drawn from the network
% seed, so the seeds of one network have different W; s2 differs between
% them by sampling only, and the Harris density uses the ensemble's
% statistics, which are the same).
names = unique(index.network, 'stable');
networks = struct('name', names(:)', 's2', [], 'mu', [], 'R_struct', [], 'harris', [], ...
    'conditions', [], 'seeds', []);
for w = 1:numel(names)
    rows = find(strcmp(index.network, names{w}));
    h = net(rows(1)).harris_in;
    H = eigdens_harris_density(h.alpha, h.f_e, h.mu_e, h.mu_i, h.sigma_e, h.sigma_i);
    H.R_class = h.R_class;
    networks(w).s2 = net(rows(1)).s2;
    networks(w).mu = net(rows(1)).mu;
    networks(w).R_struct = net(rows(1)).R_struct;
    networks(w).harris = H;
    networks(w).conditions = unique(index.condition(rows), 'stable')';
    networks(w).seeds = unique([index.seed1(rows), index.seed2(rows)], 'rows', 'stable');
end

settings = struct('grids', grids, 'measures', plan.measures, ...
    'input_index', fullfile(run_dir, plan.input_folder), 'n_runs', n_runs, ...
    'seconds', toc(t_all), 'seconds_per_run', [runs.seconds], 'created', char(datetime('now')));
out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
file = fullfile(out_dir, plan.data_file);
save(file, 'runs', 'networks', 'settings');
detail = sprintf('%s (%d runs, %d networks, %.0f s)', file, n_runs, numel(networks), settings.seconds);
vprintf(verbose, 'minimal', '  eigenvalue_density: %s\n', detail);
end
