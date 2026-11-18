function R = numerics_L_fs_check(cfg)
% NUMERICS_L_FS_CHECK Does a finer SRA1 step close check L's gap? (2026-11-17)
%
%   R = numerics_L_fs_check('out_dir', d)
%
% Exploration, not a production stage. Repeats check L of
% run_numerics_verification (noise-free Benettin LLE, full network, seeds
% [k, k+1], T_lle and lya_T_interval of the run mode, sigma_u_noise = 0) with
% SRA1 at several fs, against the ode45 Benettin values SAVED by the stage
% run in out_dir. Alongside, lambda_1 by top-K (the method of the paper's
% figures; settings of run_local_lyapunov: lya_dt 0.05, K = cfg.K, no
% auto-K) on each fs's own SRA1 trajectory, and top-K on an ode45 trajectory
% stored at each fs (ode45 RelTol = AbsTol = 1e-10, MaxStep = 1/fs as in L).
%
% The 400-Hz SRA1 Benettin is recomputed and compared with the saved value,
% which checks this driver reproduces the stage. Jobs run in a parfor over
% (condition, trial, job); the pool is the session's (pool_ensure).
%
% build_probe below is a copy of run_numerics_verification's local function
% of that name (twister pinned, preset < condition < overrides).

arguments
    cfg.out_dir    (1,:) char = ''   % the stage run's numerics_verification folder; '' -> <root>/data/mu7chain_AC_fastmedium/numerics_verification
    cfg.conditions        = {'sfa3_stdc2', 'sfa1_std1'}
    cfg.fs_list    (1,:) double = [400, 500, 800]
    cfg.K          (1,1) double = 10
    cfg.topk_lya_dt (1,1) double = 0.05
    cfg.verbose           = 'minimal'
end

setup_paths();
if isempty(cfg.out_dir)
    cfg.out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'mu7chain_AC_fastmedium', 'numerics_verification');
end
src =fullfile(cfg.out_dir, 'numerics_verification_data.mat');
S = load(src, 'results', 'settings', 'cond_names');
P = struct('preset_name', S.settings.preset_name, 'T_lle', S.settings.T_lle, ...
    'ref_tol', S.settings.ref_tol, 'K', cfg.K, 'topk_lya_dt', cfg.topk_lya_dt, 'verbose', 'near-none');
n_trials = S.settings.n_trials_lle;
conds = cellstr(cfg.conditions);

% job list: {cond, trial, solver, fs, method}
jobs = {};
for c = 1:numel(conds)
    for k = 1:n_trials
        for fs = cfg.fs_list
            jobs(end+1, :) = {conds{c}, k, 'sra1', fs, 'benettin'}; %#ok<AGROW>
            jobs(end+1, :) = {conds{c}, k, 'sra1', fs, 'topk'};     %#ok<AGROW>
            jobs(end+1, :) = {conds{c}, k, 'ode45', fs, 'topk'};    %#ok<AGROW>
        end
    end
end
n_jobs = size(jobs, 1);
pool = pool_ensure('verbose', cfg.verbose);
n_pool = 1; if ~isempty(pool); n_pool = pool.NumWorkers; end
vprintf(cfg.verbose, 'minimal', '[L_fs_check] %d jobs on %d workers\n', n_jobs, n_pool);
out = cell(n_jobs, 1);
t_all = tic;
parfor j = 1:n_jobs
    out{j} = one_job(P, jobs{j, :});
end
vprintf(cfg.verbose, 'minimal', '[L_fs_check] done in %.1f min\n', toc(t_all) / 60);

% assemble per condition
R = struct();
R.settings = struct('source', src, 'preset_name', P.preset_name, 'T_lle', P.T_lle, ...
    'lya_T_interval', [P.T_lle / 2, P.T_lle], 'fs_list', cfg.fs_list, 'K', cfg.K, ...
    'topk_lya_dt', cfg.topk_lya_dt, 'benettin_lya_dt', 0.02, 'sigma_u_noise', 0, ...
    'date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
for c = 1:numel(conds)
    i = find(strcmp(S.cond_names, conds{c}), 1);
    L = S.results(i).lle;
    C = struct('name', conds{c}, 'fs', cfg.fs_list, ...
        'ode45_benettin', arrayfun(@(t) t.ode45.LLE, L), ...
        'sra1_400_benettin_saved', arrayfun(@(t) t.sra1.LLE, L), ...
        'sra1_400_seconds_saved', arrayfun(@(t) t.sra1.seconds, L), ...
        'benettin', nan(numel(cfg.fs_list), n_trials), 'topk', nan(numel(cfg.fs_list), n_trials), ...
        'sec_benettin', nan(numel(cfg.fs_list), n_trials), 'sec_topk', nan(numel(cfg.fs_list), n_trials), ...
        'ode45_topk', nan(numel(cfg.fs_list), n_trials), 'sec_ode45_topk', nan(numel(cfg.fs_list), n_trials));
    for j = 1:n_jobs
        if ~strcmp(jobs{j, 1}, conds{c}); continue; end
        k = jobs{j, 2}; f = find(cfg.fs_list == jobs{j, 4}, 1);
        if strcmp(jobs{j, 3}, 'ode45')
            C.ode45_topk(f, k) = out{j}.LLE; C.sec_ode45_topk(f, k) = out{j}.seconds;
        elseif strcmp(jobs{j, 5}, 'benettin')
            C.benettin(f, k) = out{j}.LLE; C.sec_benettin(f, k) = out{j}.seconds;
        else
            C.topk(f, k) = out{j}.LLE; C.sec_topk(f, k) = out{j}.seconds;
        end
    end
    C.d_benettin = abs(C.benettin - C.ode45_benettin);       % fs x trial
    C.d_topk     = abs(C.topk - C.ode45_benettin);
    C.d_topk_vs_ode45topk = abs(C.topk - C.ode45_topk);      % same fs, same top-K, integrator differs
    C.d_topk_benettin = abs(C.topk - C.benettin);
    C.scatter_ode45 = std(C.ode45_benettin);
    R.(conds{c}) = C;
    print_condition(C, cfg.verbose);
end
save(fullfile(cfg.out_dir, 'L_fs_check.mat'), 'R');
vprintf(cfg.verbose, 'minimal', '[L_fs_check] -> %s\n', fullfile(cfg.out_dir, 'L_fs_check.mat'));
end

function out = one_job(P, cname, k, solver, fs, method)
args = {'rng_seeds', [k, k + 1], 'verbose', P.verbose, 'sigma_u_noise', 0, ...
    'ode_solver', solver, 'fs', fs, 'T_range', [0, P.T_lle], 'lya_method', method, ...
    'lya_T_interval', [P.T_lle / 2, P.T_lle], 'plot_deci', 4};
if strcmp(method, 'topk')
    args = [args, {'lya_K', P.K, 'lya_K_auto', false, 'lya_dt', P.topk_lya_dt}];
end
model = build_probe(P.preset_name, cname, args{:});
if strcmp(solver, 'ode45')
    model.ode_opts = odeset('RelTol', P.ref_tol, 'AbsTol', P.ref_tol, 'MaxStep', 1 / fs);
end
t0 = tic;
model.run();
out = struct('LLE', model.lya_results.LLE, 'seconds', toc(t0));
end

function print_condition(C, verbose)
vprintf(verbose, 'minimal', '\n== %s ==  ode45 Benettin scatter (std) %.4f\n', C.name, C.scatter_ode45);
vprintf(verbose, 'minimal', '  ode45 Benettin: %s\n', mat2str(C.ode45_benettin, 4));
vprintf(verbose, 'minimal', '  repro 400 Hz  : max |recomputed - saved| = %.2e\n', ...
    max(abs(C.benettin(C.fs == 400, :) - C.sra1_400_benettin_saved)));
for f = 1:numel(C.fs)
    vprintf(verbose, 'minimal', ['  fs %4d  Benettin %s | med|d| %.4f | %.1f s\n' ...
        '           top-K    %s | med|d| %.4f | %.1f s | med|topK-Ben| %.4f\n' ...
        '           top-K vs ode45 top-K at this fs %s | med %.4f   (ode45 top-K %s)\n'], C.fs(f), ...
        mat2str(C.d_benettin(f, :), 3), median(C.d_benettin(f, :)), median(C.sec_benettin(f, :)), ...
        mat2str(C.d_topk(f, :), 3), median(C.d_topk(f, :)), median(C.sec_topk(f, :)), median(C.d_topk_benettin(f, :)), ...
        mat2str(C.d_topk_vs_ode45topk(f, :), 3), median(C.d_topk_vs_ode45topk(f, :)), mat2str(C.ode45_topk(f, :), 4));
end
end

function probe = build_probe(preset_name, condition_name, varargin)
rng(0, 'twister');
[preset, ~, conditions] = srnn_param_preset(preset_name);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond  = rmfield(conditions{strcmp(names, condition_name)}, 'name');
args  = [struct2namevalue(preset), struct2namevalue(cond), varargin];
probe = SRNNNumericsProbe(args{:});
probe.build();
end
