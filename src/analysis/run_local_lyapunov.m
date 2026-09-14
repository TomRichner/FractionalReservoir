function mat_file = run_local_lyapunov(cfg)
% RUN_LOCAL_LYAPUNOV Local Lyapunov exponents and local KS entropy under a stimulus staircase.
%
%   mat_file = RUN_LOCAL_LYAPUNOV('preset_name', p, 'run_mode', m, 'out_dir', d)
%
% One network (rng_seeds [1 2]) of the preset under each of its adaptation
% conditions, run TWICE on the same Brownian path: top-K QR (exactly K, no
% retry) and Benettin (K = 1, the manuscript's estimator). T is fixed at 60 s
% for every mode because the stimulus staircase presets state their steps for
% a 60-s run (12 steps of 5 s); the exponents accumulate over [30 60] after a
% 15-s alignment. K by mode: fast 30, medium 100, production 100.
%
% What it is for (TR, 2026-09-14): how often the local rates are positive, and
% whether the local KS entropy rate sum_k max(local_k, 0) spikes at every
% change of the input and returns to zero -- the dynamical signature of a
% stable adapting network encoding a new input. Grown out of
% scripts/explorations/explore_local_lyapunov_exponents.m.
%
% Output: <out_dir>/local_lyapunov_data.mat with results(i) per condition
% (name, title, topk = lya_results, ben = lya_results, u_ex, t_ex,
% type_indices, cell_type_names, seconds) and settings. A standalone run lands
% in <root>/data/local_lyapunov.
%
% See also: fig_local_lyapunov, lyapunov_topk, SRNNCellTypePairs.lya_summary

arguments
    cfg.preset_name (1,:) char = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_steps5s_noise0p025_dualStd_3cond_mu5'
    cfg.run_mode    (1,:) char = 'production'
    cfg.out_dir     (1,:) char = ''
    cfg.verbose                 = 'minimal'
    cfg.K           (1,1) double = 0     % 0 -> per run_mode
    cfg.T           (1,1) double = 60
    cfg.n_override  (1,1) double = 0     % tests only
end

setup_paths();
switch cfg.run_mode
    case 'fast',                K = 30;
    case {'medium', 'medium2'}, K = 100;
    case 'production',          K = 100;
    otherwise
        error('run_local_lyapunov:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.K > 0; K = cfg.K; end
T = cfg.T;
if isempty(cfg.out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'local_lyapunov');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir); mkdir(out_dir); end

[~, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
title_map  = srnn_condition_titles();
titles     = cellfun(@(n) title_map(n), cond_names, 'UniformOutput', false);
seeds  = [1 2];
common = {'rng_seeds', seeds, 'fs', 400, 'T_range', [0 T], 'lya_T_interval', [T/2 T], ...
          'lya_warmup', T/4, 'verbose', cfg.verbose};
if cfg.n_override > 0
    common = [common, {'n', cfg.n_override, 'indegree', max(2, round(0.2 * cfg.n_override)), ...
        'F_tracks_network', true}];
end
vprintf(cfg.verbose, 'minimal', '[local_lyapunov] preset=%s run_mode=%s: %d conditions, T = %g s, K = %d\n', ...
    cfg.preset_name, cfg.run_mode, numel(cond_names), T, K);

t_stage = tic;
results = struct('name', cond_names, 'title', titles, 'topk', [], 'ben', [], 'u_ex', [], 't_ex', [], ...
    'type_indices', [], 'cell_type_names', [], 'seconds', []);
for i = 1:numel(cond_names)
    t0 = tic;
    rng(0, 'twister');
    m = build_from_preset(cfg.preset_name, cond_names{i}, common{:}, ...
        'lya_method', 'topk', 'lya_K', K, 'lya_K_auto', false, 'lya_dt', 0.05);
    m.run();
    results(i).topk = m.lya_results;
    results(i).u_ex = m.u_ex; results(i).t_ex = m.t_ex;
    results(i).type_indices = m.type_indices; results(i).cell_type_names = m.cell_type_names;
    rng(0, 'twister');
    m = build_from_preset(cfg.preset_name, cond_names{i}, common{:}, 'lya_method', 'benettin');
    m.run();
    results(i).ben = m.lya_results;
    results(i).seconds = toc(t0);
    vprintf(cfg.verbose, 'minimal', '  %-14s top-%d lambda_1 %+.4f | Benettin %+.4f | h_KS %.2f bit/s | %.0f s\n', ...
        cond_names{i}, K, results(i).topk.LE_spectrum(1), results(i).ben.LLE, ...
        results(i).topk.h_KS_bits, results(i).seconds);
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'seeds', seeds, ...
    'T', T, 'K', K, 'lya_T_interval', [T/2 T], 'lya_warmup', T/4, 'fs', 400, ...
    'n_override', cfg.n_override, 'minutes', toc(t_stage) / 60);
condition_titles = titles; %#ok<NASGU>
mat_file = fullfile(out_dir, 'local_lyapunov_data.mat');
save(mat_file, 'results', 'cond_names', 'condition_titles', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '[local_lyapunov] %.1f min -> %s\n', settings.minutes, mat_file);
end
