% S2B_06_PAPER_RUNS  Does the local rate that the figures of the paper show
% depend on the initial vectors? Preset size (500 neurons), input noise
% 0.025, sra1 at 500 Hz, the top-K estimator the paper runs use (K = 15, QR
% every 0.05 s), four seeds of the initial basis on ONE trajectory per case.
%  fig2G   run_paper_illustrations, 'representative_dynamics': T = [-15, 20],
%          step 0.25 over [20/3, 40/3] s, alignment [-10, 0], window [0, 20].
%  s5fast  a near-default job of the 1-D sweeps at run mode 'fast' (the data
%          the figures are drawn from now): T = [0, 10], window [5, 10],
%          warm-up 5 s.
%  s5prod  the same at 'production': T = [0, 50], window [20, 50], warm-up 5 s.
% The jobs of the sweeps differ from this one in seeds and in the swept
% parameter; this is the default network of the preset, seeds [1 2]. (The
% lya_K_auto of the sweeps can raise K to 60; the first K columns of the
% initial basis are the same for every K, so the first vector is the same.)
[ddir] = s2b_out_dirs();
pn = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
conds = {'sfa1_std1', 'sfa3_stdc2'};
seeds = [424243, 11, 12, 13];          % 424243 is the class default, the one the paper runs used
out = struct('cond', {}, 'case', {}, 'interval', {}, 'warmup', {}, 'C', {});
for c = 1:2
    for run = 1:2
        if run == 1
            input = struct('intrinsic_drive', 0, 'step_time', 20/3, 'step_off', 40/3, 'amplitude', 0.25, 'generator', @paper_midpoint_input);
            args = {'T_range', [-15 20], 'input_config', input};
            cases = {'fig2G', [0 20], 10};
        else
            args = {'T_range', [0 50]};
            cases = {'s5fast', [5 10], 5; 's5prod', [20 50], 5};
        end
        t0 = tic;
        model = build_from_preset(pn, conds{c}, args{:}, 'fs', 500, 'rng_seeds', [1 2], 'ode_solver', 'sra1', ...
            'verbose', 'near-none', 'lya_method', 'none', 'store_full_state', true);
        model.run();
        params = model.get_params();  params.u_interpolant = model.u_interpolant;
        t = model.t_out(:);  S = model.S_out;  clear model
        fprintf('%s run %d: n = %d, N = %d, simulated in %.0f s\n', conds{c}, run, params.n, size(S, 2), toc(t0));
        for q = 1:size(cases, 1)
            if strcmp(cases{q, 1}, 's5fast'); keep = t <= 10 + 1e-9; else; keep = true(size(t)); end
            C = s2b_topk_seeds(t(keep), S(keep, :), 500, params, cases{q, 2}, cases{q, 3}, seeds, 15, 0.05);
            out(end + 1) = struct('cond', conds{c}, 'case', cases{q, 1}, 'interval', cases{q, 2}, 'warmup', cases{q, 3}, 'C', C); %#ok<SAGROW>
            fprintf('   %s: top-K x %d seeds in %.0f s\n', cases{q, 1}, numel(seeds), sum([C.seconds]));
        end
        clear S
        save(fullfile(ddir, 's2b_06_paper_runs.mat'), 'out');
    end
end
s2b_06_report(out);
