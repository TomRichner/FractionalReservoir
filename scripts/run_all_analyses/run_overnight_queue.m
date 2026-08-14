%% run_overnight_queue.m
% Run several full run_all_analyses pipelines back to back, unattended.
%
% Each entry of the QUEUE below is one complete pipeline -- sensitivity, tau
% sensitivity and param space -- with its own preset, run mode and output
% directory. They run sequentially because each already saturates the parallel
% pool; running two at once would only make both slower and confuse the
% per-batch timings.
%
% WHY A SCRIPT RATHER THAN THREE CONSOLE CALLS. Two reasons that matter for an
% overnight run:
%
%   * ERROR ISOLATION. Each pipeline is wrapped in try/catch, so a failure in
%     one does not cost the others the night. The error is recorded and the
%     queue moves on; the summary at the end says which ran and which did not.
%   * STATE ISOLATION. run_all_analyses and its sub-scripts communicate through
%     base-workspace variables (master_*), and they read them with exist()
%     rather than taking arguments. A variable left behind by one entry would
%     silently apply to the next, so every master_* is cleared between entries
%     rather than overwritten.
%
% Progress is echoed with timestamps so the log can be read after the fact to
% see where the night went.
%
% See also: run_all_analyses, srnn_param_preset, analysis_run_config

setup_paths();

%% The queue
% label    : appears in the log only
% preset   : srnn_param_preset name
% run_mode : analysis_run_config mode
%
% All three are 'medium', so all three share its timings -- fs 400,
% T_range [0 20], 11 levels x 15 reps. The only differences are the preset's
% own physics: whether there is noise, and how much. The integrator follows
% from that automatically, sra1 for the two stochastic entries and rk4 for the
% deterministic one, chosen by analysis_run_config rather than named here.
queue = {
    struct('label', 'stochastic noise 0.01', ...
           'preset', 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01', ...
           'run_mode', 'medium')
    struct('label', 'deterministic', ...
           'preset', 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5', ...
           'run_mode', 'medium')
    struct('label', 'stochastic noise 0.02', ...
           'preset', 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p02', ...
           'run_mode', 'medium')
    };

%% Pre-flight
% Resolve every preset and report the integrator BEFORE committing to hours of
% compute, so a typo in a name fails in seconds rather than after entry 1.
fprintf('\n========================================\n');
fprintf('OVERNIGHT QUEUE: %d pipelines\n', numel(queue));
fprintf('Start: %s\n', datetime('now'));
fprintf('========================================\n');
for q = 1:numel(queue)
    job = queue{q};
    [d_chk, mc_chk] = srnn_param_preset(job.preset);      % errors here on a bad name
    cfg_chk = analysis_run_config('sensitivity', job.run_mode, d_chk);
    if isfield(d_chk, 'sigma_u_noise'), sig = d_chk.sigma_u_noise; else, sig = 0; end
    fprintf('  [%d] %-24s %s\n', q, job.label, job.preset);
    fprintf('      class=%s mode=%s sigma_u_noise=%g integrator=%s fs=%d\n', ...
        mc_chk, job.run_mode, sig, cfg_chk.model.ode_solver, cfg_chk.model.fs);
end
fprintf('========================================\n\n');
clear d_chk mc_chk cfg_chk sig job q;

%% Run
queue_t0      = tic;
queue_results = cell(numel(queue), 1);

for q_idx = 1:numel(queue)
    job = queue{q_idx};

    % Clear everything the pipeline communicates through, so nothing carries
    % over from the previous entry. run_all_analyses sets master_output_dir,
    % master_model_overrides, master_model_class and master_conditions itself.
    clear master_output_dir master_model_overrides master_model_class ...
          master_conditions master_save_figs ...
          run_mode preset_name cfg psa psa_tau_a psa_tau_b preset_defaults ...
          model_class conditions sens_replot_dir all_output_dirs;

    run_mode    = job.run_mode;
    preset_name = job.preset;

    fprintf('\n\n############################################################\n');
    fprintf('# QUEUE %d/%d: %s\n', q_idx, numel(queue), job.label);
    fprintf('# preset  : %s\n', job.preset);
    fprintf('# mode    : %s\n', job.run_mode);
    fprintf('# started : %s\n', datetime('now'));
    fprintf('############################################################\n\n');

    job_t0 = tic;
    try
        run_all_analyses;
        queue_results{q_idx} = struct('label', job.label, 'ok', true, ...
            'minutes', toc(job_t0), 'dir', master_output_dir, 'err', '');
        fprintf('\n### QUEUE %d/%d DONE in %.1f min -> %s\n', ...
            q_idx, numel(queue), toc(job_t0)/60, master_output_dir);
    catch job_err
        if exist('master_output_dir', 'var'), od = master_output_dir; else, od = '(none)'; end
        queue_results{q_idx} = struct('label', job.label, 'ok', false, ...
            'minutes', toc(job_t0), 'dir', od, 'err', job_err.message);
        fprintf(2, '\n### QUEUE %d/%d FAILED after %.1f min\n', q_idx, numel(queue), toc(job_t0)/60);
        fprintf(2, '### %s: %s\n', job_err.identifier, job_err.message);
        for k = 1:numel(job_err.stack)
            fprintf(2, '###   at %s (line %d)\n', job_err.stack(k).name, job_err.stack(k).line);
        end
        fprintf(2, '### continuing with the rest of the queue\n');
    end
    close all force;
end

%% Summary
fprintf('\n\n========================================\n');
fprintf('OVERNIGHT QUEUE COMPLETE\n');
fprintf('Total: %.2f hours   End: %s\n', toc(queue_t0)/3600, datetime('now'));
fprintf('========================================\n');
for q_idx = 1:numel(queue_results)
    r = queue_results{q_idx};
    if r.ok, status = 'OK    '; else, status = 'FAILED'; end
    fprintf('  [%d] %s  %-30s %6.1f min\n', q_idx, status, r.label, r.minutes);
    fprintf('        %s\n', r.dir);
    if ~r.ok, fprintf('        error: %s\n', r.err); end
end
fprintf('========================================\n');
