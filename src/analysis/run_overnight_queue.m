function results = run_overnight_queue(queue)
% RUN_OVERNIGHT_QUEUE Run several full sweep pipelines back to back, unattended.
%
%   results = RUN_OVERNIGHT_QUEUE()          % the default queue below
%   results = RUN_OVERNIGHT_QUEUE(queue)     % your own
%
% Each entry is one complete pipeline -- sensitivity, tau sensitivity and param
% space -- with its own preset and run mode. They run sequentially because each
% already saturates the parallel pool; running two at once would only make both
% slower and confuse the per-batch timings.
%
% queue is a struct array (or cell of structs) with fields:
%   label     appears in the log only
%   preset    srnn_param_preset name
%   run_mode  analysis_run_config mode
%
% WHAT THIS FILE USED TO BE. Its own header explained that it was a script
% rather than three console calls for two reasons: error isolation, and STATE
% isolation -- because run_all_analyses and its sub-scripts communicated through
% base-workspace `master_*` variables read with exist(), so a variable left
% behind by one entry would silently apply to the next. It therefore carried an
% explicit `clear master_output_dir master_model_overrides ...` of nine
% variables before every entry.
%
% That entire second reason is gone: run_all_analyses is a function, so each
% entry has its own scope and there is nothing to leak. Only error isolation is
% left, which is the try/catch below. The file shrank from 124 lines to this.
%
% See also: run_all_analyses, srnn_param_preset, analysis_run_config

arguments
    queue = default_queue()
end

setup_paths();

if isstruct(queue)
    queue = num2cell(queue);
end
n = numel(queue);

%% Pre-flight
% Resolve every preset and report the integrator BEFORE committing to hours of
% compute, so a typo in a name fails in seconds rather than after entry 1.
fprintf('\n========================================\n');
fprintf('OVERNIGHT QUEUE: %d pipelines\n', n);
fprintf('Start: %s\n', datetime('now'));
fprintf('========================================\n');
for q = 1:n
    job = queue{q};
    [d_chk, mc_chk] = srnn_param_preset(job.preset);      % errors on a bad name
    cfg_chk = analysis_run_config('sensitivity', job.run_mode, d_chk);
    if isfield(d_chk, 'sigma_u_noise'), sig = d_chk.sigma_u_noise; else, sig = 0; end
    fprintf('  [%d] %-24s %s\n', q, job.label, job.preset);
    fprintf('      class=%s mode=%s sigma_u_noise=%g integrator=%s fs=%d\n', ...
        mc_chk, job.run_mode, sig, cfg_chk.model.ode_solver, cfg_chk.model.fs);
end
fprintf('========================================\n\n');

%% Run
queue_t0 = tic;
results  = cell(n, 1);

for q = 1:n
    job = queue{q};

    fprintf('\n\n############################################################\n');
    fprintf('# QUEUE %d/%d: %s\n', q, n, job.label);
    fprintf('# preset  : %s\n', job.preset);
    fprintf('# mode    : %s\n', job.run_mode);
    fprintf('# started : %s\n', datetime('now'));
    fprintf('############################################################\n\n');

    job_t0 = tic;
    try
        run_dir = run_all_analyses(job.preset, job.run_mode);
        results{q} = struct('label', job.label, 'ok', true, ...
            'minutes', toc(job_t0)/60, 'dir', run_dir, 'err', '');
        fprintf('\n### QUEUE %d/%d DONE in %.1f min -> %s\n', ...
            q, n, toc(job_t0)/60, run_dir);
    catch job_err
        results{q} = struct('label', job.label, 'ok', false, ...
            'minutes', toc(job_t0)/60, 'dir', '(none)', 'err', job_err.message);
        fprintf(2, '\n### QUEUE %d/%d FAILED after %.1f min\n', q, n, toc(job_t0)/60);
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
for q = 1:n
    r = results{q};
    if r.ok, status = 'OK    '; else, status = 'FAILED'; end
    fprintf('  [%d] %s  %-30s %6.1f min\n', q, status, r.label, r.minutes);
    fprintf('        %s\n', r.dir);
    if ~r.ok, fprintf('        error: %s\n', r.err); end
end
fprintf('========================================\n');
end

%% ------------------------------------------------------------------------
function q = default_queue()
% The paper's operating point at two fidelities: 'fast' first as a cheap check
% that it behaves, then production for the real numbers. The integrator follows
% from the preset's noise automatically -- sra1 here, chosen by
% analysis_run_config rather than named.
%
% This was 'fast2' until that mode was removed on 2026-09-03. The check is
% therefore cheaper than it was -- 3 reps rather than 6, and half the T_range on
% the 1-D sweeps -- which is fine for "does it behave", the job of this entry.
q = { ...
    struct('label', 'dualStd (fast)', ...
           'preset', 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
           'run_mode', 'fast'), ...
    struct('label', 'dualStd (production)', ...
           'preset', 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
           'run_mode', 'production') ...
    };
end
