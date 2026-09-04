function results = run_overnight_queue(queue)
% RUN_OVERNIGHT_QUEUE Run several full paper pipelines back to back, unattended.
%
%   results = RUN_OVERNIGHT_QUEUE()          % the default queue below
%   results = RUN_OVERNIGHT_QUEUE(queue)     % your own
%
% Each entry is one COMPLETE run -- every analysis stage and then every figure,
% i.e. exactly what a *_run.m script in scripts/paper does -- driven by one
% paper_config-shaped cfg. Entries run sequentially: each saturates the parallel
% pool, so two at once would only make both slower.
%
% queue is a cell array whose entries are either
%   - a cfg struct (from paper_config or one of its wrappers), or
%   - a function handle returning one, e.g. @single_multi_TS_config.
%
% PRE-FLIGHT before any compute: every entry's preset is resolved, its run mode
% checked against run_mode_names, and its run_dir checked to be absent or empty
% (run_all_paper_analyses refuses anything else, and it is better to learn that
% in seconds than after entry 1 has run for two hours).
%
% UPDATED 2026-09-04. This used to call run_all_analyses(preset, run_mode) --
% the SWEEP pipeline only -- so a queue entry produced no memory capacity, no
% eig heatmap, no figures, and wrote to a dated directory under data/param_space
% that nothing downstream was pointed at. It now runs the same two entry points
% the scripts/paper runners do, against the named run_dir and fig_root, so an
% overnight entry and a hand-launched run are the same thing.
%
% See also: run_all_paper_analyses, make_all_paper_figures, paper_config,
%           run_mode_names

arguments
    queue = default_queue()
end

setup_paths();
project_root = fileparts(which('setup_paths'));

if ~iscell(queue); queue = num2cell(queue); end
n = numel(queue);

%% Pre-flight
cfgs = cell(n, 1);
fprintf('\n========================================\n');
fprintf('OVERNIGHT QUEUE: %d full pipelines\n', n);
fprintf('Start: %s\n', datetime('now'));
fprintf('========================================\n');
valid_modes = run_mode_names();
for q = 1:n
    entry = queue{q};
    if isa(entry, 'function_handle'); entry = entry(); end
    if ~isstruct(entry) || ~isfield(entry, 'preset_name') || ~isfield(entry, 'run_dir')
        error('run_overnight_queue:BadEntry', ...
            'Queue entry %d is not a paper_config-shaped cfg struct.', q);
    end
    cfgs{q} = entry;

    if ~ismember(entry.run_mode, valid_modes)
        error('run_overnight_queue:BadRunMode', ...
            'Entry %d: run_mode ''%s'' is not one of %s.', ...
            q, entry.run_mode, strjoin(valid_modes, ', '));
    end
    [d_chk, mc_chk] = srnn_param_preset(entry.preset_name);      % errors on a bad name
    cfg_chk = analysis_run_config('sensitivity', entry.run_mode, d_chk);

    rd = entry.run_dir;
    if ~isempty(rd) && ~is_absolute(rd); rd = fullfile(project_root, rd); end
    if ~isempty(rd) && isfolder(rd)
        listing = dir(rd); listing = listing(~ismember({listing.name}, {'.', '..'}));
        if ~isempty(listing)
            error('run_overnight_queue:TargetNotEmpty', ...
                'Entry %d: run_dir already holds %d items: %s', q, numel(listing), rd);
        end
    end

    if isfield(d_chk, 'sigma_u_noise'); sig = d_chk.sigma_u_noise; else; sig = 0; end
    fprintf('  [%d] %s\n', q, entry.preset_name);
    fprintf('      class=%s mode=%s sigma_u_noise=%g integrator=%s fs=%d\n', ...
        mc_chk, entry.run_mode, sig, cfg_chk.model.ode_solver, cfg_chk.model.fs);
    fprintf('      run_dir=%s  fig_root=%s\n', entry.run_dir, entry.fig_root);
end
fprintf('========================================\n\n');

%% Run
queue_t0 = tic;
results  = cell(n, 1);
for q = 1:n
    cfg = cfgs{q};
    label = sprintf('%s (%s)', cfg.preset_name, cfg.run_mode);

    fprintf('\n\n############################################################\n');
    fprintf('# QUEUE %d/%d: %s\n', q, n, label);
    fprintf('# run_dir : %s\n', cfg.run_dir);
    fprintf('# fig_root: %s\n', cfg.fig_root);
    fprintf('# started : %s\n', datetime('now'));
    fprintf('############################################################\n\n');

    job_t0 = tic;
    r = struct('label', label, 'ok', false, 'minutes', 0, 'run_dir', '(none)', ...
        'figs_ok', 0, 'figs_total', 0, 'err', '');
    try
        r.run_dir = run_all_paper_analyses(cfg);
        fig_results = make_all_paper_figures(cfg);
        r.figs_ok    = sum([fig_results.ok]);
        r.figs_total = numel(fig_results);
        r.ok = true;
        r.minutes = toc(job_t0) / 60;
        fprintf('\n### QUEUE %d/%d DONE in %.1f min -> %s  (%d/%d figures)\n', ...
            q, n, r.minutes, r.run_dir, r.figs_ok, r.figs_total);
    catch job_err
        r.minutes = toc(job_t0) / 60;
        r.err = job_err.message;
        fprintf(2, '\n### QUEUE %d/%d FAILED after %.1f min\n', q, n, r.minutes);
        fprintf(2, '### %s: %s\n', job_err.identifier, job_err.message);
        for k = 1:numel(job_err.stack)
            fprintf(2, '###   at %s (line %d)\n', job_err.stack(k).name, job_err.stack(k).line);
        end
        fprintf(2, '### continuing with the rest of the queue\n');
    end
    results{q} = r;
    close all force;
end

%% Summary
fprintf('\n\n========================================\n');
fprintf('OVERNIGHT QUEUE COMPLETE\n');
fprintf('Total: %.2f hours   End: %s\n', toc(queue_t0) / 3600, datetime('now'));
fprintf('========================================\n');
for q = 1:n
    r = results{q};
    if r.ok; status = 'OK    '; else; status = 'FAILED'; end
    fprintf('  [%d] %s  %-58s %6.1f min\n', q, status, r.label, r.minutes);
    fprintf('        %s\n', r.run_dir);
    if r.ok; fprintf('        figures: %d/%d\n', r.figs_ok, r.figs_total); end
    if ~r.ok; fprintf('        error: %s\n', r.err); end
end
fprintf('========================================\n');
end

%% ------------------------------------------------------------------------
function q = default_queue()
% The two mu x1.5 experiments (TR, 2026-09-04): the paper's 3-condition network
% with 50% stronger mean connectivity, then the same with three depression
% timescales. Both at 'medium', as their config files say.
q = { @sing_multi_TS_50percentStronger_config, ...
      @sing_multi_TS_50percentStronger_std3_config };
end

function tf = is_absolute(p)
tf = startsWith(p, filesep) || startsWith(p, '/') || ...
     ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end
