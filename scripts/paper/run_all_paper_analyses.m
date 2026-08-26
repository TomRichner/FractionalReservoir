function run_dir = run_all_paper_analyses(cfg)
% RUN_ALL_PAPER_ANALYSES All heavy compute for the paper, into one run directory.
%
%   run_dir = RUN_ALL_PAPER_ANALYSES()                                  % production
%   run_dir = RUN_ALL_PAPER_ANALYSES(paper_config('run_mode', 'fast'))  % smoke test
%   run_dir = RUN_ALL_PAPER_ANALYSES(paper_config())
%
% Takes ONE argument, a config STRUCT. It does not accept name-value pairs --
% those belong to paper_config, which is the single place a run is configured.
% (This header used to advertise RUN_ALL_PAPER_ANALYSES('run_mode', 'fast'),
% which errors with "Too many input arguments".)
%
% ENTRY POINT ONE OF TWO. This is the expensive half: everything that simulates
% for hours. make_all_paper_figures is the cheap half, and reads what this
% writes.
%
% Stages, all into a single dated data/param_space/run_all_<dt>/ directory:
%
%   1. the sweep pipeline   run_all_analyses -- seven 1-D sensitivity sweeps,
%                           the tau_a sweep, and the 7-D parameter-space grid
%                           (randomly subsampled), each under every adaptation
%                           condition the preset defines
%   2. memory capacity      run_memory_capacity -- paired trials across the four
%                           conditions
%   3. MC example           run_memory_capacity_example -- one network, kept
%                           per-delay reconstructions
%   4. eig heatmap          run_eig_heatmap -- pooled Jacobian eigenvalues
%
% Stage 1 writes run_manifest.mat, which is what make_all_paper_figures uses to
% find this run later. Stages 2-4 write into subfolders of it, so a run
% directory holds everything the paper was built from.
%
% ERROR ISOLATION. Each stage is wrapped: a failure is reported and the queue
% continues, so an overnight run does not lose four hours of finished sweeps to
% a crash in the last twenty minutes. The summary says what ran and what did not.
%
% WHICH NETWORK, and how much compute, come from paper_config -- edit that file,
% not this one.
%
% See also: paper_config, make_all_paper_figures, run_all_analyses,
%           run_memory_capacity, run_eig_heatmap

arguments
    cfg = paper_config()
end

setup_paths();
if ~isstruct(cfg); error('run_all_paper_analyses:BadConfig', 'cfg must be a struct.'); end

fprintf('\n========================================================\n');
fprintf('PAPER ANALYSES\n');
fprintf('  preset   : %s\n', cfg.preset_name);
fprintf('  run_mode : %s\n', cfg.run_mode);
fprintf('  start    : %s\n', datetime('now'));
fprintf('========================================================\n');
t_all = tic;

%% 1. The sweep pipeline -- creates the run directory and its manifest
run_dir = run_all_analyses(cfg.preset_name, cfg.run_mode);

results = struct('stage', {}, 'ok', {}, 'minutes', {}, 'detail', {}, 'err', {});
results = record(results, 'sweeps', true, toc(t_all)/60, run_dir, '');

%% 2-4. The figure-specific compute
% Each writes into the run directory, so the whole paper's inputs sit together.
stages = { ...
    'memory_capacity', @() run_memory_capacity( ...
        'preset_name', cfg.mc_preset, 'run_mode', cfg.run_mode, ...
        'output_dir', run_dir); ...
    'mc_example',      @() run_memory_capacity_example( ...
        'preset_name', cfg.mc_preset, 'run_mode', cfg.run_mode, ...
        'output_dir', fullfile(run_dir, 'mc_example')); ...
    'eig_heatmap',     @() run_eig_heatmap( ...
        'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
        'out_dir', fullfile(run_dir, 'eig_heatmap')) };

for k = 1:size(stages, 1)
    name = stages{k, 1};
    fprintf('\n========================================\n');
    fprintf('[%d/%d] %s\n', k, size(stages, 1), name);
    fprintf('========================================\n');
    t0 = tic;
    try
        restart_parpool();
        detail = stages{k, 2}();
        results = record(results, name, true, toc(t0)/60, detail, '');
        fprintf('  -> %s\n', detail);
    catch ME
        results = record(results, name, false, toc(t0)/60, '', ME.message);
        fprintf(2, '  FAILED: %s: %s\n', ME.identifier, ME.message);
        for s = 1:min(3, numel(ME.stack))
            fprintf(2, '    at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
        end
        fprintf(2, '  continuing with the remaining stages\n');
    end
end

%% Refresh the human-readable record
% write_run_parameters_md is called by run_all_analyses too, but stages 2-4 add
% subfolders after that; regenerating picks them up.
try
    write_run_parameters_md(run_dir);
catch ME
    warning('run_all_paper_analyses:ParametersMdFailed', ...
        'Could not refresh parameters.md: %s', ME.message);
end

%% Summary
fprintf('\n========================================================\n');
fprintf('PAPER ANALYSES COMPLETE in %.2f h\n', toc(t_all)/3600);
fprintf('  run_dir: %s\n', run_dir);
fprintf('--------------------------------------------------------\n');
for k = 1:numel(results)
    r = results(k);
    if r.ok; tag = 'OK    '; else; tag = 'FAILED'; end
    fprintf('  %s  %-18s %7.1f min\n', tag, r.stage, r.minutes);
    if ~r.ok; fprintf(2, '          %s\n', r.err); end
end
fprintf('========================================================\n');
fprintf('Next: make_all_paper_figures(paper_config(''run_dir'', run_dir))\n');
end

%% ------------------------------------------------------------------------
function results = record(results, stage, ok, minutes, detail, err)
if ~ischar(detail); detail = ''; end
results(end+1) = struct('stage', stage, 'ok', ok, 'minutes', minutes, ...
    'detail', detail, 'err', err);
end
