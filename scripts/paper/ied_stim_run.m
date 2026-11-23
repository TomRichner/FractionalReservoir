function [run_dir, fig_results] = ied_stim_run(run_mode, varargin)
% IED_STIM_RUN The stimulation bundle: the IED block run, its two stages and every figure, Figure 8 and S10 included (FD-031 items 4, 5; Plan 2).
%
%   ied_stim_run                  % 'essential': 2 networks, 1 + 1 blocks of 60 s, plus the two stages
%   ied_stim_run('fast')          % 3 networks, 2 + 2 blocks of 120 s
%   ied_stim_run('production')    % the paper's run: 12 networks, 4 + 4 five-minute blocks (Tom's to start)
%   ied_stim_run(mode, 'clinical_dir', d, ...)   % name-value overrides of ied_stim_config
%   [run_dir, fig_results] = ied_stim_run(...)
%
% setup_paths is called on the first line. Every setting comes from
% ied_stim_config(run_mode). The session's pool allocation must be set
% first (pool_set(N)); every step uses pool_ensure. In order, all into
% cfg.run_dir, so the bundle is self-contained (Plan 2, integration,
% 2026-11-21: no figure of it reads another bundle's run):
%   1. ied_stim_blocks              the modular network's blocks, summary.mat
%   2. run_ied_network_mode_stage   <run_dir>/network_mode (appendix figure)
%   3. burst_mode                   <run_dir>/burst_mode: the bursting
%                                   network of Figure 8 A to C and S10 B to
%                                   D (burst_mode_stage_plan on
%                                   cfg.burst_mode_stage, run_burst_mode_stage)
%   4. ied_stim_figures             every registered figure into cfg.fig_root:
%      fig_ied_stim_blocks/ (_rate, _supp_state)  the model's panels
%      Fig_Stim_Combined/   the old Figure 8 (appendix)
%      Fig_Main8_Stim/      Figure 8 (clinical recordings: never git add -f)
%      fig_supp10_grouped/  S10
%      fig_burst_mode/, fig_ied_network_mode/   appendix
%
% RESUMING. The blocks resume (finished blocks are skipped) and refuse a run
% directory written with other settings. Steps 2 and 3 are skipped when
% their index is there and lists every network ok (a stage starts over when
% it runs, so delete its folder to rerun it). A step that fails is reported
% and the figures are still drawn (those that read it fail on their own).
% The legacy mode 'stim_blocks_20261117' is redrawn, never simulated: no
% step 1 to 3. To redraw a finished run without simulating, call
% ied_stim_figures(ied_stim_config(run_mode)).
%
% HOW LONG (workstation, 2026-11-21; ied_stim_config and the stages' plans):
% at 'essential' the blocks 1.6 min (3 workers), network_mode a few
% minutes (1 network), burst_mode 3 networks at 3 to 4 min each; at 'fast'
% the blocks 12 jobs, network_mode 9.7 min for 3 networks (2 workers),
% burst_mode 9 networks in 5.3 min (2 workers).
arguments
    run_mode (1,:) char = 'essential'
end
arguments (Repeating)
    varargin
end
setup_paths();
cfg = ied_stim_config(run_mode, varargin{:});
legacy = strcmp(run_mode, 'stim_blocks_20261117');
assert(~legacy || isfile(fullfile(cfg.run_dir, 'summary.mat')), ...
    'ied_stim_run:Legacy', ['stim_blocks_20261117 is the finished 400-Hz run; it is redrawn, ' ...
    'not resimulated. Use ied_stim_figures(ied_stim_config(''stim_blocks_20261117'')).']);
% Both stage plans first, so a bad setting fails before any compute.
bm_plan = burst_mode_stage_plan(cfg.burst_mode_stage);
nm_plan = ied_network_mode_stage_plan(cfg);
t = tic;
steps = struct('name', {}, 'ok', {}, 'minutes', {}, 'detail', {});
if ~legacy
    sim_args = namedargs2cell(cfg.sim);
    steps = run_step(steps, 'blocks', @() blocks(cfg, sim_args));
    steps = run_step(steps, 'network_mode', @() stage_unless_done(fullfile(cfg.run_dir, nm_plan.folder), ...
        'ied_netmode_index.mat', @() run_ied_network_mode_stage(cfg, nm_plan)));
    steps = run_step(steps, 'burst_mode', @() stage_unless_done(fullfile(cfg.run_dir, bm_plan.folder), ...
        'burst_mode_index.mat', @() run_burst_mode_stage(cfg.burst_mode_stage, cfg.run_dir, bm_plan)));
end
fprintf('[ied_stim_run] simulation steps done in %.1f min\n', toc(t)/60);
fig_results = ied_stim_figures(cfg);
run_dir = cfg.run_dir;

fprintf('\n========================================================\n');
fprintf('IED STIMULATION BLOCK RUN COMPLETE (%s)\n', run_mode);
fprintf('  fs      : %g Hz\n', cfg.sim.fs);
fprintf('  run_dir : %s\n', cfg.run_dir);
for k = 1:numel(steps)
    fprintf('  %-12s: %s  %.1f min  %s\n', steps(k).name, pass_fail(steps(k).ok), steps(k).minutes, steps(k).detail);
end
fprintf('  figures : %d of %d succeeded  (%s)\n', sum([fig_results.ok]), numel(fig_results), cfg.fig_root);
fprintf('  clinical: %s\n', cfg.clinical_dir);
fprintf('  total   : %.1f min\n', toc(t)/60);
fprintf('========================================================\n');
if nargout == 0, clear run_dir fig_results; end
end

%% ------------------------------------------------------------------------
function detail = blocks(cfg, sim_args)
ied_stim_blocks('data_dir', cfg.run_dir, sim_args{:});
detail = cfg.run_dir;
end

function detail = stage_unless_done(folder, index_name, run)
% Run a stage unless its index lists every network ok.
idx = fullfile(folder, index_name);
if isfile(idx)
    L = load(idx, 'runs');
    if isfield(L, 'runs') && ~isempty(L.runs) && all(L.runs.ok)
        detail = sprintf('skipped: %s is complete', idx);
        return
    end
end
detail = run();
end

function steps = run_step(steps, name, fn)
t = tic; ok = true;
try
    detail = fn();
    if ~ischar(detail); detail = ''; end
catch ME
    ok = false; detail = sprintf('%s: %s', ME.identifier, ME.message);
    fprintf(2, '[ied_stim_run] %s FAILED: %s\n', name, detail);
end
steps(end + 1) = struct('name', name, 'ok', ok, 'minutes', toc(t) / 60, 'detail', detail);
fprintf('[ied_stim_run] %-12s %s  %.1f min\n', name, pass_fail(ok), steps(end).minutes);
end

function s = pass_fail(ok)
if ok, s = 'ok'; else, s = 'FAILED'; end
end
