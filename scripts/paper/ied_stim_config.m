function cfg = ied_stim_config(run_mode, overrides)
% IED_STIM_CONFIG The IED stimulation block run (Fig 8's model half and its supplement) as a paper bundle (FD-031 items 4, 5, 14).
%
%   cfg = ied_stim_config()                        % 'essential'
%   cfg = ied_stim_config('fast')
%   cfg = ied_stim_config('production')            % the paper's run
%   cfg = ied_stim_config('stim_blocks_20261117')  % the 400-Hz run Fig 8 was drawn from
%   cfg = ied_stim_config(mode, 'clinical_dir', d, 'run_dir', r, 'fig_root', f, ...)
%
% Run it with ied_stim_run(run_mode): ied_stim_blocks simulates into
% cfg.run_dir, then ied_stim_figures draws every entry of cfg.figures into
% cfg.fig_root/<name>/. ied_stim_figures(cfg) alone redraws from a run
% directory that is already complete (no simulation).
%
% WHAT IT RUNS. The frozen run13 modular IED-like network (ied_base_config:
% n = 350, five subnets, SFA and STD), in separate five-minute blocks, like
% the clinical recordings' no-stimulation and stimulation blocks. Per
% network seed: n_no blocks with uniform input 0, then n_stim blocks with
% uniform input cfg.sim.amplitude (0.2), constant through the block. Same
% network (connectivity, set-points from the seed) in every block, a new
% noise seed per block: base.noise_seed + (seed - 42) + 1000*block. Each
% block discards the model warmup (10 s) and a settle (10 s); the next
% block_len seconds are analysed (events, subnet-mean PSD, rate, SFA
% feedback, depression 1 - b_1).
%
% fs. 500 Hz in every mode but 'stim_blocks_20261117' (Tom, 2026-11-17,
% FD-031 item 14). The IED model does not fix its own fs: ied_base_config
% carries 400 as a plain value, and this config overrides it (cfg.sim.fs),
% so this is the one place that sets it. The 2026-11-17 data are 400 Hz and
% stay as they are; their mode states 400 so they can be redrawn.
%
% MODES (simulation jobs = networks x blocks; measured 2026-11-17 on a pool
% worker at 500 Hz: about 0.31 s per simulated second, so about 100 s per
% production block of 320 s; figures about 1 min more):
%   mode                  networks  blocks   analysed   jobs   about, 3 workers
%   essential                 2      1 + 1     60 s        4     1.6 min in all (measured)
%   fast                      3      2 + 2    120 s       12     (see ied_stim_run)
%   production               12      4 + 4    300 s       96    55 min (12 workers: 15 min)
%   stim_blocks_20261117     12      4 + 4    300 s       96    (exists; 400 Hz; redrawn only)
% 'essential' and 'fast' exist to test the pipeline end to end; their rates
% are not the paper's numbers. Agents run those two; 'production' is Tom's.
%
% WHERE THINGS LAND (relative to the repo root unless overridden):
%   cfg.run_dir   data/IED_model/<run_name>      one .mat per block, summary.mat/.csv
%   cfg.fig_root  figs/IED_model/<run_name>      one folder per figure
% with run_name = 'ied_stim_<mode>' (the dated legacy mode keeps its name).
% The simulation resumes: a block whose file exists is skipped, and the run
% refuses a run directory written with other settings (run_settings.mat).
% To rerun from scratch, delete the run directory.
%
% THE CLINICAL HALF of Fig 8 is read, never written: cfg.clinical_dir is the
% ONE setting, a results folder of cscs_stim_stability's
% run_manuscript_stim_awake (Tom chose tica_awake_best, 2026-11-17); the
% figure reads its panels/ subfolder (cfg.clinical_files). The combined
% figure holds clinical recordings: it is written under figs/ (gitignored)
% and must never be force-added to git.
%
% FIGURE 8 IN THE NEW ORDER (Plan 2, package E; Tom, 2026-11-21): entry
% Fig_Main8_Stim draws A to C the bursting network from the stage
% burst_mode, D to F the clinical panels, G to I this bundle's model.
% STIMULATION IS SELF-CONTAINED IN THIS BUNDLE (Plan 2, integration,
% 2026-11-21): ied_stim_run runs burst_mode into <run_dir>/burst_mode
% (burst_mode_stage_plan, on cfg.burst_mode_stage: this mode's size, the
% bursting preset, fs and the paper's integrator), and Fig_Main8_Stim, S10
% and fig_burst_mode read it from there. cfg.burst_mode_dir is that run
% directory, cfg.run_dir; the override 'burst_mode_dir' exists only to draw
% from another run by hand (the legacy mode's run has no burst_mode/).
% Until then it was the mu7chain_AC run of the same mode.
% Fig_Stim_Combined (the old order) is kept, in the appendix.
%
% THE NETWORK-MODE STAGE (appendix, kept separate): run_ied_network_mode_stage(cfg)
% writes <run_dir>/network_mode (ied_network_mode_stage_plan has the sizes
% by mode); ied_stim_run runs it after the blocks; entry fig_ied_network_mode
% draws it. cfg.network_mode holds its knobs (empty: the mode's sizes).
%
% GROUPS (paper_figure_groups): Fig_Main8_Stim 'main' (Figure 8),
% fig_supp10_grouped 'supplement' (S10), every other entry 'appendix'.
%
% See also: ied_stim_run, ied_stim_blocks, ied_stim_figures,
%           fig_ied_stim_blocks, ied_stim_fig_combined, fig_stim_combined,
%           ied_network_mode_stage_plan, fig_ied_network_mode
arguments
    run_mode (1,:) char {mustBeMember(run_mode, ...
        {'essential','fast','production','stim_blocks_20261117'})} = 'essential'
    overrides.run_dir (1,:) char = ''
    overrides.fig_root (1,:) char = ''
    overrides.clinical_dir (1,:) char = ''
    overrides.clinical_psd_fig (1,:) char = ''
    overrides.visible_figures (1,1) logical = false
    overrides.burst_mode_dir (1,:) char = ''
end
root = fileparts(which('setup_paths'));
cfg = struct();
cfg.bundle   = 'ied_stim';
cfg.run_mode = run_mode;

%% The experiment
seeds12 = [42 73 101 137 211 307 409 503 601 709 811 907];
sim = struct('amplitude', 0.20, 'settle', 10, 'fs', 500);
switch run_mode
    case 'essential'
        sim.seeds = seeds12(1:2);  sim.n_no = 1; sim.n_stim = 1; sim.block_len = 60;
    case 'fast'
        sim.seeds = seeds12(1:3);  sim.n_no = 2; sim.n_stim = 2; sim.block_len = 120;
    case 'production'
        sim.seeds = seeds12;       sim.n_no = 4; sim.n_stim = 4; sim.block_len = 300;
    case 'stim_blocks_20261117'
        sim.seeds = seeds12;       sim.n_no = 4; sim.n_stim = 4; sim.block_len = 300;
        sim.fs = 400;              % as simulated on 2026-11-17
end
cfg.sim = sim;

%% Where things land
if strcmp(run_mode, 'stim_blocks_20261117')
    cfg.run_name = 'stim_blocks_20261117';
else
    cfg.run_name = ['ied_stim_' run_mode];
end
cfg.run_dir  = pick(overrides.run_dir,  fullfile(root, 'data', 'IED_model', cfg.run_name));
cfg.fig_root = pick(overrides.fig_root, fullfile(root, 'figs', 'IED_model', cfg.run_name));
cfg.visible_figures = overrides.visible_figures;

%% The clinical figure folder (the one setting) and its panel files
cfg.clinical_dir = pick(overrides.clinical_dir, ...
    'C:/Users/m218089/Desktop/gitlab_repos/cscs_stim_stability/results/manuscript_stim/tica_awake_best_fig8_drawn_only');
cfg.clinical_files = struct('traces', 'Fig_Awake_Main_PanelD.fig', ...
    'rate', 'Fig_Awake_Main_PanelA_Normalized.fig', 'psd', 'Fig_Awake_Main_PanelB_PSD.fig');
% '' = panels/<clinical_files.psd>. A full path takes the PSD panel from
% another .fig (only for testing against a folder whose panels/ predate the
% no-harmonic-removal default).
cfg.clinical_psd_fig = overrides.clinical_psd_fig;

%% Figure 8 A to C: the burst_mode stage, in this bundle's run directory
% ied_stim_run runs it (burst_mode_stage_plan(cfg.burst_mode_stage), then
% run_burst_mode_stage into cfg.run_dir) after the blocks. Its config: this
% mode's size (the legacy mode, which is redrawn and never simulated, is
% given production's), the bursting figures' preset, this bundle's rate and
% the paper's integrator (SRA1, as every paper simulation; Tom,
% 2026-11-17). Knobs of the stage (cfg.burst_mode_stage.burst_mode_<name>)
% may be added to that struct; none is set here.
stage_mode = run_mode;
if strcmp(run_mode, 'stim_blocks_20261117'); stage_mode = 'production'; end
cfg.burst_mode_stage = struct('run_mode', stage_mode, 'bursting_preset', 'bursting_pairs', ...
    'fs', cfg.sim.fs, 'ode_solver', 'sra1', 'verbose', 'minimal');
cfg.burst_mode_dir = pick(overrides.burst_mode_dir, cfg.run_dir);
cfg.network_mode = struct();    % knobs of ied_network_mode_stage_plan; empty: the mode's sizes

%% The figure registry
% add(F, name, handle, group, extra-args), built by paper_figure_entry;
% ied_stim_figures calls fn('data_dir', cfg.run_dir, 'out_dir',
% <fig_root>/<name>, 'visible', ..., args{:}) in this order. The combined
% figures read the first entry's .fig, so they come after it. Functions
% with the paper's standard signature (run_dir, out_dir, save, visible)
% are called through standard_entry, which maps data_dir to run_dir.
% GROUPS (paper_figure_groups; the figure plan of 2026-11-21): Figure 8
% (Fig_Main8_Stim) 'main', S10 (fig_supp10_grouped) 'supplement', the rest
% 'appendix'. in_paper follows the group.
F = {};
% Model traces / normalised rate / PSD: the source of Fig 8's panels G-I.
F = add(F, 'fig_ied_stim_blocks',            @fig_ied_stim_blocks, 'appendix', {'which', 'main'});
F = add(F, 'fig_ied_stim_blocks_rate',       @fig_ied_stim_blocks, 'appendix', {'which', 'rate'});
% Change with stimulation of the mean rate, SFA feedback and the increase
% in depression (1 - b_1), one point per network (S10 A draws the same).
F = add(F, 'fig_ied_stim_blocks_supp_state', @fig_ied_stim_blocks, 'appendix',  {'which', 'supp_state'});
% The old Fig 8: clinical A-C over model D-F, mirror layout (appendix since
% the figure plan of 2026-11-21).
F = add(F, 'Fig_Stim_Combined',              @ied_stim_fig_combined, 'appendix', ...
        {'model_fig', fullfile(cfg.fig_root, 'fig_ied_stim_blocks', 'fig_ied_stim_blocks.fig'), ...
         'clinical_dir', cfg.clinical_dir, 'clinical_files', cfg.clinical_files, ...
         'clinical_psd_fig', cfg.clinical_psd_fig, 'run_name', cfg.run_name, 'sim', cfg.sim});
% Figure 8: the bursting network A-C (burst_mode, in cfg.burst_mode_dir),
% clinical D-F, this model G-I.
F = add(F, 'Fig_Main8_Stim',                 @ied_stim_fig_combined, 'main', ...
        {'model_fig', fullfile(cfg.fig_root, 'fig_ied_stim_blocks', 'fig_ied_stim_blocks.fig'), ...
         'clinical_dir', cfg.clinical_dir, 'clinical_files', cfg.clinical_files, ...
         'clinical_psd_fig', cfg.clinical_psd_fig, 'run_name', cfg.run_name, 'sim', cfg.sim, ...
         'burst_mode_dir', cfg.burst_mode_dir, 'tag', 'Fig_Main8_Stim'});
% S10: the modular model's state change (this run's summary.mat) and the
% bursting network's onset, crossings and synchrony (burst_mode).
F = add(F, 'fig_supp10_grouped', standard_entry(@fig_supp10_grouped), 'supplement', ...
        {'ied_run_dir', cfg.run_dir, 'burst_dir', cfg.burst_mode_dir});
% Appendix: Figure 8 A to C on their own, and the full Jacobian's network
% mode of the modular network (the network_mode stage).
F = add(F, 'fig_burst_mode', standard_entry(@fig_burst_mode), 'appendix', {'burst_dir', cfg.burst_mode_dir});
F = add(F, 'fig_ied_network_mode',           @fig_ied_network_mode, 'appendix', {});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, group, args)
F{end+1} = paper_figure_entry(name, fn, group, args);
end
function v = pick(v, default)
if isempty(v), v = default; end
end

function h = standard_entry(fn)
% A figure with the paper's standard signature (run_dir, out_dir, save,
% visible, ...) as an entry of this bundle, whose driver passes data_dir,
% out_dir and visible (ied_stim_figures). The burst_mode stage's folder is
% read from 'burst_dir' (the entry's own argument, cfg.burst_mode_dir) as
% the function's run_dir.
h = @(varargin) call_standard(fn, varargin{:});
end

function out = call_standard(fn, varargin)
names = varargin(1:2:end);
values = varargin(2:2:end);
keep = ~ismember(names, {'data_dir', 'burst_dir'});
run_dir = values{find(strcmp(names, 'data_dir'), 1)};
k = find(strcmp(names, 'burst_dir'), 1);
if ~isempty(k); run_dir = values{k}; end
rest = [names(keep); values(keep)];
out = fn('run_dir', run_dir, 'save', true, rest{:});
end
