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
% See also: ied_stim_run, ied_stim_blocks, ied_stim_figures,
%           fig_ied_stim_blocks, ied_stim_fig_combined, fig_stim_combined
arguments
    run_mode (1,:) char {mustBeMember(run_mode, ...
        {'essential','fast','production','stim_blocks_20261117'})} = 'essential'
    overrides.run_dir (1,:) char = ''
    overrides.fig_root (1,:) char = ''
    overrides.clinical_dir (1,:) char = ''
    overrides.clinical_psd_fig (1,:) char = ''
    overrides.visible_figures (1,1) logical = false
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
    'C:/Users/m218089/Desktop/gitlab_repos/cscs_stim_stability/results/manuscript_stim/tica_awake_best');
cfg.clinical_files = struct('traces', 'Fig_Awake_Main_PanelD.fig', ...
    'rate', 'Fig_Awake_Main_PanelA_Normalized.fig', 'psd', 'Fig_Awake_Main_PanelB_PSD.fig');
% '' = panels/<clinical_files.psd>. A full path takes the PSD panel from
% another .fig (only for testing against a folder whose panels/ predate the
% no-harmonic-removal default).
cfg.clinical_psd_fig = overrides.clinical_psd_fig;

%% The figure registry
% add(F, name, handle, group, extra-args), built by paper_figure_entry;
% ied_stim_figures calls fn('data_dir', cfg.run_dir, 'out_dir',
% <fig_root>/<name>, 'visible', ..., args{:}) in this order. The combined
% figure reads the first entry's .fig, so it comes last. GROUPS (Plan 2,
% package A, 2026-11-21; paper_figure_groups): Fig 8 (Fig_Stim_Combined) is
% 'main'; the other three are 'appendix' (the supplementary state-change
% figure becomes a panel of S10, built by package E; until 2026-11-21 it was
% in_paper true). in_paper follows the group.
F = {};
% Model traces / normalised rate / PSD: the source of Fig 8's panels D-F.
F = add(F, 'fig_ied_stim_blocks',            @fig_ied_stim_blocks, 'appendix', {'which', 'main'});
F = add(F, 'fig_ied_stim_blocks_rate',       @fig_ied_stim_blocks, 'appendix', {'which', 'rate'});
% Supplement: change with stimulation of the mean rate, SFA feedback and
% the increase in depression (1 - b_1), one point per network.
F = add(F, 'fig_ied_stim_blocks_supp_state', @fig_ied_stim_blocks, 'appendix',  {'which', 'supp_state'});
% Fig 8: clinical A-C over model D-F, mirror layout.
F = add(F, 'Fig_Stim_Combined',              @ied_stim_fig_combined, 'main', ...
        {'model_fig', fullfile(cfg.fig_root, 'fig_ied_stim_blocks', 'fig_ied_stim_blocks.fig'), ...
         'clinical_dir', cfg.clinical_dir, 'clinical_files', cfg.clinical_files, ...
         'clinical_psd_fig', cfg.clinical_psd_fig, 'run_name', cfg.run_name, 'sim', cfg.sim});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, group, args)
F{end+1} = paper_figure_entry(name, fn, group, args);
end
function v = pick(v, default)
if isempty(v), v = default; end
end
