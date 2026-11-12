% EXAMPLE_MU7CHAIN_MC_ESSENTIAL Memory capacity of the mu7chain network, alone, at the essential run speed.
%
% Runs ONE stage of the mu7chain bundle, the paired-trial memory-capacity
% experiment, on all four adaptation conditions, and draws the
% memory-capacity figure from it: the source of panels C and D of
% manuscript Figure 6 (cumulative memory capacity and per-delay R^2 against
% delay; memory horizon per paired trial), here with four conditions.
%
% It is the bundle run through its own entry points with a stage selection
% (FD-023), so it is the same stage, preset, integrator and figure function
% the full runs use; nothing is reimplemented here.
%
%   network    the bundle's memory-capacity preset: the mu7chain network
%              (n = 500, per-neuron SFA and depression spreads) with the
%              Wiener noise off, integrated with 'sra1'
%   readout    'synaptic': the linear readout sees each neuron's synaptic
%              output (rate times depression factor), so it sees the
%              depression state. Chain, single-pool and product routes all
%              go through it.
%   size       run mode 'essential': 3 PAIRED TRIALS (three network seeds;
%              within a trial every condition shares W, the input weights
%              and the input). The durations are those of every mode and
%              are not shortened: 10 s washout, 600 s of training, 150 s of
%              test, delays to 15 s. A shorter training window would leave
%              the readout with fewer samples than its 500 features and
%              memory capacity under-determined.
%   conditions no adaptation, single-timescale, chained multiple-timescale,
%              product multiple-timescale
%
% WHAT THREE TRIALS CAN AND CANNOT SHOW. The bands and intervals are
% bootstrap intervals over three numbers, so they are rough. No paired test
% can reach significance: with 3 trials the sign-flip test has 8 outcomes
% and its smallest two-sided p is 0.25. Read it for where the conditions
% sit, not for whether two of them differ.
%
% BEFORE RUNNING
%   setup_paths            once per session, cwd at the repository root
%   pool_set(4)            the entry points refuse without a pool
%                          allocation; three trials run at once (each needs
%                          about 16 GB of RAM), so 4 workers are enough
% and the two output folders below must not exist (delete them to rerun).
% It takes about 25 minutes (24 min 12 s when first run: the three trials
% ran one after another, because MATLAB saw room for one at a time).
%
% OUTPUT (both gitignored)
%   data/mu7chain_essential_mc              the run directory
%   figs/mu7chain_essential_mc_4cond        fig_memory_capacity/ with the
%       figure (.png, .svg, .fig) and Fig_Memory_Capacity_table.md, the
%       numbers behind it; report.md and report.pdf beside it
%
% See also: mu7chain_config, run_all_paper_analyses, run_memory_capacity,
%           fig_memory_capacity, paper_stage_plan

%% The bundle at the essential speed, with one stage selected
cfg = mu7chain_config('essential', '4cond');
cfg.stages   = {'memory_capacity'};     % nothing else runs (paper_stage_plan)
cfg.sweeps   = {};
cfg.run_dir  = 'data/mu7chain_essential_mc';
cfg.fig_root = 'figs/mu7chain_essential_mc_4cond';

% The one figure that reads this stage, taken from the bundle's full
% registry so that it is the entry every other run draws.
full_cfg = mu7chain_config('fast', '4cond');
entry_names = cellfun(@(f) f.name, full_cfg.figures, 'UniformOutput', false);
cfg.figures = full_cfg.figures(strcmp(entry_names, 'fig_memory_capacity'));
assert(isscalar(cfg.figures), 'fig_memory_capacity is not in the mu7chain registry');

%% Run the stage, then draw the figure
t_start = tic;
run_dir = run_all_paper_analyses(cfg);
t_analysis = toc(t_start);
fig_results = make_all_paper_figures(cfg);

%% Where things are
project_root = fileparts(which('setup_paths'));
table_file = fullfile(project_root, cfg.fig_root, 'fig_memory_capacity', 'Fig_Memory_Capacity_table.md');
fprintf('\nMemory capacity, essential speed: analysis %.1f min, with the figure %.1f min.\n', ...
    t_analysis / 60, toc(t_start) / 60);
fprintf('  run directory : %s\n', run_dir);
fprintf('  figure        : %s\n', fullfile(project_root, cfg.fig_root, 'fig_memory_capacity'));
if isfile(table_file)
    fprintf('  numbers       : %s\n\n', table_file);
    type(table_file);
end
