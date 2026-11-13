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
%   5. numerics             run_numerics_verification -- SRA1 reshot against a
%                           1e-10 ode45 reference (and against itself on a
%                           shared Brownian path), Benettin vs QR on a reduced
%                           network; feeds two supplemental figures
%   6. DC-LLE               run_dc_lle_analysis -- local Lyapunov exponent vs
%                           tonic DC drive, across seeds and every condition
%
% Stage 1 writes run_manifest.mat, which is what make_all_paper_figures uses to
% find this run later. Stages 2-6 write into subfolders of it, so a run
% directory holds everything the paper was built from.
%
% (The list above is from when there were six; the stages as they run today,
% in order, are paper_stage_names.)
%
% SELECTING STAGES (FD-023). By default every stage runs. A bundle may name
% the ones it wants:
%   cfg.stages   cell of stage names from paper_stage_names; absent or empty
%                = all. Order is ignored; stages run in the order below.
%   cfg.sweeps   which sweeps stage 1 runs ('sensitivity', 'tau',
%                'param_space'); absent or empty = all three.
%   cfg.eig_heatmap_examples   the eig-heatmap stage's examples; absent or
%                empty = its default three.
% OPT-IN STAGES (FD-010) run only when cfg.stages names them; "every stage"
% does not include them (paper_stage_names, third output). One so far:
%   mc_gain_sweep   block C, run_mc_gain_sweep into <run_dir>/mc_gain_sweep:
%                   memory capacity and lambda_1 over gain x bias, on
%                   cfg.mc_preset (the noise-free memory-capacity twin) with
%                   cfg.mc_ode_solver when set (else the stage's 'sra1').
%                   cfg.mc_gain_sweep (optional struct) carries, each
%                   absent or empty = the stage's own table for the mode:
%                     gain_levels, bias_levels, seeds, conditions
%                     model_args, protocol, lya_K, lya_T_sec, lya_acc_sec,
%                     segment_len   (run_mc_gain_sweep's; for tests)
%                     ei_presets    cell of E:I preset names, each run gain
%                                   only, at bias_levels = its own S_c, with
%                                   the same gains, seeds and conditions,
%                                   into <run_dir>/mc_gain_sweep_<suffix>
%                                   (suffix: the name after its last '_',
%                                   e.g. EI3to2); the result is in that
%                                   folder's mc_gain_sweep/ like the main one
%                     ei_seeds      seeds for the E:I runs (absent = seeds)
%                   An unknown field is refused in the pre-flight.
% An unknown name is an error in the pre-flight, before any compute and
% before the run directory is created (paper_stage_plan). A stage not run is
% printed as SKIPPED in the summary and recorded in run_manifest
% (stages_run, stages_skipped, sweeps_run, sweeps_skipped; absent from a run
% that selected nothing), so the directory says what it does not contain.
% The run directory, its provenance and its manifest are made whatever is
% selected. A figure whose stage was skipped errors by design
% (resolve_data_file): the bundle's figure registry must not list it.
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

%% PREFLIGHT: the run mode must be one every stage understands
% Checked HERE, before anything simulates, because the failure it guards was
% silent and expensive. Stage 1 validated run_mode itself, so a mode it accepted
% started hours of sweeps -- and stages 2-5 then hit their own switch statements,
% threw "Unknown run_mode", and were reported as FAILED rows far down a long log.
% On 2026-09-03 single_multi_TS ran that way at 'fast2': the sweeps finished,
% memory capacity / mc_example / eig_heatmap produced nothing, and
% fig_memory_capacity then fell back to a twelve-day-old .mat and reported
% SUCCESS. The stale figure was indistinguishable from a good one.
%
% run_mode_names() is now the single list every stage validates against, and
% test_run_modes asserts each stage really accepts every name in it. This check
% is the cheap front line: it costs microseconds and fails before the output
% directory is created.
valid_modes = run_mode_names();
if ~ismember(cfg.run_mode, valid_modes)
    error('run_all_paper_analyses:BadRunMode', ...
        ['run_mode ''%s'' is not one of: %s.\n\n' ...
         'Checked before any compute starts: a mode that only SOME stages ' ...
         'understand used to run the sweeps to completion and then lose the\n' ...
         'later stages one by one. (''fast2'' was removed 2026-09-03 for ' ...
         'exactly that reason -- use ''fast''.)'], ...
        cfg.run_mode, strjoin(valid_modes, ', '));
end

%% PREFLIGHT: which stages and sweeps this bundle runs (FD-023)
% cfg.stages and cfg.sweeps select; absent or empty means everything, as
% before the fields existed. paper_stage_plan validates the names against
% paper_stage_names and errors on an unknown one HERE, before any compute and
% before the run directory is created, for the same reason the run mode is
% checked here. The run mode sets sizes; this sets which stages run; the two
% are separate (run_mode_names).
plan = paper_stage_plan(cfg);

% cfg.eig_heatmap_examples (optional; absent or empty = the stage's default,
% the reference network and the two E:I-imbalanced examples) is passed to
% run_eig_heatmap as its 'examples'. A bundle whose figures read only the
% reference names that one example and saves the other two.
eig_args = {};
if isfield(cfg, 'eig_heatmap_examples') && ~isempty(cfg.eig_heatmap_examples)
    ex = cfg.eig_heatmap_examples;
    if ~(isstruct(ex) && all(isfield(ex, {'label', 'overrides'})))
        error('run_all_paper_analyses:BadEigExamples', ...
            'cfg.eig_heatmap_examples must be a struct array with fields label and overrides (see run_eig_heatmap).');
    end
    eig_args = {'examples', ex};
end

% The opt-in block C stage (FD-010), checked here when the bundle selects
% it: its options, its E:I presets and their folders. A bundle that does
% not select it is not looked at (cfg.mc_gain_sweep may be absent).
gain_sweep = [];
if ismember('mc_gain_sweep', plan.stages)
    gain_sweep = mc_gain_sweep_plan(cfg);
end

%% PREFLIGHT: this session must have a pool allocation
% Hours of sweeps on an unknown number of workers is not something to find
% out afterwards. With automatic pool creation off, a session that never
% opened a pool runs every parfor serially; and several MATLAB sessions
% share this machine, so the size is for whoever divides it to decide, not
% this function. pool_set(N) records it; pool_set(0) is an allocation too
% (run serially, on purpose) and is accepted. Checked before any compute and
% before the run directory is created.
pc = pool_config();
if ~pc.is_set
    error('run_all_paper_analyses:NoPoolAllocation', ...
        ['No parallel-pool allocation is set for this MATLAB session, so the analyses are ' ...
         'refused rather than run serially for hours.\n%s'], pc.how_to_set);
end

% KEEP THE SWEEP PLOTS OFF SCREEN unless asked otherwise (cfg.visible_figures,
% default false). The sub-analyses plot and save at the end of each stage, and
% each new window raises itself and takes keyboard focus -- so an unattended run
% repeatedly interrupts whatever else is being done on the machine. The figures
% are still drawn and still saved; only the window is suppressed.
%
% Unlike the figure entry point, there is no per-figure `visible` argument to
% pass here: run_all_analyses' sub-analyses call psa.plot_sensitivity and friends
% directly. The graphics-root default is what reaches them.
%
% with_graphics_defaults returns an onCleanup, so the session's own setting comes
% back when this returns -- including if a stage throws past the wrapper. Held in
% a variable that stays in scope for the whole function; do not clear it.
visible_figures = isfield(cfg, 'visible_figures') && cfg.visible_figures;
if ~visible_figures
    fig_visibility = with_graphics_defaults('DefaultFigureVisible', 'off'); %#ok<NASGU>
end

vprintf(cfg.verbose, 'minimal', '[paper analyses] preset=%s run_mode=%s start=%s\n', cfg.preset_name, cfg.run_mode, datetime('now'));
t_all = tic;

%% 1. The sweep pipeline -- creates the run directory and its manifest
% cfg.run_dir is the SAME field make_all_paper_figures reads from: one
% directory, two stages. Empty (the paper's case) means "make your own dated
% folder", which is what run_all_analyses does with an empty output_dir. A
% relative path resolves against the project root, so a config can say
% 'data/fast_4' without caring what the cwd is -- same rule as fig_root.
% The transcript is saved with the data (TR, 2026-09-14):
% <run_dir>/command_window.log. Since 2026-09-15 it is written by vlog, the
% sink every vprintf / vfail line is appended to, instead of diary: diary
% only sees the command window, and a run launched through the MATLAB MCP
% server has its output captured before diary does, so those logs came out
% empty. With a named run_dir the log opens before the first stage; with an
% auto-named one it can only open once run_all_analyses has created the
% directory, so the sweep stage's lines are missing from that log. One log
% at a time, process-global; the onCleanup closes it however this exits.
%
% cfg.param_space_seed (optional; absent or empty = unset) seeds the joint
% parameter-space sample, so two runs of a bundle simulate the same grid
% points (FD-016). Unset, the sample is drawn from the clock as it always
% was. Checked here, before any compute: the PSA would reject a bad value
% too, but only after the two sensitivity stages had run.
param_space_seed = [];
if isfield(cfg, 'param_space_seed'); param_space_seed = cfg.param_space_seed; end
if ~isempty(param_space_seed) && ~(isnumeric(param_space_seed) && isscalar(param_space_seed) && ...
        isreal(param_space_seed) && isfinite(param_space_seed) && param_space_seed >= 0 && ...
        param_space_seed == floor(param_space_seed) && param_space_seed <= 2^32 - 1)
    error('run_all_paper_analyses:BadParamSpaceSeed', ...
        'cfg.param_space_seed must be empty or a nonnegative integer scalar no larger than 2^32 - 1.');
end
% The sweeps to run go to run_all_analyses by argument, with the stage
% selection for run_manifest. Both are passed ONLY when something is skipped,
% so a bundle that selects nothing makes the call, and the manifest, it
% always did. run_all_analyses is called even when the 'sweeps' stage is not
% selected (with no sweeps): it is what creates the run directory, the
% provenance record and the manifest the figures find the run by.
sweep_args = {};
if plan.is_selection
    sweep_args = {'sweeps', plan.sweeps, 'manifest_fields', ...
        struct('stages_run', {plan.stages}, 'stages_skipped', {plan.stages_skipped})};
end
if isempty(cfg.run_dir)
    print_selection(plan, cfg.verbose);
    run_dir = run_all_analyses(cfg.preset_name, cfg.run_mode, 'verbose', cfg.verbose, ...
        'param_space_seed', param_space_seed, sweep_args{:});
    log_guard = start_log(fullfile(run_dir, 'command_window.log'), cfg.verbose); %#ok<NASGU>
else
    out_dir = cfg.run_dir;
    if ~is_absolute_path(out_dir)
        out_dir = fullfile(fileparts(which('setup_paths')), out_dir);
    end
    assert_empty_target(out_dir);
    if ~isfolder(out_dir); mkdir(out_dir); end
    log_guard = start_log(fullfile(out_dir, 'command_window.log'), cfg.verbose); %#ok<NASGU>
    vprintf(cfg.verbose, 'minimal', '  output   : %s\n', out_dir);
    print_selection(plan, cfg.verbose);
    run_dir = run_all_analyses(cfg.preset_name, cfg.run_mode, 'output_dir', out_dir, 'verbose', cfg.verbose, ...
        'param_space_seed', param_space_seed, sweep_args{:});
end

results = struct('stage', {}, 'ok', {}, 'minutes', {}, 'detail', {}, 'err', {});
if plan.runs.sweeps
    results = record(results, 'sweeps', true, toc(t_all)/60, run_dir, '');
end

%% 2-5. The figure-specific compute
% Each writes into the run directory, so the whole paper's inputs sit together.
%
% dc_lle runs LAST because it is the longest of the five and the only one not
% feeding an in-paper figure: if an overnight run is going to be cut short, this
% is the stage to lose. It takes cfg.preset_name -- the paper's own network --
% and sweeps tonic DC across every adaptation condition.
% Memory capacity may run on its own integrator: cfg.mc_ode_solver ('' = let
% mc_run_config pick by the preset's noise). A bundle whose MC preset is the
% paper network with sigma_u_noise = 0 names 'sra1' here so MC uses the same
% integrator as every other stage (rk4 would be the noise-free default).
% The local-Lyapunov stage runs on its own preset when the bundle names one
% (cfg.local_lyapunov_preset: the stimulus-staircase twin of the paper preset).
ll_preset = cfg.preset_name;
if isfield(cfg, 'local_lyapunov_preset') && ~isempty(cfg.local_lyapunov_preset); ll_preset = cfg.local_lyapunov_preset; end
ll_args={};
ll_fields={'K','T','step_duration_s','accumulation_start_s','warmup_s','simulation_start_s'};
for j=1:numel(ll_fields)
    key=['local_lyapunov_' ll_fields{j}];
    if isfield(cfg,key), ll_args=[ll_args,{ll_fields{j},cfg.(key)}]; end %#ok<AGROW>
end
mc_solver = '';
if isfield(cfg, 'mc_ode_solver'); mc_solver = cfg.mc_ode_solver; end
% cfg.mc_lean (absent or false = off): both memory-capacity stages take the
% lean path, which holds no state trajectory (same numbers; see run_memory_capacity).
mc_lean = isfield(cfg, 'mc_lean') && ~isempty(cfg.mc_lean) && logical(cfg.mc_lean);
gain_horizon=0; gain_duration=0;
if isfield(cfg,'transient_gain_horizon_s'), gain_horizon=cfg.transient_gain_horizon_s; end
if isfield(cfg,'transient_gain_duration_s'), gain_duration=cfg.transient_gain_duration_s; end
stages = { ...
    'memory_capacity', @() run_memory_capacity( ...
        'preset_name', cfg.mc_preset, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'ode_solver', mc_solver, 'lean', mc_lean, ...
        'output_dir', run_dir); ...
    'mc_example',      @() run_memory_capacity_example( ...
        'preset_name', cfg.mc_preset, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'ode_solver', mc_solver, 'lean', mc_lean, ...
        'output_dir', fullfile(run_dir, 'mc_example')); ...
    'eig_heatmap',     @() run_eig_heatmap( ...
        'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'out_dir', fullfile(run_dir, 'eig_heatmap'), eig_args{:}); ...
    'numerics',        @() run_numerics_verification( ...
        'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'out_dir', fullfile(run_dir, 'numerics_verification')); ...
    'lyapunov_spectrum', @() run_lyapunov_spectrum( ...
        'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'out_dir', fullfile(run_dir, 'lyapunov_spectrum')); ...
    'transient_gain', @() run_transient_gain( ...
        'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'out_dir', fullfile(run_dir, 'transient_gain'), ...
        'horizon_s',gain_horizon,'duration_s',gain_duration); ...
    'local_lyapunov', @() run_local_lyapunov( ...
        'preset_name', ll_preset, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'out_dir', fullfile(run_dir, 'local_lyapunov'),ll_args{:}); ...
    'dc_lle',          @() run_dc_lle_analysis( ...
        'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'verbose', cfg.verbose, ...
        'output_dir', run_dir) };

if isfield(cfg,'illustrations') && cfg.illustrations
    stages(end+1,:)={'paper_illustrations', @()run_paper_illustrations(cfg,run_dir)};
end
% Opt-in stages (paper_stage_names' third list) follow, in that list's order.
% They are always in the table and run only when the plan names them.
stages(end+1,:) = {'mc_gain_sweep', @() run_mc_gain_stage(cfg, run_dir, gain_sweep)};

% Keep the selected stages (FD-023). The table above must hold exactly the
% stages paper_stage_names lists after 'sweeps', then its opt-in stages, in
% that order, or a selection could name a stage that never runs.
% test_essential_mode checks that from this file's text, without running
% anything; the assert is the backstop, and it comes after the sweeps only
% because the table needs run_dir. An opt-in stage the plan does not name is
% in neither plan list, so it is compared only as a known name.
table_names = stages(:, 1)';
want_names  = [plan.stages(~strcmp(plan.stages, 'sweeps')), ...
               plan.stages_skipped(~strcmp(plan.stages_skipped, 'sweeps'))];
[all_names, ~, opt_names] = paper_stage_names();
known_names = [all_names, opt_names];
table_in    = table_names(ismember(table_names, want_names));
assert(isequal(table_in, known_names(ismember(known_names, want_names))) && ...
       isequal(table_names(ismember(table_names, opt_names)), opt_names) && ...
       all(ismember(table_names, known_names)), ...
    'run_all_paper_analyses:StageTableMismatch', ...
    'The stage table (%s) is not paper_stage_names'' list for this bundle.', strjoin(table_names, ', '));
stages = stages(ismember(table_names, plan.stages), :);

for k = 1:size(stages, 1)
    name = stages{k, 1};
    vprintf(cfg.verbose, 'minimal', '[%d/%d] %s\n', k, size(stages, 1), name);
    t0 = tic;
    lastwarn('');
    try
        pool_restart(cfg.verbose);
        detail = stages{k, 2}();
        results = record(results, name, true, toc(t0)/60, detail, '');
        vprintf(cfg.verbose, 'minimal', '  -> %s\n', detail);
    catch ME
        results = record(results, name, false, toc(t0)/60, '', ME.message);
        vfail('  FAILED: %s: %s\n', ME.identifier, ME.message);
        for s = 1:min(3, numel(ME.stack))
            vfail('    at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
        end
        vfail('  continuing with the remaining stages\n');
    end
    % MATLAB's own warnings bypass vprintf; the last one raised during the
    % stage goes to the log (not the screen, where it was already shown).
    [wmsg, wid] = lastwarn;
    if ~isempty(wmsg); vlog('append', sprintf('  [last warning] %s: %s\n', wid, wmsg)); end
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
vprintf(cfg.verbose, 'near-none', '[paper analyses] complete in %.2f h -> %s\n', toc(t_all)/3600, run_dir);
for k = 1:numel(results)
    r = results(k);
    if r.ok; tag = 'OK    '; else; tag = 'FAILED'; end
    vprintf(cfg.verbose, 'minimal', '  %s  %-18s %7.1f min\n', tag, r.stage, r.minutes);
    if ~r.ok; vfail('          %s\n', r.err); end
end
% What this run directory does NOT hold, by the bundle's own selection (also
% in run_manifest: stages_skipped, sweeps_skipped).
for k = 1:numel(plan.stages_skipped)
    vprintf(cfg.verbose, 'minimal', '  SKIPPED %-18s (not selected by cfg.stages)\n', plan.stages_skipped{k});
end
if plan.runs.sweeps && ~isempty(plan.sweeps_skipped)
    vprintf(cfg.verbose, 'minimal', '  SKIPPED within sweeps: %s (not selected by cfg.sweeps)\n', ...
        strjoin(plan.sweeps_skipped, ', '));
end
vprintf(cfg.verbose, 'verbose', 'Next: make_all_paper_figures(paper_config(''run_dir'', run_dir))\n');
end

%% ------------------------------------------------------------------------
function assert_empty_target(out_dir)
% A named run directory must be absent or empty. This is an ERROR, not a
% warning, because a reused directory ACCUMULATES rather than overwriting.
%
% Each sub-analysis appends its own timestamped folder -- ParamSpaceAnalysis2
% builds <folder_prefix>_<note>_nLevs_<N>_<dt> and descends into it -- so a
% second run adds a parallel set beside the first. What DOES get replaced is
% the top level: run_manifest, git_provenance and parameters.md then describe
% the newest run, in a directory holding several. Figures glob param_space_*,
% 1D_sensitivity_* and tau_sensitivity_* inside the run directory, and
% fig_sensitivity_medians matches those folders by PARAMETER NAME -- so it
% would find two per parameter and silently pick one.
%
% Refusing is cheap; the alternative is discovering it in a figure hours
% later. This does NOT delete anything: what to do with the old run is the
% user's call, and one of the options is "keep it".
if ~isfolder(out_dir); return; end                  % absent is fine; it gets created
listing = dir(out_dir);
listing = listing(~ismember({listing.name}, {'.', '..'}));
if isempty(listing); return; end                    % exists but empty is fine

names = {listing.name};
shown = strjoin(names(1:min(6, numel(names))), ', ');
if numel(names) > 6
    shown = sprintf('%s, ... (%d entries)', shown, numel(names));
end
error('run_all_paper_analyses:TargetNotEmpty', ...
    ['The run directory named by cfg.run_dir is not empty:\n' ...
     '  %s\n' ...
     '  contains: %s\n\n' ...
     'Analyses are refused rather than run, because a reused directory\n' ...
     'ACCUMULATES: each sweep appends its own timestamped subfolder while\n' ...
     'run_manifest and provenance at the top level are replaced, leaving a\n' ...
     'manifest that describes one run in a directory holding several.\n\n' ...
     'Delete or move it, or point cfg.run_dir somewhere else. Leave\n' ...
     'cfg.run_dir empty to get a fresh dated directory automatically.'], ...
    out_dir, shown);
end

function tf = is_absolute_path(p)
% Absolute on either platform: a leading separator, or a Windows drive letter.
% Deliberately NOT isfolder(): whether a relative path happens to exist under
% the current directory must not decide where output lands. Same helper, and
% the same reasoning, as in make_all_paper_figures.
tf = startsWith(p, filesep) || startsWith(p, '/') || ...
     ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end

function print_selection(plan, verbose)
% One line each for what a selecting bundle runs and skips; nothing for a
% bundle that runs everything.
if ~plan.is_selection; return; end
vprintf(verbose, 'minimal', '  stages   : %s\n', list_or_none(plan.stages));
vprintf(verbose, 'minimal', '  skipped  : %s (not selected by cfg.stages)\n', list_or_none(plan.stages_skipped));
end

function G = mc_gain_sweep_plan(cfg)
% The block C stage's arguments from the bundle config, validated before any
% compute. G.args: name-value pairs for run_mc_gain_sweep shared by the main
% and the E:I runs; G.ei: struct array (preset, bias, folder, seeds_args).
allowed = {'gain_levels', 'bias_levels', 'seeds', 'conditions', 'model_args', 'protocol', ...
    'lya_K', 'lya_T_sec', 'lya_acc_sec', 'segment_len', 'ei_presets', 'ei_seeds'};
if ~isfield(cfg, 'mc_preset') || isempty(cfg.mc_preset)
    error('run_all_paper_analyses:BadMcGainSweep', ...
        'The mc_gain_sweep stage runs on cfg.mc_preset, which this bundle does not set.');
end
gs = struct();
if isfield(cfg, 'mc_gain_sweep') && ~isempty(cfg.mc_gain_sweep); gs = cfg.mc_gain_sweep; end
if ~(isstruct(gs) && isscalar(gs))
    error('run_all_paper_analyses:BadMcGainSweep', 'cfg.mc_gain_sweep must be a scalar struct.');
end
unknown = setdiff(fieldnames(gs), allowed);
if ~isempty(unknown)
    error('run_all_paper_analyses:BadMcGainSweep', ...
        'cfg.mc_gain_sweep has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(allowed, ', '));
end
G = struct('args', {{}}, 'ei', struct('preset', {}, 'bias', {}, 'folder', {}, 'seeds_args', {}));
pass = setdiff(allowed, {'bias_levels', 'ei_presets', 'ei_seeds'}, 'stable');
for k = 1:numel(pass)
    if isfield(gs, pass{k}) && ~isempty(gs.(pass{k}))
        G.args = [G.args, {pass{k}, gs.(pass{k})}];
    end
end
G.bias_args = {};
if isfield(gs, 'bias_levels') && ~isempty(gs.bias_levels)
    G.bias_args = {'bias_levels', gs.bias_levels};
end
if ~isfield(gs, 'ei_presets') || isempty(gs.ei_presets); return; end
ei = cellstr(gs.ei_presets);
ei_seeds = {};
if isfield(gs, 'ei_seeds') && ~isempty(gs.ei_seeds); ei_seeds = {'seeds', gs.ei_seeds}; end
model_args = {};
if isfield(gs, 'model_args'); model_args = gs.model_args; end
for k = 1:numel(ei)
    d = srnn_param_preset(ei{k});                  % errors on an unknown name
    % Gain only, at the network's own bias: the preset's S_c, or a model_args
    % override of it (last one wins, as in the constructor).
    if ~isfield(d, 'S_c')
        error('run_all_paper_analyses:BadMcGainSweep', ...
            'E:I preset ''%s'' sets no S_c, so its own bias is not known.', ei{k});
    end
    b = d.S_c;
    i_sc = find(strcmp(model_args(1:2:end), 'S_c'), 1, 'last');
    if ~isempty(i_sc); b = model_args{2 * i_sc}; end
    parts = strsplit(ei{k}, '_');
    G.ei(k) = struct('preset', ei{k}, 'bias', b, 'folder', ['mc_gain_sweep_' parts{end}], ...
        'seeds_args', {ei_seeds});
end
folders = {G.ei.folder};
if numel(unique(folders)) < numel(folders)
    error('run_all_paper_analyses:BadMcGainSweep', ...
        'cfg.mc_gain_sweep.ei_presets end in the same suffix, so their folders would clash: %s.', ...
        strjoin(folders, ', '));
end
end

function detail = run_mc_gain_stage(cfg, run_dir, G)
% Block C: the gain-by-bias sweep into <run_dir>/mc_gain_sweep, then each E:I
% gain-only run into <run_dir>/mc_gain_sweep_<suffix>/mc_gain_sweep. Each run
% is tried in turn, so one failing does not cost the others; the stage then
% fails with every error, or returns the result folders.
common = {'run_mode', cfg.run_mode, 'verbose', cfg.verbose};
if isfield(cfg, 'mc_ode_solver') && ~isempty(cfg.mc_ode_solver)
    common = [common, {'ode_solver', cfg.mc_ode_solver}];
end
runs = struct('label', 'main', 'args', {[{'preset_name', cfg.mc_preset}, G.args, G.bias_args, ...
    {'output_dir', run_dir}]});
for k = 1:numel(G.ei)
    e = G.ei(k);
    runs(end+1) = struct('label', e.folder, 'args', {[{'preset_name', e.preset}, G.args, ...
        e.seeds_args, {'bias_levels', e.bias, 'output_dir', fullfile(run_dir, e.folder)}]}); %#ok<AGROW>
end
dirs = {}; errs = {};
for k = 1:numel(runs)
    try
        dirs{end+1} = run_mc_gain_sweep(common{:}, runs(k).args{:}); %#ok<AGROW>
    catch ME
        errs{end+1} = sprintf('%s: %s: %s', runs(k).label, ME.identifier, ME.message); %#ok<AGROW>
        vfail('  mc_gain_sweep run %s FAILED: %s\n', runs(k).label, ME.message);
    end
end
if ~isempty(errs)
    error('run_all_paper_analyses:McGainSweepFailed', '%d of %d run(s) failed:\n%s', ...
        numel(errs), numel(runs), strjoin(errs, '\n'));
end
detail = strjoin(dirs, '; ');
end

function s = list_or_none(names)
if isempty(names); s = 'none'; else; s = strjoin(names, ', '); end
end

function results = record(results, stage, ok, minutes, detail, err)
if ~ischar(detail); detail = ''; end
results(end+1) = struct('stage', stage, 'ok', ok, 'minutes', minutes, ...
    'detail', detail, 'err', err);
end

function guard = start_log(log_file, verbose)
% Append every vprintf / vfail line to log_file until the guard is destroyed
% (vlog is process-global and keeps one log; a caller's own is replaced for
% the duration and closed afterwards).
guard = vlog('open', log_file);
vprintf(verbose, 'minimal', '[log] transcript -> %s\n', log_file);
end
