function cfg = analysis_run_config(analysis, run_mode, preset_defaults)
% ANALYSIS_RUN_CONFIG Sweep size and timing settings for one analysis + run mode.
%
% Replaces the near-identical `switch run_mode` blocks that were duplicated
% across the run_all_analyses sub-scripts. run_mode controls COST/FIDELITY;
% which network is simulated is a separate axis, see srnn_param_preset.
%
% Usage:
%   cfg = analysis_run_config('sensitivity', run_mode, preset_defaults);
%   psa = ParamSpaceAnalysis2('n_levels', cfg.n_levels, ...);
%   psa.model_defaults = merge_struct(preset_defaults, cfg.model);
%   psa.add_grid_parameter('reps', 1:cfg.n_reps);
%
% Returns:
%   cfg.n_levels - levels per grid parameter
%   cfg.n_reps   - repetitions (absent for 'param_space', which has no reps axis)
%   cfg.model    - struct of SRNNModel2 settings for this mode. lya_T_interval
%                  is simply ABSENT where the class default should apply, which
%                  replaces the old `if ~isempty(lya_T_interval_mode)` guard.
%   cfg.sde_solver    - the stochastic scheme this cell would use
%   cfg.is_stochastic - whether preset_defaults asked for noise
%
% Modes:
%   'fast'       - few levels/reps, fs=200, short T_range; finishes quickly
%   'fast2'      - 'fast' with twice the reps on the 1-D sensitivity sweeps and
%                  twice the T_range on those and on param_space. The reps buy
%                  the spread WITHIN a level, which at 3 reps is too thin to
%                  read; the longer window buys an LLE that is less of a
%                  transient. tau_sensitivity is untouched, so it is still
%                  exactly 'fast'. Cheaper than the extra LEVELS that separate
%                  'fast' from 'medium'.
%   'medium'     - roughly halfway; for a usable but not publication-size run
%   'medium2'    - between medium and production, at fs=800, for overnight
%                  STOCHASTIC runs. See the note below on why its sweep sizes
%                  sit nearer medium than the true midpoint.
%   'production' - full-size sweeps at fs=400 (for real results)
%
% ON 'medium2' AND fs=800. Doubling fs is free relative to medium *per
% simulated second*: SRA1 takes two drift evaluations per step against rk4's
% four, so sra1 at fs=800 and rk4 at fs=400 both cost 1600 evaluations per
% second of simulation. The finer step buys back some of the accuracy that a
% fixed-step scheme gives up against ode45, which matters because a stochastic
% run cannot use an adaptive solver at all.
%
% Its levels and reps are NOT the medium/production midpoint. Measured against
% a completed fast2 run (19.7 min), medium alone is already ~9 h and production
% ~80 h, so the true midpoint is a multi-day job. medium2 is sized to finish
% overnight (~9 h) and therefore sits nearer medium: more levels than medium on
% the 1-D sweeps, fewer reps, a longer window, and the finer step.
%
% DETERMINISTIC AND STOCHASTIC SOLVERS. Every cell names two integrators, and
% which one is used depends on the PRESET rather than on the run mode:
%
%   deterministic   rk4 for fast / fast2 / medium, ode45 for production
%   stochastic      sra1 everywhere
%
% This exists because the two knobs stopped being orthogonal. sigma_u_noise is
% physics and so belongs to a preset, but sigma > 0 REQUIRES a stochastic
% integrator -- so a preset now constrains a value in this table. Selecting here
% is what keeps merge_struct's precedence intact: cfg.model still wins over the
% preset, it just already holds the right solver. A preset without noise picks
% the deterministic column and is bit-identical to before this existed.
%
% The stochastic column is uniformly 'sra1': it costs the same two drift
% evaluations per step as 'heun' but is strong order 1.5 rather than 1.0, and is
% CHEAPER than 'rk4' (four evaluations). pack() defaults it, so the twelve cells
% do not repeat it and cannot drift apart, while any one cell can still override.
%
% Note this loses the adaptive-tolerance fidelity of 'ode45' at 'production':
% there is no adaptive SDE solver, because step-size control is meaningless once
% the noise increments are tied to the step. A stochastic production run is
% therefore fixed-step and does NOT mean the same accuracy as its deterministic
% twin. That is accepted deliberately.
%
% See also: srnn_param_preset, merge_struct, ParamSpaceAnalysis2, sde_fixed_step

arguments
    analysis (1,:) char
    run_mode (1,:) char
    preset_defaults struct = struct()
end

valid_modes = {'fast', 'fast2', 'medium', 'medium2', 'production'};
if ~ismember(run_mode, valid_modes)
    error('analysis_run_config:badMode', ...
        'Unknown run_mode ''%s'' (expected %s).', ...
        run_mode, strjoin(cellfun(@(s) ['''' s ''''], valid_modes, 'UniformOutput', false), ', '));
end

switch analysis
    case 'sensitivity'
        % 1-D sweeps, one PSA per swept parameter, with a reps axis.
        switch run_mode
            case 'fast'
                % fs=200 keeps Benettin's lya_dt/dt guard satisfied (4>=3).
                cfg = pack(4,  3,  'rk4',   200, [0, 10], [5, 10]);
            case 'fast2'
                % fast2 differs from fast here in two ways: 6 reps rather than
                % 3, and twice the T_range. The LLE window is doubled with it,
                % so it stays the same FRACTION of the run (the second half) --
                % a longer simulation with the old [5,10] window would have
                % measured the same transient from further away.
                cfg = pack(4,  6,  'rk4',   200, [0, 20], [10, 20]);
            case 'medium'
                cfg = pack(11, 15, 'rk4',   400, [0, 20], [10, 20]);
            case 'medium2'
                % 13 levels (up from medium's 11) but 12 reps (down from 15):
                % the extra resolution along the axis is worth more than the
                % extra samples per level once the window is longer. T doubled
                % over fast2 and fs doubled over medium. ~6.5 h of the budget.
                cfg = pack(13, 12, 'rk4',   800, [0, 25], [12.5, 25]);
            case 'production'
                cfg = pack(25, 50, 'ode45', 400, [0, 50], [20, 50]);
        end

    case 'tau_sensitivity'
        % Vector-parameter sweep over tau_a_E. The swept tau reaches 60 s, so
        % medium/production use the longer T_range and the class-default LLE
        % window rather than an explicit short one.
        switch run_mode
            case {'fast', 'fast2'}
                % NOTE: this window is far shorter than the swept tau, so long-tau
                % LLE reflects a transient -- accepted for fast qualitative runs.
                % fast2 doubles reps on the 1-D sweeps above only, not here.
                cfg = pack(7,  7,  'rk4',   200, [0, 20], [10, 20]);
            case 'medium'
                % rk4 rather than ode45, so the deterministic column reads the
                % same everywhere: rk4 up to medium, ode45 only at production.
                % Safe here despite the swept tau reaching 60 s -- stiffness
                % comes from FAST modes, and the fastest is still tau_d = 0.1
                % against dt = 1/400. Also faster than ode45 with MaxStep = dt.
                cfg = pack(11, 25, 'rk4',   400, [0, 50], []);
            case 'medium2'
                % Keeps the full [0, 50] window and the class-default LLE
                % interval, as medium/production do -- the swept tau reaches
                % 60 s, so shortening the run here is what makes the long-tau
                % end meaningless. ~2.3 h of the budget.
                cfg = pack(13, 15, 'rk4',   800, [0, 50], []);
            case 'production'
                cfg = pack(25, 50, 'ode45', 400, [0, 50], []);
        end

    case 'param_space'
        % Multi-dimensional grid with no reps axis, so n_levels^n_params configs.
        switch run_mode
            case 'fast'
                cfg = pack(3, [], 'rk4',   200, [0, 20], [10, 20]);
            case 'fast2'
                % No reps axis here, so the only fast2 difference is the doubled
                % T_range (and the LLE window doubled with it).
                cfg = pack(3, [], 'rk4',   200, [0, 40], [20, 40]);
            case 'medium'
                cfg = pack(4, [], 'rk4',   400, [0, 20], [10, 20]);
            case 'medium2'
                % 5 levels, i.e. production's, because this is the cheapest of
                % the three stages: the grid is n_levels^3 but there is no reps
                % axis. ~40 min of the budget, and the n = 1000 corner is what
                % dominates it.
                cfg = pack(5, [], 'rk4',   800, [0, 25], [12.5, 25]);
            case 'production'
                cfg = pack(5, [], 'ode45', 400, [0, 50], []);
        end

    otherwise
        error('analysis_run_config:badAnalysis', ...
            ['Unknown analysis ''%s'' (expected ''sensitivity'', ' ...
            '''tau_sensitivity'', or ''param_space'').'], analysis);
end

cfg.analysis = analysis;
cfg.run_mode = run_mode;

% Swap in the stochastic scheme if the preset asked for noise. Done here, after
% the table, so the rule lives in ONE place rather than in each sub-script.
cfg.is_stochastic = preset_is_stochastic(preset_defaults);
if cfg.is_stochastic
    cfg.model.ode_solver = cfg.sde_solver;
end
end

function tf = preset_is_stochastic(p)
% The single definition of "this run is stochastic".
%
% Reads the PRESET. A sigma_u_noise made into a grid AXIS instead would not be
% seen here and would get the deterministic solver -- that fails safe, since
% ParamSpaceAnalysis2.validate_noise_settings rejects the pairing before the
% output directory is created, but it is not handled.
tf = isstruct(p) && isfield(p, 'sigma_u_noise') && any(p.sigma_u_noise(:) > 0);
end

function cfg = pack(n_levels, n_reps, ode_solver, fs, T_range, lya_T_interval, sde_solver)
% Assemble one row of the table. An empty lya_T_interval means "leave it to the
% class default", so the field is omitted rather than set empty.
%
% sde_solver defaults to 'sra1' so the twelve cells do not each repeat it. It is
% kept OUT of cfg.model: only cfg.model is merged into model_defaults, and
% sde_solver is not an SRNNModel2 property, so validate_model_defaults would
% reject it there.
if nargin < 7 || isempty(sde_solver)
    sde_solver = 'sra1';
end
cfg = struct();
cfg.n_levels = n_levels;
if ~isempty(n_reps)
    cfg.n_reps = n_reps;
end
cfg.model = struct('ode_solver', ode_solver, 'fs', fs, 'T_range', T_range);
if ~isempty(lya_T_interval)
    cfg.model.lya_T_interval = lya_T_interval;
end
cfg.sde_solver = sde_solver;
end
