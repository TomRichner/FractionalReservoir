function cfg = analysis_run_config(analysis, run_mode)
% ANALYSIS_RUN_CONFIG Sweep size and timing settings for one analysis + run mode.
%
% Replaces the near-identical `switch run_mode` blocks that were duplicated
% across the run_all_analyses sub-scripts. run_mode controls COST/FIDELITY;
% which network is simulated is a separate axis, see srnn_param_preset.
%
% Usage:
%   cfg = analysis_run_config('sensitivity', run_mode);
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
%   'production' - full-size sweeps at fs=400 (for real results)
%
% See also: srnn_param_preset, merge_struct, ParamSpaceAnalysis2

arguments
    analysis (1,:) char
    run_mode (1,:) char
end

valid_modes = {'fast', 'fast2', 'medium', 'production'};
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
                cfg = pack(4,  3,  @ode_rk4, 200, [0, 10], [5, 10]);
            case 'fast2'
                % fast2 differs from fast here in two ways: 6 reps rather than
                % 3, and twice the T_range. The LLE window is doubled with it,
                % so it stays the same FRACTION of the run (the second half) --
                % a longer simulation with the old [5,10] window would have
                % measured the same transient from further away.
                cfg = pack(4,  6,  @ode_rk4, 200, [0, 20], [10, 20]);
            case 'medium'
                cfg = pack(11, 15, @ode_rk4, 400, [0, 20], [10, 20]);
            case 'production'
                cfg = pack(25, 50, @ode45,   400, [0, 50], [20, 50]);
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
                cfg = pack(7,  7,  @ode_rk4, 200, [0, 20], [10, 20]);
            case 'medium'
                cfg = pack(11, 25, @ode45,   400, [0, 50], []);
            case 'production'
                cfg = pack(25, 50, @ode45,   400, [0, 50], []);
        end

    case 'param_space'
        % Multi-dimensional grid with no reps axis, so n_levels^n_params configs.
        switch run_mode
            case 'fast'
                cfg = pack(3, [], @ode_rk4, 200, [0, 20], [10, 20]);
            case 'fast2'
                % No reps axis here, so the only fast2 difference is the doubled
                % T_range (and the LLE window doubled with it).
                cfg = pack(3, [], @ode_rk4, 200, [0, 40], [20, 40]);
            case 'medium'
                cfg = pack(4, [], @ode_rk4, 400, [0, 20], [10, 20]);
            case 'production'
                cfg = pack(5, [], @ode45,   400, [0, 50], []);
        end

    otherwise
        error('analysis_run_config:badAnalysis', ...
            ['Unknown analysis ''%s'' (expected ''sensitivity'', ' ...
            '''tau_sensitivity'', or ''param_space'').'], analysis);
end

cfg.analysis = analysis;
cfg.run_mode = run_mode;
end

function cfg = pack(n_levels, n_reps, ode_solver, fs, T_range, lya_T_interval)
% Assemble one row of the table. An empty lya_T_interval means "leave it to the
% class default", so the field is omitted rather than set empty.
cfg = struct();
cfg.n_levels = n_levels;
if ~isempty(n_reps)
    cfg.n_reps = n_reps;
end
cfg.model = struct('ode_solver', ode_solver, 'fs', fs, 'T_range', T_range);
if ~isempty(lya_T_interval)
    cfg.model.lya_T_interval = lya_T_interval;
end
end
