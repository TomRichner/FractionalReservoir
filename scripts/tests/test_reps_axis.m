% TEST_REPS_AXIS The repetitions axis of the 1-D sweeps has exactly n_reps levels.
%
% Plan C, F2 (2026-11-12). run_sensitivity_analysis and
% run_tau_sensitivity_analysis added the reps axis with
% add_grid_parameter('reps', 1:n_reps). add_grid_parameter reads a
% two-element value as a [min max] range expanded to n_levels values, so with
% n_reps = 2 (the 'essential' mode, 4 levels x 2 reps) the axis was
% 1, 1.33, 1.67, 2: four repetitions, 16 jobs per condition where 8 were
% meant. n_reps = 1 errored; three or more were an explicit list and right.
% The drivers now call add_grid_levels, which always takes the list as given.
%
% Compute-free: grids are built with generate_grid, nothing is simulated.
% Checks that
%   * both drivers add reps through add_grid_levels (from their text);
%   * add_grid_levels gives exactly n_reps integer levels for n_reps = 1, 2,
%     3, 5, beside a range axis (sensitivity) and a vector axis (tau);
%   * at 'essential' a 1-D sweep's job table is 4 x 2 = 8 jobs per condition,
%     each (level, rep) pair once;
%   * in every other mode of run_mode_names, both drivers' job tables (configs
%     and network seeds) are identical to what add_grid_parameter gave, which
%     add_grid_parameter itself, unchanged, still shows;
%   * add_grid_parameter still reads a two-element value as a range.
%
% Asserts; silent until it throws, then prints one PASS line.
% Assumes setup_paths has run.
% See also: ParamSpaceAnalysis2.add_grid_levels, run_sensitivity_analysis,
%           run_tau_sensitivity_analysis, test_essential_mode

%% Both drivers add reps through add_grid_levels
for f = {'run_sensitivity_analysis', 'run_tau_sensitivity_analysis'}
    code = regexprep(fileread(which(f{1})), '%[^\n]*', '');   % comments off
    assert(isscalar(regexp(code, 'add_grid_levels\(''reps'',\s*1:ctx\.n_reps\)')), ...
        '%s must add reps with add_grid_levels(''reps'', 1:ctx.n_reps)', f{1});
    assert(isempty(regexp(code, 'add_grid_parameter\(''reps''', 'once')), ...
        '%s still adds reps with add_grid_parameter', f{1});
end

%% n_reps = 1, 2, 3, 5 give exactly that many integer levels
for driver = {'sensitivity', 'tau'}
    for n_reps = [1 2 3 5]
        for n_levels = [2 4 7]
            psa = sweep_psa(driver{1}, n_levels, n_reps, 'levels');
            reps = psa.param_vectors{strcmp(psa.grid_params, 'reps')};
            assert(isequal(reps, 1:n_reps), '%s, n_reps = %d, n_levels = %d: reps are %s', ...
                driver{1}, n_reps, n_levels, mat2str(reps));
            assert(psa.num_combinations == n_levels * n_reps);
        end
    end
end

%% 'essential': 4 x 2 = 8 jobs per condition
for driver = {'sensitivity', 'tau'}
    an = analysis_name(driver{1});
    cfg = analysis_run_config(an, 'essential');
    psa = sweep_psa(driver{1}, cfg.n_levels, cfg.n_reps, 'levels');
    jobs = job_table(psa);
    if strcmp(driver{1}, 'sensitivity')
        assert(cfg.n_levels == 4 && cfg.n_reps == 2);
        assert(size(jobs, 1) == 8, '%s at essential: %d jobs per condition', driver{1}, size(jobs, 1));
    else
        assert(size(jobs, 1) == cfg.n_levels * cfg.n_reps);
    end
    assert(size(unique(jobs(:, 1:2), 'rows'), 1) == size(jobs, 1), 'a (level, rep) pair repeats');
    assert(all(ismember(jobs(:, 2), 1:cfg.n_reps)), 'a rep that is not an integer in 1:n_reps');
    assert(numel(unique(jobs(:, 3))) == size(jobs, 1), 'network seeds must be distinct');
end

%% Every other mode: the same jobs and seeds as before the fix
modes = run_mode_names();
for driver = {'sensitivity', 'tau'}
    an = analysis_name(driver{1});
    for m = 1:numel(modes)
        cfg = analysis_run_config(an, modes{m});
        after = job_table(sweep_psa(driver{1}, cfg.n_levels, cfg.n_reps, 'levels'));
        before = job_table(sweep_psa(driver{1}, cfg.n_levels, cfg.n_reps, 'parameter'));
        if cfg.n_reps == 2
            % Only a two-element list was misread: n_levels repetitions.
            assert(size(before, 1) == cfg.n_levels^2 && size(after, 1) == 2 * cfg.n_levels);
        else
            assert(cfg.n_reps >= 3, '%s %s: n_reps = %d', an, modes{m}, cfg.n_reps);
            assert(isequal(after, before), '%s at %s: the job table changed', an, modes{m});
        end
        % Today only the 1-D sensitivity sweep at 'essential' has two reps
        % (the tau sweep at 'essential' runs fast's 7 x 7, test_essential_mode).
        assert((strcmp(modes{m}, 'essential') && strcmp(an, 'sensitivity')) == (cfg.n_reps == 2), ...
            '%s %s has n_reps = %d', an, modes{m}, cfg.n_reps);
    end
end

%% add_grid_parameter itself is unchanged
psa = ParamSpaceAnalysis2('n_levels', 4, 'randomize_order', false, 'verbose', 'near-none');
psa.add_grid_parameter('level_of_chaos', [1 2]);
psa.generate_grid();
assert(isequal(psa.param_vectors{1}, linspace(1, 2, 4)), 'a two-element value is still a range');
try
    psa.add_grid_parameter('reps', 1);
    error('test_reps_axis:NoThrow', 'add_grid_parameter accepted one element');
catch err
    assert(strcmp(err.identifier, 'ParamSpaceAnalysis2:InvalidInput'));
end
% add_grid_levels replaces a range set earlier under the same name, and back.
psa.add_grid_levels('level_of_chaos', [1 2]);
psa.generate_grid();
assert(isequal(psa.param_vectors{1}, [1 2]) && ~isfield(psa.param_ranges, 'level_of_chaos'));
psa.add_grid_parameter('level_of_chaos', [1 2]);
psa.generate_grid();
assert(isequal(psa.param_vectors{1}, linspace(1, 2, 4)) && ~isfield(psa.explicit_vectors, 'level_of_chaos'));

disp('PASS: test_reps_axis -- reps are exactly 1:n_reps in both drivers; only essential''s job table moved.');

%% ------------------------------------------------------------------------
function an = analysis_name(driver)
if strcmp(driver, 'sensitivity'); an = 'sensitivity'; else; an = 'tau_sensitivity'; end
end

function psa = sweep_psa(driver, n_levels, n_reps, how)
% The grid the driver builds: its swept axis, then reps, in its order and
% with randomize_order = false as both drivers set it. 'levels' is the
% drivers' call now; 'parameter' the call before the fix.
psa = ParamSpaceAnalysis2('n_levels', n_levels, 'randomize_order', false, 'verbose', 'near-none');
if strcmp(driver, 'sensitivity')
    psa.add_grid_parameter('level_of_chaos', [0.5 1.5]);   % any [min max] axis
else
    psa.add_vector_parameter('tau_a_E', 'vary_element', 'last', 'fixed_value', 0.25, ...
        'vary_range', [1, 30], 'n_elements', 3, 'spacing', 'log', 'level_spacing', 'linear');
end
if strcmp(how, 'levels')
    psa.add_grid_levels('reps', 1:n_reps);
else
    psa.add_grid_parameter('reps', 1:n_reps);
end
psa.generate_grid();
end

function jobs = job_table(psa)
% One row per job of one condition, in run order: swept level, rep, network
% seed (config_idx*100 + offset, as ParamSpaceAnalysis2.run assigns it).
idx = psa.shuffled_indices(:);
jobs = zeros(numel(idx), 3);
names = psa.grid_params;
swept = names{~strcmp(names, 'reps')};
for k = 1:numel(idx)
    c = psa.all_configs{idx(k)};
    jobs(k, :) = [c.(swept), c.reps, idx(k) * 100 + psa.network_seed_offset];
end
end
