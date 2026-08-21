% example_SRNNCellTypePairs_dualStd.m
% One SRNNCellTypePairs simulation of the DUAL-TIMESCALE STD preset, at the
% sfa_and_std operating point, with the full plotting set.
%
% The subject of this file is n_b = 2: every route depresses on two timescales,
% tau_rec = [2 4] with tau_rel = [0.25 0.5], where every other preset in the
% pipeline carries a single one. compile_synapse_config reads n_b off
% numel(tau_rec), the two b-states are integrated as columns, and the synapse
% sees their PRODUCT -- so the panel to watch is "STD (product)".
%
% THIS IS THE SCRATCH COPY -- override freely. Its sibling
% example_SRNNCellTypePairs_dualStd.m is the one that deliberately modifies
% nothing after the preset, so that it always runs exactly what the sweeps run;
% keep that property there and do the poking here.
%
% Two knobs, two different places:
%
%   the CONDITION      condition_name, in "What to simulate" below
%                      (no_adaptation | sfa_only | std_only | sfa_and_std)
%   everything else    the "Playground overrides" section, after the model is
%                      constructed -- there is a commented example there
%
% They are separate because a condition is two properties at once, n_a AND
% synapse_config. Setting one of them by hand in the overrides section would
% leave you running something other than the regime you named.
%
% Read-only Dependent properties error rather than take a value (alpha, R,
% mu_tilde, sigma_x_raw, activation_function); for connectivity use the
% _relative names.
%
% The operating point is still assembled from THREE sources, in the same
% precedence ParamSpaceAnalysis2.run_single_job uses:
%
%   srnn_param_preset            the physics -- n, connectivity, nonlinearity,
%                                noise, and the STD routes
%   analysis_run_config          the cost/fidelity knobs -- fs, T_range, solver
%   srnn_adaptation_conditions   the adaptation regime -- n_a, synapse_config
%
% See also: example_SRNNCellTypePairs_from_preset, srnn_param_preset,
%           analysis_run_config, srnn_adaptation_conditions

close all; clear; clc;

%% What to simulate
preset_name    = 'celltype_pairs_Sc0p2_noise0p025_dualStd';
condition_name = 'sfa_and_std';   % no_adaptation | sfa_only | std_only | sfa_and_std
run_mode       = 'medium';        % matches the sensitivity sweeps: fs 400, T [0 20]

% A fixed seed, because this is one illustrative network rather than a sample.
% The sweeps vary it per job, so any single run here is one draw from what they
% average.
rng_seeds = [19 20];

% Eigenvalue snapshots, in seconds. Must fall inside T_range.
eig_times = [2 10 18];

%% Assemble the operating point
[preset_defaults, model_class, conditions] = srnn_param_preset(preset_name);
cfg = analysis_run_config('sensitivity', run_mode, preset_defaults);

cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond = conditions{strcmp(cond_names, condition_name)};

% Preset first, run_mode timings second (so the integrator and fs come from the
% run mode), condition last (so n_a and synapse_config win).
model_args = namedargs2cell(merge_struct(preset_defaults, cfg.model));
model_args = [model_args, {'rng_seeds', rng_seeds}];
for f = setdiff(fieldnames(cond), {'name'})'
    model_args = [model_args, {f{1}, cond.(f{1})}]; %#ok<AGROW>
end

model = feval(model_class, model_args{:});

% Storage rather than physics: plot_eigenvalues reads the Jacobian off the
% stored trajectory, so the full state has to be kept. At n=500, fs=400 over
% [0 20] with 3250 state variables that is 8001 x 3250 doubles, about 208 MB --
% half again as much as the n_b = 1 presets, because the second timescale
% doubles the b-states. Drop it (and the plot_eigenvalues call) if you raise n
% or fs.
model.store_full_state = true;

%% Playground overrides
% Uncomment and edit. These land between the constructor and build(), which is
% the window where they still take effect -- level_of_chaos, for one, scales W
% inside build(). They are also above the report call, so whatever you set here
% gets PRINTED: what you read is what runs.
%
% To change the adaptation regime, edit condition_name at the top instead. A
% condition is n_a AND synapse_config together, and setting only one of them
% here would leave you running something you did not name.

model.level_of_chaos = 1.5;      % preset: 1.0. >1 pushes past the edge of chaos
model.S_c = 0.4;                % preset: 0.20
model.c = [0.5/3, 0];
% model.sigma_u_noise = 0.05;      % preset: 0.025 (solver is already sra1)
% model.n = 300;                   % preset: 500
% model.rng_seeds = [21 22];       % a different network draw
% model.T_range = [0 30];          % also move eig_times above to match

%% Report the configuration before spending time on it
% Pure printing -- it reads the model and the run config and writes to the
% command window, computes nothing the simulation uses. Kept as a local
% function so the interesting part of this script (what to simulate, and the
% playground overrides above) stays on one screen.
report_configuration(model, preset_name, model_class, cond.name, run_mode, cfg);

%% Build and run
model.build();
model.run();

%% Plots
% Compact summary: one panel per quantity. Where types or routes share an axes
% they are drawn back-to-front, so E sits on top of I.
model.plot();

% Per-cell-type view: one COLUMN per type, every neuron's trace. The
% "STD (product)" row is prod(b) across BOTH timescales and every route out of
% that type -- the quantity the synapse actually multiplies by.
model.plot_celltypes();

% Effective-Jacobian spectrum at three times. Watch whether the deeper
% dual-timescale depression pulls eigenvalues back inside the unit circle
% faster than the single-timescale presets do.
model.plot_eigenvalues(eig_times);

% W on a diverging colormap -- the SCALED W, level_of_chaos included.
model.plot_W();
model.plot_W_spectrum();
model.plot_weight_histogram();

%% Result
fprintf('\n=== Result ===\n');
fprintf('LLE                   : %.6f\n', model.lya_results.LLE);
fprintf('mean rate per type    :');
for q = 1:model.n_cellTypes
    fprintf('  %s = %.4f', model.cell_type_names{q}, ...
        mean(model.plot_data.r.(model.cell_type_names{q}), 'all'));
end
fprintf('\ndead-end states       : %d\n', model.dead_state_count);

%% ==================== Local functions ====================
function report_configuration(model, preset_name, model_class, cond_name, ...
    run_mode, cfg)
%REPORT_CONFIGURATION Print the assembled operating point. No side effects.
%
% Everything here is read off the model AFTER construction, so any playground
% override placed above the call shows up in what is printed -- which is the
% point: what you read is what will run. Move an override below this call and
% the report quietly describes a network you did not simulate.
fprintf('\n========================================================\n');
fprintf('PRESET: %s\n', preset_name);
fprintf('  model class : %s\n', model_class);
fprintf('  condition   : %s\n', cond_name);
fprintf('  run mode    : %s  (%d levels x %d reps in the sweeps)\n', ...
    run_mode, cfg.n_levels, cfg.n_reps);
fprintf('--------------------------------------------------------\n');
fprintf('NETWORK\n');
fprintf('  n / indegree        : %d / %d   (alpha = %.4f)\n', ...
    model.n, model.indegree, model.alpha);
fprintf('  cell types          : %s, f = %s\n', ...
    strjoin(model.cell_type_names, '/'), mat2str(model.f));
fprintf('  mu_tilde_relative   : %s   (post <- pre)\n', mat2str(model.mu_tilde_relative));
fprintf('  sigma_tilde_relative: %s\n', mat2str(model.sigma_tilde_relative));
fprintf('  level_of_chaos      : %g\n', model.level_of_chaos);
fprintf('  F_tracks_network    : %d', model.F_tracks_network);
if ~model.F_tracks_network
    fprintf('   (pinned at n=%d, indegree=%d)', model.F_ref_n, model.F_ref_indegree);
end
fprintf('\n');
fprintf('  F (default_val)     : %.6g\n', model.default_val);
fprintf('  bulk radius R       : %.4f\n', model.R);
fprintf('  outlier eigenvalues : %s\n', mat2str(round(model.lambda_O(:)', 4)));
fprintf('--------------------------------------------------------\n');
fprintf('NONLINEARITY AND ADAPTATION\n');
fprintf('  activation          : %s  (S_a = %g, S_c = %g)\n', ...
    model.activation, model.S_a, model.S_c);
fprintf('  tau_d               : %g\n', model.tau_d);
fprintf('  n_a (SFA per type)  : %s\n', mat2str(model.n_a));
for q = 1:model.n_cellTypes
    if model.n_a(q) > 0
        fprintf('    %-3s tau_a = %s,  c = %g\n', model.cell_type_names{q}, ...
            mat2str(model.tau_a{q}, 4), model.c(q));
    end
end
% mat2str, not %g: the whole point of this preset is that tau_rec and tau_rel
% are ROWS, and %g would silently print only their first element.
fprintf('  STD routes (pre->post)\n');
sc = model.synapse_config;
if isempty(fieldnames(sc))
    fprintf('    (none)\n');
else
    for pre = fieldnames(sc)'
        for post = fieldnames(sc.(pre{1}))'
            s = sc.(pre{1}).(post{1}).std;
            fprintf('    %s->%-3s tau_rec = %-9s tau_rel = %-12s (n_b = %d)\n', ...
                pre{1}, post{1}, mat2str(s.tau_rec), mat2str(s.tau_rel), ...
                numel(s.tau_rec));
        end
    end
end
fprintf('  n_b (pre x post)    : %s\n', mat2str(model.n_b_pairs));

% Steady-state depression, as a standing reminder that the b-states MULTIPLY.
% Each settles at 1/(1 + r*tau_rec/tau_rel); the synapse sees their product. Read
% off the first route's config rather than hardcoded, so this stays honest if the
% preset's taus change. Evaluated at a nominal rate, for orientation only.
pre_names = fieldnames(sc);
if ~isempty(pre_names)
    post_names = fieldnames(sc.(pre_names{1}));
    s_ref = sc.(pre_names{1}).(post_names{1}).std;
    r_ref = 0.3;
    b_ref = 1 ./ (1 + r_ref * s_ref.tau_rec ./ s_ref.tau_rel);
    fprintf('  steady-state b at r = %.2g: %s  ->  product %.4f\n', ...
        r_ref, mat2str(round(b_ref, 4)), prod(b_ref));
end
fprintf('--------------------------------------------------------\n');
fprintf('INTEGRATION\n');
fprintf('  T_range / fs        : %s / %d Hz\n', mat2str(model.T_range), model.fs);
fprintf('  solver              : %s', model.ode_solver);
if cfg.is_stochastic
    fprintf('   (stochastic: chosen because the preset carries noise)\n');
else
    fprintf('   (deterministic)\n');
end
fprintf('  sigma_u_noise       : %g   (input-referred)\n', model.sigma_u_noise);
if model.sigma_u_noise > 0
    fprintf('    diffusion coeff   : %.6g  (sigma_u_noise / tau_d)\n', model.sigma_x_raw);
    fprintf('    x_noise_std       : %.6g  (sigma_u_noise / sqrt(2 tau_d))\n', model.x_noise_std);
end
fprintf('  lya_method / window : %s / %s\n', model.lya_method, mat2str(model.lya_T_interval));
fprintf('  intrinsic_drive     : %g\n', model.input_config.intrinsic_drive);
fprintf('  state dimension     : %d\n', model.N_sys_eqs);
fprintf('========================================================\n\n');
end
