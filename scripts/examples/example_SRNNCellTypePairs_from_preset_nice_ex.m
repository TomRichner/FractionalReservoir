% example_SRNNCellTypePairs_from_preset.m
% One SRNNCellTypePairs simulation at exactly the operating point the
% sensitivity sweeps use, with the full plotting set.
%
% The SRNNCellTypePairs analogue of test_SRNN2_defaults.m, which is just
%
%   model = SRNNModel2('n_a_E', 3, 'n_b_E', 1, 'f', 0.5);
%
% This needs more than that, and the extra is worth understanding rather than
% copying: a sweep's operating point is assembled from THREE sources, and an
% example that hardcodes the numbers instead would drift away from the sweeps
% the moment any of them changed.
%
%   srnn_param_preset      the physics -- n, connectivity, nonlinearity, noise --
%                          AND the adaptation conditions, its third output. It
%                          used to delegate those to srnn_adaptation_conditions;
%                          each preset now states its own regimes outright.
%   analysis_run_config    the cost/fidelity knobs -- fs, T_range, integrator
%
% They are combined here in the same precedence ParamSpaceAnalysis2 uses in
% run_single_job: preset defaults first, run-mode timings over them, condition
% fields last. So this script simulates one grid point of the sensitivity sweep,
% and changing the sweep changes this too.
%
% See also: test_SRNN2_defaults, srnn_param_preset, analysis_run_config,
%           run_sensitivity_analysis

close all; clear; clc;

%% What to simulate
% Edit these three lines to move to another preset, regime or fidelity.
preset_name    = 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01';
condition_name = 'sfa_and_std';   % SFA on E only + STD on all four routes
run_mode       = 'medium';        % matches the sensitivity sweeps: fs 400, T [0 20]

% A fixed seed, because this is one illustrative network rather than a sample.
% The sweeps vary it per job -- ParamSpaceAnalysis2 gives every grid point its
% own network seed -- so any single run here is one draw from what they average.
rng_seeds = [19 20];

% Eigenvalue snapshots, in seconds. Must fall inside T_range.
eig_times = [2 10 18];

%% Assemble the operating point
[preset_defaults, model_class, conditions] = srnn_param_preset(preset_name);
cfg = analysis_run_config('sensitivity', run_mode, preset_defaults);

cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond = conditions{strcmp(cond_names, condition_name)};

% Preset first, run_mode timings second (so the integrator and fs come from the
% run mode), condition last (so n_a and synapse_config win) -- run_single_job's
% order exactly.
model_args = namedargs2cell(merge_struct(preset_defaults, cfg.model));
model_args = [model_args, {'rng_seeds', rng_seeds}];
for f = setdiff(fieldnames(cond), {'name'})'
    model_args = [model_args, {f{1}, cond.(f{1})}]; %#ok<AGROW>
end

model = feval(model_class, model_args{:});

% plot_eigenvalues reads the Jacobian off the stored trajectory, so the full
% state has to be kept. At n=500, fs=400 over [0 20] that is ~8001 x 2250
% doubles, about 140 MB -- fine here, but the sweeps leave it false for a
% reason. Drop it if you raise n or fs.
model.store_full_state = true;

%% Report the configuration before spending time on it
fprintf('\n========================================================\n');
fprintf('PRESET: %s\n', preset_name);
fprintf('  model class : %s\n', model_class);
fprintf('  condition   : %s\n', cond.name);
fprintf('  run mode    : %s  (%d levels x %d reps in the sweeps)\n', ...
    run_mode, cfg.n_levels, cfg.n_reps);
fprintf('--------------------------------------------------------\n');
fprintf('NETWORK\n');
fprintf('  n / indegree        : %d / %d   (alpha = %.4f)\n', model.n, model.indegree, model.alpha);
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
fprintf('  mu_tilde (absolute) : %s\n', mat2str(model.mu_tilde, 4));
fprintf('  bulk radius R       : %.4f\n', model.R);
fprintf('  outlier eigenvalues : %s\n', mat2str(round(model.lambda_O(:)', 4)));
fprintf('--------------------------------------------------------\n');
fprintf('NONLINEARITY AND ADAPTATION\n');
fprintf('  activation          : %s  (S_a = %g, S_c = %g)\n', ...
    model.activation, model.S_a, model.S_c);
if isempty(model.mu_S_c) && ~any(model.sigma_S_c > 0)
    fprintf('  setpoints           : homogeneous (S_c_vec stays empty)\n');
else
    fprintf('  setpoints           : mu_S_c = %s, sigma_S_c = %s\n', ...
        mat2str(model.mu_S_c), mat2str(model.sigma_S_c));
end
fprintf('  tau_d               : %g\n', model.tau_d);
fprintf('  n_a (SFA per type)  : %s\n', mat2str(model.n_a));
for q = 1:model.n_cellTypes
    if model.n_a(q) > 0
        fprintf('    %-3s tau_a = %s,  c = %g\n', model.cell_type_names{q}, ...
            mat2str(model.tau_a{q}, 4), model.c(q));
    end
end
fprintf('  STD routes (pre->post)\n');
sc = model.synapse_config;
if isempty(fieldnames(sc))
    fprintf('    (none)\n');
else
    for pre = fieldnames(sc)'
        for post = fieldnames(sc.(pre{1}))'
            s = sc.(pre{1}).(post{1}).std;
            fprintf('    %s->%-3s tau_rec = %-5g tau_rel = %g\n', ...
                pre{1}, post{1}, s.tau_rec, s.tau_rel);
        end
    end
end
fprintf('  n_b (pre x post)    : %s\n', mat2str(model.n_b_pairs));
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
    % Read the derived values off the model: both are Dependent, so they cannot
    % disagree with sigma_u_noise the way hand-computed numbers could.
    fprintf('    diffusion coeff   : %.6g  (sigma_u_noise / tau_d)\n', model.sigma_x_raw);
    fprintf('    x_noise_std       : %.6g  (sigma_u_noise / sqrt(2 tau_d))\n', model.x_noise_std);
    if strcmp(model.activation, 'piecewise')
        % piecewiseSigmoid saturates at c +/- (1 - S_a/2): with a = S_a/2 its
        % outer breakpoints are x1 = c + a - 1 and x4 = c + 1 - a. So the
        % half-span is 1 - S_a/2, NOT anything linear in S_a -- at S_a = 0.8
        % that is 0.6, which is the number the preset comments quote against.
        half_span = 1 - model.S_a / 2;
        fprintf('    vs the piecewise half-span (+/- %.3g): %.1f%%\n', ...
            half_span, 100 * model.x_noise_std / half_span);
    end
end
fprintf('  lya_method / window : %s / %s\n', model.lya_method, mat2str(model.lya_T_interval));
fprintf('  lya_warmup          : %g s  (iteration starts %g s before the window)\n', ...
    model.lya_warmup, model.lya_warmup);
fprintf('  intrinsic_drive     : %g\n', model.input_config.intrinsic_drive);
fprintf('  state dimension     : %d\n', model.N_sys_eqs);
fprintf('========================================================\n\n');

%% Build and run
model.S_c = [0.20];
model.level_of_chaos = 1;
model.sigma_u_noise = 0.03;
model.plot_deci = model.fs/20; % go from 400 Hz to 20

sc = model.synapse_config;
for pre = fieldnames(sc)'
    for post = fieldnames(sc.(pre{1}))'
        sc.(pre{1}).(post{1}).std.tau_rec = 2;
        sc.(pre{1}).(post{1}).std.tau_rel = 0.25;
    end
end
model.synapse_config = sc;

model.build();
model.run();

%% Plots
model.plot();                    % compact summary: one panel per quantity

% Per-cell-type view: one COLUMN per type, every neuron's trace, with b and g
% collapsed across routes as prod(b) and coloured by target type.
model.plot_celltypes();

% Effective-Jacobian spectrum at three times: early (still settling), middle and
% late. Watch whether adaptation pulls eigenvalues back inside the unit circle
% as the network engages.
model.plot_eigenvalues(eig_times);

% W itself on a diverging colormap, so zero is white and the sign of each
% synapse reads at a glance. Cell-type boundaries mark the (post, pre) blocks.
% This is the SCALED W -- what was simulated, level_of_chaos included.
model.plot_W();

% Eigenvalues of W. This is where the block structure shows up directly: the
% outlier leaving the bulk disk of radius R, rather than only as the printed
% lambda_O.
model.plot_W_spectrum();

% Weight distribution per cell type. The negative tail on the E columns is the
% Gaussian sampler, not a bug -- Dale's law here is statistical, so a fraction
% ~normcdf(-mu_tilde/sigma_tilde) of E synapses come out negative.
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
