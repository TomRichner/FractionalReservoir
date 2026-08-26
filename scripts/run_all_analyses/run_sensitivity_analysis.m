function out_dirs = run_sensitivity_analysis(ctx)
% RUN_SENSITIVITY_ANALYSIS 1-D parameter sweeps, one ParamSpaceAnalysis2 per axis.
%
%   out_dirs = RUN_SENSITIVITY_ANALYSIS(ctx)
%   out_dirs = RUN_SENSITIVITY_ANALYSIS()          % standalone, class defaults
%
% Sweeps individual parameters across multiple levels and repetitions, comparing
% the four adaptation conditions. Each parameter gets its own PSA run, so each
% is a clean 1-D sweep with the other axes held at the preset's operating point.
%
% ctx comes from resolve_run_context('sensitivity', ...). Called with no
% argument it resolves its own at SRNNModel2 class defaults, so the function is
% still runnable by hand.
%
% Returns a cellstr of the output directory per swept parameter, in the order
% they were swept.
%
% WAS A SCRIPT. It read master_output_dir / master_save_figs /
% master_model_overrides / master_model_class / master_conditions / run_mode /
% save_figs out of the base workspace with exist(). See resolve_run_context for
% why that had to go.
%
% See also: resolve_run_context, ParamSpaceAnalysis2, run_all_analyses

arguments
    ctx struct = resolve_run_context('sensitivity')
end

setup_paths();

note = 'sensitivity';

% LLE histogram y-axis range for the sensitivity heatmaps (plot_sensitivity).
lle_hist_range = [-2, 2];

%% Which parameters to sweep
% {param_name, [min, max]}. The fraction-excitatory axis is named per class --
% see ctx.f_param.
% level_of_chaos [0.25, 2.5] rather than [0.5, 1.5]: measured, not guessed. A
% medium sweep of all seven regimes puts their edge-of-chaos crossings at
%
%   no_adaptation 0.54   sfa_only_oneTS 0.59   sfa_only 0.48
%   std_only_oneTS 1.18  std_only 2.27
%   sfa3_std1 0.98       sfa_and_std 2.43
%
% so [0.5, 1.5] bracketed only three of the seven and reported the depressing
% regimes as flat and stable throughout -- true, but carrying no information
% about WHERE they stop being stable. [0.25, 2.5] contains all seven, confirmed
% by two independent sweeps at different resolutions agreeing within ~0.1 on six
% of them. Reproduce with scripts/examples/explore_sensitivity_range.m.
params_to_sweep = {
    'n',              [100, 1000];
    ctx.f_param,      [0.2, 0.8];
    'level_of_chaos', [0.25, 2.5];
    };

% The four connectivity blocks, swept relative to the preset's own operating
% point. mu_block_from_preset decides the ranges and is shared with
% run_param_space_analysis, so the 1-D levels and the grid axes cover exactly
% the same span.
if strcmp(ctx.model_class, 'SRNNCellTypePairs')
    mu_ranges = mu_block_from_preset(ctx.preset_defaults);
    for b_idx = 1:size(mu_ranges, 1)
        params_to_sweep(end+1, :) = mu_ranges(b_idx, :); %#ok<AGROW>
        if ctx.verbose
            fprintf('[sensitivity] %s sweeping [%+g %+g]\n', ...
                mu_ranges{b_idx, 1}, mu_ranges{b_idx, 2}(1), mu_ranges{b_idx, 2}(2));
        end
    end
end

%% Run one sweep per parameter
n_params = size(params_to_sweep, 1);
out_dirs = cell(1, n_params);

for p_idx = 1:n_params
    param_name  = params_to_sweep{p_idx, 1};
    param_range = params_to_sweep{p_idx, 2};

    fprintf('\n========================================\n');
    fprintf('=== Sensitivity: %s [%.3g, %.3g] (%d/%d) ===\n', ...
        param_name, param_range(1), param_range(2), p_idx, n_params);
    fprintf('========================================\n');

    psa = ParamSpaceAnalysis2( ...
        'n_levels', ctx.n_levels, ...
        'batch_size', 25, ...
        'note', sprintf('%s_%s', note, param_name), ...
        'randomize_order', false, ...   % ordered axes matter for a sensitivity sweep
        'verbose', ctx.verbose);
    psa.folder_prefix   = '1D_sensitivity';
    psa.model_class     = ctx.model_class;
    psa.integer_params  = ctx.integer_params;
    psa.model_defaults  = ctx.model_defaults;
    if ~isempty(ctx.output_dir)
        psa.output_dir = ctx.output_dir;
    end

    psa.set_conditions(ctx.conditions);

    if strcmp(ctx.model_class, 'SRNNModel2')
        % STD with two recovery timescales. Give the STD conditions n_b_E = 2
        % (two E depression timescales, product of the two b columns) with a
        % 1x2 recovery-time-constant vector; the release constant tau_b_E_rel
        % stays scalar/shared. n_b_E must come from the CONDITION
        % (ParamSpaceAnalysis2 ignores it in model_defaults), while tau_b_E_rec
        % flows through model_defaults. STD is on E neurons only.
        %
        % SRNNCellTypePairs says the same thing per route inside synapse_config,
        % so there is nothing to override for that class -- and tau_b_E_rec is
        % not one of its properties, so setting it would be a hard
        % validate_model_defaults error rather than a no-op.
        for ci = 1:numel(psa.conditions)
            if psa.conditions{ci}.n_b_E > 0
                psa.conditions{ci}.n_b_E = 2;
            end
        end
        psa.model_defaults.tau_b_E_rec = [0.5, 5];   % two E recovery timescales (s)
    end

    psa.add_grid_parameter(param_name, param_range);
    psa.add_grid_parameter('reps', 1:ctx.n_reps);

    psa.run();

    % Copy this file for reproducibility. mfilename is the function name here,
    % which is also the file name, so this still resolves.
    copyfile([mfilename('fullpath') '.m'], psa.output_dir);

    psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_hist_range);
    psa.plot_sensitivity('metric', 'mean_rate');

    % psa_object.mat is written by run() itself -- once before batching so a
    % crashed run stays recoverable, and again on completion.

    out_dirs{p_idx} = psa.output_dir;

    if ctx.save_figs
        fig_dir = fullfile(psa.output_dir, 'figures');
        save_some_figs_to_folder_2(fig_dir, ...
            sprintf('sensitivity_%s', param_name), [], {'fig', 'png'});
        fprintf('Figures saved to %s\n', fig_dir);
    end
    close all;
end

%% Summary
fprintf('\n=== Sensitivity Analysis Complete ===\n');
fprintf('Parameters analyzed:\n');
for p_idx = 1:n_params
    fprintf('  %s: %s\n', params_to_sweep{p_idx, 1}, out_dirs{p_idx});
end
end
