function out_dir = run_tau_sensitivity_analysis(ctx)
% RUN_TAU_SENSITIVITY_ANALYSIS Sweep the maximum SFA timescale tau_a_E(end).
%
%   out_dir = RUN_TAU_SENSITIVITY_ANALYSIS(ctx)
%   out_dir = RUN_TAU_SENSITIVITY_ANALYSIS()       % standalone, class defaults
%
% Vector-parameter sweep: tau_a_E is a length-3 logspaced vector and the LAST
% (largest) element is swept from 1 to 30 s, under the SFA+STD condition only.
%
% ctx comes from resolve_run_context('tau_sensitivity', ...).
%
% WAS A SCRIPT reading the master_* base-workspace protocol; see
% resolve_run_context.
%
% See also: resolve_run_context, ParamSpaceAnalysis2, run_all_analyses

arguments
    ctx struct = resolve_run_context('tau_sensitivity')
end

setup_paths();

note = 'tau_timescales';

% Condition: SFA + STD, taken from the shared four-regime set rather than
% respelled here, so this sweep uses exactly the regime the other analyses call
% 'sfa_and_std'. Whichever class is in play, it gives E three SFA timescales --
% which is what n_elements = 3 below is coupled to.
is_sfa_and_std = cellfun(@(c) strcmp(c.name, 'sfa_and_std'), ctx.conditions);
if ~any(is_sfa_and_std)
    error('run_tau_sensitivity_analysis:NoSfaAndStd', ...
        ['No condition named ''sfa_and_std'' in the preset''s condition set ' ...
         '(found: %s). This sweep is defined on that regime.'], ...
        strjoin(cellfun(@(c) c.name, ctx.conditions, 'UniformOutput', false), ', '));
end
condition = ctx.conditions(is_sfa_and_std);

%% tau_a_E(end) sweep -- vector parameter
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity: tau_a_E(end) [1, 30] ===\n');
fprintf('========================================\n');

psa = ParamSpaceAnalysis2( ...
    'n_levels', ctx.n_levels, ...
    'batch_size', 25, ...
    'note', sprintf('%s_tau_a_E_max', note), ...
    'randomize_order', false, ...
    'verbose', ctx.verbose);
psa.folder_prefix  = 'tau_sensitivity';
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = ctx.model_defaults;
if ~isempty(ctx.output_dir)
    psa.output_dir = ctx.output_dir;
end

psa.set_conditions(condition);

% tau_a_E is a vector of length 3, logspaced from 0.25 to max, and we sweep the
% max (last element) from 1 to 30 s. It is a real property on SRNNModel2 and a
% scalar-row alias onto tau_a{1} on SRNNCellTypePairs, so the same axis name
% works for both. n_elements stays coupled BY HAND to the condition's three SFA
% timescales above.
%
% WHAT THIS SWEEP IS FOR. When the network is stable, the slowest linear mode of
% the closed-loop system sets the largest Lyapunov exponent -- and with SFA
% present that mode is usually the slowest adaptation state, giving LLE close to
% -1/tau_a_E(end). This axis therefore demonstrates directly that the adaptation
% timescale CONTROLS the exponent in the stable regime. That is a property of
% the coupled system, not an artefact: a is a state variable like any other, and
% the measured exponents sit slightly below -1/tau because the a<->x coupling
% perturbs the bare eigenvalue.
%
% The range was [5, 60] and is now [1, 30]. The wider relative span (a factor of
% 30 in -1/tau rather than 12) is the point, and the fast end matters most: with
% STD also active its own slowest mode is around -(1/tau_rec + r/tau_rel) ~ -0.6,
% so once -1/tau_a_E(end) drops below that the STD mode becomes the slowest and
% takes over. The prediction to check is therefore NOT a bare -1/tau line but
% max(-1/tau_a_E(end), STD mode) -- a knee around tau_a_E(end) ~ 1.6 s, which
% [1, 30] brackets and [5, 60] did not reach.
psa.add_vector_parameter('tau_a_E', ...
    'vary_element', 'last', ...
    'fixed_value', 0.25, ...
    'vary_range', [1, 30], ...
    'n_elements', 3, ...
    'spacing', 'log', ...
    'level_spacing', 'linear');

psa.add_grid_parameter('reps', 1:ctx.n_reps);

psa.run();

copyfile([mfilename('fullpath') '.m'], psa.output_dir);

psa.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1]);
psa.plot_sensitivity('metric', 'mean_rate');

% run() writes psa_object.mat itself, always under the name `psa`. This script
% used to save it again as `psa_tau_a`, which is why readers had to guess the
% variable name -- and why Fig_sfa_EOC's original loader broke on newer runs.

if ctx.save_figs
    fig_dir = fullfile(psa.output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'tau_sensitivity_tau_a', [], {'fig', 'png'});
    fprintf('Figures saved to %s\n', fig_dir);
end
close all;

out_dir = psa.output_dir;

% NOTE: a tau_b_E_rec sweep used to sit here, commented out. It has been removed
% rather than carried as dead code -- it was SRNNModel2-only (tau_b_E_rec is not
% a SRNNCellTypePairs property, so it would be a hard validate_model_defaults
% error on the current preset), and git history has it at run_tau_sensitivity_
% analysis.m before this commit if it is ever wanted back.

%% Summary
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity Analysis Complete ===\n');
fprintf('tau_a_E results: %s\n', out_dir);
fprintf('========================================\n');
end
