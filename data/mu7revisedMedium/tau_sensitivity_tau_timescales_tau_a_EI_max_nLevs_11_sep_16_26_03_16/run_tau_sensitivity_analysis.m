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
% WHY ONE CONDITION, AND WHAT THAT COSTS -- three limitations that are all the
% same underlying fact: a VECTOR grid parameter is configured once per PSA,
% while the number of SFA timescales is a property of each CONDITION.
%
%   1. n_elements IS PER-SWEEP, NOT PER-CONDITION. add_vector_parameter takes a
%      single n_elements for the whole run, so a sweep spanning regimes with
%      DIFFERENT timescale counts -- comparing one-timescale against
%      three-timescale adaptation while varying tau, say -- is not expressible.
%      That is why this filters to the full-adaptation regime rather than
%      sweeping the whole condition set. The regime is matched by NAME, and
%      presets use two: 'sfa_and_std' (4- and 7-condition) and 'sfa3_std2'
%      (3-condition), identical physics under names chosen for each set's own
%      comparison. Both are accepted; see the lookup below.
%
%   2. THE GRID OVERRIDES THE CONDITION, silently. Precedence in run_single_job
%      is model_defaults < condition < grid, and tau_a_E's setter writes
%      tau_a{1} outright. So running this against a multi-condition set would
%      hand no_adaptation a three-element tau_a and collapse every regime into
%      the same one, with nothing reporting it. The single-condition filter is
%      what prevents that, not a check.
%
%   3. n_elements = 3 IS COUPLED BY HAND to the condition's ladder (see the
%      comment at the add_vector_parameter call). Nothing verifies the two
%      agree. Since each preset now states its conditions in full -- writing
%      log_ladder(lo, hi, K) into its own case -- a preset retuned to K = 2
%      would leave this sweep imposing three timescales on a two-timescale
%      regime, and c/K keeps the TOTAL adaptation right, so even the firing
%      rates would look plausible. All three dualStd presets currently use
%      K = 3, so the value is correct today. If you change a preset's ladder,
%      change n_elements below to match. (Deriving it from the condition --
%      numel(condition{1}.tau_a{1}) -- would remove the coupling entirely and
%      is the obvious fix if this bites.)
%
% A FOURTH, currently unreachable: with vary_element = 'first' AND
% n_elements = 1, ParamSpaceAnalysis2 builds the vector with
% logspace(start, end, 1), which returns END -- the fixed value -- so the swept
% axis would have no effect and every grid point would be identical. This sweep
% uses 'last', where logspace returns the varied value and the behaviour is
% correct. log_ladder is the fix if that configuration is ever wanted.
%
% WAS A SCRIPT reading the master_* base-workspace protocol; see
% resolve_run_context.
%
% See also: resolve_run_context, ParamSpaceAnalysis2, run_all_analyses,
%           log_ladder, srnn_param_preset

arguments
    ctx struct = resolve_run_context('tau_sensitivity')
end

setup_paths();
% Figures this stage draws are created INVISIBLE (TR, 2026-09-14): a new figure
% window raises itself and takes keyboard and mouse focus, and a sweep draws
% many. They still save; with_graphics_defaults restores the root default when
% the guard goes out of scope.
fig_guard = with_graphics_defaults('DefaultFigureVisible', 'off'); %#ok<NASGU>

note = 'tau_timescales';

% Condition: the FULL-ADAPTATION regime -- all SFA timescales and all depression
% timescales -- taken from the preset rather than respelled here, so this sweep
% runs exactly the regime the other analyses do. Whichever class is in play it
% gives E three SFA timescales, which is what n_elements = 3 below is coupled to.
%
% RESOLVED, NOT NAMED. Condition names state their adaptation structure, so the
% full-adaptation regime is called sfa3_std2 in the paper's presets, sfa3_std1
% where depression has a single timescale, and sfa1_std1_stf1 in the facilitation
% preset. No literal covers those, which is what made single_multi_TS_run fail
% here with "No condition named 'sfa_and_std'" -- twice, since a two-name list
% was only a wider guess. full_adaptation_condition errors if nothing conforms
% rather than falling back to a positional guess.
cond_names = cellfun(@(c) c.name, ctx.conditions, 'UniformOutput', false);
full_name  = full_adaptation_condition(ctx.conditions);
idx        = find(strcmp(cond_names, full_name), 1);
condition  = ctx.conditions(idx);
vprintf(ctx.verbose, 'verbose', 'Condition: %s (the full-adaptation regime)\n', full_name);

%% tau_a_E(end) sweep -- vector parameter
vprintf(ctx.verbose, 'verbose', '\n========================================\n');
% WHICH TYPES ARE SWEPT (2026-09-14). The reference preset adapts E AND I with
% the same ladder, and the manuscript's Methods say the slowest time constant
% was varied in both. The earlier sweep moved tau_a_E only -- E's ladder -- so
% it did not match that description. The axis is now chosen from the
% condition: tau_a_EI (one ladder written to every type) when the condition
% gives >= 2 types identical non-empty ladders, else tau_a_E as before. Which
% one ran is recorded in tau_levels.md next to the resolved ladders per level.
if strcmp(ctx.model_class, 'SRNNCellTypePairs') && iscell(condition{1}.tau_a) && ...
        numel(condition{1}.tau_a) >= 2 && ~isempty(condition{1}.tau_a{1}) && ...
        all(cellfun(@(r) isequal(reshape(r, 1, []), reshape(condition{1}.tau_a{1}, 1, [])), condition{1}.tau_a))
    tau_axis = 'tau_a_EI';
else
    tau_axis = 'tau_a_E';
end
n_elements = numel(condition{1}.tau_a{1});   % derived from the condition, not hand-coupled
vprintf(ctx.verbose, 'verbose', '=== Tau Sensitivity: %s(end) [1, 30], %d elements ===\n', tau_axis, n_elements);

vprintf(ctx.verbose, 'verbose', '========================================\n');

psa = ParamSpaceAnalysis2( ...
    'n_levels', ctx.n_levels, ...
    'batch_size', 25, ...
    'note', sprintf('%s_%s_max', note, tau_axis), ...   % folder names the axis that ran
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
psa.add_vector_parameter(tau_axis, ...
    'vary_element', 'last', ...
    'fixed_value', 0.25, ...
    'vary_range', [1, 30], ...
    'n_elements', n_elements, ...
    'spacing', 'log', ...
    'level_spacing', 'linear');

psa.add_grid_parameter('reps', 1:ctx.n_reps);

psa.run();

copyfile([mfilename('fullpath') '.m'], psa.output_dir);

% The manifest that proves which types were swept: the resolved E and I
% ladders at every level, read back through effective_param from one
% successful job per level, into tau_levels.mat and tau_levels.md.
write_tau_levels(psa, tau_axis, ctx.model_class);

psa.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1]);
psa.plot_sensitivity('metric', 'mean_rate');

% run() writes psa_object.mat itself, always under the name `psa`. This script
% used to save it again as `psa_tau_a`, which is why readers had to guess the
% variable name -- and why Fig_sfa_EOC's original loader broke on newer runs.

if ctx.save_figs
    fig_dir = fullfile(psa.output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'tau_sensitivity_tau_a', [], {'fig', 'png'});
    vprintf(ctx.verbose, 'verbose', 'Figures saved to %s\n', fig_dir);
end
close all;

out_dir = psa.output_dir;

% NOTE: a tau_b_E_rec sweep used to sit here, commented out. It has been removed
% rather than carried as dead code -- it was SRNNModel2-only (tau_b_E_rec is not
% a SRNNCellTypePairs property, so it would be a hard validate_model_defaults
% error on the current preset), and git history has it at run_tau_sensitivity_
% analysis.m before this commit if it is ever wanted back.

%% Summary
vprintf(ctx.verbose, 'minimal', '[tau_sensitivity] axis %s, %d elements\n', tau_axis, n_elements);
vprintf(ctx.verbose, 'verbose', '\n========================================\n');
vprintf(ctx.verbose, 'verbose', '=== Tau Sensitivity Analysis Complete ===\n');
vprintf(ctx.verbose, 'minimal', '[tau_sensitivity] complete: %s\n', out_dir);
vprintf(ctx.verbose, 'verbose', '========================================\n');
end

%% ------------------------------------------------------------------------
function write_tau_levels(psa, tau_axis, model_class)
% Resolve the SFA ladders that actually ran at every level, for both cell
% types, by constructing the model each job would have built (effective_param
% gives the swept vector; the other type's ladder comes from the condition),
% and record them. A row per level: index, the swept vector, tau_a{1},
% tau_a{2} (SRNNCellTypePairs) or tau_a_E (SRNNModel2).
lookup = psa.vector_param_lookup.(tau_axis);
n_lev  = numel(lookup);
levels = struct('level', num2cell(1:n_lev), 'swept', reshape(lookup, 1, []), 'tau_a_E', cell(1, n_lev), 'tau_a_I', cell(1, n_lev));
cond = psa.conditions{1};
for k = 1:n_lev
    args = [struct2namevalue(psa.model_defaults), struct2namevalue(rmfield(cond, 'name')), ...
        {tau_axis, lookup{k}}];
    m = feval(model_class, args{:});     % construction only: no build, no run
    if strcmp(model_class, 'SRNNCellTypePairs')
        levels(k).tau_a_E = m.tau_a{1};
        if numel(m.tau_a) >= 2; levels(k).tau_a_I = m.tau_a{2}; end
    else
        levels(k).tau_a_E = m.tau_a_E;
        levels(k).tau_a_I = [];
    end
end
save(fullfile(psa.output_dir, 'tau_levels.mat'), 'levels', 'tau_axis', 'model_class');
fid = fopen(fullfile(psa.output_dir, 'tau_levels.md'), 'w');
if fid > 0
    fprintf(fid, '# tau sweep levels\n\nAxis: `%s` (model %s). Condition: %s.\n\n', tau_axis, model_class, cond.name);
    fprintf(fid, '| level | swept vector | tau_a E | tau_a I |\n|---|---|---|---|\n');
    for k = 1:n_lev
        fprintf(fid, '| %d | %s | %s | %s |\n', k, mat2str(levels(k).swept, 4), ...
            mat2str(levels(k).tau_a_E, 4), mat2str(levels(k).tau_a_I, 4));
    end
    fclose(fid);
end
end
