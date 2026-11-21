function test_paper_figure_registry()
% TEST_PAPER_FIGURE_REGISTRY Every figure the manuscript embeds is a registered entry of a run bundle (FD-031 item 11).
%
% Tom (2026-11-17): every figure fix must come back on a new run from a
% clean main, with no hand step. That holds only while each figure the
% manuscript embeds is an entry of some bundle's cfg.figures, so the run
% draws it. This test lists the paper's figures (PAPER_FIGURES below: the
% PNG the manuscript embeds, the bundle that draws it, the registry entry)
% and checks each against that bundle's cfg.figures, in the run modes the
% paper is built from ('fast' and 'fastmedium', figure set '3cond'). Removing
% a registration fails it.
%
% EXPECTED AFTER MERGE. A row whose fifth column names a job (e.g. 'J3') is
% registered by work on another branch. While that row's bundle config
% function does not exist on the path, the row is skipped with a SKIP line
% naming the job; as soon as the function exists (the branch is merged), the
% row is checked like any other, so the test becomes strict by itself, with
% no edit. A row whose bundle exists but lacks the entry is a failure, merged
% or not. When every job has merged, the fifth column can be emptied
% ('' everywhere); nothing else changes.
%
% THE LIST AGAINST THE PAPER REPO. When the manuscript repo is found beside
% this repo's main clone (StochasticPlasticDynamicalSystemPaper, whose
% figure_sources.json lists every embedded PNG under "figures"), every PNG
% named there must have a row here, so a new embed without a registration
% fails too. When the paper repo is not found the cross-check is skipped
% with a SKIP line (the registry checks still run).
%
% in_paper. An entry that the paper embeds but marks in_paper false is
% reported as a NOTE line, not a failure: the flag only labels the manifest,
% and setting it is the owning config's decision.
%
% Quick: builds four bundle configs (no simulation), about 10 s.

%% The paper's figures (manuscript embeds, 2026-11-17)
% {png embedded, bundle, registry entry, what it is, expected-after-merge job}
% Bundles: 'AC' = mu7chain_AC_config (blocks A + C; Figs 1-6 and the
% supplement), 'C' = mu7chain_gainbias_config (block C; Fig 7), 'IED' =
% ied_stim_config (the IED stimulation run; Fig 8 and its supplement, FD-031
% items 4-5, job J3).
PAPER_FIGURES = {
    'Fig_Main1_Grouped.png',              'AC',  'fig_main1_grouped',               'Fig 1',   ''
    'Fig_Main2_Grouped.png',              'AC',  'fig_main2_grouped',               'Fig 2',   ''
    % Fig 3: the new layout (no lambda row) is the entry's default output
    % since J2 (FD-031 item 7), fig_main3_grouped/Fig_Main3_Grouped.png. The
    % old name Fig_Main3_Grouped_noC.png is retired (RETIRED_NAMES below).
    'Fig_Main3_Grouped.png',              'AC',  'fig_main3_grouped',               'Fig 3',   ''
    'Fig_Main4_Grouped.png',              'AC',  'fig_main4_grouped',               'Fig 4',   ''
    'Fig_Main5_Grouped.png',              'AC',  'fig_main5_grouped',               'Fig 5',   ''
    'Fig_Main6_Grouped.png',              'AC',  'fig_main6_grouped',               'Fig 6',   ''
    'Fig_Main_MC_Grouped.png',            'C',   'fig_main_mc_grouped',             'Fig 7',   ''
    'Fig_Main_MC_Grouped.png',            'AC',  'fig_main_mc_grouped',             'Fig 7 (A + C bundle)', ''
    'Fig_Stim_Combined.png',              'IED', 'Fig_Stim_Combined',               'Fig 8',   ''
    'fig_ied_stim_blocks_supp_state.png', 'IED', 'fig_ied_stim_blocks_supp_state',  'supplement', ''
    'Fig_EI_Weights_ParamSpace_LLE.png',  'AC',  'fig_EI_weights_param_space',      'supplement', ''
    'Fig_EI_Weights_ParamSpace_mean_rate.png', 'AC', 'fig_EI_weights_param_space',  'supplement', ''
    'Fig_Eig_Heatmap_Imbalance.png',      'AC',  'fig_eig_heatmap_imbalance',       'supplement', ''
    'Fig_Local_vs_Finite_LLE.png',        'AC',  'fig_local_vs_finite_lle',         'supplement', ''
    'Fig_MC_Example.png',                 'AC',  'fig_memory_capacity_example',     'supplement', ''
    'Fig_Numerics_Jacobian.png',          'AC',  'fig_numerics_jacobian',           'supplement', ''
    'Fig_numerics_lya_method.png',        'AC',  'fig_numerics_lya_method',         'supplement', ''
    'Fig_numerics_solver.png',            'AC',  'fig_numerics_solver',             'supplement', ''
    'Fig_SFA_steady_state.png',           'AC',  'fig_SFA_steady_state',            'supplement', ''
    'Fig_STD_steady_state.png',           'AC',  'fig_STD_steady_state',            'supplement', ''
    'Fig_Sensitivity_LLE_core.png',       'AC',  'fig_sensitivity_analysis_allStd', 'supplement', ''
    'Fig_Sensitivity_LLE_mu.png',         'AC',  'fig_sensitivity_analysis_allStd', 'supplement', ''
    'Fig_Sensitivity_mean_rate_core.png', 'AC',  'fig_sensitivity_analysis_allStd', 'supplement', ''
    'Fig_Sensitivity_mean_rate_mu.png',   'AC',  'fig_sensitivity_analysis_allStd', 'supplement', ''
    'Fig_Sensitivity_mean_rate_medians.png', 'AC', 'fig_sensitivity_medians',       'supplement', ''
    'Fig_single_neuron_revised.png',      'AC',  'fig_single_neuron_revised',       'supplement', ''
    'bursting_synchrony.png',             'AC',  'fig_stim_engages_adaptation',     'supplement', ''
    'bursting_timeseries.png',            'AC',  'fig_stim_engages_adaptation',     'supplement', ''
    };
% PNG names the paper repo may still list that a figure step no longer
% writes, each mapped to the name that replaced it: {old, current}. The
% cross-check reads figure_sources.json through this table, so the paper
% repo can be updated after this repo without the test failing in between
% (a NOTE line says the old name is still listed there).
RETIRED_NAMES = {
    'Fig_Main3_Grouped_noC.png',          'Fig_Main3_Grouped.png'     % Fig 3, FD-031 item 7
    };
BUNDLE_FN = struct('AC', 'mu7chain_AC_config', 'C', 'mu7chain_gainbias_config', 'IED', 'ied_stim_config');
MODES = {'fast', 'fastmedium'};

n_bad = 0; n_ok = 0; n_skip = 0;

%% The registries, per bundle and mode
regs = struct();
for b = fieldnames(BUNDLE_FN)'
    fn = BUNDLE_FN.(b{1});
    regs.(b{1}) = {};
    if ~exist(fn, 'file'); continue; end
    if strcmp(b{1}, 'IED')
        cfgs = {feval(fn)};                       % its own default mode
    else
        cfgs = cellfun(@(m) feval(fn, m, '3cond'), MODES, 'UniformOutput', false);
    end
    regs.(b{1}) = cfgs;
end

%% Each paper figure against its bundle's registry
for r = 1:size(PAPER_FIGURES, 1)
    [png, b, entry, what, job] = PAPER_FIGURES{r, :};
    fn = BUNDLE_FN.(b);
    if ~exist(fn, 'file')
        if ~isempty(job)
            fprintf('SKIP %-40s %s: bundle %s not on this branch yet; expected after merge of %s\n', ...
                png, what, fn, job);
            n_skip = n_skip + 1;
        else
            fprintf('FAIL %-40s %s: bundle %s does not exist\n', png, what, fn);
            n_bad = n_bad + 1;
        end
        continue
    end
    for k = 1:numel(regs.(b))
        cfg = regs.(b){k};
        names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
        hit = find(strcmp(names, entry), 1);
        label = sprintf('%s(%s)', fn, mode_of(cfg));
        if isempty(hit)
            fprintf('FAIL %-40s %s: no entry %s in %s\n', png, what, entry, label);
            n_bad = n_bad + 1;
        else
            n_ok = n_ok + 1;
            if isfield(cfg.figures{hit}, 'in_paper') && ~cfg.figures{hit}.in_paper
                fprintf('NOTE %-40s %s: entry %s in %s is marked in_paper false\n', png, what, entry, label);
            end
        end
    end
end

%% The list against the paper repo's figure_sources.json
src = paper_figure_sources();
if isempty(src)
    fprintf('SKIP cross-check: the paper repo''s figure_sources.json was not found beside this clone\n');
else
    listed = embedded_pngs(src);
    [is_old, k_old] = ismember(listed, RETIRED_NAMES(:, 1));
    for k = find(is_old(:))'
        fprintf('NOTE %-40s is retired; figure_sources.json still lists it (read as %s)\n', ...
            listed{k}, RETIRED_NAMES{k_old(k), 2});
        listed{k} = RETIRED_NAMES{k_old(k), 2};
    end
    listed = unique(listed);
    missing = setdiff(listed, PAPER_FIGURES(:, 1));
    for k = 1:numel(missing)
        fprintf('FAIL %-40s embedded by the manuscript (%s) but has no row in PAPER_FIGURES\n', missing{k}, src);
        n_bad = n_bad + 1;
    end
    extra = setdiff(PAPER_FIGURES(:, 1), listed);
    for k = 1:numel(extra)
        fprintf('NOTE %-40s has a row here but is not in figure_sources.json\n', extra{k});
    end
end

fprintf('\ntest_paper_figure_registry: %d registrations found, %d skipped (expected after merge), %d missing\n', ...
    n_ok, n_skip, n_bad);
assert(n_bad == 0, 'test_paper_figure_registry:Missing', ...
    '%d paper figure(s) have no registered entry in their bundle (see the lines above).', n_bad);
end

%% ------------------------------------------------------------------------
function m = mode_of(cfg)
m = '';
if isfield(cfg, 'run_mode'); m = cfg.run_mode; end
if isfield(cfg, 'figure_set'); m = [m ', ' cfg.figure_set]; end
end

function f = paper_figure_sources()
% <dir of the main clone>/../StochasticPlasticDynamicalSystemPaper/figure_sources.json,
% found through git's common dir so a worktree finds it too; '' if absent.
f = '';
root = fileparts(which('setup_paths'));
[st, out] = system(sprintf('git -C "%s" rev-parse --path-format=absolute --git-common-dir', root));
cands = {fileparts(root)};
if st == 0
    common = strtrim(out);
    cands = [{fileparts(fileparts(common))}, cands];   % <main clone>/.git -> parent of main clone
end
for k = 1:numel(cands)
    c = fullfile(cands{k}, 'StochasticPlasticDynamicalSystemPaper', 'figure_sources.json');
    if isfile(c); f = c; return; end
end
end

function keys = embedded_pngs(file)
% The keys (PNG names) of the "figures" object of figure_sources.json,
% verbatim: jsondecode would mangle them into valid field names, so they are
% read from the text, and their count checked against the decoded object.
txt = fileread(file);
J = jsondecode(txt);
assert(isfield(J, 'figures'), 'test_paper_figure_registry:NoFigures', '%s has no "figures".', file);
keys = regexp(txt, '"([^"]+\.png)"\s*:\s*\{', 'tokens');
keys = unique(cellfun(@(c) c{1}, keys, 'UniformOutput', false));
assert(numel(keys) == numel(fieldnames(J.figures)), 'test_paper_figure_registry:Keys', ...
    'Read %d PNG keys from %s but it has %d figures.', numel(keys), file, numel(fieldnames(J.figures)));
end
