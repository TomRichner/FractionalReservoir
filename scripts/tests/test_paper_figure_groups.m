function test_paper_figure_groups()
% TEST_PAPER_FIGURE_GROUPS Every entry of every paper bundle is in a group, and the paper keeps to its limits (Plan 2, package A).
%
% The figure plan (Tom, 2026-11-21; the manuscript repository's
% docs/Advancing_the_narrative/figures_main_supp_appen.md) sorts the figures
% into three groups: main (eight, no more), supplement (about ten, held at
% ten) and appendix (the rest). paper_figure_groups names them and their
% limits; paper_figure_entry builds a registry row with its group and derives
% in_paper from it. This test fails when
%   1. an entry of a paper bundle (mu7chain_config, mu7chain_AC_config,
%      mu7chain_gainbias_config, ied_stim_config), in any run mode and
%      figure set, has no group or one paper_figure_groups does not name;
%   2. an entry's in_paper disagrees with its group (main and supplement are
%      in the paper, appendix is not);
%   3. two entries of one bundle share a name, or the registry's rows do not
%      share one field order (they must concatenate into a struct array);
%   4. the bundles that make the paper, mu7chain_AC_config (Figures 1 to 7
%      and the supplement) and ied_stim_config (Figure 8), together register
%      more main or supplementary figures than the limits allow, or name one
%      entry in both;
%   5. a name has one group in mu7chain_config or mu7chain_gainbias_config
%      and another in mu7chain_AC_config, which composes them;
%   6. paper_figure_entry accepts a group it does not know.
%
% Quick: builds configs only, no simulation (about 7 s).
%
% See also: paper_figure_groups, paper_figure_entry, test_paper_figure_registry

G = paper_figure_groups();
group_names = {G.name};
modes = run_mode_names();
n_bad = 0;

%% 1 to 3, per bundle, run mode and figure set
% The mu7chain bundles at every run mode, with the mode's own set and the
% other set asked for by name; the IED bundle at each of its modes.
cfgs = {};
for m = modes
    for set_ = {'', '3cond', '4cond'}
        cfgs{end + 1} = {'mu7chain_config', mu7chain_config(m{1}, set_{1})}; %#ok<AGROW>
        cfgs{end + 1} = {'mu7chain_AC_config', mu7chain_AC_config(m{1}, set_{1})}; %#ok<AGROW>
        cfgs{end + 1} = {'mu7chain_gainbias_config', mu7chain_gainbias_config(m{1}, set_{1})}; %#ok<AGROW>
    end
end
for m = {'essential', 'fast', 'production', 'stim_blocks_20261117'}
    cfgs{end + 1} = {'ied_stim_config', ied_stim_config(m{1})}; %#ok<AGROW>
end
for k = 1:numel(cfgs)
    [bundle, cfg] = cfgs{k}{:};
    label = sprintf('%s(%s%s)', bundle, cfg.run_mode, set_text(cfg));
    n_bad = n_bad + check_registry(label, cfg.figures, G);
end
fprintf('%d bundle configs checked: every entry grouped, in_paper derived, names unique\n', numel(cfgs));

%% 4. The paper: mu7chain_AC_config with ied_stim_config, against the limits
for m = modes
    AC = mu7chain_AC_config(m{1});
    ied_mode = m{1};
    if ~ismember(ied_mode, {'essential', 'fast', 'production'}); ied_mode = 'production'; end
    IED = ied_stim_config(ied_mode);
    names = [entry_names(AC.figures), entry_names(IED.figures)];
    groups = [entry_groups(AC.figures), entry_groups(IED.figures)];
    both = intersect(entry_names(AC.figures), entry_names(IED.figures));
    if ~isempty(both)
        fprintf('BAD  paper at %s: entries named in both bundles: %s\n', m{1}, strjoin(both, ', '));
        n_bad = n_bad + 1;
    end
    counts = cell(1, numel(G));
    for j = 1:numel(G)
        n = numel(unique(names(strcmp(groups, G(j).name))));
        counts{j} = sprintf('%s %d', G(j).name, n);
        if n > G(j).limit
            fprintf('BAD  paper at %s: %d %s figures, the limit is %d: %s\n', m{1}, n, G(j).name, ...
                G(j).limit, strjoin(unique(names(strcmp(groups, G(j).name))), ', '));
            n_bad = n_bad + 1;
        end
    end
    fprintf('paper at %-10s (AC + IED %s): %s\n', m{1}, ied_mode, strjoin(counts, ', '));
end

%% 5. The parts agree with the composition
for m = modes
    AC = mu7chain_AC_config(m{1});
    ac_names = entry_names(AC.figures); ac_groups = entry_groups(AC.figures);
    for part = {mu7chain_config(m{1}), mu7chain_gainbias_config(m{1})}
        p_names = entry_names(part{1}.figures); p_groups = entry_groups(part{1}.figures);
        [shared, ia, ib] = intersect(p_names, ac_names);
        differ = ~strcmp(p_groups(ia), ac_groups(ib));
        if any(differ)
            fprintf('BAD  %s: group differs between a part and mu7chain_AC_config: %s\n', m{1}, ...
                strjoin(shared(differ), ', '));
            n_bad = n_bad + 1;
        end
    end
end

%% 6. An unknown group is refused where the row is written
try
    paper_figure_entry('x', @disp, 'supplementary');
    fprintf('BAD  paper_figure_entry accepted the group ''supplementary''\n');
    n_bad = n_bad + 1;
catch ME
    if ~strcmp(ME.identifier, 'paper_figure_entry:BadGroup')
        fprintf('BAD  paper_figure_entry raised %s, not :BadGroup\n', ME.identifier);
        n_bad = n_bad + 1;
    end
end
e = paper_figure_entry('x', @disp, 'supplement', {'a', 1});
if ~(e.in_paper && isequal(fieldnames(e)', {'name', 'fn', 'group', 'in_paper', 'args'}) && isequal(e.args, {'a', 1}))
    fprintf('BAD  paper_figure_entry built %s\n', strjoin(fieldnames(e)', ', '));
    n_bad = n_bad + 1;
end
if ~isequal(group_names, {'main', 'supplement', 'appendix'}) || ~isequal([G.limit], [8 10 Inf])
    fprintf('BAD  paper_figure_groups changed: %s, limits %s\n', strjoin(group_names, ', '), mat2str([G.limit]));
    n_bad = n_bad + 1;
end

assert(n_bad == 0, 'test_paper_figure_groups:Bad', '%d problem(s); see the BAD lines above.', n_bad);
fprintf('test_paper_figure_groups: all checks passed\n');
end

%% ------------------------------------------------------------------------
function n_bad = check_registry(label, figures, G)
n_bad = 0;
names = entry_names(figures);
[~, first] = unique(names, 'stable');
dup = unique(names(setdiff(1:numel(names), first)));
if ~isempty(dup)
    fprintf('BAD  %s: names used twice: %s\n', label, strjoin(dup, ', '));
    n_bad = n_bad + 1;
end
for k = 1:numel(figures)
    f = figures{k};
    if ~isfield(f, 'group')
        fprintf('BAD  %s: %s has no group\n', label, f.name);
        n_bad = n_bad + 1;
        continue
    end
    j = find(strcmp({G.name}, f.group), 1);
    if isempty(j)
        fprintf('BAD  %s: %s has the unknown group ''%s''\n', label, f.name, f.group);
        n_bad = n_bad + 1;
    elseif f.in_paper ~= G(j).in_paper
        fprintf('BAD  %s: %s is %s but in_paper %d\n', label, f.name, f.group, f.in_paper);
        n_bad = n_bad + 1;
    end
end
try
    F = [figures{:}]; %#ok<NASGU>
catch ME
    fprintf('BAD  %s: the rows do not concatenate (%s)\n', label, ME.message);
    n_bad = n_bad + 1;
end
end

function n = entry_names(figures)
n = cellfun(@(f) f.name, figures, 'UniformOutput', false);
end

function g = entry_groups(figures)
g = cellfun(@(f) f.group, figures, 'UniformOutput', false);
end

function s = set_text(cfg)
s = '';
if isfield(cfg, 'figure_set') && ~isempty(cfg.figure_set); s = [', ' cfg.figure_set]; end
end
