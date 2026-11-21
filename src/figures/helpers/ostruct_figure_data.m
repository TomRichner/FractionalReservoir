function F = ostruct_figure_data(cfg, caller)
% OSTRUCT_FIGURE_DATA What the outlier_structures figures share: the stage's measures, the networks in their groups and order, the conditions, and the statistics over seeds.
%
%   F = OSTRUCT_FIGURE_DATA(cfg, caller)
%
%   cfg     the figure's cfg: data_file (the stage's ostruct_index.mat), run_dir, conditions
%   caller  the figure's name, for the error identifiers
%
% Finds <run_dir>/outlier_structures/ostruct_index.mat through
% resolve_data_file (ONE run directory or the explicit file, no fallback)
% and reads the stage through ostruct_load, its one reader.
%
% F
%   file, folder, S (ostruct_load's output: runs, checks, settings, T, ...)
%   T            the measures, one row per finished run (all seeds)
%   cond_names   select_conditions over the conditions the runs hold, in
%                the plan's order; colors (manuscript_style), titles
%                ('sentence'), st (manuscript_style)
%   keys         n_net x 1, the networks with at least one finished run, in
%                the figures' order: by group (rotation, local, blockmean),
%                then as the exploration sorted them (added structure by
%                strength, then rotation; local outliers before the
%                reference; real block-mean outliers before the loops)
%   group        n_net x 1, 1 rotation, 2 local, 3 blockmean
%   group_names  {'rotation', 'local', 'blockmean'}
%   short        n_net x 1 cellstr, a tick label (TeX), with the structural
%                eigenvalue in bulk radii (median over seeds)
%   med, lo, hi  structs: for every numeric column of T, n_net x n_cond
%                arrays of the median, minimum and maximum over the seeds
%                (NaN where a network has no finished run of a condition)
%   n_seed       n_net x n_cond, the finished runs behind each statistic
%   row_of       @(w, c) the rows of T of network w and condition c
%   seeds_text   '1 seed(s)', or '3 to 5 seeds' when networks differ
%
% Errors <caller>:NoRuns when the stage finished no run.
%
% See also: ostruct_load, ostruct_write_table, fig_ostruct_summary,
%           fig_ostruct_plane, fig_ostruct_rotation, fig_ostruct_h1

F = struct();
F.file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'outlier_structures')}, 'ostruct_index.mat', ...
    'Run the opt-in stage outlier_structures (run_outlier_structures_stage) first');
F.folder = fileparts(F.file);
F.S = ostruct_load(F.folder);
T = F.S.T;
if isempty(T) || height(T) == 0
    error([caller ':NoRuns'], 'The outlier_structures stage in %s finished no run (its index lists why).', F.folder);
end
F.T = T;

plan_conds = F.S.settings.conditions;
available = plan_conds(ismember(plan_conds, T.condition));
F.cond_names = select_conditions(available, cfg.conditions);
F.st = manuscript_style();
F.colors = cellfun(@(c) F.st.condition_color(c), F.cond_names, 'UniformOutput', false);
F.titles = grouped_condition_label(F.cond_names, 'sentence');

%% The networks: group and order
F.group_names = {'rotation', 'local', 'blockmean'};
nets = F.S.networks;
keys = unique(T.key, 'stable');
n_net = numel(keys);
group = zeros(n_net, 1); sort_by = zeros(n_net, 1); short = cell(n_net, 1);
for w = 1:n_net
    nw = nets(strcmp({nets.key}, keys{w}));
    rows = strcmp(T.key, keys{w});
    lre = median(T.lam0_re_R(rows)); lim = median(T.lam0_im_R(rows)); labs = median(T.lam0_abs_R(rows));
    lam = sprintf('%.2f', lre);
    if abs(lim) > 1e-6; lam = sprintf('%.2f \\pm %.2fi', lre, lim); end
    group(w) = find(strcmp(F.group_names, nw.group));
    switch nw.group
        case 'rotation'
            sort_by(w) = 100 * (lre > 2) + lim;               % the 1.4 R sweep by rotation, then the stronger ones
            short{w} = sprintf('added, (%s) R', lam);
        case 'local'
            no_outlier = isnumeric(nw.expect) && isscalar(nw.expect) && isnan(nw.expect);
            if no_outlier
                sort_by(w) = Inf;                             % the reference comparison last
                short{w} = nw.label;
            else
                sort_by(w) = -labs;                           % the clearest local outlier first
                short{w} = sprintf('local, (%s) R', lam);
            end
        otherwise
            if contains(nw.name, '-')
                sort_by(w) = 1 + lim;                         % the loops after the single routes
                short{w} = sprintf('E-I loop, (%s) R', lam);
            else
                sort_by(w) = 0;
                short{w} = nw.label;
            end
    end
end
[~, order] = sortrows([group, sort_by, (1:n_net)']);
F.keys = keys(order); F.group = group(order); F.short = short(order);

%% Statistics over seeds
n_c = numel(F.cond_names);
F.row_of = @(w, c) find(strcmp(T.key, F.keys{w}) & strcmp(T.condition, F.cond_names{c}));
cols = T.Properties.VariableNames;
cols = cols(cellfun(@(v) isnumeric(T.(v)) || islogical(T.(v)), cols));
F.med = struct(); F.lo = struct(); F.hi = struct();
F.n_seed = zeros(n_net, n_c);
for q = 1:numel(cols)
    F.med.(cols{q}) = nan(n_net, n_c); F.lo.(cols{q}) = nan(n_net, n_c); F.hi.(cols{q}) = nan(n_net, n_c);
end
for w = 1:n_net
    for c = 1:n_c
        k = F.row_of(w, c);
        F.n_seed(w, c) = numel(k);
        if isempty(k); continue; end
        for q = 1:numel(cols)
            v = double(T.(cols{q})(k));
            v = v(isfinite(v));
            if isempty(v); continue; end
            F.med.(cols{q})(w, c) = median(v);
            F.lo.(cols{q})(w, c) = min(v);
            F.hi.(cols{q})(w, c) = max(v);
        end
    end
end
n = unique(F.n_seed(F.n_seed > 0));
if isscalar(n)
    F.seeds_text = sprintf('%d seed(s)', n);
else
    F.seeds_text = sprintf('%d to %d seeds', min(n), max(n));
end
end
