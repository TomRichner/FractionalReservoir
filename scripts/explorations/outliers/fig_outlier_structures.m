function files = fig_outlier_structures(data_dir, fig_dir, opts)
%FIG_OUTLIER_STRUCTURES The figures of explore_outlier_structures, each with a sidecar .md that carries the numbers and a caption.
%
%   files = fig_outlier_structures()
%   files = fig_outlier_structures(data_dir, fig_dir, 'keep_open', false)
%
% Reads every run in data_dir (outlier_mode_table) and writes, with nothing
% on the panels but axis labels and one title line:
%
%   structures_plane_rotation, structures_plane_local,
%   structures_plane_blockmean
%                    the complex plane, rows: networks, columns:
%                    conditions. The whole frozen-gain spectrum at the
%                    calibration states as rates, the null-calibrated edge
%                    of the bulk (gray circle), the INTENDED mode ringed
%                    (followed by its eigenvector; with its conjugate when
%                    complex), where the weights alone put it (open
%                    diamond: uniform gains), zero (dashed).
%   structures_rotation   the rotation sweep: against the mode's rotation
%                    rate, the fraction of its structural detachment the
%                    gains removed, its real part frozen and in the full
%                    Jacobian, and its rotation rate in the full Jacobian.
%   structures_h1    every network and condition: the observed relative
%                    shift of the mode's eigenvalue from the uniform-gain
%                    picture against the first-order prediction from the
%                    gains of its carriers (real and imaginary parts).
%   structures_summary    every network: detachment of the intended mode
%                    (with the uniform-gain value), the responsive
%                    fractions, and the mode's rate frozen and full.
%
% The networks are grouped by what made the intended eigenvalue
% (outlier_network's target and the name): added low-rank structure,
% sorted by rotation; local outliers and the reference; block means.
% Figures are made invisible and closed unless keep_open is true (then
% they stay open, tagged 'fig_outlier_structures'). Colors:
% outlier_cond_color.
%
% See also: explore_outlier_structures, outlier_mode_table, outlier_mode_measures, fig_outliers

arguments
    data_dir (1,:) char = ''
    fig_dir  (1,:) char = ''
    opts.keep_open (1,1) logical = false
    opts.resolution (1,1) double = 200
end

root = fileparts(which('setup_paths'));
if isempty(data_dir); data_dir = fullfile(root, 'data', 'explorations', 'outlier_structures'); end
if isempty(fig_dir);  fig_dir  = fullfile(root, 'figs', 'explorations', 'outlier_structures'); end
if ~isfolder(fig_dir); mkdir(fig_dir); end
delete(findall(groot, 'Type', 'figure', 'Tag', 'fig_outlier_structures'));   % this function's own, nothing else

[T, M, MM, Rs] = outlier_mode_table(data_dir, 'print', false, 'keep_runs', true);
st = manuscript_style();
cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
conds = cond_order(ismember(cond_order, T.condition));
n_cond = numel(conds);
gray = [0.55 0.55 0.55];
fs_t = st.tick_fs - 2; fs_l = st.label_fs - 2;
vis = 'off'; if opts.keep_open; vis = 'on'; end
files = {};

% The networks (a name and a seed), their group and their order.
T.key = strcat(T.network, '_s', arrayfun(@(s) sprintf('%d', s), T.seed, 'UniformOutput', false));
[keys, first] = unique(T.key, 'stable');
group = zeros(numel(keys), 1); sort_by = zeros(numel(keys), 1);
for w = 1:numel(keys)
    r = T(first(w), :);
    if strcmp(r.target{1}, 'lowrank')
        group(w) = 1; sort_by(w) = 100 * (r.lam0_re_R > 2) + r.lam0_im_R;   % the weaker sweep, then the stronger pair
    elseif strcmp(r.target{1}, 'rightmost')
        group(w) = 2; sort_by(w) = -r.lam0_abs_R;          % the clearest local outlier first, the reference last
    else
        group(w) = 3; sort_by(w) = abs(r.lam0_im_R) > 1e-6;   % real block-mean outliers, then the loops
    end
end
[~, order] = sortrows([group, sort_by]);
keys = keys(order); group = group(order); first = first(order);
n_net = numel(keys);
row_of = @(w, c) find(strcmp(T.key, keys{w}) & strcmp(T.condition, conds{c}), 1);
label = cell(n_net, 1); short = cell(n_net, 1);
for w = 1:n_net
    [label{w}, short{w}] = net_names(T(first(w), :));
end
ccol = cellfun(@(c) outlier_cond_color(st, c), conds, 'UniformOutput', false);

table_cols = {'lam0_re_R', 'lam0_im_R', 'coh0', 'cond0', 'lambda_1', 'mode_rho', 'mode_rho_uniform', 'mode_reined', ...
    'mode_ov', 'mode_frac_rightmost', 'rate_struct', 'rate_frozen', 'rate_full', 'rot_struct', 'rot_frozen', 'rot_full', ...
    'frac_on_E', 'frac_on_I', 'on_top', 'on_rest'};

%% 1. The complex plane, one figure per group
group_tag = {'rotation', 'local', 'blockmean'};
group_caption = { ...
    ['Structure ADDED to the reference weight matrix on random patterns of neurons (dense, ignoring Dale''s law): a real ' ...
    'eigenvalue (rank one, top row) and conjugate pairs of increasing imaginary part (rank two) with a real part ' ...
    'of about 1.4 structural bulk radii, then a stronger real one and a stronger pair at about 2.5 bulk radii (the row ' ...
    'titles give the structural eigenvalue in bulk radii).'], ...
    ['Networks with NO structure put in. Top rows: seeds chosen for a local outlier of the balanced sparse Dale matrix ' ...
    '(real in the reference network on seed 27; complex at indegree 200 on seed 25). Last row: the reference network on ' ...
    'seed 7, where the intended mode is simply the rightmost eigenvalue of W, a bulk-edge mode.'], ...
    ['Outliers made by the block means: single routes changed (real outliers), and the E-I loop strengthened together ' ...
    'with E-to-E (a complex pair).']};
for g = 1:3
    ws_ = find(group == g);
    if isempty(ws_); continue; end
    fig = new_fig([1 1 12 3.5 * numel(ws_)], vis);
    tl = tiledlayout(fig, numel(ws_), n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
    for q = 1:numel(ws_)
        w = ws_(q);
        axs = gobjects(1, n_cond); xr = [-22 5]; rad_max = 0;
        for c = 1:n_cond
            ax = nexttile(tl); hold(ax, 'on'); axs(c) = ax;
            k = row_of(w, c);
            if ~isempty(k)
                R = Rs{k}; m = M{k}; tau = R.tau_d;
                ev = ([R.cal.ev_M] - 1) ./ tau;                % n x n_cal, rates
                th = linspace(0, 2 * pi, 200);
                rad = m.kappa * median(R.R_eff) / tau;
                rad_max = max(rad_max, rad);
                plot(ax, -1 / tau + rad * cos(th), rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
                scatter(ax, real(ev(:)), imag(ev(:)), 5, ccol{c}, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
                md = outlier_with_conjugates(([R.dense.mode_lam_M] - 1) ./ tau);
                plot(ax, real(md), imag(md), 'o', 'MarkerSize', 6, 'MarkerFaceColor', ccol{c}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
                u = outlier_with_conjugates((R.ref.lam0 * T.d_rms(k) - 1) / tau);
                plot(ax, real(u), imag(u), 'd', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
                xline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
                xr = [min(xr(1), -1 / tau - 1.15 * rad), max([xr(2), max(real(u)) + 3, max(real(md)) + 2])];
            end
            set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
            if q == 1; title(ax, st.condition_title(conds{c}), 'FontSize', st.title_fs, 'FontWeight', 'normal'); end
            if c == 1; ylabel(ax, sprintf('%s\nIm (s^{-1})', label{w}), 'FontSize', fs_l); end
            if q == numel(ws_); xlabel(ax, 'Re (s^{-1})', 'FontSize', fs_l); end
        end
        for c = 1:n_cond
            axis(axs(c), 'equal'); xlim(axs(c), xr); ylim(axs(c), [-1 1] * max(0.40 * diff(xr), 1.08 * rad_max));
        end
    end
    name = ['structures_plane_' group_tag{g}];
    files{end + 1} = save_fig(fig, fig_dir, name, opts); %#ok<AGROW>
    write_sidecar(fig_dir, name, [group_caption{g} ' Rows are the networks, columns the conditions. Dots: every eigenvalue ' ...
        'of the frozen-gain effective connectivity (each neuron''s outgoing weights times its current gain, slow variables ' ...
        'held), drawn as a rate (lambda - 1)/tau_d, at 8 states spread over 20 s after 60 s of settling, overlaid. Gray ' ...
        'circle: the null-calibrated edge of the bulk (median over the window), centered on -1/tau_d = -10 1/s. Ringed ' ...
        'dots: the INTENDED mode at each of the 8 states, followed by its eigenvector (the eigenvector of the frozen-gain ' ...
        'matrix that overlaps most with the structural one; mode_ov is the median overlap, and where it is low the mode has ' ...
        'not survived as one eigenvector and the ring marks only the best match), with its conjugate when it is complex. ' ...
        'Open diamonds: where the same mode would be if every neuron had the same gain (the rms of the observed gains, so ' ...
        'the bulk keeps its size): what the weights alone give. Dashed: zero. mode_frac_rightmost is the fraction of states ' ...
        'at which the intended mode is also the rightmost eigenvalue. The view is not extended to the left: where the ' ...
        'followed mode has gone far left of the bulk (a population mode whose excitatory carriers are pinned becomes the ' ...
        'inhibitory, negative eigenvalue; rate_frozen of -25 to -60 1/s) its rings are outside the panel. One seed per network.'], ...
        T(ismember(T.key, keys(ws_)), :), table_cols);
end

%% 2. The rotation sweep
ws_ = find(group == 1);
if ~isempty(ws_)
    fig = new_fig([1 1 13 4.2], vis);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    ycol = {'mode_reined', 'rate_full', 'rot_full'};
    ylab = {'fraction of detachment removed', 'real part of the mode (s^{-1})', 'rotation rate, full Jacobian (rad s^{-1})'};
    for p = 1:3
        ax = nexttile(tl); hold(ax, 'on');
        for strong = [false true]                              % the sweep at 1.4 R (circles), the stronger pair (squares)
            mk = 'o'; if strong; mk = 's'; end
            for c = 1:n_cond
                ks = arrayfun(@(w) row_of(w, c), ws_, 'UniformOutput', false);
                ks = [ks{:}];
                ks = ks((T.lam0_re_R(ks) > 2) == strong);
                if isempty(ks); continue; end
                x = T.rot_struct(ks);
                [x, o] = sort(x); ks = ks(o);
                if p == 2
                    plot(ax, x, T.rate_frozen(ks), ['--' mk], 'Color', ccol{c}, 'LineWidth', 1.25, 'MarkerSize', 7, ...
                        'MarkerFaceColor', 'w', 'MarkerEdgeColor', ccol{c});
                end
                plot(ax, x, T.(ycol{p})(ks), ['-' mk], 'Color', ccol{c}, 'LineWidth', 2.0, 'MarkerSize', 8 + strong, ...
                    'MarkerFaceColor', ccol{c}, 'MarkerEdgeColor', 'w');
            end
        end
        if p == 1; yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5); end
        if p == 2; yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0); end
        if p == 3
            lim = [0, max([xlim(ax), ylim(ax)])];
            plot(ax, lim, lim, '-', 'Color', gray, 'LineWidth', 1.5);
        end
        set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        xlabel(ax, 'rotation rate of the mode, uniform gain (rad s^{-1})', 'FontSize', fs_l);
        ylabel(ax, ylab{p}, 'FontSize', fs_l);
    end
    files{end + 1} = save_fig(fig, fig_dir, 'structures_rotation', opts);
    write_sidecar(fig_dir, 'structures_rotation', ...
        ['The rotation sweep: added low-rank structure with a real part of about 1.4 structural bulk radii and an imaginary ' ...
        'part from 0 to about one bulk radius (circles; the rank-one, real case is at 0), and a stronger pair of networks ' ...
        'with a real part of about 2.5 bulk radii, real and with an imaginary part of one bulk radius (squares). Horizontal axis: the rotation rate of the ' ...
        'mode in the uniform-gain picture, Im(lambda_0) d_rms / tau_d, which differs between conditions because the rms ' ...
        'gain does. One line per condition. LEFT: the fraction of the mode''s structural detachment from the bulk edge ' ...
        'that the actual gains removed (1, gray line: back at the edge; more: inside the bulk). MIDDLE: the real part of the ' ...
        'mode with the slow variables frozen (open, dashed) and free, in the full Jacobian (filled, solid); dashed black: ' ...
        'zero. RIGHT: the rotation rate of the mode''s continuation in the full Jacobian; gray: equal to the uniform-gain ' ...
        'rate. The full-Jacobian values are medians over 8 states; the frozen ones over 201. One seed, [7 8].'], ...
        T(ismember(T.key, keys(ws_)), :), [table_cols, {'slow_shift', 'cont_ov', 'gain_ratio', 'mean_over_rms', 'shift_obs_re', 'shift_pred_re'}]);
end

%% 3. H1: observed shift against the first-order prediction
fig = new_fig([1 1 9.5 4.6], vis);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
marks = {'o', 's', '^'};
parts = {'re', 'im'};
plab = {'real part', 'imaginary part'};
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    xo = T.(['shift_pred_' parts{p}]); yo = T.(['shift_obs_' parts{p}]);
    lim = [min([xo; yo; 0]) - 0.05, max([xo; yo; 0]) + 0.05];
    plot(ax, lim, lim, '-', 'Color', gray, 'LineWidth', 1.5);
    xline(ax, 0, ':', 'Color', gray); yline(ax, 0, ':', 'Color', gray);
    for w = 1:n_net
        for c = 1:n_cond
            k = row_of(w, c);
            if isempty(k); continue; end
            plot(ax, xo(k), yo(k), marks{group(w)}, 'MarkerSize', 8, 'MarkerFaceColor', ccol{c}, ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.75);
        end
    end
    axis(ax, 'equal'); xlim(ax, lim); ylim(ax, lim);
    set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    xlabel(ax, sprintf('predicted shift, %s', plab{p}), 'FontSize', fs_l);
    ylabel(ax, sprintf('observed shift, %s', plab{p}), 'FontSize', fs_l);
end
files{end + 1} = save_fig(fig, fig_dir, 'structures_h1', opts);
write_sidecar(fig_dir, 'structures_h1', ...
    ['H1 across every network and condition. The shift of the intended mode''s eigenvalue from the uniform-gain picture, ' ...
    'lambda / (d_rms lambda_0) - 1 (complex; 0: the mode moved exactly as the bulk did; negative real part: it lost more ' ...
    'than the bulk), observed (vertical; the eigenvalue of the frozen-gain matrix followed by its eigenvector) against the ' ...
    'first-order prediction from the gains of the neurons that carry it, sum_j p_j d_j / d_rms - 1 with p_j the ' ...
    'participation of neuron j from the left and right structural eigenvectors (horizontal). LEFT: real parts. RIGHT: ' ...
    'imaginary parts. One point per network and condition (median over 201 states); color: condition; circles: added ' ...
    'low-rank structure; squares: local outliers and the reference (no structure put in); triangles: block-mean ' ...
    'outliers. Gray: equality. mean_over_rms is the mean gain over the rms gain (what a mode spread evenly over all ' ...
    'neurons would lose against the bulk from the spread of the gains alone); gain_ratio is the participation-weighted ' ...
    'gain over the mean gain (what the carriers lose beyond the average neuron). mode_ov is the overlap of the followed ' ...
    'eigenvector with the structural one: where it is low the observed value is only the best match.'], ...
    T, {'lam0_re_R', 'lam0_im_R', 'cond0', 'p_E', 'mode_ov', 'shift_obs_re', 'shift_pred_re', 'shift_obs_im', 'shift_pred_im', ...
    'mean_over_rms', 'gain_ratio', 'gain_top_rest', 'on_top', 'on_rest', 'frac_on_E', 'frac_on_I', 'gain_ratio_EI'});

%% 4. The summary across networks
fig = new_fig([1 1 max(9, 1.15 * n_net + 2) 9.5], vis);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
off = linspace(-0.25, 0.25, n_cond);
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    for w = 1:n_net
        for c = 1:n_cond
            k = row_of(w, c);
            if isempty(k); continue; end
            x = w + off(c);
            if p == 1
                rho = sort(MM{k}.rho(isfinite(MM{k}.rho)));
                if ~isempty(rho)
                    lo = max(rho(max(1, round(0.05 * numel(rho)))), 0.2);
                    plot(ax, [x x], [lo, rho(round(0.95 * numel(rho)))], '-', 'Color', ccol{c}, 'LineWidth', 2.0);
                end
                plot(ax, x, max(T.mode_rho(k), 0.2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', ccol{c}, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                plot(ax, x, T.mode_rho_uniform(k), 'd', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
            elseif p == 2
                plot(ax, x, T.frac_on_E(k), 'o', 'MarkerSize', 8, 'MarkerFaceColor', ccol{c}, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                plot(ax, x, T.frac_on_I(k), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ccol{c}, 'LineWidth', 1.5);
            else
                plot(ax, [x x], [T.rate_frozen(k) T.rate_full(k)], '-', 'Color', ccol{c}, 'LineWidth', 1.0);
                plot(ax, x, T.rate_frozen(k), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ccol{c}, 'LineWidth', 1.5);
                plot(ax, x, T.rate_full(k), 'o', 'MarkerSize', 8, 'MarkerFaceColor', ccol{c}, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
            end
        end
    end
    xline(ax, find(diff(group) ~= 0) + 0.5, ':', 'Color', gray, 'LineWidth', 1.0);
    if p == 1
        yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
        set(ax, 'YScale', 'log', 'YLim', [0.2 8], 'YTick', [0.25 0.5 1 2 4 8]);
        ylabel(ax, 'detachment  Re \lambda_{mode} / bulk edge', 'FontSize', fs_l);
    elseif p == 2
        ylim(ax, [0 1.05]);
        ylabel(ax, 'responsive fraction', 'FontSize', fs_l);
    else
        yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
        ylabel(ax, 'rate of the mode (s^{-1})', 'FontSize', fs_l);
    end
    set(ax, 'XTick', 1:n_net, 'XTickLabel', {}, 'XLim', [0.4, n_net + 0.6], 'FontSize', fs_t, ...
        'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    if p == 3; set(ax, 'XTickLabel', short, 'XTickLabelRotation', 35); end
end
files{end + 1} = save_fig(fig, fig_dir, 'structures_summary', opts);
write_sidecar(fig_dir, 'structures_summary', ...
    ['Every network, in three groups separated by dotted lines: added low-rank structure (by increasing rotation), local ' ...
    'outliers and the reference, block-mean outliers. Within a network the three conditions. TOP: detachment of the ' ...
    'INTENDED mode, the real part of its eigenvalue in the frozen-gain connectivity over the null-calibrated edge of the ' ...
    'bulk (dot: median over 20 s; bar: 5th to 95th percentile; log axis, values below 0.2 drawn at 0.2). Gray line: the ' ...
    'edge. Open diamonds: the same with every gain made equal, what the weights alone give. MIDDLE: the fraction of ' ...
    'neurons off the floor and ceiling of the nonlinearity, excitatory (filled) and inhibitory (open): detachment at or ' ...
    'below 1 with a low responsive fraction is a pinned network, not a controlled one. BOTTOM: the real part of the mode ' ...
    'with the slow variables frozen (open) and free, its continuation in the full Jacobian (filled; median over 8 ' ...
    'states); dashed: zero. right_rho and right_reined are the detachment and reined fraction of the RIGHTMOST eigenvalue ' ...
    '(outlier_measures), for the states where that is another mode. One seed per network.'], ...
    T, [table_cols, {'right_rho', 'right_rho_uniform', 'right_reined', 'full_abscissa', 'full_frac_unstable', 'slow_shift', ...
    'gain_ratio', 'mean_over_rms'}]);

if ~opts.keep_open
    delete(findall(groot, 'Type', 'figure', 'Tag', 'fig_outlier_structures'));
end
end

%% ------------------------------------------------------------------------
function fig = new_fig(pos, vis)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', pos, 'Tag', 'fig_outlier_structures', 'Visible', vis);
end

function [label, short] = net_names(r)
% A row title (two lines at most, TeX) and a tick label for one network.
net = outlier_network(r.network{1});
label = net.label;
lam = sprintf('%.2f', r.lam0_re_R);
if abs(r.lam0_im_R) > 1e-6; lam = sprintf('%.2f \\pm %.2fi', r.lam0_re_R, r.lam0_im_R); end
if strcmp(r.target{1}, 'lowrank')
    short = sprintf('added, (%s) R', lam);
    label = short;
elseif strcmp(r.target{1}, 'rightmost')
    if r.lam0_abs_R > 1.2
        short = sprintf('local, (%s) R', lam);
        label = sprintf('%s, seed %d', label, r.seed);
    else
        short = sprintf('%s, seed %d', label, r.seed);
        label = short;
    end
else
    tokens = strsplit(r.network{1}, '-');
    if numel(tokens) > 1
        short = sprintf('E-I loop, (%s) R', lam);
        label = short;
    else
        short = label;
    end
end
end

function file = save_fig(fig, fig_dir, name, opts)
file = fullfile(fig_dir, [name '.png']);
exportgraphics(fig, file, 'Resolution', opts.resolution);
fprintf('[fig_outlier_structures] wrote %s\n', file);
end

function write_sidecar(fig_dir, name, caption, T, cols)
% The caption paragraph and the numbers of the figure's runs, as markdown.
fid = fopen(fullfile(fig_dir, [name '.md']), 'w');
closer = onCleanup(@() fclose(fid));
fprintf(fid, '# %s\n\n%s\n\n', name, caption);
fprintf(fid, ['Rates are in 1/s, rotation rates in rad/s, structural eigenvalues in units of the structural bulk radius R. ' ...
    'Definitions: `scripts/explorations/outliers/outlier_mode_measures.m`; the networks: ' ...
    '`scripts/explorations/outliers/explore_outlier_structures.m`.\n\n']);
head = [{'network', 'seed', 'condition'}, cols];
fprintf(fid, '| %s |\n', strjoin(head, ' | '));
fprintf(fid, '|%s\n', repmat('---|', 1, numel(head)));
for k = 1:height(T)
    cells = cell(1, numel(cols));
    for q = 1:numel(cols)
        cells{q} = sprintf('%.3g', T.(cols{q})(k));
    end
    fprintf(fid, '| %s | %d | %s | %s |\n', T.network{k}, T.seed(k), T.condition{k}, strjoin(cells, ' | '));
end
end
