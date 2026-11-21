function out = fig_eigenvalue_density(opts)
%FIG_EIGENVALUE_DENSITY The figures and the table of the eigenvalue-density exploration.
%
%   out = fig_eigenvalue_density()
%   out = fig_eigenvalue_density('data_file', f, 'out_dir', d, 'network', 'reference')
%
% Reads eigenvalue_density.mat (explore_eigenvalue_density) and writes, to
% figs/explorations/eigenvalue_density, four figures of ONE network
% (default the reference network; its two seeds pooled) and one table of
% every run. Each PNG has a .md sidecar that says what is drawn; nothing
% is written on the panels beyond axis labels and one title line.
%
%   eigenvalue_density_radial    the answer: how the eigenvalues of the
%                                frozen-gain connectivity are distributed
%                                in radius, by condition, with the theory
%   eigenvalue_density_gains     the gain distributions behind it
%   eigenvalue_density_plane     one spectrum per condition in the plane
%   eigenvalue_density_measures  the expressivity measures by condition
%   eigenvalue_density_table.md  every run, every measure
%
% Colors are manuscript_style's condition colors, with no adaptation in
% dark red [0.75 0.10 0.15] (FR-004).
%
% See also: explore_eigenvalue_density, edens_measures

arguments
    opts.data_file (1,:) char = ''
    opts.out_dir   (1,:) char = ''
    opts.network   (1,:) char = 'reference'
    opts.visible   (1,1) logical = false
end
root = fileparts(which('setup_paths'));
if isempty(opts.data_file)
    opts.data_file = fullfile(root, 'data', 'explorations', 'eigenvalue_density', 'eigenvalue_density.mat');
end
if isempty(opts.out_dir)
    opts.out_dir = fullfile(root, 'figs', 'explorations', 'eigenvalue_density');
end
if ~isfolder(opts.out_dir); mkdir(opts.out_dir); end
L = load(opts.data_file, 'E');
E = L.E;
st = manuscript_style();
conds = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
colors = cellfun(@(c) st.condition_color(c), conds, 'UniformOutput', false);
colors{1} = [0.75 0.10 0.15];
titles = cellfun(@(c) st.condition_title(c), conds, 'UniformOutput', false);
vis = matlab.lang.OnOffSwitchState(opts.visible);
u = E.settings.u; u2e = E.settings.u2_edges; u2c = (u2e(1:end - 1) + u2e(2:end)) / 2;
S = E.settings.measures;
uniform_shell = 1 - S.shell_r^2;

runs = E.runs(strcmp({E.runs.network}, opts.network));
if isempty(runs)
    error('fig_eigenvalue_density:NoRuns', 'No runs of network ''%s'' in %s.', opts.network, opts.data_file);
end
by_cond = cellfun(@(c) runs(strcmp({runs.condition}, c)), conds, 'UniformOutput', false);
pool = @(rs, f) mean(cat(1, rs.(f)), 1);
suffix = '';
if ~strcmp(opts.network, 'reference'); suffix = ['_' opts.network]; end
ax_style = {'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out'};

%% 1. Radial distribution
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 26 11]);
c1 = onCleanup(@() close(fig));
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, u, min(u.^2, 1), ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8);
for c = 1:3
    plot(ax, u, pool(by_cond{c}, 'cdf_M'), '-', 'Color', colors{c}, 'LineWidth', 3.5);
end
for c = 1:3
    plot(ax, u, pool(by_cond{c}, 'cdf_theory'), '--', 'Color', 'k', 'LineWidth', 1.1);
end
set(ax, ax_style{:}, 'XLim', [0 1.2], 'YLim', [0 1], 'XTick', 0:0.2:1.2, 'YTick', 0:0.2:1);
xlabel(ax, 'r / R_{eff}', 'FontSize', st.label_fs);
ylabel(ax, 'fraction of eigenvalues within r', 'FontSize', st.label_fs);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, [0 1.44], [1 1], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8);
inside = u2c < 1;
for c = 1:3
    stairs(ax, u2e, [pool(by_cond{c}, 'dens_M'), NaN], '-', 'Color', colors{c}, 'LineWidth', 2.5);
end
for c = 1:3
    th = pool(by_cond{c}, 'dens_theory');
    plot(ax, u2c(inside), th(inside), '--', 'Color', 'k', 'LineWidth', 1.1);
end
set(ax, ax_style{:}, 'XLim', [0 1.2], 'YScale', 'log', 'YLim', [0.03 5], 'XTick', 0:0.2:1.2, ...
    'YTick', [0.03 0.1 0.3 1 3], 'YTickLabel', {'0.03', '0.1', '0.3', '1', '3'});
xlabel(ax, '(r / R_{eff})^2', 'FontSize', st.label_fs);
ylabel(ax, 'density / uniform disk', 'FontSize', st.label_fs);
tag1 = ['eigenvalue_density_radial' suffix];
exportgraphics(fig, fullfile(opts.out_dir, [tag1 '.png']), 'Resolution', 200);

%% 2. Gains
fig2 = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 26 11]);
c2 = onCleanup(@() close(fig2));
tl = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
for c = 1:3
    g = sort(reshape(cat(2, by_cond{c}.gains_cal), [], 1));
    stairs(ax, [0; g; 1], [0; (1:numel(g))' / numel(g); 1], '-', 'Color', colors{c}, 'LineWidth', 2.5);
end
set(ax, ax_style{:}, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
xlabel(ax, 'gain d', 'FontSize', st.label_fs);
ylabel(ax, 'fraction of neurons with gain \leq d', 'FontSize', st.label_fs);
ax = nexttile(tl); hold(ax, 'on');
for c = 1:3
    for k = 1:numel(by_cond{c})
        r = by_cond{c}(k);
        plot(ax, 1 - r.p0_t, r.pr_t, '.', 'Color', colors{c}, 'MarkerSize', 9);
    end
end
plot(ax, [0 1], [0 1], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8);
set(ax, ax_style{:}, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
xlabel(ax, 'fraction of neurons with nonzero gain', 'FontSize', st.label_fs);
ylabel(ax, '\langle d^2 \rangle^2 / \langle d^4 \rangle', 'FontSize', st.label_fs);
tag2 = ['eigenvalue_density_gains' suffix];
exportgraphics(fig2, fullfile(opts.out_dir, [tag2 '.png']), 'Resolution', 200);

%% 3. One spectrum per condition
fig3 = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 30 10.5]);
c3 = onCleanup(@() close(fig3));
tl = tiledlayout(fig3, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% Wide enough for the bulk and for the boundary at Re = 1.
lim = max(1.15 * max(arrayfun(@(c) by_cond{c}(end).R_eff_example, 1:3)), 1.15);
tick_step = 1; if lim < 1.4; tick_step = 0.5; end
th = linspace(0, 2 * pi, 361);
for c = 1:3
    r = by_cond{c}(end);
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, [1 1], [-lim lim], '-', 'Color', st.zeroline_color, 'LineWidth', 1.2);
    plot(ax, r.R_eff_example * cos(th), r.R_eff_example * sin(th), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
    plot(ax, real(r.ev_example), imag(r.ev_example), '.', 'Color', colors{c}, 'MarkerSize', 8);
    axis(ax, 'equal');
    set(ax, ax_style{:}, 'XLim', [-lim lim], 'YLim', [-lim lim], 'XTick', -2:tick_step:2, 'YTick', -2:tick_step:2);
    title(ax, titles{c}, 'FontSize', st.title_fs, 'FontWeight', 'normal');
    xlabel(ax, 'Re \lambda(W D)', 'FontSize', st.label_fs);
    if c == 1; ylabel(ax, 'Im \lambda(W D)', 'FontSize', st.label_fs); end
end
tag3 = ['eigenvalue_density_plane' suffix];
exportgraphics(fig3, fullfile(opts.out_dir, [tag3 '.png']), 'Resolution', 200);

%% 4. Measures
fig4 = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 34 10.5]);
c4 = onCleanup(@() close(fig4));
tl = tiledlayout(fig4, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
short = cellfun(@(c) st.condition_short(c), conds, 'UniformOutput', false);
gray = [0.45 0.45 0.45];
panels = { ...
    @(r) mean(r.m_M(1, :)) / uniform_shell,  @(r) mean(r.m_theory(1, :)) / uniform_shell, 'rim occupancy / uniform disk'; ...
    @(r) mean(r.m_M(3, :)),                  @(r) mean(r.m_uniform(3, :)),                 'eigenvalues, -2 to 0 s^{-1}'; ...
    @(r) mean(r.m_M(4, :)),                  @(r) mean(r.full.n_unstable),                 'eigenvalues above 0'; ...
    @(r) mean(r.full.n_band_x),              @(r) NaN,                                     'full Jacobian, -2 to 0 s^{-1}'};
for q = 1:size(panels, 1)
    ax = nexttile(tl); hold(ax, 'on');
    if q == 1; plot(ax, [0.4 3.6], [1 1], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8); end
    for c = 1:3
        rs = by_cond{c};
        x = c + linspace(-0.12, 0.12, numel(rs));
        if numel(rs) == 1; x = c; end
        for k = 1:numel(rs)
            plot(ax, x(k), panels{q, 2}(rs(k)), 'd', 'Color', gray, 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.3);
            plot(ax, x(k), panels{q, 1}(rs(k)), 'o', 'Color', colors{c}, 'MarkerFaceColor', colors{c}, 'MarkerSize', 8);
        end
    end
    yl = ylim(ax);
    set(ax, ax_style{:}, 'XLim', [0.4 3.6], 'XTick', 1:3, 'XTickLabel', short, 'YLim', [0 yl(2) * 1.05], ...
        'XTickLabelRotation', 30);
    ylabel(ax, panels{q, 3}, 'FontSize', st.label_fs);
end
tag4 = ['eigenvalue_density_measures' suffix];
exportgraphics(fig4, fullfile(opts.out_dir, [tag4 '.png']), 'Resolution', 200);

%% Sidecars and the table
net_note = sprintf('Network `%s`, %d runs (seeds pooled), 12 states per run over 20 s.', opts.network, numel(runs));
write_md(opts.out_dir, tag1, [ ...
    'How the eigenvalues of the frozen-gain effective connectivity M = W D are spread in radius. ' net_note ' ' ...
    'Radii are in units of the predicted bulk radius R_eff at each state. ' ...
    'LEFT: fraction of eigenvalues within radius r. Thick colored lines: measured (dark red: no adaptation; orange: one ' ...
    'timescale; blue: chained multiple timescales). Thin dashed black: the random-matrix theory evaluated with the ' ...
    'gains measured at the same states (`edens_radial_theory_blocks`). Dotted gray: a uniform disk, (r/R)^2, which is ' ...
    'what the weights alone give in this network (Harris et al. 2023, Eq. 20: g = 0, equal E and I variances). The ' ...
    'value at r = 0 is the fraction of neurons with zero gain: their eigenvalues are exactly at the origin. A curve ' ...
    'above the dotted line is a spectrum concentrated toward the center. ' ...
    'RIGHT: the density relative to a uniform disk of the same radius, in equal-area bins of (r/R_eff)^2, the ' ...
    'eigenvalues exactly at the origin left out (log axis). Steps: measured; dashed black: theory; dotted gray: uniform. ' ...
    'Beyond 1 are the eigenvalues outside the predicted bulk (finite size, and the local outliers of a matrix with ' ...
    'E and I means). Numbers: `eigenvalue_density_table.md`.']);
write_md(opts.out_dir, tag2, [ ...
    'The gains behind the densities. ' net_note ' ' ...
    'LEFT: cumulative distribution of the per-neuron gain d (slope of the nonlinearity times the available vesicle ' ...
    'fraction) over neurons and states; the step at d = 0 is the fraction of neurons pinned at the floor or the ' ...
    'ceiling of the nonlinearity. Colors as in the radial figure. Without adaptation the gains are nearly binary (0 or ' ...
    'close to 1); with adaptation few are zero and the rest are graded by depression. ' ...
    'RIGHT: one dot per state (every 0.2 s): the participation ratio of the squared gains, <d^2>^2 / <d^4>, against the ' ...
    'fraction of neurons with nonzero gain. The participation ratio is the density at the rim of the bulk relative to a ' ...
    'uniform disk (Ahmadian et al. 2015, Eq. 2.18). Dotted: the two are equal when the nonzero gains are all the same; ' ...
    'the distance below it is what the spread of the nonzero gains costs.']);
write_md(opts.out_dir, tag3, [ ...
    'One spectrum of M = W D per condition (the middle of the 12 states of the last seed of network `' opts.network '`), ' ...
    'in units of W. Gray circle: the predicted bulk radius R_eff from the gains at that state. Green line: Re = 1, the ' ...
    'stability boundary of the dendritic block ((lambda - 1)/tau_d = 0). Zero-gain neurons contribute eigenvalues ' ...
    'exactly at the origin, drawn on top of each other.']);
write_md(opts.out_dir, tag4, [ ...
    'The expressivity measures by condition; one pair of markers per seed. ' net_note ' Definitions: `edens_measures`. ' ...
    'PANEL 1: fraction of the eigenvalues of M beyond 0.8 R_eff, divided by the 0.36 of a uniform disk. Filled: ' ...
    'measured; open diamond: theory with the measured gains; dotted: uniform. ' ...
    'PANEL 2: number of eigenvalues of the frozen-gain dendritic block with rate between -2 and 0 1/s (of 500). ' ...
    'Filled: measured; open diamond: a uniform disk of the same radius. ' ...
    'PANEL 3: number with positive rate. Filled: frozen gains; open diamond: the full Jacobian at the same states (slow ' ...
    'variables free). ' ...
    'PANEL 4: number of eigenvalues of the FULL Jacobian between -2 and 0 1/s whose eigenvector lies mostly (at least ' ...
    'half its squared norm) in the dendritic variables.']);

names = runs(1).measure_names;
fid = fopen(fullfile(opts.out_dir, 'eigenvalue_density_table.md'), 'w');
fprintf(fid, '# eigenvalue_density_table\n\n');
fprintf(fid, ['Every saved run of the outliers exploration (500 neurons, 12 states over 20 s each), means over the states. ' ...
    'Definitions in `scripts/explorations/eigenvalue_density/edens_measures.m`.\n\n' ...
    '- `p0`: fraction of neurons with zero gain. `cv`: coefficient of variation of the nonzero gains. `PR`: <d^2>^2/<d^4>.\n' ...
    '- `R_eff`: predicted bulk radius of M = W D (units of W; the boundary is Re = 1).\n' ...
    '- `rim`: fraction of eigenvalues beyond %.1f R_eff divided by the uniform disk''s %.2f: measured (M), the null draws ' ...
    '(same gains, shuffled mean-free W), the theory with the measured gains, and the actual W under a uniform gain.\n' ...
    '- `r50`: median modulus / R_eff (uniform disk: 0.707).\n' ...
    '- `band`: eigenvalues of the frozen-gain block with rate in (-%g, 0] 1/s; `unst`: with rate > 0. M: measured; th: theory with the ' ...
    'measured gains; uni: a uniform disk of radius R_eff.\n' ...
    '- `full`: the full Jacobian at the same states: all eigenvalues in the band, those mostly in the dendritic variables (x), and ' ...
    'the unstable ones.\n' ...
    '- `ev`, `mom`, `unseen`: checks of the gain recovery (largest error in the reproduced spectrum and in the gain moments by ' ...
    'type; neurons whose gain could not be recovered).\n\n'], S.shell_r, uniform_shell, S.band / S.tau_d);
fprintf(fid, ['| network | condition | seed | p0 | cv | PR | R_eff | rim M | rim null | rim th | rim unigain | r50 M | r50 th | ' ...
    'band M | band th | band uni | unst M | unst th | unst uni | full band | full band x | full unst | ev | mom | unseen |\n']);
fprintf(fid, '|%s\n', repmat('---|', 1, 25));
all_runs = E.runs;
for k = 1:numel(all_runs)
    r = all_runs(k);
    m = @(M, name) mean(M(strcmp(names, name), :));
    fprintf(fid, ['| %s | %s | %d | %.3f | %.2f | %.3f | %.2f | %.3f | %.3f | %.3f | %.2f | %.3f | %.3f | %.1f | %.1f | %.1f | ' ...
        '%.1f | %.1f | %.1f | %.0f | %.1f | %.1f | %.0e | %.0e | %d |\n'], ...
        r.network, r.condition, r.seed, mean(r.p0), mean(r.cv_nz), mean(r.pr), mean(r.R_eff), ...
        m(r.m_M, 'shell') / uniform_shell, m(r.m_null, 'shell') / uniform_shell, m(r.m_theory, 'shell') / uniform_shell, ...
        m(r.m_unigain, 'shell') / uniform_shell, m(r.m_M, 'r50'), m(r.m_theory, 'r50'), ...
        m(r.m_M, 'n_band'), m(r.m_theory, 'n_band'), m(r.m_uniform, 'n_band'), ...
        m(r.m_M, 'n_unstable'), m(r.m_theory, 'n_unstable'), m(r.m_uniform, 'n_unstable'), ...
        mean(r.full.n_band), mean(r.full.n_band_x), mean(r.full.n_unstable), ...
        r.check_ev, r.check_moments_cal, r.n_unseen);
end
if ~isempty(E.switching)
    sw = E.switching;
    fprintf(fid, ['\n## The switching run of Figure 3 (full Jacobian only)\n\n' ...
        'Eigenvalues of the full Jacobian at 40 states per condition, per state. The per-neuron gains of this run were not ' ...
        'saved, so its frozen-gain density cannot be formed; `uniform` is what a uniform disk of the rms-gain radius would ' ...
        'put in the same places.\n\n' ...
        '| condition | state dimension | rms gain | in (-%g, 0] 1/s | uniform disk | above 0 | uniform disk |\n|---|---|---|---|---|---|---|\n'], S.band / S.tau_d);
    for c = 1:numel(sw.cond_names)
        fprintf(fid, '| %s | %d | %.3f | %.1f | %.1f | %.2f | %.1f |\n', sw.cond_names{c}, sw.N(c), sw.rms_gain(c), ...
            sw.n_band(c), sw.uniform_band(c), sw.n_unstable(c), sw.uniform_unstable(c));
    end
end
if ~isempty(E.harris)
    fprintf(fid, ['\n## W alone (Harris et al. 2023)\n\nN sigma_se^2 = %.4f, N sigma_si^2 = %.4f, g = %.3g (Eq. 20), R = %.4f ' ...
        '(Eq. 18; the class''s R is %.4f): the density of W itself is uniform.\n'], ...
        E.harris.v_e, E.harris.v_i, E.harris.g, E.harris.R, E.harris.R_class);
end
fclose(fid);

out = struct('out_dir', opts.out_dir, 'tags', {{tag1, tag2, tag3, tag4}}, ...
    'table', fullfile(opts.out_dir, 'eigenvalue_density_table.md'));
end

function write_md(out_dir, tag, text)
fid = fopen(fullfile(out_dir, [tag '.md']), 'w');
fprintf(fid, '# %s\n\n%s\n', tag, text);
fclose(fid);
end
