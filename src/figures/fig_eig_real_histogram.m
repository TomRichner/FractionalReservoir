function out = fig_eig_real_histogram(cfg)
% FIG_EIG_REAL_HISTOGRAM Jacobian eigenvalues projected on the real axis, against a uniform disk.
%
%   out = FIG_EIG_REAL_HISTOGRAM('run_dir', d)
%   out = FIG_EIG_REAL_HISTOGRAM('run_dir', d, 'disk_radius', 'unit', 'count_scale', 'clipped')
%
% One panel per condition: a histogram of Re(lambda) of the eigenvalues of
% the full Jacobian at the sampled states (the data fig_eig_heatmap's
% density pools; the reference example), in counts per state, and ONE
% dotted curve: what the histogram would be if the network's eigenvalues
% filled a disk uniformly, as those of a random matrix do (Tom, 2026-11-19).
%
% THE REFERENCE. The dendritic block of the Jacobian is (-I + W G) / tau_d,
% with G the diagonal matrix of frozen gains. For a random W the eigenvalues
% of W G fill a disk uniformly, so those of the block fill a disk of center
% c = -1/tau_d and some radius rho, n of them (the number of neurons). The
% projection of a uniform disk on the real axis is a semicircle: the
% expected count in a bin [x1, x2] is n (F(u2) - F(u1)), u = (x - c) / rho
% clipped to [-1, 1], F(u) = 1/2 + (u sqrt(1 - u^2) + asin(u)) / pi; for a
% narrow bin of width w that is n w (2 / (pi rho^2)) sqrt(rho^2 - (x - c)^2).
% It is integrated per bin here, so it sums to n over the disk exactly.
% The adaptation and depression variables add the other eigenvalues of the
% full Jacobian (N - n per state), which the curve does not describe: that
% contrast is part of what the figure shows.
%
% THE RADIUS, 'disk_radius':
%   'gain'  (default) rho = R * rms(d) / tau_d, with R the theoretical bulk
%           radius of W (the class's R, level_of_chaos included) and d_j the
%           frozen gain of neuron j at a sampled state (the slope of its
%           nonlinearity times its route's output factor;
%           SRNNCellTypePairs.frozen_gain), root mean square over neurons,
%           averaged over the sampled states. The bulk radius follows the
%           rms gain, so the comparison isolates the SHAPE of the
%           distribution from its size. Needs the gains in the run
%           (run_eig_heatmap with 'disk_gains' or the Ritz pass); without
%           them the figure says so and uses 'unit'.
%   'unit'  rho = R / tau_d: the same connectivity with every neuron at
%           unit slope and no depression. Needs nothing beyond the preset.
%   a number  that radius (1/s).
% Both radii are in the table whichever is drawn.
%
% THE COUNT AXIS, 'count_scale'. The slow variables put hundreds of
% eigenvalues per state into a few bins, which flattens the disk on a
% linear axis. 'log' (default): logarithmic counts per state, where the
% disk part, the dotted curve and the spikes can be read together.
% 'clipped': linear, cut at 1.5 times the reference curve's peak; bars that
% exceed it are drawn to the top and marked with a small triangle.
%
% x limits are fig_eig_heatmap's window (the same quantile rule over the
% drawn conditions), so the two figures read together. No text on the
% panels: the numbers are in <tag>_table.md, with a caption paragraph.
%
% See also: fig_eig_heatmap, run_eig_heatmap, SRNNCellTypePairs.frozen_gain

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.conditions                = {}
    cfg.disk_radius               = 'gain'   % 'gain' | 'unit' | a positive number (1/s)
    cfg.count_scale (1,:) char {mustBeMember(cfg.count_scale, {'log', 'clipped'})} = 'log'
    cfg.n_bins      (1,1) double {mustBeInteger, mustBePositive} = 200
end
setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st = manuscript_style();
data_file = resolve_data_file(cfg.data_file, cfg.run_dir, {fullfile(cfg.run_dir, 'eig_heatmap')}, ...
    'eig_heatmap_data.mat', 'Run run_eig_heatmap first');
D = load(data_file);
if isfield(D, 'cond_names'); held = D.cond_names; else; held = D.condition_titles; end
[cond_names, idx] = select_conditions(held, cfg.conditions);
titles = titles_by_name(cond_names, D.condition_titles(idx));
evals = D.evals_by_cond(idx);
n_cond = numel(cond_names);
radius_is_name = ischar(cfg.disk_radius) || isstring(cfg.disk_radius);
if radius_is_name
    mode = char(cfg.disk_radius);
    if ~ismember(mode, {'gain', 'unit'})
        error('fig_eig_real_histogram:BadRadius', 'disk_radius must be ''gain'', ''unit'' or a positive number.');
    end
elseif ~(isnumeric(cfg.disk_radius) && isscalar(cfg.disk_radius) && cfg.disk_radius > 0)
    error('fig_eig_real_histogram:BadRadius', 'disk_radius must be ''gain'', ''unit'' or a positive number.');
else
    mode = 'given';
end

% fig_eig_heatmap's window: square, holding 99.9% of the drawn eigenvalues.
keep_frac = 0.999;
all_evals = vertcat(evals{:});
tail = (1 - keep_frac) / 2;
re_lo = quantile(real(all_evals), tail); re_hi = quantile(real(all_evals), 1 - tail);
span = max(re_hi - re_lo, 2 * quantile(abs(imag(all_evals)), keep_frac));
re_lim = (re_lo + re_hi) / 2 + [-0.5, 0.5] * span;
edges = linspace(re_lim(1), re_lim(2), cfg.n_bins + 1);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
width = edges(2) - edges(1);

% The reference disk of each condition.
have_gains = isfield(D, 'disk_by_cond');
if strcmp(mode, 'gain') && ~have_gains
    vprintf(cfg.verbose, 'minimal', ['[fig_eig_real_histogram] %s holds no frozen gains (run_eig_heatmap with ' ...
        '''disk_gains'' or the Ritz pass); drawing the unit-gain disk instead.\n'], data_file);
    mode = 'unit';
end
ref = struct('n', cell(1, n_cond), 'N', [], 'n_states', [], 'c', [], 'rho', [], 'rho_unit', [], ...
    'rho_gain', [], 'rms_gain', [], 'expected', [], 'counts', [], 'mode', []);
for i = 1:n_cond
    if have_gains
        G = D.disk_by_cond{idx(i)};
        n = G.n; tau_d = G.tau_d; R = G.R;
        ref(i).rms_gain = mean(G.rms_gain);
        ref(i).rho_gain = R * ref(i).rms_gain / tau_d;
    else
        % The preset alone: R does not depend on the network size.
        m = build_from_preset(D.settings.preset_name, cond_names{i}, 'verbose', 'near-none');
        n = D.settings.n; tau_d = m.tau_d; R = m.R;
        ref(i).rms_gain = NaN; ref(i).rho_gain = NaN;
    end
    ref(i).n = n; ref(i).c = -1 / tau_d; ref(i).rho_unit = R / tau_d;
    switch mode
        case 'gain';  ref(i).rho = ref(i).rho_gain;
        case 'unit';  ref(i).rho = ref(i).rho_unit;
        case 'given'; ref(i).rho = cfg.disk_radius;
    end
    ref(i).mode = mode;
    ref(i).n_states = numel(D.J_times_by_cond{idx(i)});
    ref(i).N = numel(evals{i}) / ref(i).n_states;
    ref(i).expected = n * diff(semicircle_cdf((edges - ref(i).c) / ref(i).rho));
    ref(i).counts = histcounts(real(evals{i}), edges) / ref(i).n_states;
end

%% Draw
fig = figure('Color', 'w', 'Position', [150, 150, 430 * n_cond, 330]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 1, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
floor_count = 0.5 / max([ref.n_states]);            % half an eigenvalue in one state: the log axis's floor
top_log = 1.5 * max(arrayfun(@(r) max(r.counts), ref));
ax = gobjects(1, n_cond);
for i = 1:n_cond
    ax(i) = nexttile(tl); hold(ax(i), 'on');
    col = [0 0 0];
    if st.condition_color.isKey(cond_names{i}); col = st.condition_color(cond_names{i}); end
    c = ref(i).counts; e = ref(i).expected;
    xs = [edges(1:end - 1); edges(2:end)];
    switch cfg.count_scale
        case 'log'
            base = floor_count;
            cc = max(c, base);
            patch(ax(i), [xs(1); xs(:); xs(end)], [base; reshape([cc; cc], [], 1); base], col, ...
                'FaceAlpha', 0.55, 'EdgeColor', 'none', 'Tag', 'eig_real_histogram');
            ee = e; ee(ee < base) = NaN;
            plot(ax(i), centers, ee, ':', 'Color', [0 0 0], 'LineWidth', 1.6, 'Tag', 'uniform_disk_reference');
            set(ax(i), 'YScale', 'log'); ylim(ax(i), [base, top_log]);
        case 'clipped'
            top = 1.5 * max(e);
            cc = min(c, top);
            patch(ax(i), [xs(1); xs(:); xs(end)], [0; reshape([cc; cc], [], 1); 0], col, ...
                'FaceAlpha', 0.55, 'EdgeColor', 'none', 'Tag', 'eig_real_histogram');
            plot(ax(i), centers, e, ':', 'Color', [0 0 0], 'LineWidth', 1.6, 'Tag', 'uniform_disk_reference');
            over = c > top;
            plot(ax(i), centers(over), top * ones(1, nnz(over)), '^', 'Color', col, 'MarkerFaceColor', col, ...
                'MarkerSize', 4, 'Tag', 'clipped_bar');
            ylim(ax(i), [0, 1.06 * top]);
    end
    plot(ax(i), [0 0], ylim(ax(i)), '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw, ...
        'Tag', 'zero_eigenvalue_reference');
    xlim(ax(i), re_lim); box(ax(i), 'off'); set(ax(i), 'TickDir', 'out');
    title(ax(i), titles{i}, 'FontWeight', 'normal', 'FontSize', 14);
    xlabel(ax(i), 'Re(\lambda) (1/s)');
    if i == 1; ylabel(ax(i), 'eigenvalues per state'); end
end

%% Table
fig_tag = 'fig_eig_real_histogram';
if strcmp(cfg.count_scale, 'clipped'); fig_tag = [fig_tag '_clipped']; end
L = {};
L{end + 1} = '# Jacobian eigenvalues projected on the real axis, against a uniform disk';
L{end + 1} = '';
L{end + 1} = sprintf(['Histogram: Re(lambda) of the eigenvalues of the full Jacobian at the sampled states, counts per state, ' ...
    '%d bins of width %.4g 1/s over [%.4g, %.4g] 1/s (the density figure''s window). Dotted curve: the projection of a ' ...
    'uniform disk holding n eigenvalues, radius drawn: `%s`. Count axis: `%s`.'], cfg.n_bins, width, re_lim(1), re_lim(2), ...
    mode, cfg.count_scale);
L{end + 1} = '';
L{end + 1} = ['| Condition | eigenvalues per state | n in the disk | center c (1/s) | radius drawn (1/s) | radius, unit gain | ' ...
    'radius, rms gain | rms gain | outside the disk''s real-axis extent | Re > 0 | Re(lambda) 0.1 and 99.9 percentiles | ' ...
    'within the extent: histogram / n | peak bin: histogram / curve at the center | right of zero: histogram / curve |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|';
for i = 1:n_cond
    re = real(evals{i}); r = ref(i);
    inside = re >= r.c - r.rho & re <= r.c + r.rho;
    [~, ic] = min(abs(centers - r.c));
    exp_pos = r.n * (1 - semicircle_cdf((0 - r.c) / r.rho));
    L{end + 1} = sprintf(['| %s | %g | %d | %.3g | %.3g | %.3g | %.3g | %.3g | %.1f%% | %.2f%% (%.2f per state) | %+.3g, %+.3g | ' ...
        '%.2f | %.3g / %.3g | %.2f / %.2f per state |'], titles{i}, r.N, r.n, r.c, r.rho, r.rho_unit, r.rho_gain, r.rms_gain, ...
        100 * mean(~inside), 100 * mean(re > 0), sum(re > 0) / r.n_states, quantile(re, 0.001), quantile(re, 0.999), ...
        (sum(inside) / r.n_states) / r.n, r.counts(ic), r.expected(ic), sum(re > 0) / r.n_states, exp_pos); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = ['Columns. **eigenvalues per state**: the state dimension N. **n in the disk**: the number of neurons, the ' ...
    'eigenvalues the reference describes. **radius, unit gain**: R / tau_d, R the theoretical bulk radius of W. **radius, ' ...
    'rms gain**: R rms(d) / tau_d, d the frozen gains (slope times the route''s output factor), rms over neurons, mean over ' ...
    'states (NaN when the run holds no gains). **outside the extent**: share of all eigenvalues with Re outside ' ...
    '[c - rho, c + rho]. **within the extent: histogram / n**: eigenvalues per state inside that range over n; above 1 ' ...
    'because the slow variables'' eigenvalues lie there too. **peak bin**: counts per state in the bin at the center ' ...
    'against the curve there. **right of zero**: eigenvalues per state with Re > 0 against the n the uniform disk puts ' ...
    'there.'];
L{end + 1} = '';
L{end + 1} = '## Caption paragraph, to adapt';
L{end + 1} = '';
L{end + 1} = ['Histograms: real parts of the eigenvalues of the Jacobian at the sampled states (the eigenvalues of the ' ...
    'density figure, projected on the real axis; counts per state). Dotted: the projection expected if the network''s n ' ...
    'eigenvalues filled a disk uniformly, as the eigenvalues of a random matrix do: a semicircle centered at -1/tau_d ' ...
    'whose radius is the random-matrix radius of the connectivity times the root-mean-square gain of the neurons at the ' ...
    'sampled states. The curve describes only the n eigenvalues of the dendritic block; the adaptation and depression ' ...
    'variables add the others, which appear as the narrow peaks near their own rates. Dashed: zero.'];
table_text = strjoin(L, newline);
vprintf(cfg.verbose, 'verbose', '%s\n', table_text);

out = struct('figs', fig, 'files', {{}}, 'source', data_file, 'reference', ref, 'edges', edges);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n', table_text);
        fclose(fid);
    end
end
end

function F = semicircle_cdf(u)
% Share of a uniform unit disk left of the vertical line at u.
u = min(1, max(-1, u));
F = 0.5 + (u .* sqrt(1 - u .^ 2) + asin(u)) / pi;
end
