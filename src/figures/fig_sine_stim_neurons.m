function out = fig_sine_stim_neurons(cfg)
% FIG_SINE_STIM_NEURONS Single neurons of the adapting network under the 2 s sine train (a review figure).
%
%   out = FIG_SINE_STIM_NEURONS('run_dir', d)                % d holds stim_neurons/
%   out = FIG_SINE_STIM_NEURONS('run_dir', d, 'out_dir', o, 'window', [30 36])
%
% For Tom's review, not the poster (TR, 2026-10-20): what the individual
% neurons of the multiple-timescale network do over a few modulation cycles
% of the 2 s sine train, from run_stim_neurons (a rerun of one
% run_stim_evoked run that keeps single neurons; it reproduces the saved
% population mean to 4e-9). Rows, sharing the time axis:
%
%   (a) external input u_i (the 1 ms pulses; zero on undriven neurons)
%   (b) the Wiener term as an equivalent input, sigma_u xi_1 / sqrt(dt):
%       what the integrator adds to x in each 0.5 ms step, expressed in the
%       units of u (run_stim_neurons)
%   (c) dendritic potential x_i
%   (d) firing rate r_i
%   (e) synaptic output r_i prod_m b_im (route 1; all routes are identical)
%   (f) SFA feedback (c/K) sum_k a_ik
%   (g) STD product prod_m b_im
%
% Thin lines: the 25 neurons sampled with a fixed seed, coloured by class
% (driven E, undriven E, I). Thick lines: the mean of each class over the
% whole network. Dashed lines: the envelope peaks. Equations:
% docs/EquationsParametersDocs/Equations_stability_paper.md.
%
% 16 x 10.8 in, 14 pt; layout by hand (see fig_sine_paired_stim for why).
% The SVG's text is OUTLINED all the same: R2026a's vector export outlines
% every text object once the figure holds this many points (about 2 x 10^6:
% 25 neurons x 7 rows x 12,000 samples), as it does for a TeX Greek command or
% the Unicode square root (which is why (b) writes (Delta t)^{1/2}). The PNG is
% the reference for this figure.
%
% See also: run_stim_neurons, fig_sine_paired_stim, fig_example_timeseries

arguments
    cfg.run_dir  (1,:) char    = ''
    cfg.out_dir  (1,:) char    = ''
    cfg.mat_file (1,:) char    = ''
    cfg.preset_name (1,:) char = ''      % unused; accepted for a uniform call
    cfg.window   (1,2) double  = [30 36]  % s of stimulation time shown (3 cycles)
    cfg.save     (1,1) logical = true
    cfg.visible  (1,1) logical = true
    cfg.verbose  (1,:) char    = 'minimal'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
if isempty(cfg.mat_file)
    mat_file = resolve_data_file('', cfg.run_dir, {fullfile(cfg.run_dir, 'stim_neurons')}, ...
        'stim_neurons_data.mat', 'run_stim_neurons(sine_paired_stim_config())');
else
    mat_file = cfg.mat_file;
end
Z = load(mat_file, 'data', 'settings'); D = Z.data;
FN = 'Arial'; fsz = 14;
PI = char(928); SIG = char(931); XI = char(958); SG = char(963); DLT = char(916);

tcol = SRNNCellTypePairs.type_colors(2);
ccol = [tcol(1, :) * 0.75; tcol(1, :) + 0.45 * (1 - tcol(1, :)); tcol(2, :) + 0.35 * (1 - tcol(2, :))];   % driven E, undriven E, I (thin)
mcol = [0.35 0.1 0.0; 0.95 0.55 0.2; 0.0 0.15 0.45];   % the class means: dark brown, orange, navy

keep = D.t >= cfg.window(1) - 1e-9 & D.t <= cfg.window(2) + 1e-9;
t = D.t(keep);
rows = { ...
    'u',     sprintf('(a)  External input (u_i)'); ...
    'noise', sprintf('(b)  Noise input as an equivalent input (%s_u %s_{1,i} / (%st)^{1/2})', SG, XI, DLT); ...
    'x',     '(c)  Dendritic potential (x_i)'; ...
    'r',     '(d)  Firing rate (r_i)'; ...
    'syn',   sprintf('(e)  Synaptic output (r_i %s_m b_{im})', PI); ...
    'sfa',   sprintf('(f)  SFA feedback ((c/K) %s_k a_{ik})', SIG); ...
    'std',   sprintf('(g)  Synaptic resources, STD product (%s_m b_{im})', PI)};

W = 16; H = 10.8;
fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'inches', 'Position', [0.2 0.2 W H]);
lx = 1.0; lw = W - lx - 0.3;
n_r = size(rows, 1); rh = 1.0; slot = 1.33; b0 = 1.1;
peaks = (ceil(cfg.window(1) / D.mod_period - 0.5):floor(cfg.window(2) / D.mod_period - 0.5)) * D.mod_period ...
    + D.mod_period / 2;
ax_all = gobjects(1, n_r);
for i = 1:n_r
    b = b0 + (n_r - i) * slot;
    ax = axes(fig, 'Units', 'inches', 'Position', [lx b lw rh], 'FontName', FN); hold(ax, 'on');
    ax_all(i) = ax;
    v = rows{i, 1};
    Y = double(D.traces.(v)(:, keep));
    % Thin: the sample, drawn I first so E lands on top (draw_order).
    for c = [3 2 1]
        for j = find(D.sample_class == c)
            plot(ax, t, Y(j, :), '-', 'Color', ccol(c, :), 'LineWidth', 0.6);
        end
    end
    % Thick: the class means over the whole network.
    for c = [3 2 1]
        plot(ax, t, double(D.means.(v)(c, keep)), '-', 'Color', mcol(c, :), 'LineWidth', 2.6);
    end
    for pk = peaks, xline(ax, pk, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1); end
    xlim(ax, cfg.window);
    switch v
        case 'u',     ylim(ax, [-3 43]);
        case 'noise', ylim(ax, [-4.5 4.5]);
    end
    text(ax, 0, 1.03, rows{i, 2}, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
        'FontSize', fsz, 'FontName', FN, 'Interpreter', 'tex');
    set(ax, 'FontSize', fsz - 1, 'FontName', FN, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');
    if i < n_r, set(ax, 'XTickLabel', []); end
end
xlabel(ax_all(end), 'Time since the train started (s)', 'FontSize', fsz, 'FontName', FN);
% Legend: classes, how many of each are drawn, and what the thick lines are.
% Column-major with NumColumns 3: each class's thin swatch above its mean.
hl = gobjects(1, 6); lab = cell(1, 6);
for c = 1:3
    hl(2 * c - 1) = plot(ax_all(1), NaN, NaN, '-', 'Color', ccol(c, :), 'LineWidth', 1.2);
    lab{2 * c - 1} = sprintf('%s: %d of the %d neurons (thin)', D.class_names{c}, D.sample_counts(c), D.counts(c));
    hl(2 * c) = plot(ax_all(1), NaN, NaN, '-', 'Color', mcol(c, :), 'LineWidth', 2.6);
    lab{2 * c} = sprintf('%s: mean of all %d (thick)', D.class_names{c}, D.counts(c));
end
lg = legend(ax_all(1), hl, lab, 'NumColumns', 3, 'FontSize', fsz - 1, 'FontName', FN, 'Box', 'off');
lg.Units = 'inches'; lg.Position(1:2) = [(W - lg.Position(3)) / 2, 0.05];
ttl = sprintf('Multiple-timescale adaptation network, %g s sine train: 25 random neurons (seed %d) and the class means', ...
    D.mod_period, Z.settings.opts.sample_seed);
% In two pieces: a single text wider than ~6.8 in vanishes from the 600-dpi PNG.
th = text(ax_all(1), 0, 1.03 + 0.42 / rh, ttl(1:strfind(ttl, ':')), 'Units', 'normalized', ...
    'VerticalAlignment', 'bottom', 'FontSize', fsz + 1, 'FontName', FN);
text(ax_all(1), th.Extent(3) + 0.005, 1.03 + 0.42 / rh, ttl(strfind(ttl, ':') + 1:end), 'Units', 'normalized', ...
    'VerticalAlignment', 'bottom', 'FontSize', fsz + 1, 'FontName', FN);

out = struct('figs', fig, 'files', {{}}, 'source', sprintf('%s; window %s s; git %s', mat_file, ...
    mat2str(cfg.window), Z.settings.git.commit_short), 'counts', D.counts, 'sample_counts', D.sample_counts);
if cfg.save
    save_figure_stable(out_dir, 'fig_sine_stim_neurons', fig);
    out.files = existing_outputs(out_dir, 'fig_sine_stim_neurons');
end
vprintf(cfg.verbose, 'minimal', '  [sine_stim_neurons] %d files in %s\n', numel(out.files), out_dir);
if cfg.visible, set(fig, 'Visible', 'on'); end
end
