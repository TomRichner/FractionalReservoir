function out = fig_outlier_slow_variables(cfg)
% FIG_OUTLIER_SLOW_VARIABLES What the slow variables do to the rightmost frozen-gain mode: the full Jacobian's spectrum near zero.
%
%   out = FIG_OUTLIER_SLOW_VARIABLES('run_dir', d)
%   out = FIG_OUTLIER_SLOW_VARIABLES('run_dir', d, 'seed', [1 2], 're_lim', [-14 10], 'im_lim', [-9 9])
%
% Draws what the opt-in stage outlier_spectra saved (FD-040), through
% outlier_fig_data. Reads one run directory, or the index file given as
% data_file; errors when neither holds the stage, or when the networks
% asked for have no dense full-Jacobian spectra; simulates nothing.
% FD-038 (Tom, 2026-11-20: "should go in the supplement"). Supplementary.
%
% ROWS are the networks (the reference and E-to-E x 1.5), COLUMNS the
% conditions. Each panel is a close-up, near the origin, of the spectrum
% of the FULL Jacobian (adaptation and depression free to respond) at the
% calibration states of ONE seed (cfg.seed; the first of the index),
% overlaid, in 1/s:
%   dots            the network modes: eigenvalues whose eigenvector has at
%                   least x_thr (0.1) of its squared norm in the dendritic
%                   variables. The single-neuron slow modes, thousands on
%                   the negative real axis, are left out.
%   open circles    the rightmost eigenvalue of the FROZEN-gain matrix at
%                   each state, as a rate, with its conjugate when complex
%   filled circles  its CONTINUATION: the eigenvalue of the full Jacobian
%                   whose dendritic pattern overlaps most with that frozen
%                   mode, with its conjugate (condition color, black edge)
% Without adaptation there are no slow variables and the two coincide (the
% open circle is under the filled one).
%
% AXES (complex_plane_axes; docs/Figure_style_guide.md): through the
% origin, 1:1, IDENTICAL in every panel, no break. The window is given
% (re_lim, im_lim) or fitted once to the marks of the drawn panels: it
% holds every open and filled circle, with room above and below, and
% extends one half-height left of the leftmost mark. The imaginary axis is
% the stability boundary. Ticks only where the spectrum is not: on the
% real axis right of zero, on the imaginary axis one tick near each end
% (the window is a cut through the bulk, which fills its left side, so
% there are no ticks on the negative real axis).
%
% Text on the figure: the condition names, in their colors, the network
% names, the axis names and ticks. Fig_Outlier_Slow_Variables_table.md has
% a caption paragraph and, per panel, median [min, max] over ALL seeds of
% the per-run medians: the frozen rightmost rate, the continuation and its
% OVERLAP with the frozen mode (how safely it is identified), the number
% of eigenvalues beyond zero, the full Jacobian's abscissa, and how often
% the continuation is the same eigenvalue when the dendritic fraction
% that admits a candidate is 0.02 or 0.3 in place of 0.1.
%
% See also: outlier_fig_data, complex_plane_axes, ospec_dense_spectrum,
%           fig_outlier_plane, fig_outlier_dwell

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.networks                  = {'reference', 'EEx1p5'}
    cfg.conditions                = {}
    cfg.seed double               = []    % the seed drawn; [] = the first of the index
    cfg.x_thr (1,1) double {mustBeGreaterThanOrEqual(cfg.x_thr, 0), mustBeLessThanOrEqual(cfg.x_thr, 1)} = 0.1
    cfg.re_lim double             = []    % [lo hi] 1/s; [] = fitted to the marks of the drawn panels
    cfg.im_lim double             = []
    cfg.dot_size (1,1) double     = 5
    cfg.dot_alpha (1,1) double    = 0.2
    cfg.mark_size (1,1) double    = 7
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = outlier_fig_data(cfg, 'fig_outlier_slow_variables', 'slow', 'networks', cfg.networks, 'x_thr', cfg.x_thr);
st = F.st;
[n_net, n_cond, n_seed] = size(F.D);

%% The marks of the drawn panels, and one window for all
P = F.D(:, :, F.i_seed);
drawn = [P{~cellfun(@isempty, P)}];
if isempty(drawn)
    error('fig_outlier_slow_variables:NoSeed', 'No finished run of seed [%d %d] among the panels.', F.seeds(F.i_seed, 1), F.seeds(F.i_seed, 2));
end
marks = [drawn.frozen, drawn.cont];
marks = marks(isfinite(marks));
half = 1.45 * max([abs(imag(marks)), 1]);
im_lim = [-half, half];
re_lim = [min([real(marks), 0]) - half, max([real(marks), 0]) + 0.15 * half];
if ~isempty(cfg.im_lim); im_lim = cfg.im_lim(:)'; end
if ~isempty(cfg.re_lim); re_lim = cfg.re_lim(:)'; end
step = nice_step(im_lim(2) / 4);
im_tick = step * floor(0.92 * im_lim(2) / step);
re_step = nice_step(max(re_lim(2), eps) / 2);
re_ticks = re_step:re_step:re_lim(2);
im_ticks = [-im_tick, im_tick]; im_ticks = im_ticks(im_ticks ~= 0 & im_ticks >= im_lim(1) & im_ticks <= im_lim(2));
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 're_ticks', re_ticks, 'im_ticks', unique(im_ticks), ...
    'font_size', st.tick_fs, 'line_width', st.axis_lw};
g = complex_plane_axes([], axes_args{:});

%% Panels placed by hand (inches): identical sizes, 1:1
fig_w = 13; left = 0.6; right = 0.15; col_gap = 0.35; top_m = 0.5; row_gap = 0.45; bottom = 0.5;
pw = (fig_w - left - right - (n_cond - 1) * col_gap) / n_cond;
ph = pw * diff(g.y_lim) / diff(g.x_lim);
fig_h = top_m + n_net * ph + (n_net - 1) * row_gap + bottom;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', vis);
n_marks_out = zeros(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        ax = axes(fig, 'Units', 'inches', 'Position', ...
            [left + (c - 1) * (pw + col_gap), bottom + (n_net - w) * (ph + row_gap), pw, ph]);
        names = {'re_label', '', 'im_label', ''};
        if w == n_net; names{2} = 'Re (s^{-1})'; end
        if c == 1;     names{4} = 'Im (s^{-1})'; end
        cp = complex_plane_axes(ax, axes_args{:}, names{:});
        d = P{w, c};
        if ~isempty(d)
            col = F.colors{c};
            z = d.net(~isnan(d.net));
            z = z(real(z) >= re_lim(1) & real(z) <= re_lim(2) & abs(imag(z)) <= im_lim(2));
            scatter(ax, cp.x(real(z)), imag(z), cfg.dot_size, col, 'filled', 'MarkerFaceAlpha', cfg.dot_alpha, 'MarkerEdgeColor', 'none');
            fr = with_conjugates(d.frozen);
            ct = with_conjugates(d.cont);
            plot(ax, cp.x(real(fr)), imag(fr), 'o', 'MarkerSize', cfg.mark_size, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
            plot(ax, cp.x(real(ct)), imag(ct), 'o', 'MarkerSize', cfg.mark_size, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
            m = [fr, ct];
            n_marks_out(w, c) = nnz(real(m) < re_lim(1) | real(m) > re_lim(2) | abs(imag(m)) > im_lim(2));
        end
        set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);
        if w == 1
            title(ax, F.titles{c}, 'FontSize', st.title_fs, 'FontWeight', 'normal', 'Color', F.colors{c});
        end
        if c == 1
            set(ax.YLabel, 'String', F.net_labels{w}, 'FontSize', st.label_fs, 'Visible', 'on', 'Color', 'k');
        end
    end
end

%% The numbers, per panel over every seed (one value per run: its median over the states)
S = cell(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        runs = squeeze(F.D(w, c, :)); runs = runs(~cellfun(@isempty, runs));
        if isempty(runs); continue; end
        runs = [runs{:}];
        s = struct();
        s.n_seeds = numel(runs);
        s.frozen_re = arrayfun(@(d) median(real(d.frozen)), runs);
        s.cont_re = arrayfun(@(d) median(real(d.cont)), runs);
        s.cont_im = arrayfun(@(d) median(imag(d.cont)), runs);
        s.cont_ov = arrayfun(@(d) median(d.cont_ov), runs);
        s.cont_ov_all = [runs.cont_ov];                               % every state of every seed
        s.cont_xfrac = arrayfun(@(d) median(d.cont_xfrac), runs);
        s.n_pos = arrayfun(@(d) median(d.n_pos), runs);
        s.frac_unstable = arrayfun(@(d) mean(d.n_pos > 0), runs);
        s.abscissa = arrayfun(@(d) median(d.abscissa), runs);
        s.same_lo = arrayfun(@(d) mean(d.same_alt(1, :)), runs);
        s.same_hi = arrayfun(@(d) mean(d.same_alt(2, :)), runs);
        S{w, c} = s;
    end
end

fig_tag = 'Fig_Outlier_Slow_Variables';
out = struct('figs', fig, 'files', {{}}, 'source', F.file, 'window', g, 'numbers', {S});
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        first = drawn(1);
        fprintf(fid, '# %s\n\n', fig_tag);
        fprintf(fid, ['Rows: the networks (%s). Columns: the conditions. A close-up near the origin of the spectrum of ' ...
            'the FULL Jacobian (adaptation and depression variables free to respond), in 1/s, at the %d calibration ' ...
            'states of one run, network seed [%d %d] (%d neurons), overlaid. DOTS: eigenvalues whose eigenvector has ' ...
            'at least %g of its squared norm in the dendritic variables (network modes; the single-neuron slow modes ' ...
            'on the negative real axis are left out). OPEN CIRCLES: the rightmost eigenvalue of the frozen-gain ' ...
            'effective connectivity at each state, as a rate (lambda - 1)/tau_d (slow variables held), with its ' ...
            'conjugate when it is complex. FILLED CIRCLES: the eigenvalue of the full Jacobian whose dendritic ' ...
            'pattern overlaps most with that frozen mode (its continuation), with its conjugate. Without adaptation ' ...
            'there are no slow variables, the two coincide, and the open circle lies under the filled one. With ' ...
            'adaptation the full spectrum solves z = (lambda(M(z)) - 1)/tau_d with the gains filtered by the slow ' ...
            'variables, so a frozen mode that is real, or nearly so, can continue as a complex pair, a slow ' ...
            'oscillation in the band between the adaptation filter and the dendritic low-pass; the table says ' ...
            'where each panel''s continuation lies.\n\n'], ...
            strjoin(F.networks, ', '), first.n_cal, F.seeds(F.i_seed, :), first.n, cfg.x_thr);
        fprintf(fid, ['AXES: through the origin (`complex_plane_axes`), one unit of Re as long as one unit of Im, ' ...
            'identical in every panel, no break: real axis [%.3g, %.3g] 1/s, imaginary axis [%.3g, %.3g] 1/s. The ' ...
            'imaginary axis is the stability boundary. The view is a cut: the spectrum continues to the left, above ' ...
            'and below. Ticks only right of zero on the real axis and near the ends of the imaginary axis, clear of ' ...
            'the bulk. Marks (open or filled) outside the view, per panel in reading order: %s. Data: `%s`.\n\n'], ...
            g.re_lim, g.y_lim, strjoin(arrayfun(@(v) sprintf('%d', v), reshape(n_marks_out.', 1, []), 'UniformOutput', false), ', '), F.file);
        fprintf(fid, ['Per panel, over the %d seed(s) of the index: each run gives the median over its calibration ' ...
            'states, and the table the median of those over seeds with [min, max]. `frozen` is the real part of the ' ...
            'rightmost frozen-gain rate; `continuation` the real and imaginary part of the full-Jacobian eigenvalue ' ...
            'that carries its pattern; `overlap` the modulus of the inner product of the two dendritic patterns, ' ...
            'each normalized (1: the same pattern), per run and, in the next column, the smallest and largest over ' ...
            'every state of every seed; `x fraction` the share of the continuation''s eigenvector in the dendritic ' ...
            'variables; `beyond zero` the number of eigenvalues of the full Jacobian with positive real part and ' ...
            '`unstable states` the fraction of states with any; `abscissa` its rightmost real part (with adaptation, ' ...
            'in a stable network, a single-neuron slow mode); `same at 0.02`, `same at 0.3` the fraction of states ' ...
            'at which the continuation is the same eigenvalue when a candidate needs that dendritic fraction ' ...
            'instead of 0.1. Rates in 1/s.\n\n'], n_seed);
        fprintf(fid, ['| network | condition | seeds | frozen Re | continuation Re | continuation Im | overlap | ' ...
            'overlap, all states | x fraction | beyond zero | unstable states | abscissa | same at 0.02 | same at 0.3 |\n']);
        fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
        for w = 1:n_net
            for c = 1:n_cond
                s = S{w, c};
                if isempty(s); continue; end
                fprintf(fid, '| %s | %s | %d | %s | %s | %s | %s | [%.2f, %.2f] | %s | %s | %s | %s | %s | %s |\n', ...
                    F.networks{w}, F.cond_names{c}, s.n_seeds, mmm(s.frozen_re, '%.3g'), mmm(s.cont_re, '%.3g'), ...
                    mmm(s.cont_im, '%.3g'), mmm(s.cont_ov, '%.2f'), min(s.cont_ov_all), max(s.cont_ov_all), ...
                    mmm(s.cont_xfrac, '%.2f'), mmm(s.n_pos, '%.3g'), mmm(s.frac_unstable, '%.2f'), ...
                    mmm(s.abscissa, '%.3g'), mmm(s.same_lo, '%.2f'), mmm(s.same_hi, '%.2f'));
            end
        end
        if F.n_failed > 0
            fprintf(fid, '\n%d run(s) of these networks did not finish and are not in the numbers.\n', F.n_failed);
        end
        clear closer
        out.files{end + 1} = table_name;
    end
end
for w = 1:n_net
    for c = 1:n_cond
        if isempty(S{w, c}); continue; end
        vprintf(cfg.verbose, 'verbose', '%s | %s: continuation %s %+si 1/s, overlap %s\n', F.networks{w}, F.cond_names{c}, ...
            mmm(S{w, c}.cont_re, '%.3g'), mmm(S{w, c}.cont_im, '%.3g'), mmm(S{w, c}.cont_ov, '%.2f'));
    end
end
end

%% ------------------------------------------------------------------------
function z = with_conjugates(z)
% The conjugate of each complex eigenvalue appended; a real one (|Im| below
% 1e-9 of the larger of its modulus and 1) is kept once, on the real axis.
z = z(:).';
z = z(isfinite(z));
is_real = abs(imag(z)) <= 1e-9 * max(abs(z), 1);
z(is_real) = real(z(is_real));
z = [z, conj(z(~is_real))];
end

function s = mmm(v, fmt)
% median [min, max] over seeds, collapsed when they coincide; NaN ignored.
v = v(~isnan(v));
if isempty(v)
    s = '-';
elseif isscalar(v) || max(v) == min(v)
    s = sprintf(fmt, median(v));
else
    s = sprintf([fmt ' [' fmt ', ' fmt ']'], median(v), min(v), max(v));
end
end

function s = nice_step(x)
% 1, 2 or 5 times a power of ten, the largest not above x.
e = 10^floor(log10(x));
m = [1 2 5 10];
s = e * m(find(m * e <= x, 1, 'last'));
end
