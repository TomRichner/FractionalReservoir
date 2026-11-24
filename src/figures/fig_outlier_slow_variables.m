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
%                   the negative real axis, are not network modes.
%   gray dots       those left-out eigenvalues (dendritic share below
%                   x_thr), faint and under everything else, so that
%                   the picture is the whole spectrum of the window.
%                   Until 2026-11-22 they were not drawn, which left
%                   white lens-shaped voids about the negative real
%                   axis in the adapting conditions: eigenvalues there
%                   are mostly slow single-neuron modes with a small
%                   imaginary part, and a reader took the hole for a
%                   feature of the spectrum. cfg.show_filtered = false
%                   leaves them out again.
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
% the stability boundary. Ticks only where the spectrum is not: one on the
% real axis right of zero, one near each end of the imaginary axis (the
% window is a cut through the bulk, which fills its left side, so there
% are no ticks on the negative real axis). The panels share their window
% and their ticks, so the tick LABELS are on the outer panels only, as for
% any shared axes: the real ones on the bottom row, the imaginary ones on
% the first column (the reference network's outliers fill the right side
% of its panels, where a label would lie on them whatever the window). The
% real tick is the largest multiple of a round step whose label no open or
% filled circle of the bottom row comes within 1.5 pt of; the imaginary
% one likewise for the first column (the sidecar names a tick left out).
%
% THE AXIS NAMES, the same in every panel that has one. "Re" (bottom row)
% is outside the window, from the right end of the real axis rightward,
% under its tick labels; the gaps between the columns are made wide
% enough for it. "Im" (first column) is right of the top of the imaginary
% axis, as high in the window as it goes, and raised over the highest
% open or filled circle that would be under it (1.5 pt clear); what then
% reaches above the window gets a strip over each row, the titles raised
% by it. So they read at any window: set inside the window by fixed
% rules, they met the outliers and the "-8" once the window shrank to
% +/-8 (the fastmedium run, 2026-11-24).
%
% SIZE: the figure is manuscript_style canvas_w wide (or fig_w, which a
% composite passes: its block's width), so it is drawn at the size it is
% printed; three columns of that width hold the condition titles side by
% side. The height follows from the window (1:1).
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
    cfg.show_filtered (1,1) logical = true    % the eigenvalues below x_thr, faint gray, under the network modes
    cfg.filtered_color (1,3) double = [0.62 0.62 0.62]
    cfg.filtered_alpha (1,1) double = 0.18
    cfg.fig_w double = []                   % inches; [] = manuscript_style canvas_w (a composite passes its block's width)
    cfg.panel_w double = []                 % inches; [] = the panels fill the width. Smaller: the same columns (titles), wider gaps, shorter rows (S7, FD-057)
    cfg.row_gap (1,1) double = 0.25         % inches between the network rows
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
im_tick = step * (1:floor(0.92 * im_lim(2) / step));          % the candidates; one pair is drawn (below)
re_step = nice_step(max(re_lim(2), eps) / 2);
re_ticks = re_step:re_step:re_lim(2);
im_ticks = [-im_tick, im_tick]; im_ticks = im_ticks(im_ticks ~= 0 & im_ticks >= im_lim(1) & im_ticks <= im_lim(2));
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 're_ticks', re_ticks, 'im_ticks', unique(im_ticks), ...
    'font_size', st.tick_fs, 'line_width', st.axis_lw};
g = complex_plane_axes([], axes_args{:});

%% Panels placed by hand (inches): identical sizes, 1:1
% The width is the common canvas (manuscript_style canvas_w) or the block a
% composite gives it. The margins are what the type needs, in inches: the
% network name left of the first column, the condition titles above, the
% axis names outside each window (Re right of it, Im over it); the columns
% as close as that allows.
fig_w = cfg.fig_w; if isempty(fig_w); fig_w = st.canvas_w; end
pt_in = 1 / 72;
RE_NAME = 'Re (s^{-1})'; IM_NAME = 'Im (s^{-1})';
name_sz = text_inches({RE_NAME, IM_NAME}, st.tick_fs);
re_room = name_sz(1, 1) + 0.12;                     % right of every plane: the Re name and its gap to the next
left = 0.5; right = 0.05; bottom = 0.08;
col_gap = max(0.22, re_room); trail = max(0, re_room - right);
pw = (fig_w - left - right - trail - (n_cond - 1) * col_gap) / n_cond;
if ~isempty(cfg.panel_w) && cfg.panel_w < pw
    % Smaller planes in the same columns: each panel centered in its column,
    % so the condition titles keep their room and the rows get shorter.
    pw_col = pw; pw = cfg.panel_w;
    left = left + (pw_col - pw) / 2; col_gap = col_gap + (pw_col - pw);
end
ph = pw * diff(g.y_lim) / diff(g.x_lim);
upi = diff(g.x_lim) / pw;                           % axis units per inch (1:1)
marks_of = @(D) cellfun(@(d) [with_conjugates(d.frozen), with_conjugates(d.cont)], D(~cellfun(@isempty, D)), 'UniformOutput', false);
m_re = marks_of(P(n_net, :)); m_im = marks_of(P(:, 1));
% The Im name (first column): right of the top of the imaginary axis, as
% high in the window as it goes, raised over the highest open or filled
% circle under it (1.5 pt clear) where one would be under it; the part
% above the window gets a strip of its own over each row (the titles go
% up by it).
tl = 0.018 * diff(g.y_lim);
r_mark = (cfg.mark_size / 2 + 1 + 1.5) / 72 * upi;
z = [m_im{:}]; z = z(isfinite(z));
under = real(z) >= 1.2 * tl - r_mark & real(z) <= 1.2 * tl + name_sz(2, 1) * upi + r_mark;
im_y = g.y_lim(2) - name_sz(2, 2) * upi;            % hanging from the top of the window
if any(under); im_y = max(im_y, max(imag(z(under))) + r_mark); end
im_room = max(0, im_y + name_sz(2, 2) * upi - g.y_lim(2)) / upi;
row_gap = max(cfg.row_gap, im_room + 0.04);
top_m = 0.1 + 1.25 * st.title_fs * pt_in + im_room;
% The tick labels are on the outer panels (real: bottom row; imaginary:
% first column). On each axis the one tick drawn is the largest of the
% candidates whose label none of those panels' marks comes near; when
% every candidate is covered, the largest (the sidecar says so).
hidden_re = covered_ticks(g, [m_re{:}], cfg.mark_size, st.tick_fs, upi);
[~, hidden_im] = covered_ticks(g, [m_im{:}], cfg.mark_size, st.tick_fs, upi);
re_keep = setdiff(g.re_ticks, hidden_re); im_keep = setdiff(g.im_ticks, hidden_im);
if isempty(re_keep); re_keep = g.re_ticks; hidden_re = []; end
if isempty(im_keep); im_keep = g.im_ticks; hidden_im = []; end
re_keep = max(re_keep);
im_keep = im_keep(abs(im_keep) == max(abs(im_keep)));
hidden_re = hidden_re(hidden_re > re_keep);                  % the larger ticks passed over
hidden_im = hidden_im(abs(hidden_im) > max(abs(im_keep)));
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 're_ticks', re_keep, 'im_ticks', im_keep, ...
    'font_size', st.tick_fs, 'line_width', st.axis_lw};
g = complex_plane_axes([], axes_args{:});
g.hidden_re_ticks = hidden_re; g.hidden_im_ticks = hidden_im;
fig_h = top_m + n_net * ph + (n_net - 1) * row_gap + bottom;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', vis);
n_marks_out = zeros(n_net, n_cond);
n_filtered = zeros(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        ax = axes(fig, 'Units', 'inches', 'Position', ...
            [left + (c - 1) * (pw + col_gap), bottom + (n_net - w) * (ph + row_gap), pw, ph]);
        % Names and tick labels on the outer panels: Re on the bottom row,
        % Im on the first column; the names outside the window.
        names = {'re_label', '', 'im_label', '', 're_label_side', 'right', 'im_label_y', im_y, ...
            're_tick_labels', w == n_net, 'im_tick_labels', c == 1};
        if w == n_net; names{2} = RE_NAME; end
        if c == 1;     names{4} = IM_NAME; end
        cp = complex_plane_axes(ax, axes_args{:}, names{:});
        d = P{w, c};
        if ~isempty(d)
            col = F.colors{c};
            if cfg.show_filtered
                zf = filtered_modes(F.folder, d);
                zf = zf(real(zf) >= re_lim(1) & real(zf) <= re_lim(2) & abs(imag(zf)) <= im_lim(2));
                n_filtered(w, c) = numel(zf);
                scatter(ax, cp.x(real(zf)), imag(zf), cfg.dot_size, cfg.filtered_color, 'filled', ...
                    'MarkerFaceAlpha', cfg.filtered_alpha, 'MarkerEdgeColor', 'none', 'Tag', 'slow_variables_filtered');
            end
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
            % Raised by the strip that holds the Im name over the window.
            tt = title(ax, F.titles{c}, 'FontSize', st.title_fs, 'FontWeight', 'normal', 'Color', F.colors{c});
            tt.Units = 'inches'; tt.Position(2) = tt.Position(2) + im_room; tt.Units = 'data';
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
        if cfg.show_filtered
            fprintf(fid, ['GRAY DOTS, faint and under the others: the eigenvalues that filter leaves out (dendritic share ' ...
                'below %g), mostly the slow modes of single neurons, on and close to the negative real axis. They are ' ...
                'drawn so that the window shows the whole spectrum: without them the adapting conditions have white ' ...
                'lens-shaped holes about the real axis, which are the filter and no property of the network modes. ' ...
                'Drawn in the window, per panel in reading order: %s.\n\n'], cfg.x_thr, ...
                strjoin(arrayfun(@(v) sprintf('%d', v), reshape(n_filtered.', 1, []), 'UniformOutput', false), ', '));
        else
            fprintf(fid, ['The eigenvalues that filter leaves out are not drawn (show_filtered false): white holes ' ...
                'about the negative real axis in the adapting conditions are the filter, not the spectrum.\n\n']);
        end
        fprintf(fid, ['AXES: through the origin (`complex_plane_axes`), one unit of Re as long as one unit of Im, ' ...
            'identical in every panel, no break: real axis [%.3g, %.3g] 1/s, imaginary axis [%.3g, %.3g] 1/s. The ' ...
            'imaginary axis is the stability boundary. The view is a cut: the spectrum continues to the left, above ' ...
            'and below. One tick right of zero on the real axis and one near each end of the imaginary axis, clear of ' ...
            'the bulk, labelled on the outer panels only (real: the bottom row; imaginary: the first column). %s ' ...
            'Marks (open or filled) outside the view, per panel in reading order: %s. Data: `%s`.\n\n'], ...
            g.re_lim, g.y_lim, hidden_ticks_text(g), strjoin(arrayfun(@(v) sprintf('%d', v), reshape(n_marks_out.', 1, []), 'UniformOutput', false), ', '), F.file);
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
function z = filtered_modes(folder, d)
% The eigenvalues of the dense full-Jacobian spectra of d's run that the
% dendritic-share filter leaves out (xfrac below d.x_thr), every
% calibration state, as a row. outlier_fig_data keeps the network modes
% only, so the run's file is read once more here.
runs = ospec_read_index(folder);
k = find(strcmp(runs.network, d.network) & strcmp(runs.condition, d.condition) & ...
    runs.seed1 == d.seed(1) & runs.seed2 == d.seed(2) & runs.ok, 1);
z = zeros(1, 0);
if isempty(k); return; end
R = ospec_load_run(folder, runs(k, :));
if ~R.has_dense || isempty(R.dense); return; end
lam = double([R.dense.lam]); xf = double([R.dense.xfrac]);
z = reshape(lam(xf < d.x_thr & isfinite(lam)), 1, []);
end

function [hidden_re, hidden_im] = covered_ticks(g, z, mark_size, fs, units_per_inch)
% The ticks whose LABEL an open or filled circle would lie on. The label
% boxes are those of complex_plane_axes: under the real axis, centered on
% the tick; left of the imaginary axis, centered on the tick's height.
% Sizes from the font (a digit is about 0.56 em wide), a little generous,
% plus the marker's radius. (fig_outlier_plane does the same for its dots
% and diamonds.)
pt = units_per_inch / 72;                            % one point, in axis units
tl = 0.018 * diff(g.y_lim);                          % complex_plane_axes' tick length
z = z(isfinite(z));
x = real(z); y = imag(z);                            % no break in this figure: the plotted x is the real part
r = (mark_size / 2 + 1 + 1.5) * pt;                  % the marker, its edge and 1.5 pt of clearance
label_w = @(t) 0.56 * fs * pt * (numel(sprintf('%g', abs(t))) + (t < 0));
label_h = 0.85 * fs * pt;
top = -1.5 * tl - 0.1 * fs * pt;                     % the glyphs begin a little under the top of the text box
hidden_re = [];
for t = g.re_ticks
    if any(abs(x - t) <= label_w(t) / 2 + r & y <= top + r & y >= top - label_h - r)
        hidden_re(end + 1) = t; %#ok<AGROW>
    end
end
hidden_im = [];
for t = g.im_ticks
    if any(x <= -1.6 * tl + r & x >= -1.6 * tl - label_w(t) - r & abs(y - t) <= label_h / 2 + r)
        hidden_im(end + 1) = t; %#ok<AGROW>
    end
end
end

function s = hidden_ticks_text(g)
% The sidecar's sentence on the ticks left out.
if isempty(g.hidden_re_ticks) && isempty(g.hidden_im_ticks)
    s = 'No tick was left out under a mark.';
    return;
end
parts = {};
if ~isempty(g.hidden_re_ticks); parts{end + 1} = ['on the real axis at ' strjoin(compose('%g', g.hidden_re_ticks), ', ') ' 1/s']; end
if ~isempty(g.hidden_im_ticks); parts{end + 1} = ['on the imaginary axis at ' strjoin(compose('%g', g.hidden_im_ticks), ', ') ' 1/s']; end
s = ['TICKS PASSED OVER, in every panel, because an open or filled circle of a panel that carries the tick labels ' ...
    'comes near where the label would be: ' strjoin(parts, ' and ') '; the largest clear tick is drawn. The window is unchanged.'];
end

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

function sz = text_inches(strs, fs)
% The width and height, in inches, of each string set at fs points (tex), one row each.
f = figure('Visible', 'off', 'Units', 'inches', 'Position', [0 0 4 2]);
cleaner = onCleanup(@() close(f));
ax = axes(f, 'Units', 'inches', 'Position', [0 0 4 2], 'Visible', 'off');
sz = zeros(numel(strs), 2);
for k = 1:numel(strs)
    h = text(ax, 0, 0, strs{k}, 'FontSize', fs, 'Interpreter', 'tex', 'Units', 'inches');
    e = h.Extent; sz(k, :) = e(3:4);
end
end
