function out = fig_freqresp_network(cfg)
% FIG_FREQRESP_NETWORK The network frequency response: bulk radius, rightmost real part and loop gain against frequency, per network.
%
%   out = FIG_FREQRESP_NETWORK('run_dir', d)
%   out = FIG_FREQRESP_NETWORK('data_file', f, 'out_dir', o)
%
% Draws what the opt-in stage frequency_response saved (FD-041;
% run_frequency_response_stage). It reads one run directory,
% <run_dir>/frequency_response/frequency_response_data.mat, or the file
% given, errors when neither holds it, and computes nothing beyond medians.
% Supplementary.
%
% THE OBJECT. With each neuron's adaptation and depression variables
% eliminated, a mode going as exp(z t) sees M(z) = W .* d(z), the weights
% with every presynaptic neuron's column scaled by its complex gain at z,
% and z is an eigenvalue of the full Jacobian exactly when 1 + tau_d z is
% an eigenvalue of M(z). The curves are for z = 2 pi i f.
%
% ONE ROW PER NETWORK (the stage's networks, in its order), three panels:
%   left    the predicted bulk radius of M(i omega); gray: |1 + i omega tau_d|,
%           the modulus an eigenvalue of M(i omega) must have for a mode
%           at that frequency to be marginal
%   middle  the largest real part among the eigenvalues of M(i omega); gray
%           line at 1
%   right   the loop gain, bulk radius over |1 + i omega tau_d|; gray line
%           at 1: above it the bulk holds growing modes rotating at that
%           frequency
% One curve per condition, in the condition's color (manuscript_style):
% the MEDIAN over every seed and calibration state, with the quartiles as
% a band. In a network whose ensemble mean matrix has an eigenvalue
% outside the structural bulk (E-to-E x 1.5), dashed: the leading
% eigenvalue of the mean part K(i omega), its real part (middle) and its
% modulus over |1 + i omega tau_d| (right), medians.
% MARKS on the loop-gain panel: a triangle at the frequency of the leftover
% mode (the eigenvalue of the full Jacobian the frozen-gain outlier
% continues into, from the outlier_spectra stage; median imaginary part
% over seeds and states, in Hz), for the conditions with slow variables;
% bars at the bottom: the range of rotation rates of the frozen bulk, from
% the lowest frequency drawn to R_frozen / (2 pi tau_d), per condition.
% (bulk_bars = false leaves the bars out; the table keeps their numbers.)
% The three columns each share one y range over the rows.
%
% NOTHING BUT THE DRAWING ON THE PANELS. Fig_Freqresp_Network_table.md
% holds the caption and, per network and condition, the numbers.
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions, in the order given; empty
% draws every condition the stage saved (select_conditions).
%
% See also: run_frequency_response_stage, fresp_run_one, fig_freqresp_plane,
%           fig_freqresp_phase, fig_freqresp_single_neuron, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the networks are recorded in the .mat
    cfg.conditions                = {}    % names to draw, in order; empty = all the run holds (select_conditions)
    cfg.networks                  = {}    % networks (rows) to draw, in order; empty = all the stage ran
    cfg.panels      (1,:) double {mustBeMember(cfg.panels, [1 2 3])} = [1 2 3]   % columns: 1 bulk radius, 2 rightmost real part, 3 loop gain
    cfg.bulk_bars   (1,1) logical = true    % the bars at the foot of the loop-gain panel (false: Main Figure 6, where nothing on the panel explains them)
    cfg.fig_size    double        = []      % [width height] of the figure in pixels; [] = 1250 x 350 per row for three panels, narrower in proportion for fewer
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
[D, data_file] = fresp_fig_load(cfg.data_file, cfg.run_dir);
networks = fresp_pick_networks(D.settings.networks, cfg.networks);
conditions = select_conditions(D.settings.conditions, cfg.conditions);
n_net = numel(networks); n_cond = numel(conditions);
gray = [0.6 0.6 0.6];

vis = 'on'; if ~cfg.visible; vis = 'off'; end
panels = unique(cfg.panels, 'stable'); n_p = numel(panels);
fig_size = cfg.fig_size;
if isempty(fig_size); fig_size = [round(1250 * n_p / 3), 350 * n_net]; end
fig = figure('Color', 'w', 'Position', [80 80 fig_size(1:2)], 'Visible', vis);
% A column that is not asked for is drawn on a scrap figure and thrown
% away: the table needs all three, and the drawing code stays one piece.
scrap = figure('Visible', 'off'); scrap_closer = onCleanup(@() close(scrap));
tl = tiledlayout(fig, n_net, n_p, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(n_net, 3);
rows = {};
y_hi = [2.2 2.2 1.6]; y_lo = [0 0 0];        % the ranges at 500 neurons; widened below when the data need it
marks = cell(n_net, n_cond);
for w = 1:n_net
    for c = 1:3
        k = find(panels == c);
        if isempty(k)
            ax(w, c) = axes(scrap); %#ok<LAXES> creates one, does not raise one
        else
            ax(w, c) = nexttile(tl, (w - 1) * n_p + k);
        end
        hold(ax(w, c), 'on');
        ax(w, c).Tag = sprintf('freqresp_network_%s_c%d', networks{w}, c);
    end
    for c = 1:n_cond
        runs = D.runs(strcmp({D.runs.network}, networks{w}) & strcmp({D.runs.condition}, conditions{c}));
        if isempty(runs)
            error('fig_freqresp_network:NoRun', 'The stage saved no run of network %s under condition %s.', ...
                networks{w}, conditions{c});
        end
        col = st.condition_color(conditions{c});
        f = runs(1).f; tau_d = runs(1).tau_d;
        den = abs(1 + 2i * pi * f * tau_d);
        P = pooled(runs);                       % every seed and state, one row each
        if c == 1
            plot(ax(w, 1), f, den, '-', 'Color', gray, 'LineWidth', 1.2);
        end
        Rm = median(P.R_both, 1);
        band(ax(w, 1), f, P.R_both, col);
        plot(ax(w, 1), f, Rm, '-', 'Color', col, 'LineWidth', st.line_lw);
        band(ax(w, 2), f, P.max_re, col);
        plot(ax(w, 2), f, median(P.max_re, 1), '-', 'Color', col, 'LineWidth', st.line_lw);
        band(ax(w, 3), f, P.R_both ./ den, col);
        plot(ax(w, 3), f, Rm ./ den, '-', 'Color', col, 'LineWidth', st.line_lw);
        K1_gain = median(abs(P.K1), 1) ./ den;
        if runs(1).structural.has_outlier
            plot(ax(w, 2), f, median(real(P.K1), 1), '--', 'Color', col, 'LineWidth', 1.4);
            plot(ax(w, 3), f, K1_gain, '--', 'Color', col, 'LineWidth', 1.4);
        end
        y_hi = max(y_hi, 1.05 * [max(quart(P.R_both, 75)), max(quart(P.max_re, 75)), max(quart(P.R_both ./ den, 75))]);
        y_lo = min(y_lo, [0, 1.05 * min(quart(P.max_re, 25)), 0]);

        % The leftover mode's frequency and the frozen bulk's range of
        % rotation rates, drawn once the y range is known.
        cont = P.cont(isfinite(P.cont));
        f_cont = median(imag(cont)) / (2 * pi);
        f_bulk = median(P.R_inf) / tau_d / (2 * pi);
        marks{w, c} = struct('f_cont', f_cont, 'f_bulk', f_bulk, 'has_slow', runs(1).has_slow, 'col', col, 'f1', f(1));

        % Numbers.
        Lg = Rm ./ den; [pk, ipk] = max(Lg);
        pk_seed = arrayfun(@(r) max(median(r.R.both, 1) ./ den), runs);
        i90 = find(Rm < 0.9 * median(P.R_inf), 1, 'last');
        f90 = NaN;
        if ~isempty(i90) && i90 < numel(f); f90 = f(i90 + 1); end
        planes = [runs.plane];
        rows(end + 1, :) = {networks{w}, conditions{c}, numel(runs), size(P.R_both, 1), ...
            median(P.R_inf), median(P.R_zero_both), median(P.R_zero_sfa), median(P.R_zero_std), f90, ...
            pk, f(ipk), sprintf('%.3g to %.3g', min(pk_seed), max(pk_seed)), ...
            max(median(P.R_sfa, 1) ./ den), max(median(P.R_std, 1) ./ den), ...
            median(P.edge, 'omitnan'), median([planes.abscissa]), ...
            median(real(cont)), median(imag(cont)), f_cont, f_bulk, ...
            median(real(P.K_inf1)), median(real(P.K_zero1)), max(K1_gain)}; %#ok<AGROW>
    end
end
for w = 1:n_net
    for c = 1:n_cond
        m = marks{w, c};
        y_bar = (0.025 + 0.022 * (c - 1)) * y_hi(3);
        if cfg.bulk_bars
            plot(ax(w, 3), [m.f1, m.f_bulk], [y_bar y_bar], '-', 'Color', m.col, 'LineWidth', 3);
        end
        if m.has_slow && isfinite(m.f_cont) && m.f_cont > 0
            plot(ax(w, 3), m.f_cont, 0.95 * y_hi(3), 'v', 'Color', m.col, 'MarkerFaceColor', m.col, 'MarkerSize', 8);
        end
    end
    yline(ax(w, 2), 1, '-', 'Color', gray, 'LineWidth', 1.0);
    yline(ax(w, 3), 1, '-', 'Color', gray, 'LineWidth', 1.0);
    for c = 1:3
        set(ax(w, c), 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', ...
            'Layer', 'top', 'XScale', 'log', 'XLim', D.settings.f([1 end]), 'XTick', decade_ticks(D.settings.f), ...
            'XMinorTick', 'off', 'YLim', [y_lo(c), y_hi(c)]);
        if w == n_net
            xlabel(ax(w, c), 'Frequency (Hz)', 'FontSize', st.label_fs);
        else
            set(ax(w, c), 'XTickLabel', []);
        end
    end
    ylabel(ax(w, 1), 'Bulk radius of M(i\omega)', 'FontSize', st.label_fs);
    ylabel(ax(w, 2), 'Rightmost Re eig M(i\omega)', 'FontSize', st.label_fs);
    ylabel(ax(w, 3), 'Loop gain', 'FontSize', st.label_fs);
end

clear scrap_closer
r1 = D.runs(1);
head = {'network', 'condition', 'seeds', 'states pooled', 'R frozen', 'R(0) both', 'R(0) adaptation alone', ...
    'R(0) depression alone', 'f where R reaches 90% of frozen (Hz)', 'peak loop gain', 'at Hz', ...
    'peak loop gain, range over seeds', 'peak, adaptation alone', 'peak, depression alone', ...
    'predicted bulk abscissa (1/s)', 'full Jacobian abscissa at the check states (1/s)', ...
    'leftover mode Re (1/s)', 'leftover mode Im (rad/s)', 'leftover mode Hz', 'frozen bulk rotates up to (Hz)', ...
    'K1 frozen', 'K1(0)', 'peak mean-mode loop gain'};
caption = [ ...
    'The network frequency response. With the slow variables eliminated neuron by neuron, a mode going as exp(z t) sees the matrix M(z) = W .* d(z): ' ...
    'the weights with every presynaptic neuron''s column scaled by that neuron''s complex gain at z. z is an eigenvalue of the full Jacobian exactly when ' ...
    '1 + tau_d z is an eigenvalue of M(z). The curves are for z = i 2 pi f: the MEDIAN over every seed and calibration state of the stage `outlier_spectra` ' ...
    sprintf('(%d seed(s), %d states each; %d neurons), with the quartiles as a band. ', size(D.settings.seeds, 1), numel(r1.t_cal), r1.n) ...
    sprintf('ROWS, top to bottom: the networks %s (`reference` is the paper''s network, `EEx1p5` the E-to-E mean times 1.5). ', strjoin(networks, ', ')) ...
    sprintf('COLOR: the condition (%s), in the manuscript''s condition colors; a condition without slow variables does not depend on frequency. ', strjoin(conditions, ', ')) ...
    'LEFT: the bulk radius of M(i omega) predicted from the gains (Perron root of the block variances times the summed squared gain moduli). ' ...
    'Gray: |1 + i omega tau_d|, the modulus an eigenvalue of M(i omega) must have for a mode at that frequency to be marginal (the dendritic low-pass). ' ...
    'MIDDLE: the largest real part among the eigenvalues of M(i omega); gray line at 1. Dashed, in a network whose block means put an eigenvalue of the mean matrix outside the bulk before any gain acts (judged on the ensemble the weights were drawn from): ' ...
    'the real part of the leading eigenvalue of the mean part K(i omega), where the structural outlier sits when it is outside the bulk. ' ...
    'RIGHT: the loop gain, bulk radius over |1 + i omega tau_d|; above 1 (gray line) the bulk holds growing modes rotating at that frequency. Dashed: the same for the mean-part eigenvalue. ' ...
    'Triangles: the frequency of the leftover mode, the eigenvalue of the full Jacobian that the frozen-gain outlier continues into (stage `outlier_spectra`; median imaginary part over seeds and states, in Hz). ' ...
    'Bars at the bottom: the range of rotation rates of the frozen bulk, up to R_frozen / (2 pi tau_d), per condition. ' ...
    'Table: medians over seeds and states. R in units of W (1 is marginal at zero frequency). "peak loop gain, range over seeds": the peak of each seed''s own median curve. ' ...
    '"predicted bulk abscissa" is the rightmost point of the region R(z) >= |1 + tau_d z| of the complex plane ' ...
    sprintf('(grid of %.3g 1/s by %.3g rad/s, real part from -8 to %.3g 1/s; a value at the grid''s right end means the region runs past it); ', ...
        D.settings.sigma(2) - D.settings.sigma(1), D.settings.w(2) - D.settings.w(1), D.settings.sigma(end)) ...
    'the full Jacobian''s abscissa is from the dense spectrum at the check states only, and with adaptation it can be the edge of the slow single-neuron band, not a network mode. ' ...
    'K1: the leading eigenvalue of the mean part (real part), frozen and at zero frequency.'];
out = fresp_fig_save(fig, cfg, out_dir, 'Fig_Freqresp_Network', data_file, ...
    'Network frequency response', caption, head, rows);
end

%% ------------------------------------------------------------------------
function P = pooled(runs)
% Every seed's calibration states stacked: one row per (seed, state).
R = [runs.R]; Z = [runs.R_zero];
P.R_both = cat(1, R.both); P.R_sfa = cat(1, R.sfa); P.R_std = cat(1, R.std);
P.max_re = cat(1, runs.max_re);
P.R_inf = cat(1, runs.R_inf);
P.R_zero_both = cat(1, Z.both); P.R_zero_sfa = cat(1, Z.sfa); P.R_zero_std = cat(1, Z.std);
K1 = arrayfun(@(r) permute(r.K_eig(1, :, :), [3 2 1]), runs, 'UniformOutput', false);   % states x nf each
P.K1 = cat(1, K1{:});
P.K_inf1 = cell2mat(arrayfun(@(r) r.K_inf(1, :), runs, 'UniformOutput', false));
P.K_zero1 = cell2mat(arrayfun(@(r) r.K_zero(1, :), runs, 'UniformOutput', false));
P.cont = [runs.cont];
P.edge = [runs.edge];
end

function q = quart(M, p)
% The p-th percentile down the rows (the samples); a single sample is itself.
if size(M, 1) == 1; q = M; else; q = prctile(M, p, 1); end
end

function band(ax, x, M, col)
if size(M, 1) < 2; return; end
lo = quart(M, 25); hi = quart(M, 75);
patch(ax, [x, fliplr(x)], [lo, fliplr(hi)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

function t = decade_ticks(f)
t = 10 .^ (ceil(log10(f(1)) - 1e-9):floor(log10(f(end)) + 1e-9));
end
