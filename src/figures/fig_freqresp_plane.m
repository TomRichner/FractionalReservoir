function out = fig_freqresp_plane(cfg)
% FIG_FREQRESP_PLANE The full Jacobian's eigenvalues in the complex plane, with the frozen-gain circle and the frequency-dependent edge of the bulk.
%
%   out = FIG_FREQRESP_PLANE('run_dir', d)
%   out = FIG_FREQRESP_PLANE('data_file', f, 'out_dir', o)
%
% Draws what the opt-in stage frequency_response saved (FD-041;
% run_frequency_response_stage). It reads one run directory,
% <run_dir>/frequency_response/frequency_response_data.mat, or the file
% given, errors when neither holds it, and computes nothing. Supplementary.
%
% ONE ROW PER NETWORK, ONE COLUMN PER CONDITION (titled by the condition,
% in its color), upper half plane (the spectrum is symmetric about the real
% axis), ONE SEED AND ONE STATE: the first seed of the stage (cfg.seed
% picks another) at the stage's shown calibration state. Per panel:
%   gray dots     every eigenvalue of the dense full Jacobian
%   dashed        the frozen-gain edge, the circle |1 + tau_d z| = R_frozen
%                 about -1/tau_d: where the bulk would end with the slow
%                 variables held
%   colored line  the edge with the slow variables free, the curve
%                 R(z) = |1 + tau_d z|, R(z) the bulk radius predicted from
%                 the neurons' gains at the complex frequency z
%   filled dot    (conditions with slow variables) the leftover mode: the
%                 eigenvalue the frozen-gain outlier continues into
%                 (stage outlier_spectra)
%   cross         where the frozen mean part has an eigenvalue outside the
%                 frozen bulk: the root of kappa(z) = 1 + tau_d z, kappa
%                 the leading eigenvalue of the mean part K(z)
%                 (fresp_outlier_root): the population mode from the mean
%                 gains alone
% Axes by complex_plane_axes: through the origin (the imaginary axis is the
% stability boundary), 1:1, one scale in every panel, ticks outside the
% largest frozen circle only, the axis names once (Re under the bottom
% row, Im on the first column). Eigenvalues outside the window (the far slow
% band, a far-left real eigenvalue) are not drawn; the table counts them.
%
% NOTHING BUT THE DRAWING ON THE PANELS. Fig_Freqresp_Plane_table.md holds
% the caption, the seed and state drawn, and the numbers per panel.
%
% Options beyond the standard ones
%   seed      [] (the first seed of the stage) or a 1 x 2 seed pair
%   re_lim, im_lim   the window, [-12 8] 1/s and [0 20] rad/s (the stage's
%             complex grid must cover it for the edge to be drawn across it)
%   tick_step 4
%
% See also: run_frequency_response_stage, fresp_run_one, complex_plane_axes,
%           fig_freqresp_network, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the networks are recorded in the .mat
    cfg.conditions                = {}    % names to draw, in order; empty = all the run holds (select_conditions)
    cfg.seed        double        = []
    cfg.re_lim      (1,2) double  = [-12 8]
    cfg.im_lim      (1,2) double  = [0 20]
    cfg.tick_step   (1,1) double  = 4
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
[D, data_file] = fresp_fig_load(cfg.data_file, cfg.run_dir);
networks = D.settings.networks;
conditions = select_conditions(D.settings.conditions, cfg.conditions);
titles = titles_by_name(conditions);
n_net = numel(networks); n_cond = numel(conditions);
seed = cfg.seed;
if isempty(seed); seed = D.settings.seeds(1, :); end
if ~(numel(seed) == 2 && ismember(seed(:)', D.settings.seeds, 'rows'))
    error('fig_freqresp_plane:NoSuchSeed', 'The stage saved no seed %s (it has %s).', ...
        mat2str(seed), mat2str(D.settings.seeds));
end

% The panels' data first: the ticks keep clear of the largest frozen circle.
A = cell(n_net, n_cond); P = A;
r_max = 0;
for w = 1:n_net
    for c = 1:n_cond
        k = find(strcmp({D.runs.network}, networks{w}) & strcmp({D.runs.condition}, conditions{c}) & ...
            arrayfun(@(r) isequal(r.seed(:)', seed(:)'), D.runs), 1);
        if isempty(k)
            error('fig_freqresp_plane:NoRun', 'The stage saved no run of network %s, condition %s, seed %s.', ...
                networks{w}, conditions{c}, mat2str(seed));
        end
        A{w, c} = D.runs(k);
        P{w, c} = A{w, c}.plane([A{w, c}.plane.state] == A{w, c}.shown_state);
        r_max = max(r_max, P{w, c}.R_inf / A{w, c}.tau_d);
    end
end
tau_d = A{1, 1}.tau_d;
axes_args = plane_axes_args(cfg, -1 / tau_d, r_max, st);

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Position', [80 80 420 * n_cond, 400 * n_net], 'Visible', vis);
tl = tiledlayout(fig, n_net, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
for w = 1:n_net
    for c = 1:n_cond
        a = A{w, c}; p = P{w, c};
        col = st.condition_color(conditions{c});
        ax = nexttile(tl, (w - 1) * n_cond + c);
        % The axis names once: Re under the bottom row, Im on the first column.
        re_label = ''; im_label = '';
        if w == n_net; re_label = 'Re (s^{-1})'; end
        if c == 1; im_label = 'Im (rad/s)'; end
        cp = complex_plane_axes(ax, axes_args{:}, 're_label', re_label, 'im_label', im_label);

        lam = double(p.lam);
        lam = lam(imag(lam) >= 0);
        shown = real(lam) >= cfg.re_lim(1) & real(lam) <= cfg.re_lim(2) & imag(lam) <= cfg.im_lim(2);
        plot(ax, cp.x(real(lam(shown))), imag(lam(shown)), '.', 'Color', [0.45 0.45 0.45], 'MarkerSize', 5);

        % Frozen edge: |1 + tau_d z| = R_inf, a circle about -1/tau_d.
        th = linspace(0, pi, 200);
        zc = (p.R_inf * exp(1i * th) - 1) / a.tau_d;
        plot(ax, cp.x(real(zc)), imag(zc), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);

        % Edge with the slow variables free: R(z) = |1 + tau_d z|.
        [SG, WW] = ndgrid(D.settings.sigma, D.settings.w);
        h = double(p.R) - abs(1 + a.tau_d * complex(SG, WW));
        if any(h(:) > 0) && any(h(:) < 0)
            contour(ax, cp.x(D.settings.sigma), D.settings.w, h', [0 0], 'Color', col, 'LineWidth', st.line_lw);
        end

        if a.has_slow && isfinite(p.cont)
            plot(ax, cp.x(real(p.cont)), imag(p.cont), 'o', 'Color', 'k', 'MarkerFaceColor', col, ...
                'MarkerSize', 7, 'LineWidth', 1.0);
        end
        if isfinite(p.root)
            plot(ax, cp.x(real(p.root)), imag(p.root), 'x', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
        end
        if w == 1
            title(ax, titles{c}, 'FontWeight', 'normal', 'FontSize', st.title_fs, 'Color', col);
        end

        edges = a.edge(isfinite(a.edge));
        planes = a.plane;
        rows(end + 1, :) = {networks{w}, conditions{c}, (p.R_inf - 1) / a.tau_d, p.edge, p.abscissa, ...
            real(p.cont), imag(p.cont), real(p.root), imag(p.root) + 0 * real(p.root), ...
            numel(lam), nnz(~shown), median(edges), ...
            mean(abs([planes.edge] - [planes.abscissa]))}; %#ok<AGROW>
    end
end

a = A{1, 1}; p = P{1, 1};
head = {'network', 'condition', 'frozen edge (1/s)', 'predicted edge (1/s)', 'full Jacobian abscissa (1/s)', ...
    'leftover mode Re (1/s)', 'leftover mode Im (rad/s)', 'mean-part root Re', 'mean-part root Im (rad/s)', ...
    'eigenvalues, upper half plane', 'of them outside the window', ...
    'predicted edge, median over this seed''s states', 'mean |predicted edge - abscissa| over the check states'};
caption = [ ...
    'Where the frequency response puts the full Jacobian''s network eigenvalues. ' ...
    sprintf('ONE SEED AND ONE STATE ARE DRAWN: seed [%d %d] (the first of the stage unless another was asked for), calibration state %d (t = %.4g s); %d neurons. ', ...
        seed(1), seed(2), p.state, p.t, a.n) ...
    sprintf('ROWS, top to bottom: the networks %s (`reference` is the paper''s network, `EEx1p5` the E-to-E mean times 1.5). COLUMNS: %s. ', ...
        strjoin(networks, ', '), strjoin(conditions, ', ')) ...
    'Gray dots: all eigenvalues of the dense full Jacobian (upper half plane; the spectrum is symmetric about the real axis). ' ...
    sprintf('Dashed: the frozen-gain edge, the circle |1 + tau_d z| = R_frozen about -1/tau_d = %.3g 1/s: where the bulk would end if the slow variables were held. ', -1 / a.tau_d) ...
    'Colored line: the edge with the slow variables free, the curve R(z) = |1 + tau_d z|, with R(z) the bulk radius predicted from the neurons'' gains at the complex frequency z (`fresp_radius`); ' ...
    'network eigenvalues are expected inside it. Without slow variables the two edges are the same circle. ' ...
    'Filled dot (conditions with slow variables): the leftover mode, the eigenvalue the frozen-gain outlier continues into (stage `outlier_spectra`). ' ...
    'Cross (only where the frozen mean part has an eigenvalue outside the frozen bulk): the root of kappa(z) = 1 + tau_d z, ' ...
    'kappa the leading eigenvalue of the mean part K(z) (`fresp_outlier_root`): where the population mode should be from the mean gains alone. ' ...
    'The axes are drawn through the origin, so the imaginary axis is the stability boundary; one scale in every panel; ticks outside the largest frozen circle only. ' ...
    sprintf('The window is Re in [%g, %g] 1/s, Im in [%g, %g] rad/s; eigenvalues outside it (the band of slow single-neuron modes on the far negative real axis, any far-left real eigenvalue) are not drawn and are counted in the table. ', ...
        cfg.re_lim, cfg.im_lim) ...
    'The structure hugging the negative real axis is the band of slow single-neuron modes (poles of the gains), not network modes; the prediction is a bulk edge for large networks and individual eigenvalues scatter about it. ' ...
    sprintf('Table: "predicted edge" is the rightmost point of the region R(z) >= |1 + tau_d z| (grid of %.3g 1/s by %.3g rad/s, real part from -8 to %.3g 1/s; a value at the grid''s right end means the region runs past it).', ...
        D.settings.sigma(2) - D.settings.sigma(1), D.settings.w(2) - D.settings.w(1), D.settings.sigma(end))];
out = fresp_fig_save(fig, cfg, out_dir, 'Fig_Freqresp_Plane', data_file, ...
    'Full-Jacobian eigenvalues and the frequency-dependent edge', caption, head, rows);
end

%% ------------------------------------------------------------------------
function args = plane_axes_args(cfg, center, radius, st)
% One geometry for every panel. Ticks at multiples of tick_step outside the
% largest frozen circle; when the circle covers the whole window on an
% axis, that axis gets the default ticks instead of none.
args = {'re_lim', cfg.re_lim, 'im_lim', cfg.im_lim, 'tick_step', cfg.tick_step, ...
    'font_size', st.tick_fs, 'line_width', st.axis_lw, ...
    're_label_side', 'above'};          % below would fall off the bottom row's tile
mult = @(lim) (ceil(lim(1) / cfg.tick_step):floor(lim(2) / cfg.tick_step)) * cfg.tick_step;
re_t = mult(cfg.re_lim); re_t = re_t(re_t ~= 0 & abs(re_t - center) >= radius);
im_t = mult(cfg.im_lim); im_t = im_t(im_t ~= 0 & abs(im_t) >= radius);
if numel(im_t) > 1; im_t = im_t(im_t < cfg.im_lim(2)); end     % the top end is where a column title sits
if ~isempty(re_t); args = [args, {'re_ticks', re_t, 're_clear', center + [-1 1] * radius}]; end
if ~isempty(im_t); args = [args, {'im_ticks', im_t, 'im_clear', [-1 1] * radius}]; end
end
