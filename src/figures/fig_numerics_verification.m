function out = fig_numerics_verification(cfg)
% FIG_NUMERICS_VERIFICATION Supplemental figures on the precision of the numerics.
%
%   out = FIG_NUMERICS_VERIFICATION('run_dir', d, 'variant', 'solver')
%   out = FIG_NUMERICS_VERIFICATION('run_dir', d, 'variant', 'lya_method')
%
% The PLOT half of run_numerics_verification; the .mat it reads is produced by
% that stage, inside the run directory, and there is no fallback. Two
% variants, registered as two entries so each gets its own output folder:
%
%   'solver'      One column per adaptation regime. Row 1: x of two neurons
%                 from the free-running ode45 (1e-10) and SRA1 (fs 400) runs,
%                 with both Benettin LLEs -- where the LLE is positive the
%                 traces must part, whatever the precision. Row 2: the
%                 reshooting error of SRA1 at 400 Hz against the ode45
%                 reference, per state family, over one 0.02 s segment,
%                 relative to that family's RMS -- the error free of chaotic
%                 amplification. Row 3: that error vs step size, noise-free
%                 (SRA1's drift is second-order RK: slope 2) and with noise on
%                 a shared Brownian path (strong order 1.5), with the fitted
%                 slopes. See the stage's header for why a noise-free slope
%                 between 1 and 2 is expected on a piecewise activation.
%   'lya_method'  One column per regime, on the REDUCED network the stage
%                 states in settings.n_small. Row 1: the top-K local QR
%                 exponents as a band with Benettin's local exponent on top
%                 (the overlay from test_benettin_vs_qr). Row 2: the sorted QR
%                 spectrum with Benettin's LLE as a line through lambda_1.
%
% See also: run_numerics_verification, resolve_data_file, test_benettin_vs_qr

arguments
    cfg.variant     (1,:) char {mustBeMember(cfg.variant, {'solver', 'lya_method'})} = 'solver'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.n_qr_show   (1,1) double  = 20    % local QR exponents drawn in the band
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'numerics_verification')}, ...
    'numerics_verification_data.mat', ...
    'Run run_numerics_verification first');
D = load(data_file);
R = D.results;
S = D.settings;
n_cond = numel(R);

switch cfg.variant
    case 'solver'
        fig_tag = 'Fig_numerics_solver';
        fig = plot_solver(R, S, st, n_cond, cfg.visible);
    case 'lya_method'
        fig_tag = 'Fig_numerics_lya_method';
        fig = plot_lya_method(R, S, st, n_cond, cfg.visible, cfg.n_qr_show);
end

out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

%% ------------------------------------------------------------------------
function fig = plot_solver(R, S, st, n_cond, visible)
fig = figure('Position', [100, 80, 430 * n_cond, 900], 'Visible', onoff(visible));
tl  = tiledlayout(fig, 3, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
blocks = S.block_names;
j400 = find([R(1).free.fs] == S.fs_lle, 1);
if isempty(j400); j400 = 1; end

for i = 1:n_cond
    col = cond_color(st, R(i).name, i);
    L   = R(i).lle;

    % Row 1: free-running overlay
    ax = nexttile(tl, i);
    hold(ax, 'on');
    n_show = min(2, size(L.ode45.x_examples, 1));
    offs = 0;
    for k = 1:n_show
        xo = L.ode45.x_examples(k, :);
        xs = L.sra1.x_examples(k, :);
        plot(ax, L.ode45.t, xo + offs, '-',  'Color', col, 'LineWidth', 1.4, 'HandleVisibility', vis(k));
        plot(ax, L.sra1.t,  xs + offs, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', vis(k));
        offs = offs + 1.2 * (max(xo) - min(xo));
    end
    hold(ax, 'off');
    xlim(ax, [L.ode45.t(1), L.ode45.t(end)]);
    set(ax, 'YTick', [], 'FontSize', st.tick_fs);
    xlabel(ax, 'time (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'x, two neurons', 'FontSize', st.label_fs); end
    title(ax, {R(i).title, sprintf('%s: ode45 %+.3f, SRA1 %+.3f', ...
        st.label_lle, L.ode45.LLE, L.sra1.LLE)}, 'FontWeight', 'normal', ...
        'FontSize', st.title_fs);
    if i == 1
        legend(ax, {sprintf('ode45, tol %g', S.ref_tol), ...
            sprintf('SRA1, %d Hz', S.fs_lle)}, 'Location', 'best', 'FontSize', st.tick_fs - 2);
    end

    % Row 2: reshooting error per state family at the paper's rate
    ax = nexttile(tl, n_cond + i);
    F  = R(i).free(j400);
    hold(ax, 'on');
    styles = {'-', '--', ':'};
    labels = {};
    for b = 1:numel(blocks)
        rel = F.err_long(:, b) / F.block_rms(b);
        if all(~isfinite(rel)); continue; end
        plot(ax, F.t_long, rel, styles{b}, 'Color', col, 'LineWidth', 1.2);
        labels{end+1} = sprintf('%s (rms %.2g)', blocks{b}, F.block_rms(b)); %#ok<AGROW>
    end
    hold(ax, 'off');
    set(ax, 'YScale', 'log', 'FontSize', st.tick_fs);
    xlim(ax, [F.t_long(1), F.t_long(end)]);
    xlabel(ax, 'restart time (s)', 'FontSize', st.label_fs);
    if i == 1
        ylabel(ax, sprintf('|error| / rms over %g s', S.seg_long), 'FontSize', st.label_fs);
    end
    title(ax, sprintf('SRA1 %d Hz reshot from ode45', F.fs), 'FontWeight', 'normal', ...
        'FontSize', st.title_fs);
    legend(ax, labels, 'Location', 'best', 'FontSize', st.tick_fs - 2);

    % Row 3: convergence with step size
    ax = nexttile(tl, 2 * n_cond + i);
    hold(ax, 'on');
    dt  = 1 ./ [R(i).free.fs];
    e0  = [R(i).free.err_long_total_rms];
    p0  = polyfit(log(dt), log(e0), 1);
    loglog(ax, dt, e0, 'o-', 'Color', col, 'MarkerFaceColor', col, 'LineWidth', st.line_lw, ...
        'DisplayName', sprintf('noise-free, slope %.2f', p0(1)));
    loglog(ax, dt, e0(end) * (dt / dt(end)).^2, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
        'DisplayName', 'slope 2 (RK2 drift)');
    if ~isempty(R(i).noisy)
        e1 = [R(i).noisy.err_long_total_rms];
        p1 = polyfit(log(dt), log(e1), 1);
        loglog(ax, dt, e1, 's--', 'Color', col, 'MarkerFaceColor', 'w', 'LineWidth', st.line_lw, ...
            'DisplayName', sprintf('noise \\sigma_u = %g, slope %.2f', S.sigma_u_noise, p1(1)));
        loglog(ax, dt, e1(end) * (dt / dt(end)).^1.5, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'DisplayName', 'slope 1.5 (SRA1 strong)');
    end
    hold(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', st.tick_fs, 'XDir', 'reverse');
    xticks(ax, sort(dt)); xticklabels(ax, arrayfun(@(f) sprintf('1/%d', f), sort([R(i).free.fs], 'descend'), 'UniformOutput', false));
    xlabel(ax, 'step (s)', 'FontSize', st.label_fs);
    if i == 1
        ylabel(ax, sprintf('rms |error| over %g s', S.seg_long), 'FontSize', st.label_fs);
    end
    title(ax, 'convergence of the reshooting error', 'FontWeight', 'normal', 'FontSize', st.title_fs);
    legend(ax, 'Location', 'best', 'FontSize', st.tick_fs - 2);
end

title(tl, {sprintf('Numerical precision on the n = %d network (%s)', S.n, S.preset_name), ...
    sprintf(['reference: ode45 RelTol = AbsTol = %g; reshooting resets SRA1 to the reference ' ...
             'every segment so chaos cannot amplify the error'], S.ref_tol)}, ...
    'FontWeight', 'bold', 'Interpreter', 'none');
end

function fig = plot_lya_method(R, S, st, n_cond, visible, K_show)
fig = figure('Position', [100, 80, 430 * n_cond, 720], 'Visible', onoff(visible));
tl  = tiledlayout(fig, 2, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_cond
    col = cond_color(st, R(i).name, i);
    Q   = R(i).qr;
    LLE_b = Q.benettin.LLE;
    LLE_q = Q.qr.LE_spectrum(1);
    K = min(K_show, size(Q.qr.local_LE_spectrum_t, 2));

    % Row 1: local exponents overlay
    ax = nexttile(tl, i);
    hold(ax, 'on');
    cmap = parula(K);
    for j = K:-1:1
        h = plot(ax, Q.qr.t_lya, Q.qr.local_LE_spectrum_t(:, j), 'Color', cmap(j, :), 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        if j == K; set(h, 'DisplayName', sprintf('QR local \\lambda_{%d}', K), 'HandleVisibility', 'on'); end
        if j == 1; set(h, 'DisplayName', 'QR local \lambda_1', 'HandleVisibility', 'on'); end
    end
    plot(ax, Q.benettin.t_lya, Q.benettin.local_lya, 'Color', col, 'LineWidth', 1.6, ...
        'DisplayName', 'Benettin local');
    yline(ax, LLE_b, '--', 'Color', col, 'LineWidth', 1.2, 'DisplayName', 'Benettin LLE');
    yline(ax, LLE_q, ':', 'Color', [0 0 0], 'LineWidth', 1.2, 'DisplayName', 'QR \lambda_1');
    hold(ax, 'off');
    set(ax, 'FontSize', st.tick_fs);
    xlim(ax, [min(Q.qr.t_lya(1), Q.benettin.t_lya(1)), max(Q.qr.t_lya(end), Q.benettin.t_lya(end))]);
    xlabel(ax, 'time (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'local Lyapunov exponent', 'FontSize', st.label_fs); end
    title(ax, {R(i).title, sprintf('Benettin %+.4f   QR \\lambda_1 %+.4f', LLE_b, LLE_q)}, ...
        'FontWeight', 'normal', 'FontSize', st.title_fs);
    if i == 1; legend(ax, 'Location', 'best', 'FontSize', st.tick_fs - 2); end

    % Row 2: the sorted spectrum
    ax = nexttile(tl, n_cond + i);
    spec = Q.qr.LE_spectrum(:);
    hold(ax, 'on');
    plot(ax, 1:numel(spec), spec, '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 6);
    plot(ax, 1, spec(1), 'o', 'Color', [0 0 0], 'MarkerFaceColor', 'w', 'MarkerSize', 7);
    yline(ax, LLE_b, '--', 'Color', col, 'LineWidth', 1.4);
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    hold(ax, 'off');
    set(ax, 'FontSize', st.tick_fs);
    xlim(ax, [0, numel(spec) + 1]);
    xlabel(ax, 'exponent index', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'Lyapunov exponent (1/s)', 'FontSize', st.label_fs); end
    title(ax, sprintf('QR spectrum, %d states; |\\Delta\\lambda_1| = %.4f', ...
        numel(spec), abs(LLE_b - LLE_q)), 'FontWeight', 'normal', 'FontSize', st.title_fs);
    if i == 1
        legend(ax, {'QR spectrum', 'QR \lambda_1', 'Benettin LLE'}, ...
            'Location', 'best', 'FontSize', st.tick_fs - 2);
    end
end

title(tl, {sprintf('Benettin vs QR on a reduced network: n = %d, same preset physics (%s)', ...
    Q.n, S.preset_name), ...
    sprintf('noise-free, ode45; T = [%g, %g] s with %g s warmup; QR is O(N^2) so the full n = %d network is out of reach', ...
    Q.T_range, S.lya_small_warmup, S.n)}, 'FontWeight', 'bold', 'Interpreter', 'none');
end

%% ------------------------------------------------------------------------
function c = cond_color(st, name, i)
if st.condition_color.isKey(name)
    c = st.condition_color(name);
else
    cm = lines(max(i, 7));
    c = cm(i, :);
end
end

function s = onoff(tf)
if tf; s = 'on'; else; s = 'off'; end
end

function v = vis(k)
% Only the first neuron's pair carries legend entries.
if k == 1; v = 'on'; else; v = 'off'; end
end
