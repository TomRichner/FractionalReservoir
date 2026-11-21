function files = fig_frequency_response(opts)
%FIG_FREQUENCY_RESPONSE The figures of the frequency-response exploration, each with a markdown sidecar.
%
%   files = fig_frequency_response()
%   files = fig_frequency_response('which', {'single', 'network'})
%
% Reads data/explorations/frequency_response/ (freqresp_single,
% freqresp_network) and, for the eigenvalues of the full Jacobian and the
% continuation modes, yesterday's outlier runs in
% data/explorations/outliers/ (read only). Writes into
% figs/explorations/frequency_response/:
%
%   freqresp_single_neuron   magnitude and phase of one neuron's response,
%                            at three operating points taken from the runs
%   freqresp_network         the network curves: bulk radius, rightmost
%                            real part and loop gain against frequency
%   freqresp_one_vs_multiple the ratio of the one- and multiple-timescale
%                            attenuation, matched and at their own
%                            operating points, and the phase
%   freqresp_plane           the full Jacobian's eigenvalues with the edge
%                            the frequency-dependent radius predicts
%
% CONVENTIONS. Panels carry axis labels and nothing else: no legends,
% counts or annotations. The condition is the color (manuscript_style;
% no adaptation in dark red, FR-004), the mechanism the line style:
% solid = adaptation and depression free, dashed = adaptation alone,
% dotted = depression alone. Everything else is in the sidecar.
%
% Name-value options
%   which       cellstr of the figures to draw (default: all four)
%   seed        [7 8]
%   state       the calibration state the plane figure shows (7)
%   resolution  300
%   visible     'off'
%
% See also: explore_frequency_response, freqresp_network, freqresp_single

arguments
    opts.which (1,:) cell = {'single', 'network', 'one_vs_multiple', 'plane'}
    opts.seed (1,2) double = [7 8]
    opts.state (1,1) double = 7
    opts.resolution (1,1) double = 300
    opts.visible (1,:) char = 'off'
end

root = fileparts(which('setup_paths'));
D.data = fullfile(root, 'data', 'explorations', 'frequency_response');
D.outl = fullfile(root, 'data', 'explorations', 'outliers');
D.figs = fullfile(root, 'figs', 'explorations', 'frequency_response');
if ~isfolder(D.figs); mkdir(D.figs); end
D.st = manuscript_style();
D.opts = opts;
D.networks = {'reference', 'EEx1p5'};
D.conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};

files = {};
for w = opts.which
    switch w{1}
        case 'single';          files{end + 1} = fig_single(D);          %#ok<AGROW>
        case 'network';         files{end + 1} = fig_network(D);         %#ok<AGROW>
        case 'one_vs_multiple'; files{end + 1} = fig_one_vs_multiple(D); %#ok<AGROW>
        case 'plane';           files{end + 1} = fig_plane(D);           %#ok<AGROW>
        otherwise
            error('fig_frequency_response:Which', 'Unknown figure ''%s''.', w{1});
    end
end
end

%% ------------------------------------------------------------------------
function file = fig_single(D)
[~, B] = freqresp_single(D.opts.seed);
f = B.f;
labels = {'E median', 'I median', 'E driven'};
modes = {'both', 'sfa', 'std'};
styles = {'-', '--', ':'};

fig = new_fig([1 1 11 6], D);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:3
    ax1 = nexttile(tl, p); hold(ax1, 'on');
    ax2 = nexttile(tl, 3 + p); hold(ax2, 'on');
    for o = B.own(strcmp({B.own.label}, labels{p}))
        col = cond_color(D.st, o.condition);
        for m = 1:3
            G = o.G.(modes{m});
            plot(ax1, f, abs(G), styles{m}, 'Color', col, 'LineWidth', D.st.line_lw);
            plot(ax2, f, angle(G) * 180 / pi, styles{m}, 'Color', col, 'LineWidth', D.st.line_lw);
        end
    end
    style_axis(ax1, D); style_axis(ax2, D);
    set(ax1, 'XScale', 'log', 'XLim', f([1 end]), 'YLim', [0 1.05], 'XTickLabel', []);
    set(ax2, 'XScale', 'log', 'XLim', f([1 end]), 'YLim', [0 80], 'XTick', 10.^(-3:2));
    set(ax1, 'XTick', 10.^(-3:2));
    xlabel(ax2, 'Frequency (Hz)', 'FontSize', D.st.label_fs);
    if p == 1
        ylabel(ax1, 'Gain / frozen gain', 'FontSize', D.st.label_fs);
        ylabel(ax2, 'Phase lead (deg)', 'FontSize', D.st.label_fs);
    end
end
file = save_fig(fig, D, 'freqresp_single_neuron');

rows = {};
for o = B.own
    for m = 1:3
        s = o.summary.(modes{m});
        rows(end + 1, :) = {o.condition, o.network, o.label, modes{m}, o.neuron, o.slope, o.rate, o.F, ...
            s.G0, s.dB0, s.f_half, s.f_3dB, s.f_90, s.phase_max, s.f_phase, s.peak, o.max_diff}; %#ok<AGROW>
    end
end
head = {'condition', 'network', 'neuron picked', 'free', 'neuron', 'slope', 'rate', 'output factor F', ...
    'G(0)', 'G(0) dB', 'f_half Hz', 'f_3dB Hz', 'f_90 Hz', 'max phase lead deg', 'at Hz', 'peak', 'closed vs numeric'};
caption = [ ...
    'One neuron''s frequency response: the gain from its dendritic input to its synaptic output with its slow variables free, ' ...
    'divided by the frozen gain (slope times available vesicle fraction, the limit of high frequency). TOP: magnitude; BOTTOM: phase lead. ' ...
    'COLUMNS, left to right: the excitatory neuron of median rate in the reference network; the inhibitory neuron of median rate there; ' ...
    'a hard-driven excitatory neuron (90th percentile of rate) in the network with the E-to-E mean times 1.5; each among the neurons with nonzero slope, ' ...
    sprintf('at calibration state %d (t = 70.9 s) of seed [%d %d], so each condition has its own neuron and operating point (table). ', B.state, B.seed(1), B.seed(2)) ...
    'COLOR: orange one timescale (`sfa1_std1`), blue multiple timescales (`sfa3_stdc2`). LINE: solid adaptation and depression both free; ' ...
    'dashed adaptation alone (depression held); dotted depression alone (adaptation held). A value of 1 means the slow variables do nothing at that frequency. ' ...
    'The curves are the closed form (`freqresp_gain`); the same response computed from the model''s own Jacobian blocks for that neuron agrees to the last column of the table (relative). ' ...
    'G(0): gain at zero frequency over the frozen gain. f_half: where half the attenuation in dB is gone. f_3dB, f_90: above these the gain stays within 3 dB, within 10%, of the frozen gain. ' ...
    'peak: the largest magnitude (1 = no resonance of the neuron itself).'];
write_sidecar(D, 'freqresp_single_neuron', caption, head, rows);
end

%% ------------------------------------------------------------------------
function file = fig_network(D)
n_net = numel(D.networks); n_cond = numel(D.conditions);
fig = new_fig([1 1 12.5 7], D);
tl = tiledlayout(fig, n_net, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
for w = 1:n_net
    ax = gobjects(1, 3);
    for c = 1:3
        ax(c) = nexttile(tl, (w - 1) * 3 + c); hold(ax(c), 'on');
    end
    for c = 1:n_cond
        A = load_network(D, D.networks{w}, D.conditions{c});
        col = cond_color(D.st, D.conditions{c});
        f = A.f; om = 2 * pi * f;
        den = abs(1 + 1i * om * A.tau_d);
        if c == 1
            plot(ax(1), f, den, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
        end
        Rm = median(A.R.both, 1);
        band(ax(1), f, prctile(A.R.both, 25, 1), prctile(A.R.both, 75, 1), col);
        plot(ax(1), f, Rm, '-', 'Color', col, 'LineWidth', D.st.line_lw);

        max_re = squeeze(max(real(A.ev), [], 1));           % nf_ev x n_states
        if size(A.ev, 2) == 1; max_re = repmat(max_re(:)', numel(f), 1); end
        plot(ax(2), f, median(max_re, 2), '-', 'Color', col, 'LineWidth', D.st.line_lw);

        plot(ax(3), f, Rm ./ den, '-', 'Color', col, 'LineWidth', D.st.line_lw);
        K1 = squeeze(A.K_eig(1, :, :));                    % nf x n_states
        if strcmp(D.networks{w}, 'EEx1p5')
            plot(ax(2), f, median(real(K1), 2), '--', 'Color', col, 'LineWidth', 1.4);
            plot(ax(3), f, median(abs(K1), 2)' ./ den, '--', 'Color', col, 'LineWidth', 1.4);
        end

        % The frequencies of the modes found in the outlier runs, and the
        % range of rotation rates of the frozen bulk (up to R_inf / tau_d).
        O = load_outlier(D, D.networks{w}, D.conditions{c});
        cont = [O.dense.cont];
        f_cont = median(imag(cont)) / (2 * pi);
        f_bulk = median(A.R_inf) / A.tau_d / (2 * pi);
        y_bar = 0.04 + 0.035 * (c - 1);
        plot(ax(3), [f(1) f_bulk], [y_bar y_bar], '-', 'Color', col, 'LineWidth', 3);
        if c > 1
            plot(ax(3), f_cont, 1.52, 'v', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 8);
        end

        % Numbers.
        Lg = Rm ./ den; [pk, ipk] = max(Lg);
        Ls = median(A.R.sfa, 1) ./ den; Ld = median(A.R.std, 1) ./ den;
        f90 = f(find(Rm < 0.9 * median(A.R_inf), 1, 'last') + 1);
        if isempty(f90); f90 = NaN; end
        edge = plane_abscissa(A);
        rows(end + 1, :) = {D.networks{w}, D.conditions{c}, median(A.R_inf), median(A.R_zero.both), ...
            median(A.R_zero.sfa), median(A.R_zero.std), f90, pk, f(ipk), max(Ls), max(Ld), ...
            median(edge), median([O.dense.abscissa]), corr_(edge, [O.dense.abscissa]), ...
            median(real(cont)), median(imag(cont)), f_cont, f_bulk, ...
            median(real(A.K_inf(1, :))), median(real(A.K_zero(1, :))), max(median(abs(K1), 2)' ./ den)}; %#ok<AGROW>
    end
    for c = 1:3
        style_axis(ax(c), D);
        set(ax(c), 'XScale', 'log', 'XLim', [1e-3 1e2], 'XTick', 10.^(-3:2));
        if w == n_net
            xlabel(ax(c), 'Frequency (Hz)', 'FontSize', D.st.label_fs);
        else
            set(ax(c), 'XTickLabel', []);
        end
    end
    yline(ax(2), 1, '-', 'Color', [0.6 0.6 0.6]); yline(ax(3), 1, '-', 'Color', [0.6 0.6 0.6]);
    set(ax(1), 'YLim', [0 2.2]); set(ax(2), 'YLim', [0 2.2]); set(ax(3), 'YLim', [0 1.6]);
    ylabel(ax(1), 'Bulk radius of M(i\omega)', 'FontSize', D.st.label_fs);
    ylabel(ax(2), 'Rightmost Re eig M(i\omega)', 'FontSize', D.st.label_fs);
    ylabel(ax(3), 'Loop gain', 'FontSize', D.st.label_fs);
end
file = save_fig(fig, D, 'freqresp_network');

head = {'network', 'condition', 'R frozen', 'R(0) both', 'R(0) adaptation alone', 'R(0) depression alone', ...
    'f where R reaches 90% of frozen (Hz)', 'peak loop gain', 'at Hz', 'peak, adaptation alone', 'peak, depression alone', ...
    'predicted bulk abscissa (1/s)', 'full Jacobian abscissa (1/s)', 'their correlation over 12 states', ...
    'continuation Re (1/s)', 'continuation Im (rad/s)', 'continuation Hz', 'frozen bulk rotates up to (Hz)', ...
    'K1 frozen', 'K1(0)', 'peak mean-mode loop gain'};
caption = [ ...
    'The network frequency response. With the slow variables eliminated neuron by neuron, a mode going as exp(z t) sees the matrix M(z) = W .* d(z): ' ...
    'the weights with every presynaptic neuron''s column scaled by that neuron''s complex gain at z. z is an eigenvalue of the full Jacobian exactly when ' ...
    '1 + tau_d z is an eigenvalue of M(z). The curves are for z = i 2 pi f, medians over the 12 calibration states of the outlier runs (seed [7 8]; 500 neurons). ' ...
    'ROWS: the reference network; E-to-E mean times 1.5. COLOR: dark red no adaptation (no slow variables, so no dependence on frequency), orange one timescale, blue multiple timescales. ' ...
    'LEFT: the bulk radius of M(i omega) predicted from the gains (Perron root of the block variances times the summed squared gain moduli; band: quartiles over states). ' ...
    'Gray: |1 + i omega tau_d|, the modulus an eigenvalue of M(i omega) must have for a mode at that frequency to be marginal (the dendritic low-pass). ' ...
    'MIDDLE: the largest real part among the eigenvalues of M(i omega); gray line at 1. Bottom row, dashed: the real part of the leading eigenvalue of the mean part K(i omega), where the structural outlier sits when it is outside the bulk. ' ...
    'RIGHT: the loop gain, bulk radius over |1 + i omega tau_d|; above 1 (gray line) the bulk holds growing modes rotating at that frequency. Bottom row, dashed: the same for the mean-part eigenvalue. ' ...
    'Triangles: the frequency of the continuation mode found in the outlier runs (median imaginary part over the states, in Hz). ' ...
    'Bars at the bottom: the range of rotation rates of the frozen bulk, 0 to R_frozen / tau_d (in Hz), per condition. ' ...
    'Read left to right: adaptation and depression lower the radius below a few tenths of a hertz and not above about 1 to 3 Hz; the dendrite cuts above 1/(2 pi tau_d) = 1.6 Hz; ' ...
    'what is left is a band near 0.5 to 1 Hz, which is where the continuation modes are. ' ...
    'Table: R in units of W (1 is marginal at zero frequency). "predicted bulk abscissa" is the rightmost point of the region R(z) >= |1 + tau_d z| of the complex plane, ' ...
    'per state, and is compared with the rightmost eigenvalue of the dense full Jacobian at the same state (median of each, and their correlation over the 12 states; ' ...
    'in the E-to-E times 1.5 network the full abscissa is the edge of the slow single-neuron band, not a network mode, so the two are not the same thing there).'];
write_sidecar(D, 'freqresp_network', caption, head, rows);
end

%% ------------------------------------------------------------------------
function file = fig_one_vs_multiple(D)
[~, B] = freqresp_single(D.opts.seed);
f = B.f;
modes = {'both', 'sfa', 'std'};
styles = {'-', '--', ':'};
r_show = 0.25;
gray = [0.55 0.55 0.55];

fig = new_fig([1 1 12.5 4], D);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) matched operating point
ax = nexttile(tl, 1); hold(ax, 'on');
m1 = B.matched(strcmp({B.matched.condition}, 'sfa1_std1') & [B.matched.r0] == r_show);
m3 = B.matched(strcmp({B.matched.condition}, 'sfa3_stdc2') & [B.matched.r0] == r_show);
ratio_abs = (m1.F .* abs(m1.G.both)) ./ (m3.F .* abs(m3.G.both));
plot(ax, f, median(ratio_abs, 1), '-', 'Color', gray, 'LineWidth', 4);
for m = 1:3
    ratio = abs(m1.G.(modes{m})) ./ abs(m3.G.(modes{m}));
    if m == 1
        band(ax, f, prctile(ratio, 25, 1), prctile(ratio, 75, 1), [0 0 0]);
    end
    plot(ax, f, median(ratio, 1), styles{m}, 'Color', 'k', 'LineWidth', D.st.line_lw);
end
finish_ratio(ax, D, f);
ylabel(ax, 'Gain, one / multiple timescales', 'FontSize', D.st.label_fs);

% (b) each condition's own operating point: the network's radius
ax = nexttile(tl, 2); hold(ax, 'on');
A1 = load_network(D, 'reference', 'sfa1_std1');
A3 = load_network(D, 'reference', 'sfa3_stdc2');
fn = A1.f;
plot(ax, fn, median(A1.R.both, 1) ./ median(A3.R.both, 1), '-', 'Color', gray, 'LineWidth', 4);
for m = 1:3
    n1 = median(A1.R.(modes{m}), 1) ./ median(A1.R_inf);
    n3 = median(A3.R.(modes{m}), 1) ./ median(A3.R_inf);
    plot(ax, fn, n1 ./ n3, styles{m}, 'Color', 'k', 'LineWidth', D.st.line_lw);
end
finish_ratio(ax, D, f);
ylabel(ax, 'Bulk radius, one / multiple timescales', 'FontSize', D.st.label_fs);

% (c) phase at the matched operating point
ax = nexttile(tl, 3); hold(ax, 'on');
for mt = {m1, m3}
    col = cond_color(D.st, mt{1}.condition);
    for m = 1:3
        plot(ax, f, median(angle(mt{1}.G.(modes{m})), 1) * 180 / pi, styles{m}, 'Color', col, 'LineWidth', D.st.line_lw);
    end
end
plot(ax, f, atan(2 * pi * f * A1.tau_d) * 180 / pi, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
style_axis(ax, D);
set(ax, 'XScale', 'log', 'XLim', f([1 end]), 'XTick', 10.^(-3:2), 'YLim', [0 80]);
xlabel(ax, 'Frequency (Hz)', 'FontSize', D.st.label_fs);
ylabel(ax, 'Phase lead (deg)', 'FontSize', D.st.label_fs);
file = save_fig(fig, D, 'freqresp_one_vs_multiple');

% Table: the ratios at chosen frequencies.
f_tab = [0.001 0.01 0.03 0.1 0.3 0.7 1 3 10];
rows = {};
for r0 = unique([B.matched.r0])
    a1 = B.matched(strcmp({B.matched.condition}, 'sfa1_std1') & [B.matched.r0] == r0);
    a3 = B.matched(strcmp({B.matched.condition}, 'sfa3_stdc2') & [B.matched.r0] == r0);
    for m = 1:3
        ratio = abs(a1.G.(modes{m})) ./ abs(a3.G.(modes{m}));
        rows(end + 1, :) = [{sprintf('matched, rate %.2f', r0), modes{m}, 'normalized', max(ratio(:))}, ...
            num2cell(interp1(log(f), median(ratio, 1), log(f_tab)))]; %#ok<AGROW>
    end
    ratio = (a1.F .* abs(a1.G.both)) ./ (a3.F .* abs(a3.G.both));
    rows(end + 1, :) = [{sprintf('matched, rate %.2f', r0), 'both', 'absolute', max(ratio(:))}, ...
        num2cell(interp1(log(f), median(ratio, 1), log(f_tab)))]; %#ok<AGROW>
end
for w = 1:numel(D.networks)
    A1 = load_network(D, D.networks{w}, 'sfa1_std1');
    A3 = load_network(D, D.networks{w}, 'sfa3_stdc2');
    for m = 1:3
        ratio = (median(A1.R.(modes{m}), 1) ./ median(A1.R_inf)) ./ (median(A3.R.(modes{m}), 1) ./ median(A3.R_inf));
        rows(end + 1, :) = [{['own, ' D.networks{w}], modes{m}, 'normalized', max(ratio)}, ...
            num2cell(interp1(log(A1.f), ratio, log(f_tab)))]; %#ok<AGROW>
    end
    ratio = median(A1.R.both, 1) ./ median(A3.R.both, 1);
    rows(end + 1, :) = [{['own, ' D.networks{w}], 'both', 'absolute', max(ratio)}, ...
        num2cell(interp1(log(A1.f), ratio, log(f_tab)))]; %#ok<AGROW>
end
head = [{'operating point', 'free', 'ratio of', 'largest over frequency (and neurons)'}, ...
    arrayfun(@(v) sprintf('%g Hz', v), f_tab, 'UniformOutput', false)];
caption = [ ...
    'One timescale against multiple timescales: the ratio of the two gains against frequency. Below 1 the single-timescale condition attenuates MORE than the multiple-timescale one at that frequency. ' ...
    'Black lines are ratios of NORMALIZED gains (each divided by its own frozen gain): solid adaptation and depression free, dashed adaptation alone, dotted depression alone. ' ...
    'The thick gray line is the ratio of ABSOLUTE gains (frozen gain included), which is what the network sees. ' ...
    sprintf('LEFT: matched operating point. Every excitatory neuron of the reference network at rate %.2f and slope 1, slow variables at their steady state for that rate ', r_show) ...
    '(`freqresp_matched_state`); median over the 250 neurons, each compared with itself (same fast adaptation constant and releasable-pool times in both conditions); band: quartiles, both free. ' ...
    'MIDDLE: each condition at its own operating point: the predicted bulk radius of the reference network (median over 12 states), normalized by its frozen radius (black) and absolute (gray). ' ...
    'RIGHT: the phase lead of the normalized gain at the matched operating point (median over neurons), orange one timescale, blue multiple; gray: the phase lag of the dendrite, atan(omega tau_d), which a mode''s gain must cancel. ' ...
    'WHAT IT SHOWS. Dashed black line, left: for adaptation alone the ratio never exceeds 1 (largest value over all neurons and frequencies in the table): one timescale carrying the whole budget on the fastest rung attenuates at least as much as the ladder at every frequency. ' ...
    'Depression goes the other way at low frequency (the chain''s reserve deepens the steady state) and the same way above about 0.1 Hz. ' ...
    'Normalized, the multiple-timescale gain is the larger one from about 0.1 Hz up (the single-timescale filter is the stronger one there, by up to about 25%). ' ...
    'The absolute gain at the matched operating point is nevertheless smaller with multiple timescales at every frequency, because the frozen gain itself (the vesicle fraction the slow reserve has set) is smaller. ' ...
    'At the conditions'' own operating points (middle) the two effects nearly cancel between 0.1 and 0.7 Hz (absolute ratio 0.95 to 1.1 in the reference network) and the frozen gain wins above.'];
write_sidecar(D, 'freqresp_one_vs_multiple', caption, head, rows);
end

function finish_ratio(ax, D, f)
yline(ax, 1, '-', 'Color', [0.6 0.6 0.6]);
style_axis(ax, D);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', f([1 end]), 'XTick', 10.^(-3:2), ...
    'YLim', [0.5 5], 'YTick', [0.5 1 2 4]);
xlabel(ax, 'Frequency (Hz)', 'FontSize', D.st.label_fs);
end

%% ------------------------------------------------------------------------
function file = fig_plane(D)
n_net = numel(D.networks); n_cond = numel(D.conditions);
k = D.opts.state;
fig = new_fig([1 1 12.5 7.5], D);
tl = tiledlayout(fig, n_net, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
for w = 1:n_net
    for c = 1:n_cond
        A = load_network(D, D.networks{w}, D.conditions{c});
        O = load_outlier(D, D.networks{w}, D.conditions{c});
        col = cond_color(D.st, D.conditions{c});
        ax = nexttile(tl, (w - 1) * n_cond + c); hold(ax, 'on');
        lam = double(O.dense(k).lam);
        lam = lam(imag(lam) >= 0);
        plot(ax, real(lam), imag(lam), '.', 'Color', [0.45 0.45 0.45], 'MarkerSize', 5);
        xline(ax, 0, '-', 'Color', [0.75 0.75 0.75]);

        % Frozen edge: |1 + tau_d z| = R_inf, a circle about -1/tau_d.
        th = linspace(0, pi, 200);
        zc = (A.R_inf(k) * exp(1i * th) - 1) / A.tau_d;
        plot(ax, real(zc), imag(zc), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);

        % Edge with the slow variables free: R(z) = |1 + tau_d z|.
        P = A.plane;
        [SG, WW] = ndgrid(P.sigma, P.w);
        h = double(squeeze(P.R(k, :, :))) - abs(1 + A.tau_d * complex(SG, WW));
        contour(ax, P.sigma, P.w, h', [0 0], 'Color', col, 'LineWidth', D.st.line_lw);

        if c > 1
            cont = O.dense(k).cont;
            plot(ax, real(cont), imag(cont), 'o', 'Color', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);
        end
        style_axis(ax, D);
        set(ax, 'XLim', [-12 8], 'YLim', [0 20]);
        if w == n_net; xlabel(ax, 'Re (1/s)', 'FontSize', D.st.label_fs); else; set(ax, 'XTickLabel', []); end
        if c == 1; ylabel(ax, 'Im (rad/s)', 'FontSize', D.st.label_fs); else; set(ax, 'YTickLabel', []); end

        % The mean-part root (the population mode from the 2 x 2 mean
        % matrix alone), where the frozen mean part has an eigenvalue
        % outside the frozen bulk.
        root = NaN;
        if c > 1 && real(A.K_inf(1, k)) > A.R_inf(k)
            [~, params, ws] = freqresp_model(D.networks{w}, D.conditions{c}, D.opts.seed);
            [~, Q] = freqresp_states(D.networks{w}, D.conditions{c}, D.opts.seed);
            root = freqresp_outlier_root(Q.S_cal(:, k), params, ws, O.dense(k).cont);
            plot(ax, real(root), imag(root), 'x', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
        end

        edge = plane_abscissa(A);
        ab = [O.dense.abscissa];
        rows(end + 1, :) = {D.networks{w}, D.conditions{c}, (A.R_inf(k) - 1) / A.tau_d, edge(k), ab(k), ...
            real(O.dense(k).cont), imag(O.dense(k).cont), real(root), imag(root) + 0 * real(root), median(edge), median(ab), corr_(edge, ab), ...
            mean(abs(edge - ab))}; %#ok<AGROW>
    end
end
file = save_fig(fig, D, 'freqresp_plane');

head = {'network', 'condition', 'frozen edge (1/s), this state', 'predicted edge (1/s), this state', ...
    'full Jacobian abscissa (1/s), this state', 'continuation Re', 'continuation Im (rad/s)', ...
    'mean-part root Re', 'mean-part root Im (rad/s)', ...
    'predicted edge, median of 12 states', 'full abscissa, median of 12 states', 'correlation over states', 'mean |difference| over states'};
caption = [ ...
    sprintf('Where the frequency response puts the full Jacobian''s network eigenvalues. Calibration state %d (t = 70.9 s), seed [7 8]; rows: reference network, E-to-E mean times 1.5; columns: no adaptation, one timescale, multiple timescales. ', k) ...
    'Gray dots: all eigenvalues of the dense full Jacobian from the outlier run (upper half plane; the spectrum is symmetric about the real axis). ' ...
    'Dashed: the frozen-gain edge, the circle |1 + tau_d z| = R_frozen about -1/tau_d = -10 1/s: where the bulk would end if the slow variables were held. ' ...
    'Colored line: the edge with the slow variables free, the curve R(z) = |1 + tau_d z|, with R(z) the bulk radius predicted from the neurons'' gains at the complex frequency z (`freqresp_radius`); ' ...
    'network eigenvalues are expected inside it. Black ring: the continuation mode the outlier run picked. Vertical line: zero. ' ...
    'Black cross (only where the frozen mean part has an eigenvalue outside the frozen bulk: E-to-E times 1.5 with multiple timescales): the root of kappa(z) = 1 + tau_d z, ' ...
    'kappa the leading eigenvalue of the 2 x 2 mean part K(z) (`freqresp_outlier_root`): where the population mode should be from the mean gains alone. There the ring is outside the bulk edge because it is the structural outlier, not a bulk mode. ' ...
    'Without adaptation the two edges are the same circle. With adaptation the edge is pulled left where the imaginary part is small (slowly rotating modes, which the slow variables follow) and returns to the frozen circle as the rotation rate grows: ' ...
    'the dent is the reach of adaptation. The rightmost point of the colored curve is off the real axis, at a few rad/s. ' ...
    'The structure hugging the negative real axis is the band of slow single-neuron modes (poles of the gains), not network modes; the prediction is a bulk edge for large networks and individual eigenvalues scatter about it (500 neurons). ' ...
    'Table: "predicted edge" is the rightmost point of the region R(z) >= |1 + tau_d z| (grid of 0.1 1/s by 0.25 rad/s, real part above -8).'];
write_sidecar(D, 'freqresp_plane', caption, head, rows);
end

%% ------------------------------------------------------------------------
function edge = plane_abscissa(A)
% Per state: the rightmost point of the region R(z) >= |1 + tau_d z|
% (real part above -8 1/s, to stay off the far slow band).
P = A.plane;
[SG, WW] = ndgrid(P.sigma, P.w);
target = abs(1 + A.tau_d * complex(SG, WW));
n_states = size(P.R, 1);
edge = nan(1, n_states);
for k = 1:n_states
    inside = squeeze(P.R(k, :, :)) >= target;
    inside(P.sigma < -8, :) = false;
    rows = find(any(inside, 2));
    if ~isempty(rows); edge(k) = P.sigma(rows(end)); end
end
end

function r = corr_(a, b)
a = a(:); b = b(:);
ok = isfinite(a) & isfinite(b);
if nnz(ok) < 3 || std(a(ok)) == 0 || std(b(ok)) == 0
    r = NaN;
else
    cc = corrcoef(a(ok), b(ok));
    r = cc(1, 2);
end
end

function A = load_network(D, network, condition)
L = load(fullfile(D.data, sprintf('network_%s_%s_s%d-%d.mat', network, condition, D.opts.seed(1), D.opts.seed(2))));
A = L.A;
end

function R = load_outlier(D, network, condition)
L = load(fullfile(D.outl, outlier_run_name(network, condition, D.opts.seed, '')), 'R');
R = L.R;
end

function fig = new_fig(pos, D)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', pos, 'Visible', D.opts.visible, ...
    'Tag', 'fig_frequency_response');
end

function style_axis(ax, D)
set(ax, 'FontSize', D.st.tick_fs, 'LineWidth', D.st.axis_lw, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top');
end

function band(ax, x, lo, hi, col)
patch(ax, [x, fliplr(x)], [lo, fliplr(hi)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

function col = cond_color(st, name)
% No adaptation in dark red (FR-004, Tom 2026-11-20), as fig_outliers
% does; the other conditions from manuscript_style.
if strcmp(name, 'no_adaptation')
    col = [0.75 0.10 0.15];
else
    col = st.condition_color(name);
end
end

function file = save_fig(fig, D, name)
file = fullfile(D.figs, [name '.png']);
exportgraphics(fig, file, 'Resolution', D.opts.resolution);
close(fig);
fprintf('[fig_frequency_response] wrote %s\n', file);
end

function write_sidecar(D, name, caption, head, rows)
fid = fopen(fullfile(D.figs, [name '.md']), 'w');
closer = onCleanup(@() fclose(fid));
fprintf(fid, '# %s\n\n%s\n\n', name, caption);
fprintf(fid, ['Made by `scripts/explorations/frequency_response/fig_frequency_response.m` from ' ...
    '`data/explorations/frequency_response/`.\n\n']);
fprintf(fid, '| %s |\n', strjoin(head, ' | '));
fprintf(fid, '|%s\n', repmat('---|', 1, numel(head)));
for r = 1:size(rows, 1)
    cells = cell(1, size(rows, 2));
    for q = 1:size(rows, 2)
        v = rows{r, q};
        if ischar(v)
            cells{q} = v;
        else
            cells{q} = sprintf('%.3g', v);
        end
    end
    fprintf(fid, '| %s |\n', strjoin(cells, ' | '));
end
end
