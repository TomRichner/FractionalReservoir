function out = p5_dc_vs_pulses_compute(S)
% P5_DC_VS_PULSES_COMPUTE The work behind p5_dc_vs_pulses.m: simulate, measure, draw, save.
%
%   out = P5_DC_VS_PULSES_COMPUTE(S)
%
% S is the settings struct the script builds (its header says what each
% field means and what to expect). Returns
%   out.cases    struct array, one per simulation: condition, form, mean
%                level, the measures before and during, the input check
%   out.table    the same as a table (one row per case)
%   out.figs     figure handles (one per mean level, the original's, the summary)
%   out.files    files written (empty unless S.save_figs)
%
% The conditions come from bursting_condition_model (the one definition);
% the input from pulse_train_input. The burst definition is pilot 0b's
% (bursting_chain_mts_figure), restated in detect_bursts below because the
% one there is a local function and cannot be called.
%
% See also: p5_dc_vs_pulses, pulse_train_input, bursting_condition_model

%% The cases
cases = struct([]);
for c = 1:numel(S.conditions)
    cond = S.conditions{c};
    for m = 1:numel(S.mean_levels)
        for f = 1:numel(S.forms_Hz)
            fp = S.forms_Hz(f);
            if strcmp(cond, 'original') && ~ismember(fp, S.original_forms_Hz), continue; end
            cases(end+1).condition = cond; %#ok<AGROW>
            cases(end).pulse_freq = fp;
            cases(end).pulse_mean = S.mean_levels(m);
        end
    end
end
nK = numel(cases);

%% Simulate and measure
t0 = tic;
for k = 1:nK
    ck = cases(k);
    ic = struct('n_steps', 1, 'step_density', struct('E', 0, 'I', 0), 'amp', 0, ...
        'no_stim_pattern', true, 'intrinsic_drive', 0, 'positive_only', false, ...
        'generator', @pulse_train_input, 'pulse_mean', ck.pulse_mean, ...
        'pulse_freq', ck.pulse_freq, 'pulse_duty', S.duty, ...
        't_on', S.t_on, 't_off', S.T_show);
    [model, V] = bursting_condition_model(ck.condition, 'rng_seeds', S.rng_seeds, ...
        'T_settle', S.T_settle, 'T_show', S.T_show, 'input_config', ic, ...
        'plot_deci', S.plot_deci, 'model_args', {'fs', S.fs});
    model.run();
    D = extract(model, S);
    cases(k).tag    = case_tag(ck);
    cases(k).label  = case_label(ck);
    cases(k).values = V;
    cases(k).D      = D;
    cases(k).input_check = input_check(model, ck, S);
    cases(k).bursts = detect_bursts(D.t, D.mE, S);
    cases(k).before = window_measures(D, cases(k).bursts, S.win_before, ck, S);
    cases(k).during = window_measures(D, cases(k).bursts, S.win_during, ck, S);
    vprintf('minimal', 'minimal', '%2d/%d %-28s bursts/s %.3f -> %.3f, mean r_E %.4f -> %.4f (%.0f s)\n', ...
        k, nK, cases(k).label, cases(k).before.burst_rate, cases(k).during.burst_rate, ...
        cases(k).before.mean_rE, cases(k).during.mean_rE, toc(t0));
    clear model
end
run_time = toc(t0);

%% The before windows of one condition must be identical (input is zero there)
for cond = unique({cases.condition})
    idx = find(strcmp({cases.condition}, cond{1}));
    ref = cases(idx(1)).D;
    w = ref.t < S.t_on;
    dmax = 0;
    for i = idx(2:end)
        dmax = max(dmax, max(abs(cases(i).D.mE(w) - ref.mE(w))));
    end
    vprintf('minimal', 'minimal', 'before-window max |diff| of mean r_E across %s runs: %g\n', cond{1}, dmax);
    for i = idx, cases(i).before_identity_maxdiff = dmax; end
end

tbl = cases_table(cases);
out = struct('cases', cases, 'table', tbl, 'run_time_s', run_time, 'figs', gobjects(0), 'files', {{}});

%% Figures
figs = gobjects(0); tags = {};
for m = 1:numel(S.mean_levels)
    sel = strcmp({cases.condition}, 'chain_mts') & [cases.pulse_mean] == S.mean_levels(m);
    figs(end+1) = draw_level(cases(sel), S, sprintf(['chain MTS, mean input %.2g: DC versus pulses ' ...
        '(%g %% duty), stimulus on at %g s'], S.mean_levels(m), 100*S.duty, S.t_on)); %#ok<AGROW>
    tags{end+1} = sprintf('p5_chain_mts_mean%s', strrep(sprintf('%.2g', S.mean_levels(m)), '.', 'p')); %#ok<AGROW>
end
sel = strcmp({cases.condition}, 'original');
if any(sel)
    of = S.original_forms_Hz;
    figs(end+1) = draw_level(cases(sel), S, sprintf(['original (sfa3\\_std1), reference: %s ' ...
        'at both mean levels, stimulus on at %g s'], ...
        strjoin(arrayfun(@form_label, of, 'UniformOutput', false), ' and '), S.t_on));
    tags{end+1} = 'p5_original_reference';
end
figs(end+1) = draw_summary(cases, S);
tags{end+1} = 'p5_summary';
out.figs = figs;

%% Save
if S.save_figs
    if ~isfolder(S.out_dir), mkdir(S.out_dir); end
    files = {};
    for i = 1:numel(figs)
        f_png = fullfile(S.out_dir, [tags{i} '.png']);
        f_fig = fullfile(S.out_dir, [tags{i} '.fig']);
        exportgraphics(figs(i), f_png, 'Resolution', 110);
        old = figs(i).CreateFcn;
        figs(i).CreateFcn = 'set(gcbo, ''Visible'', ''on'')';
        savefig(figs(i), f_fig);
        figs(i).CreateFcn = old;
        files = [files, {f_png, f_fig}]; %#ok<AGROW>
    end
    traces = struct();
    for k = 1:nK
        c = cases(k); D = c.D;
        traces.(c.tag) = struct('label', c.label, 'condition', c.condition, ...
            'pulse_mean', c.pulse_mean, 'pulse_freq', c.pulse_freq, 'pulse_duty', S.duty, ...
            't', D.t, 'mean_r_E', D.mE, 'mean_r_I', D.mI, ...
            'mean_adaptation', D.adapt, 'mean_depression', D.dep, ...
            'depression_routes', {D.dep_routes}, 'input', D.u, ...
            't_input', D.t_u, 'input_full_rate', D.u_full, ...
            't_rates_100Hz', D.t100, 'r_E_100Hz', D.rE100, 'r_I_100Hz', D.rI100, ...
            'bursts', c.bursts, 'before', c.before, 'during', c.during, ...
            'input_check', c.input_check, 'values', c.values);
    end
    traces.settings = S;
    traces.measures_table = tbl;
    traces.run_time_s = run_time;
    f_mat = fullfile(S.out_dir, 'traces.mat');
    save(f_mat, '-struct', 'traces', '-v7.3');
    writetable(tbl, fullfile(S.out_dir, 'measures.csv'));
    files = [files, {f_mat, fullfile(S.out_dir, 'measures.csv')}];
    out.files = files;
end
end

%% ------------------------------------------------------------------------
function D = extract(model, S)
p = model.plot_data;
params = model.get_params();
keep = p.t >= 0 & p.t <= S.T_show + 1e-9;
D.t  = p.t(keep); D.t = D.t(:)';
rE = p.r.E(:, keep); rI = p.r.I(:, keep);
D.mE = mean(rE, 1); D.mI = mean(rI, 1);
step100 = round(0.01 / (D.t(2) - D.t(1)));
D.t100  = D.t(1:step100:end);
D.rE100 = single(rE(:, 1:step100:end));
D.rI100 = single(rI(:, 1:step100:end));
K = size(p.a.E, 2);
D.adapt = squeeze(mean(model.c(1) / K * sum(p.a.E(:, :, keep), 2), 1))';
[factor, ~] = SRNNCellTypePairs.std_factors(p.b, params);
names = model.cell_type_names; acc = []; rn = {};
for i = 1:numel(names)
    for j = 1:numel(names)
        v = factor.(names{i}).(names{j});
        if isempty(v), continue; end
        acc = [acc; v(:, keep)]; %#ok<AGROW>
        rn{end+1} = sprintf('%s->%s', names{i}, names{j}); %#ok<AGROW>
    end
end
D.dep = mean(acc, 1); D.dep_routes = rn;
% the input on neuron 1 (identical on every neuron), at the plot rate and at
% the integration rate over [0, T_show]
D.u = p.u.E(1, keep);
wk = model.t_ex >= -1e-9;
D.t_u = model.t_ex(wk)'; D.u_full = model.u_ex(1, wk);
end

function chk = input_check(model, ck, S)
t = model.t_ex(:)'; U = model.u_ex;
w = t >= S.t_on - 1e-9 & t < S.T_show - 1e-9;
chk.same_on_every_neuron = all(all(U == U(1, :)));
chk.mean_samples_on   = mean(U(1, w));
% the integrator sees the linear interpolant; its integral from one sample
% before t_on (where the first ramp starts) to T_show, over T_show - t_on
ww = t >= S.t_on - 1/S.fs - 1e-9 & t <= S.T_show + 1e-9;
chk.mean_interpolant_on = trapz(t(ww), U(1, ww)) / (S.T_show - S.t_on);
chk.mean_before = mean(U(1, t < S.t_on - 1e-9));
chk.max_input = max(U(1, :));
if ck.pulse_freq > 0
    chk.period_samples = S.fs / ck.pulse_freq;
    chk.width_samples  = S.duty * S.fs / ck.pulse_freq;
    chk.width_s        = S.duty / ck.pulse_freq;
    chk.n_pulses       = (S.T_show - S.t_on) * ck.pulse_freq;
else
    chk.period_samples = NaN; chk.width_samples = NaN; chk.width_s = NaN; chk.n_pulses = 0;
end
chk.rel_err_mean = abs(chk.mean_samples_on - ck.pulse_mean) / ck.pulse_mean;
end

function s = detect_bursts(t, mE, S)
% Pilot 0b's burst definition (bursting_chain_mts_figure>detect_bursts):
% the population mean E rate above S.burst_thresh for at least
% S.burst_min_dur s; crossings separated by less than S.burst_merge s are one.
dt = t(2) - t(1);
above = mE > S.burst_thresh;
d = diff([false, above, false]);
on = find(d == 1); off = find(d == -1) - 1;
k = 1;
while k < numel(on)
    if (on(k+1) - off(k) - 1) * dt < S.burst_merge
        off(k) = off(k+1); on(k+1) = []; off(k+1) = [];
    else
        k = k + 1;
    end
end
dur = (off - on + 1) * dt;
ok = dur >= S.burst_min_dur;
on = on(ok); off = off(ok);
s.t_on = t(on); s.t_off = t(off); s.duration = dur(ok);
s.peak_t = nan(1, numel(on)); s.peak_mE = nan(1, numel(on));
for b = 1:numel(on)
    [s.peak_mE(b), i] = max(mE(on(b):off(b)));
    s.peak_t(b) = t(on(b) + i - 1);
end
end

function M = window_measures(D, B, win, ck, S)
w = D.t >= win(1) - 1e-9 & D.t < win(2) - 1e-9;
in = B.t_on >= win(1) & B.t_on < win(2);
L = win(2) - win(1);
M.window = win;
M.n_bursts = nnz(in);
M.burst_rate = M.n_bursts / L;
M.median_burst_dur = median(B.duration(in));
M.median_ibi = median(diff(B.t_on(in)));
% share of the window spent inside a detected burst, whatever its onset: a
% burst merged from blips that starts before the window (the rate grazing
% the threshold at every pulse) is not counted in n_bursts but shows here
inb = false(size(D.t));
for b = 1:numel(B.t_on)
    inb = inb | (D.t >= B.t_on(b) & D.t <= B.t_off(b));
end
M.frac_time_in_burst = mean(inb(w));
M.mean_rE = mean(D.mE(w));
M.mean_rI = mean(D.mI(w));
M.max_mean_rE = max(D.mE(w));
M.mean_adapt = mean(D.adapt(w));
M.mean_dep = mean(D.dep(w));
M.min_dep = min(D.dep(w));
% locking to the pulse frequency (computed in the before window too, at the
% same frequency, as the reference for "no pulses")
fp = ck.pulse_freq;
if fp > 0
    x = D.mE(w); fsp = 1 / (D.t(2) - D.t(1));
    x = x - mean(x);
    N = numel(x);
    X = abs(fft(x)).^2;
    fr = (0:N-1) * fsp / N;
    pos = fr > 0 & fr < fsp / 2;
    [~, ib] = min(abs(fr - fp));
    M.power_frac_at_fp = X(ib) / sum(X(pos));
    harm = false(size(fr));
    for h = 1:floor((fsp/2 - 1e-9) / fp)
        [~, ih] = min(abs(fr - h * fp)); harm(ih) = true;
    end
    M.power_frac_harmonics = sum(X(harm & pos)) / sum(X(pos));
    % pulse-triggered average over whole cycles, phase from t_on
    P = round(fsp / fp);
    tw = D.t(w);
    first = find(mod(round((tw - S.t_on) * fsp), P) == 0, 1);
    nc = floor((numel(tw) - first + 1) / P);
    seg = reshape(D.mE(find(w, 1) + first - 1 + (0:nc*P-1)), P, nc);
    pta = mean(seg, 2)';
    M.pta = pta;
    M.pta_swing = max(pta) - min(pta);
    % burst onsets against the pulses
    if M.n_bursts > 0
        ph = mod(B.t_on(in) - S.t_on, 1 / fp) * fp;   % phase in cycles, 0 = pulse onset
        M.onset_phase = ph;
        M.vector_strength = abs(mean(exp(2i * pi * ph)));
        M.median_onset_latency_s = median(ph) / fp;
        M.bursts_per_pulse = M.n_bursts / (L * fp);
    else
        M.onset_phase = []; M.vector_strength = NaN; M.median_onset_latency_s = NaN;
        M.bursts_per_pulse = 0;
    end
else
    M.power_frac_at_fp = NaN; M.power_frac_harmonics = NaN; M.pta = [];
    M.pta_swing = NaN; M.onset_phase = []; M.vector_strength = NaN;
    M.median_onset_latency_s = NaN; M.bursts_per_pulse = NaN;
end
end

function tag = case_tag(ck)
if ck.pulse_freq == 0, form = 'DC'; else, form = sprintf('p%gHz', ck.pulse_freq); end
tag = sprintf('%s_%s_mean%s', ck.condition, form, strrep(sprintf('%.2g', ck.pulse_mean), '.', 'p'));
end

function s = form_label(fp)
if fp == 0, s = 'DC'; else, s = sprintf('%g Hz pulses', fp); end
end

function lab = case_label(ck)
lab = sprintf('%s, %s, mean %.2g', strrep(ck.condition, '_', ' '), form_label(ck.pulse_freq), ck.pulse_mean);
end

function T = cases_table(cases)
n = numel(cases);
col = @(f) arrayfun(f, cases(:));
T = table({cases.condition}', col(@(c) c.pulse_mean), col(@(c) c.pulse_freq), ...
    col(@(c) c.input_check.mean_samples_on), col(@(c) c.input_check.mean_interpolant_on), ...
    col(@(c) c.input_check.width_samples), ...
    col(@(c) c.before.burst_rate), col(@(c) c.during.burst_rate), ...
    col(@(c) c.before.frac_time_in_burst), col(@(c) c.during.frac_time_in_burst), ...
    col(@(c) c.before.mean_rE), col(@(c) c.during.mean_rE), ...
    col(@(c) c.before.mean_rI), col(@(c) c.during.mean_rI), ...
    col(@(c) c.before.max_mean_rE), col(@(c) c.during.max_mean_rE), ...
    col(@(c) c.before.mean_adapt), col(@(c) c.during.mean_adapt), ...
    col(@(c) c.before.mean_dep), col(@(c) c.during.mean_dep), ...
    col(@(c) c.before.power_frac_at_fp), col(@(c) c.during.power_frac_at_fp), ...
    col(@(c) c.during.power_frac_harmonics), ...
    col(@(c) c.before.pta_swing), col(@(c) c.during.pta_swing), ...
    col(@(c) c.during.vector_strength), col(@(c) c.during.bursts_per_pulse), ...
    col(@(c) c.during.median_onset_latency_s), ...
    'VariableNames', {'condition', 'mean_input', 'pulse_Hz', 'input_mean_samples', ...
    'input_mean_interp', 'pulse_width_samples', ...
    'bursts_per_s_before', 'bursts_per_s_during', 'frac_in_burst_before', ...
    'frac_in_burst_during', 'mean_rE_before', 'mean_rE_during', ...
    'mean_rI_before', 'mean_rI_during', 'max_mean_rE_before', 'max_mean_rE_during', ...
    'mean_adapt_before', 'mean_adapt_during', 'mean_dep_before', 'mean_dep_during', ...
    'pow_frac_fp_before', 'pow_frac_fp_during', 'pow_frac_harm_during', ...
    'pta_swing_before', 'pta_swing_during', 'vector_strength_during', ...
    'bursts_per_pulse_during', 'onset_latency_s_during'});
assert(height(T) == n);
end

function fig = draw_level(cs, S, ttl)
nC = numel(cs);
fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'pixels', ...
    'Position', [60 40 1500 1080], 'NumberTitle', 'off', 'Name', ttl);
tl = tiledlayout(fig, 5, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
cE = [0.85 0.33 0.10]; cI = [0.00 0.45 0.74]; cU = [0.55 0.55 0.55];
ax = gobjects(5, nC);
ya = [0, max(arrayfun(@(c) max(c.D.adapt), cs)) * 1.05 + eps];
yd = [min(arrayfun(@(c) min(c.D.dep), cs)) * 0.97, 1];
nE = size(cs(1).D.rE100, 1);
for k = 1:nC
    D = cs(k).D; B = cs(k).bursts; Mb = cs(k).before; Md = cs(k).during;
    ax(1,k) = nexttile(tl, k);
    imagesc(ax(1,k), D.t100, 1:(nE + size(D.rI100, 1)), double([D.rE100; D.rI100]), [0 1]);
    colormap(ax(1,k), gray); set(ax(1,k), 'YDir', 'reverse'); hold(ax(1,k), 'on');
    yline(ax(1,k), nE + 0.5, '-', 'Color', cI);
    xline(ax(1,k), S.t_on, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    hold(ax(1,k), 'off');
    title(ax(1,k), sprintf(['%s\nbursts/s %.2f \\rightarrow %.2f, mean r_E %.3f \\rightarrow %.3f'], ...
        cs(k).label, Mb.burst_rate, Md.burst_rate, Mb.mean_rE, Md.mean_rE), 'FontSize', 9);
    if k == 1, ylabel(ax(1,k), sprintf('neuron (E 1-%d, I below)', nE)); end

    ax(2,k) = nexttile(tl, nC + k); hold(ax(2,k), 'on');
    for b = 1:numel(B.t_on)
        patch(ax(2,k), [B.t_on(b) B.t_off(b) B.t_off(b) B.t_on(b)], [0 0 1 1], ...
            [0.88 0.88 0.88], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hI = plot(ax(2,k), D.t, D.mI, 'Color', cI);
    hE = plot(ax(2,k), D.t, D.mE, 'Color', cE);
    yline(ax(2,k), S.burst_thresh, ':k');
    xline(ax(2,k), S.t_on, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    hold(ax(2,k), 'off'); ylim(ax(2,k), [0 1]);
    if k == 1
        ylabel(ax(2,k), 'population mean rate');
        legend(ax(2,k), [hE hI], {'E', 'I'}, 'Location', 'northwest', 'Box', 'off');
    end

    ax(3,k) = nexttile(tl, 2*nC + k);
    plot(ax(3,k), D.t, D.adapt, 'Color', cE); ylim(ax(3,k), ya);
    xline(ax(3,k), S.t_on, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    if k == 1, ylabel(ax(3,k), 'mean (c/K)\Sigma_k a_k, E'); end

    ax(4,k) = nexttile(tl, 3*nC + k);
    plot(ax(4,k), D.t, D.dep, 'Color', [0.3 0.3 0.3]); ylim(ax(4,k), yd);
    xline(ax(4,k), S.t_on, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    if k == 1, ylabel(ax(4,k), sprintf('mean STD factor\n(%s)', strjoin(D.dep_routes, ', '))); end
    xlabel(ax(4,k), 'time (s)');

    ax(5,k) = nexttile(tl, 4*nC + k); hold(ax(5,k), 'on');
    wz = D.t_u >= S.zoom(1) & D.t_u <= S.zoom(2);
    area(ax(5,k), D.t_u(wz), D.u_full(wz), 'FaceColor', cU, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    wz2 = D.t >= S.zoom(1) & D.t <= S.zoom(2);
    plot(ax(5,k), D.t(wz2), D.mE(wz2), 'Color', cE, 'LineWidth', 1);
    yline(ax(5,k), S.burst_thresh, ':k');
    hold(ax(5,k), 'off'); xlim(ax(5,k), S.zoom);
    ylim(ax(5,k), [0, max(1, max(arrayfun(@(c) c.input_check.max_input, cs))) * 1.05]);
    xlabel(ax(5,k), sprintf('time (s), zoom %g-%g s', S.zoom));
    if k == 1, ylabel(ax(5,k), sprintf('input u (grey)\nmean r_E (orange)')); end
    if Md.pta_swing >= 0
        title(ax(5,k), sprintf('pulse-triggered swing of mean r_E %.3f', Md.pta_swing), ...
            'FontWeight', 'normal', 'FontSize', 9);
    end
end
linkaxes(ax(1:4, :), 'x'); xlim(ax(1,1), [0 S.T_show]);
for r = 1:5, linkaxes(ax(r, :), 'y'); end
for a = ax(:)', set(a, 'Box', 'off', 'TickDir', 'out'); end
title(tl, sprintf(['%s. bursting\\_pairs, seed %s, no noise, rk4 at %g Hz. Burst (shaded): mean r_E > %.2g ' ...
    'for >= %.2g s. Green line: onset. Measures: before %g-%g s, during %g-%g s.'], ttl, ...
    mat2str(S.rng_seeds), S.fs, S.burst_thresh, S.burst_min_dur, S.win_before, S.win_during), 'FontSize', 10);
end

function fig = draw_summary(cases, S)
fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'pixels', ...
    'Position', [60 60 1400 800], 'NumberTitle', 'off', 'Name', 'p5 summary');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
forms = S.forms_Hz; xl = arrayfun(@form_label, forms, 'UniformOutput', false);
xl = strrep(xl, ' pulses', ' Hz'); xl = strrep(xl, 'Hz Hz', 'Hz');
meas = {'burst_rate', 'bursts per second'; 'mean_rE', 'mean E rate'; ...
        'mean_dep', 'mean depression factor'; 'mean_adapt', 'mean adaptation term (c/K)\Sigma a, E'};
cols = [0.00 0.45 0.74; 0.85 0.33 0.10];
for i = 1:size(meas, 1)
    ax = nexttile(tl); hold(ax, 'on');
    h = gobjects(0); lab = {};
    for m = 1:numel(S.mean_levels)
        for cond = {'chain_mts', 'original'}
            sel = strcmp({cases.condition}, cond{1}) & [cases.pulse_mean] == S.mean_levels(m);
            if ~any(sel), continue; end
            cs = cases(sel);
            [~, xi] = ismember([cs.pulse_freq], forms);
            y = arrayfun(@(c) c.during.(meas{i,1}), cs);
            off = (m - 1.5) * 0.12 + strcmp(cond{1}, 'original') * 0.05;
            if strcmp(cond{1}, 'chain_mts')
                h(end+1) = plot(ax, xi + off, y, '-o', 'Color', cols(m, :), ...
                    'MarkerFaceColor', cols(m, :), 'LineWidth', 1.2); %#ok<AGROW>
                lab{end+1} = sprintf('chain MTS, mean %.2g', S.mean_levels(m)); %#ok<AGROW>
            else
                h(end+1) = plot(ax, xi + off, y, 's', 'Color', cols(m, :), ...
                    'MarkerSize', 8, 'LineWidth', 1.2); %#ok<AGROW>
                lab{end+1} = sprintf('original, mean %.2g', S.mean_levels(m)); %#ok<AGROW>
            end
        end
    end
    cb = cases(find(strcmp({cases.condition}, 'chain_mts'), 1)).before.(meas{i,1});
    h(end+1) = yline(ax, cb, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1); lab{end+1} = 'chain MTS before'; %#ok<AGROW>
    io = find(strcmp({cases.condition}, 'original'), 1);
    if ~isempty(io)
        h(end+1) = yline(ax, cases(io).before.(meas{i,1}), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1); %#ok<AGROW>
        lab{end+1} = 'original before'; %#ok<AGROW>
    end
    hold(ax, 'off');
    set(ax, 'XTick', 1:numel(forms), 'XTickLabel', xl, 'Box', 'off', 'TickDir', 'out');
    xlim(ax, [0.5 numel(forms) + 0.5]);
    yl = ylim(ax); ylim(ax, [0, yl(2) * 1.08]);
    if strcmp(meas{i,1}, 'mean_dep'), ylim(ax, [min(yl(1), 0.6) 1]); end
    ylabel(ax, meas{i,2}); title(ax, [meas{i,2} ', during stimulation'], 'FontWeight', 'normal');
    if i == 1, legend(ax, h, lab, 'Location', 'best', 'Box', 'off', 'FontSize', 8); end
end
title(tl, sprintf(['Pilot 5 summary: during %g-%g s (stimulus on at %g s; pulses %g %% duty, amplitude = mean / %g). ' ...
    'Lines: before %g-%g s.'], S.win_during, S.t_on, 100*S.duty, S.duty, S.win_before), 'FontSize', 10);
end
