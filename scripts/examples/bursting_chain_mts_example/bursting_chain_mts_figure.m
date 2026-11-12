function out = bursting_chain_mts_figure(opts)
% BURSTING_CHAIN_MTS_FIGURE Simulate the bursting network under original, chain MTS and product MTS; one figure.
%
%   out = BURSTING_CHAIN_MTS_FIGURE(name, value, ...)
%
% The work behind example_bursting_chain_mts.m, which states the physics,
% the parameter values and the numbers to expect; read that header first.
% Kept as a function so the example's settings reach it as arguments and so
% it can be called with the figure invisible (a MATLAB with no desktop).
%
% Each condition is built by bursting_condition_model (the one definition
% of the three conditions) on the same preset, seed, input and noise path.
%
% out.figure       the figure handle
% out.stats(k)     per condition: name, burst count, rate, durations,
%                  fraction of neurons active at the peak, mean rates
% out.traces       the struct saved as traces.mat (see save_traces)
% out.files        files written (empty unless opts.save_figs)
%
% With save_figs, out_dir receives <fig_tag>.fig, <fig_tag>.png and
% traces.mat (variables listed in example_bursting_chain_mts.m).
%
% See also: example_bursting_chain_mts, bursting_condition_model,
%           SRNNCellTypePairs.std_route_factor

arguments
    opts.preset_name   (1,:) char    = 'bursting_pairs'
    opts.rng_seeds     (1,2) double  = [42 42]
    opts.T_settle      (1,1) double  = 10      % s simulated before t = 0, not shown
    opts.T_show        (1,1) double  = 100     % s shown and analysed, [0, T_show]
    opts.plot_deci     (1,1) double  = 4       % 400 Hz integration -> 100 Hz traces
    opts.input_config                = []      % [] = no external input at all
    opts.model_args    (1,:) cell    = {}      % extra constructor arguments, applied last
    opts.burst_thresh  (1,1) double  = 0.2     % population mean E rate
    opts.burst_min_dur (1,1) double  = 0.1     % s above threshold
    opts.burst_merge   (1,1) double  = 0.2     % s; shorter dips do not split a burst
    opts.active_thresh (1,1) double  = 0.1     % a neuron with r_i above this is active
    opts.visible       (1,1) logical = true
    opts.save_figs     (1,1) logical = false
    opts.out_dir       (1,:) char    = ''
    opts.fig_tag       (1,:) char    = 'bursting_chain_mts'
end

input_desc = 'input as given';
if isempty(opts.input_config), input_desc = 'no input'; end

%% The three conditions, from the one definition
cond_names = {'original', 'chain_mts', 'product_mts'};
nC = numel(cond_names);
D  = cell(1, nC);
V  = cell(1, nC);
for k = 1:nC
    [model, V{k}] = bursting_condition_model(cond_names{k}, ...
        'preset_name', opts.preset_name, 'rng_seeds', opts.rng_seeds, ...
        'T_settle', opts.T_settle, 'T_show', opts.T_show, ...
        'input_config', opts.input_config, 'plot_deci', opts.plot_deci, ...
        'model_args', opts.model_args);
    model.run();
    D{k} = condition_traces(model, opts);
end
labels = {sprintf('original (%s)', V{1}.preset_condition), 'chain MTS', 'product MTS'};
for k = 1:nC, D{k}.label = labels{k}; end

%% Bursts
stats = struct([]);
for k = 1:nC
    s = detect_bursts(D{k}, opts);
    s.name = labels{k};
    stats = [stats, s]; %#ok<AGROW>
end

%% The traces, as saved
traces = struct();
for k = 1:nC
    traces.(cond_names{k}) = struct( ...
        'label',          labels{k}, ...
        't',              D{k}.t, ...
        'r_E',            D{k}.rE, ...
        'r_I',            D{k}.rI, ...
        'mean_r_E',       mean(D{k}.rE, 1), ...
        'mean_r_I',       mean(D{k}.rI, 1), ...
        'mean_adaptation', D{k}.adapt, ...
        'mean_depression', D{k}.dep, ...
        'depression_routes', {D{k}.dep_routes}, ...
        'bursts',         stats(k), ...
        'values',         V{k});
end
traces.burst_definition = struct('thresh_mean_r_E', opts.burst_thresh, ...
    'min_duration_s', opts.burst_min_dur, 'merge_gap_s', opts.burst_merge, ...
    'active_thresh_r', opts.active_thresh);

%% The figure
fig = draw(D, stats, opts, input_desc);
out = struct('figure', fig, 'stats', stats, 'traces', traces, 'files', {{}});

if opts.save_figs
    if isempty(opts.out_dir)
        error('bursting_chain_mts_figure:NoOutDir', 'save_figs is set but out_dir is empty.');
    end
    if ~isfolder(opts.out_dir), mkdir(opts.out_dir); end
    f_fig = fullfile(opts.out_dir, [opts.fig_tag '.fig']);
    f_png = fullfile(opts.out_dir, [opts.fig_tag '.png']);
    f_mat = fullfile(opts.out_dir, 'traces.mat');
    % 1500 px wide on screen at 96 px/in; 110 dpi gives about 1720 px.
    exportgraphics(fig, f_png, 'Resolution', 110);
    % A .fig saved from an invisible figure reopens invisible; the
    % CreateFcn makes it show itself however it is opened.
    old_fcn = fig.CreateFcn;
    fig.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';
    savefig(fig, f_fig);
    fig.CreateFcn = old_fcn;
    save(f_mat, '-struct', 'traces');
    out.files = {f_fig, f_png, f_mat};
end
end

%% ------------------------------------------------------------------------
function T = condition_traces(model, opts)
% Every trace the figure and the statistics need, over [0, T_show].
p      = model.plot_data;
params = model.get_params();
keep   = p.t >= 0 & p.t <= opts.T_show;
T.t    = p.t(keep);
T.t    = T.t(:)';
T.sigma_u_noise = model.sigma_u_noise;
T.rE   = p.r.E(:, keep);
T.rI   = p.r.I(:, keep);

% Mean adaptation term (c/K) sum_k a_k over the E neurons (I has no SFA here).
K = size(p.a.E, 2);
if K > 0
    T.adapt = squeeze(mean(model.c(1) / K * sum(p.a.E(:, :, keep), 2), 1))';
else
    T.adapt = zeros(size(T.t));
end

% Mean depression factor over the depressed routes, from the model.
[factor, ~] = SRNNCellTypePairs.std_factors(p.b, params);
names = model.cell_type_names;
acc = []; route_names = {};
for i = 1:numel(names)
    for j = 1:numel(names)
        v = factor.(names{i}).(names{j});
        if isempty(v), continue; end
        acc = [acc; v(:, keep)]; %#ok<AGROW>
        route_names{end+1} = sprintf('%s->%s', names{i}, names{j}); %#ok<AGROW>
    end
end
T.dep = mean(acc, 1);
T.dep_routes = route_names;
end

function s = detect_bursts(T, opts)
% A burst: the population mean E rate above burst_thresh for at least
% burst_min_dur s; crossings separated by less than burst_merge s are one.
dt   = T.t(2) - T.t(1);
mE   = mean(T.rE, 1);
above = mE > opts.burst_thresh;
d = diff([false, above, false]);
on = find(d == 1); off = find(d == -1) - 1;
% merge short gaps
k = 1;
while k < numel(on)
    if (on(k+1) - off(k) - 1) * dt < opts.burst_merge
        off(k) = off(k+1); on(k+1) = []; off(k+1) = [];
    else
        k = k + 1;
    end
end
dur = (off - on + 1) * dt;
ok  = dur >= opts.burst_min_dur;
on = on(ok); off = off(ok); dur = dur(ok);

nB = numel(on);
peak_t = nan(1, nB); frac_peak = nan(1, nB); frac_any = nan(1, nB);
Rall = [T.rE; T.rI];
for b = 1:nB
    [~, i] = max(mE(on(b):off(b)));
    ip = on(b) + i - 1;
    peak_t(b)    = T.t(ip);
    frac_peak(b) = mean(Rall(:, ip) > opts.active_thresh);
    frac_any(b)  = mean(any(Rall(:, on(b):off(b)) > opts.active_thresh, 2));
end
s.n_bursts      = nB;
s.burst_rate    = nB / (T.t(end) - T.t(1));
s.t_on          = T.t(on);
s.t_off         = T.t(off);
s.peak_t        = peak_t;
s.duration      = dur;
s.median_dur    = median(dur);
s.median_ibi    = median(diff(peak_t));
s.frac_active_peak  = frac_peak;
s.median_frac_peak  = median(frac_peak);
s.median_frac_any   = median(frac_any);
s.mean_rE       = mean(T.rE(:));
s.mean_rI       = mean(T.rI(:));
end

function fig = draw(D, stats, opts, input_desc)
nC = numel(D);
vis = 'off'; if opts.visible, vis = 'on'; end
fig = figure('Color', 'w', 'Visible', vis, 'Units', 'pixels', ...
    'Position', [60 60 1500 900], 'NumberTitle', 'off', ...
    'Name', 'bursting network: original, chain MTS, product MTS');
tl = tiledlayout(fig, 4, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
cE = [0.85 0.33 0.10]; cI = [0.00 0.45 0.74];
ax = gobjects(4, nC);
ylim_adapt = [0, max(cellfun(@(T) max(T.adapt), D)) * 1.05 + eps];
ylim_dep   = [min(cellfun(@(T) min(T.dep), D)) * 0.95, 1];
nE = size(D{1}.rE, 1);
for k = 1:nC
    T = D{k}; s = stats(k);
    % rates image
    ax(1,k) = nexttile(tl, k);
    imagesc(ax(1,k), T.t, 1:(nE + size(T.rI,1)), [T.rE; T.rI], [0 1]);
    colormap(ax(1,k), gray); set(ax(1,k), 'YDir', 'reverse'); hold(ax(1,k), 'on');
    yline(ax(1,k), nE + 0.5, '-', 'Color', cI, 'LineWidth', 1);
    hold(ax(1,k), 'off');
    title(ax(1,k), sprintf('%s\n%d bursts, mean r_E %.3f, r_I %.3f', ...
        strrep(T.label, '_', '\_'), s.n_bursts, s.mean_rE, s.mean_rI), 'Interpreter', 'tex');
    if k == 1, ylabel(ax(1,k), sprintf('neuron (E 1-%d, I below)', nE)); end
    % population mean rates, bursts marked
    ax(2,k) = nexttile(tl, nC + k); hold(ax(2,k), 'on');
    for b = 1:s.n_bursts
        patch(ax(2,k), [s.t_on(b) s.t_off(b) s.t_off(b) s.t_on(b)], [0 0 1 1], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    % I first, so E (type 1) is drawn on top; legend in natural order.
    hI = plot(ax(2,k), T.t, mean(T.rI, 1), 'Color', cI, 'DisplayName', 'I');
    hE = plot(ax(2,k), T.t, mean(T.rE, 1), 'Color', cE, 'DisplayName', 'E');
    yline(ax(2,k), opts.burst_thresh, ':k', 'HandleVisibility', 'off');
    hold(ax(2,k), 'off'); ylim(ax(2,k), [0 1]);
    if k == 1
        ylabel(ax(2,k), 'mean rate');
        legend(ax(2,k), [hE hI], 'Location', 'northeast', 'Box', 'off');
    end
    % adaptation
    ax(3,k) = nexttile(tl, 2*nC + k);
    plot(ax(3,k), T.t, T.adapt, 'Color', cE);
    ylim(ax(3,k), ylim_adapt);
    if k == 1, ylabel(ax(3,k), 'mean (c/K)\Sigma_k a_k, E'); end
    % depression
    ax(4,k) = nexttile(tl, 3*nC + k);
    plot(ax(4,k), T.t, T.dep, 'Color', [0.3 0.3 0.3]);
    ylim(ax(4,k), ylim_dep);
    if k == 1
        ylabel(ax(4,k), sprintf('mean STD factor\n(%s)', strjoin(T.dep_routes, ', ')));
    end
    xlabel(ax(4,k), 'time (s)');
end
linkaxes(ax(:), 'x'); xlim(ax(1,1), [0 opts.T_show]);
for r = 1:4, linkaxes(ax(r,:), 'y'); end
for a = ax(:)', set(a, 'Box', 'off', 'TickDir', 'out'); end
title(tl, sprintf(['%s, seed %s, %s, sigma_u_noise %g. Burst: mean E rate > %.2f for >= %.2f s ' ...
    '(gaps < %.2f s merged); shaded.'], opts.preset_name, mat2str(opts.rng_seeds), ...
    input_desc, D{1}.sigma_u_noise, ...
    opts.burst_thresh, opts.burst_min_dur, opts.burst_merge), 'Interpreter', 'none', 'FontSize', 10);
end
