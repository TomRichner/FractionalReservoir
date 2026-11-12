% RUN_P4_ABERRANT_MEASURE Block B pilot 4: one measure of aberrant activity for bursts and IED-like events.
%
% Evaluates candidate measures of aberrant population activity, each
% computed from per-neuron rates alone (no modules, no per-network
% thresholds), against the two networks' own detectors:
%   bursting network   bursts of the mean E rate (pilot 0b), 33 / 34 / 0
%   run13 network      module-recruitment events (IED detector v1, pilot
%                      0a), 30 / 12 / 1
% on 16 traces:
%   6 saved           original, chain MTS, product MTS of each network
%                     (figs/explore_B/bursting_chain_mts, run13_chain_mts)
%   4 bursting arms   pilot 2's intact / frozen x stimulus / none, re-run
%                     here with per-neuron rates (p4_simulate_extra)
%   4 run13 windows   original and chain MTS, 50 s without and 50 s with a
%                     uniform input of 0.20 (p4_simulate_extra)
%   2 paper network   mu7chain, no adaptation and chained MTS, no input
% The candidates are defined in p4_candidates; the recommended one is
% coactivation_excess.
%
% Run setup_paths first (cwd at the repo root). Serial, no pool. About
% 45 s of simulation (cached in p4_extra_traces.mat) and 30 s of analysis.
% Output, in figs/explore_B/p4_aberrant_measure/ of the MAIN checkout:
%   p4_candidates_table.csv / .mat   candidate by trace, with the detector
%                                    counts
%   p4_robustness.csv                the two halves of each window, k = 2
%                                    and 4, and for the recommended measure
%                                    the floor and z over a factor of two
%   p4_recommended_unstimulated.png / .fig, p4_recommended_bursting_arms,
%   p4_recommended_stim_and_paper, p4_candidates_summary
%
% See also: coactivation_excess, p4_candidates, p4_simulate_extra

main_root = 'C:\Users\m218089\Desktop\github_repos\FractionalReservoir';
out_dir = fullfile(main_root, 'figs', 'explore_B', 'p4_aberrant_measure');
if ~isfolder(out_dir), mkdir(out_dir); end
t_start = tic;

%% The traces
B = load(fullfile(main_root, 'figs', 'explore_B', 'bursting_chain_mts', 'traces.mat'));
R = load(fullfile(main_root, 'figs', 'explore_B', 'run13_chain_mts', 'traces.mat'));
extra_file = fullfile(out_dir, 'p4_extra_traces.mat');
if isfile(extra_file)
    S = load(extra_file);  X = S.X;
else
    X = p4_simulate_extra(extra_file);
end

conds = {'original', 'chain_mts', 'product_mts'};
cond_labels = {'original', 'chain MTS', 'product MTS'};
T = struct('name', {}, 'label', {}, 'group', {}, 't', {}, 'r', {}, 'is_E', {}, ...
    'det_count', {}, 'det_times', {}, 'det_note', {});
for k = 1:3
    s = B.(conds{k});
    T(end + 1) = struct('name', ['bursting_' conds{k}], 'label', ['bursting, ' cond_labels{k}], ...
        'group', 'bursting saved', 't', s.t, 'r', [s.r_E; s.r_I], ...
        'is_E', [true(size(s.r_E, 1), 1); false(size(s.r_I, 1), 1)], ...
        'det_count', s.bursts.n_bursts, 'det_times', s.bursts.peak_t, 'det_note', 'burst detector'); %#ok<SAGROW>
end
for k = 1:3
    s = R.(['traces_' conds{k}]);
    isE = false(size(s.rate, 1), 1);  isE(s.E_idx) = true;
    T(end + 1) = struct('name', ['run13_' conds{k}], 'label', ['run13, ' cond_labels{k}], ...
        'group', 'run13 saved', 't', s.t_rate, 'r', double(s.rate), 'is_E', isE, ...
        'det_count', height(s.events), 'det_times', s.events.Time', 'det_note', 'IED detector v1'); %#ok<SAGROW>
end
arm_labels = containers.Map( ...
    {'bursting_intact_nostim', 'bursting_intact_stim', 'bursting_frozen_stim', 'bursting_frozen_nostim', ...
     'run13_original_nostim', 'run13_original_stim', 'run13_chain_mts_nostim', 'run13_chain_mts_stim', ...
     'mu7chain_no_adaptation', 'mu7chain_sfa3_stdc2'}, ...
    {'bursting, intact, no stimulus', 'bursting, intact, stimulated', 'bursting, frozen, stimulated', ...
     'bursting, frozen, no stimulus', 'run13 original, no input', 'run13 original, input 0.20', ...
     'run13 chain MTS, no input', 'run13 chain MTS, input 0.20', ...
     'paper network, no adaptation', 'paper network, chained MTS'});
for i = 1:numel(X)
    w = X(i).t >= X(i).window(1) & X(i).t <= X(i).window(2);
    tt = X(i).t(w);  rr = double(X(i).r(:, w));
    switch X(i).network
        case 'bursting arms'
            [n_det, t_det] = burst_detector(tt, mean(rr(X(i).is_E, :), 1));
            note = 'burst detector';
            if n_det == 1 && min(mean(rr(X(i).is_E, :), 1)) > 0.2
                note = 'burst detector: the mean E rate never falls below 0.2, one unbroken crossing';
            end
        case 'run13 stimulated'
            d = ied_detect_events(tt, double(X(i).x(:, w)), rr, X(i).groups, struct());
            n_det = d.metrics.event_count;  t_det = [d.events.time];
            note = 'IED detector v1, this window''s own thresholds';
        otherwise
            n_det = NaN;  t_det = [];  note = 'no detector (no modules, no tuned threshold)';
    end
    T(end + 1) = struct('name', X(i).name, 'label', arm_labels(X(i).name), 'group', X(i).network, ...
        't', tt, 'r', rr, 'is_E', X(i).is_E, 'det_count', n_det, 'det_times', t_det, 'det_note', note); %#ok<SAGROW>
end
nT = numel(T);

%% Candidate by trace
rows = cell(nT, 1);
for i = 1:nT
    v = p4_candidates(T(i).r, T(i).is_E, 3);
    v.class = string(v.class);
    head = struct('trace', string(T(i).name), 'group', string(T(i).group), 'n', size(T(i).r, 1), ...
        'window_s', T(i).t(end) - T(i).t(1), 'detector_count', T(i).det_count, ...
        'detector_per_min', T(i).det_count / (T(i).t(end) - T(i).t(1)) * 60);
    % Agreement in time with the detector: the share of its events that have
    % c(t) above the chance level within 0.15 s, and the share of the
    % measure's own episodes (c > L) that have a detector event within 0.15 s.
    [~, ~, det] = coactivation_excess(T(i).r, 'dt', T(i).t(2) - T(i).t(1));
    t_above = T(i).t(det.c > det.level);
    hit = NaN;  matched = NaN;
    if ~isempty(T(i).det_times) && ~isnan(T(i).det_count)
        hit = mean(arrayfun(@(te) any(abs(t_above - te) <= 0.15), T(i).det_times));
        if ~isempty(det.event_on)
            matched = mean(arrayfun(@(a, b) any(T(i).det_times >= T(i).t(a) - 0.15 & ...
                T(i).det_times <= T(i).t(b) + 0.15), det.event_on, det.event_off));
        end
    end
    agree = table(det.events_per_min, hit, matched, 'VariableNames', ...
        {'coact_episodes_per_min', 'detector_events_with_coactivation', 'coact_episodes_with_detector_event'});
    rows{i} = [struct2table(head), struct2table(v), agree, ...
        table(string(T(i).det_note), 'VariableNames', {'detector_note'})];
end
tbl = vertcat(rows{:});
writetable(tbl, fullfile(out_dir, 'p4_candidates_table.csv'));
disp(tbl(:, {'trace', 'detector_count', 'coact_excess', 'coact_time', 'coact_q99', 'coact_time10', ...
    'pop_time', 'pop_mass', 'synchrony', 'skewness', 'mean_rate', 'pinned_fraction', 'class'}));
disp(tbl(:, {'trace', 'detector_per_min', 'coact_episodes_per_min', 'detector_events_with_coactivation', ...
    'coact_episodes_with_detector_event', 'gated_mass', 'p99_over_median', 'mean_rate_E'}));

%% Robustness: halves of the window, k over a factor of two; floor and z for the recommended measure
cand = {'pop_time', 'pop_mass', 'coact_time10', 'coact_q99', 'coact_excess', 'coact_time', ...
    'gated_mass', 'synchrony', 'skewness', 'p99_over_median'};
rob = table();
for i = 1:nT
    nt = size(T(i).r, 2);  h = floor(nt / 2);
    v = {p4_candidates(T(i).r, T(i).is_E, 3), ...
         p4_candidates(T(i).r(:, 1:h), T(i).is_E, 3), p4_candidates(T(i).r(:, h + 1:end), T(i).is_E, 3), ...
         p4_candidates(T(i).r, T(i).is_E, 2), p4_candidates(T(i).r, T(i).is_E, 4)};
    for c = 1:numel(cand)
        vals = cellfun(@(s) s.(cand{c}), v);
        rob = [rob; table(string(T(i).name), string(cand{c}), vals(1), vals(2), vals(3), vals(4), vals(5), ...
            'VariableNames', {'trace', 'candidate', 'full_k3', 'first_half', 'second_half', 'k2', 'k4'})]; %#ok<AGROW>
    end
end
writetable(rob, fullfile(out_dir, 'p4_robustness.csv'));
rec = table();
for i = 1:nT
    a = @(varargin) coactivation_excess(T(i).r, varargin{:});
    rec = [rec; table(string(T(i).name), a(), a('floor', 0.05), a('floor', 0.2), a('z', 3), a('z', 7), ...
        a('k', 2), a('k', 4), ...
        'VariableNames', {'trace', 'A', 'floor_0p05', 'floor_0p2', 'z_3', 'z_7', 'k_2', 'k_4'})]; %#ok<AGROW>
end
writetable(rec, fullfile(out_dir, 'p4_recommended_constants.csv'));
disp(rec);
tic; coactivation_excess(T(end).r); t_one = toc;
fprintf('coactivation_excess on %d x %d: %.2f s\n', size(T(end).r), t_one);

%% Figures
files = {};
% 1. the six unstimulated traces
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1650 1050]);
tl = tiledlayout(fig, 6, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
order = [1 4 2 5 3 6];           % rows: condition; columns: network
for j = 1:6
    i = order(j);  col = 2 - mod(j, 2);  row = ceil(j / 2);
    ax1 = nexttile(tl, (2 * row - 2) * 2 + col);
    ax2 = nexttile(tl, (2 * row - 1) * 2 + col);
    trace_panels(ax1, ax2, T(i), tbl(i, :));
    if row == 3, xlabel(ax2, 'time (s)'); end
end
title(tl, ['Recommended measure on the six unstimulated traces: mean E rate (top of each pair) and co-activation c(t) (bottom); ' ...
    'dashed = chance level L, triangles = the network''s own detector'], 'FontSize', 10);
files = [files, save_fig(fig, out_dir, 'p4_recommended_unstimulated')];

% 2. the four bursting arms
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1650 520]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'tight', 'Padding', 'compact');
for j = 1:4
    i = 6 + j;
    ax1 = nexttile(tl, j);  ax2 = nexttile(tl, 4 + j);
    trace_panels(ax1, ax2, T(i), tbl(i, :));
    ylim(ax1, [0 1]);  ylim(ax2, [0 1]);  xlabel(ax2, 'time (s)');
end
title(tl, ['The four arms of the bursting network (pilot 2), 25 to 60 s: mean E rate (top) and co-activation c(t) (bottom). ' ...
    'A is 0 for the steady, the saturated and the silent arm; the class comes from the companion'], 'FontSize', 10);
files = [files, save_fig(fig, out_dir, 'p4_recommended_bursting_arms')];

% 3. run13 with and without input, and the paper's network
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1650 780]);
tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'tight', 'Padding', 'compact');
order = [11 13 15; 12 14 16];
for rr = 1:2
    for cc = 1:3
        i = order(rr, cc);
        ax1 = nexttile(tl, (2 * rr - 2) * 3 + cc);  ax2 = nexttile(tl, (2 * rr - 1) * 3 + cc);
        trace_panels(ax1, ax2, T(i), tbl(i, :));
        if rr == 2, xlabel(ax2, 'time (s)'); end
    end
end
title(tl, ['Beyond the saved traces: run13 original and chain MTS without (top pair) and with (bottom pair) a uniform input of 0.20; ' ...
    'right column: the paper''s network with no adaptation (top) and chained MTS (bottom)'], 'FontSize', 10);
files = [files, save_fig(fig, out_dir, 'p4_recommended_stim_and_paper')];

% 4. summary: candidates against the detector counts
show = {'coact_excess', 'coact_time', 'coact_q99', 'coact_time10', 'pop_time', 'pop_mass', 'synchrony', 'skewness'};
show_lab = {'3c excess co-activation A (recommended)', '3d time above chance level', '3b 99th pct of c', ...
    '3a time c > 0.10', '1 time in excursion (pop. mean)', '2 excursion mass (pop. mean)', '4 synchrony', '5 skewness (pop. mean)'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1650 560]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
sets = {1:3, 4:6};  set_names = {'bursting network', 'run13 network'};
cols = [0 0 0; lines(7); 0.45 0.45 0.45];
for sidx = 1:2
    ax = nexttile(tl);
    idx = sets{sidx};
    M = zeros(3, 1 + numel(show));
    M(:, 1) = tbl.detector_count(idx) / tbl.detector_count(idx(1));
    for c = 1:numel(show)
        M(:, c + 1) = tbl.(show{c})(idx) / tbl.(show{c})(idx(1));
    end
    b = bar(ax, M, 'grouped');
    for c = 1:numel(b), b(c).FaceColor = cols(c, :); b(c).EdgeColor = 'none'; end
    set(ax, 'XTickLabel', cond_labels);  ylim(ax, [0 1.4]);  yline(ax, 1, ':k');
    ylabel(ax, 'reading / reading in original');
    title(ax, sprintf('%s: detector %d : %d : %d', set_names{sidx}, tbl.detector_count(idx)));
    if sidx == 1
        legend(ax, [{'detector count'}, show_lab], 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
    end
end
ax = nexttile(tl);
pairs = [7 8; 11 12; 13 14];  pair_names = {'bursting, intact', 'run13 original', 'run13 chain MTS'};
M = nan(3, 1 + numel(show));
for p = 1:3
    M(p, 1) = tbl.detector_count(pairs(p, 2)) / tbl.detector_count(pairs(p, 1));
    for c = 1:numel(show)
        M(p, c + 1) = tbl.(show{c})(pairs(p, 2)) / tbl.(show{c})(pairs(p, 1));
    end
end
b = bar(ax, M, 'grouped');
for c = 1:numel(b), b(c).FaceColor = cols(c, :); b(c).EdgeColor = 'none'; end
set(ax, 'XTickLabel', pair_names);  ylim(ax, [0 1.4]);  yline(ax, 1, ':k');
ylabel(ax, 'stimulated / unstimulated');
title(ax, {'under stimulation (a bar at 1 or above = no suppression seen)', ...
    'bursting: the burst detector reads the steady rate as 1 burst'}, 'FontSize', 9);
title(tl, 'Candidates against the detector counts, each normalised to its own reading in the reference trace (bars cut at 1.4)', 'FontSize', 10);
files = [files, save_fig(fig, out_dir, 'p4_candidates_summary')];

run_time_s = toc(t_start);
save(fullfile(out_dir, 'p4_candidates_table.mat'), 'tbl', 'rob', 'rec', 'run_time_s');
fprintf('run time %.0f s\n', run_time_s);
cellfun(@(f) fprintf('  %s\n', f), files);

%% Local functions
function [n, t_peak] = burst_detector(t, mE)
% Pilot 0b's burst definition (bursting_chain_mts_figure, detect_bursts):
% the mean E rate above 0.2 for at least 0.1 s; dips shorter than 0.2 s do
% not split a burst.
dt = t(2) - t(1);
d = diff([false, mE > 0.2, false]);
on = find(d == 1);  off = find(d == -1) - 1;
k = 1;
while k < numel(on)
    if (on(k + 1) - off(k) - 1) * dt < 0.2
        off(k) = off(k + 1);  on(k + 1) = [];  off(k + 1) = [];
    else
        k = k + 1;
    end
end
ok = (off - on + 1) * dt >= 0.1;
on = on(ok);  off = off(ok);
n = numel(on);
t_peak = zeros(1, n);
for b = 1:n
    [~, j] = max(mE(on(b):off(b)));
    t_peak(b) = t(on(b) + j - 1);
end
end

function trace_panels(ax1, ax2, Tr, row)
% Top: mean E rate with the detector's events. Bottom: c(t) and the level.
[A, comp, det] = coactivation_excess(Tr.r);
m = mean(Tr.r(Tr.is_E, :), 1);
plot(ax1, Tr.t, m, 'Color', [0.80 0.30 0.10], 'LineWidth', 0.6);
hold(ax1, 'on');
if ~isempty(Tr.det_times)
    plot(ax1, Tr.det_times, max(m) * 1.08 * ones(size(Tr.det_times)), 'kv', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
end
xlim(ax1, [Tr.t(1) Tr.t(end)]);  ylim(ax1, [0, max(max(m) * 1.18, 0.05)]);
set(ax1, 'XTickLabel', [], 'FontSize', 8);  ylabel(ax1, 'mean r_E');
if isnan(row.detector_count), dtxt = 'no detector'; else, dtxt = sprintf('detector %d', row.detector_count); end
title(ax1, {Tr.label, sprintf('%s; A = %.2f x 10^{-3}; mean rate %.3f; class %s', ...
    dtxt, 1e3 * A, comp.mean_rate, comp.class)}, 'FontSize', 8, 'FontWeight', 'normal');
plot(ax2, Tr.t, det.c, 'Color', [0.10 0.35 0.70], 'LineWidth', 0.6);
hold(ax2, 'on');
yline(ax2, det.level, '--k');
xlim(ax2, [Tr.t(1) Tr.t(end)]);  ylim(ax2, [0, max(max(det.c) * 1.1, 2 * det.level)]);
set(ax2, 'FontSize', 8);  ylabel(ax2, 'c(t)');
end

function files = save_fig(fig, out_dir, tag)
png = fullfile(out_dir, [tag '.png']);
exportgraphics(fig, png, 'Resolution', 100);
fig.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';
savefig(fig, fullfile(out_dir, [tag '.fig']));
close(fig);
files = {png};
end
