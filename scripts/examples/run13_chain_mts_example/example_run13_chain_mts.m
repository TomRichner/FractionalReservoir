% EXAMPLE_RUN13_CHAIN_MTS The run13 IED-like network with multiple-timescale adaptation, chained and product depression.
%
% Pilot 0a, an exploratory look requested by Tom on 2026-11-12 at the start
% of planning block B: the frozen run13 network of the IED pilot
% (ied_base_config: 350 neurons, 245 E / 105 I, five weakly coupled
% 70-neuron modules, Wiener noise 0.025, no input) simulated for 100 s
% after the IED code's 10-s settling interval, under three conditions on
% the SAME network, setpoints, initial x and noise path. The conditions are
% defined in ONE place, run13_condition_model (this folder), which states
% every constant; in short:
%
%   original     run13 exactly as ied_build_network builds it: one SFA
%                timescale per neuron (E 0.40 s, I 0.75 s), one depression
%                pool on E outputs (tau_rec 2 s, tau_rel 0.25 s), none on I
%   chain MTS    SFA ladders E [0.40 1.58 10] s, I [0.75 1.58 10] s (same
%                total c), E-output depression as the chained vesicle pool,
%                tau_rec [2 5] s, tau_rel 0.25 s, n2 2
%   product MTS  the same SFA ladders, E-output depression as the product of
%                two pools, tau_rec [2 4] s, tau_rel [0.25 0.5] s
% (run13's value where run13 has the parameter, the mu7chain bundle's
% chained-MTS default where only the MTS form has it; Tom, 2026-11-12.)
% Nothing is tuned.
%
% One figure opens, one column per condition: x of every neuron grouped by
% module (E before I within a module), the module-mean x, the mean rate of E
% and I, the mean adaptation term (c/K) sum_k a_k of E and I, the mean
% depression factor of the E->E route (the model's own std_route_factor:
% the releasable pool b1 on a chain route, b1 b2 on a product route), and
% the detector's group recruitment with the events (triangles). Events are
% detector v1 (ied_detect_events) on the 100-s window.
%
% A regression check re-runs 'original' at run13's own length (10 + 50 s)
% and compares it with the documented run13. It needs its own run because
% the SRA1 integrator draws its two Wiener streams as randn(n, n_steps) one
% after the other, so the second stream, and with it the trajectory,
% depends on the run's length: the 100-s runs do not repeat run13's first
% 50 s (x differs by up to about 0.014).
%
% Run setup_paths once in the session first (cwd at the repository root);
% like every script in scripts/examples this file has no path code
% (setup_paths also puts scripts/explorations/IED_model on the path). It
% takes about 70 s, serial. Nothing is saved unless save_figs is set.
%
% WHAT TO EXPECT with the settings below (main at eddd8ee or later):
%   original     30 events (18.0/min), mean rate E 0.1248 I 0.1379,
%                SFA term E 0.0624 I 0.0690, depression factor 0.6074
%   chain MTS    12 events (7.2/min),  mean rate E 0.0913 I 0.0952,
%                SFA term E 0.0449 I 0.0471, depression factor 0.5845
%   product MTS   1 event  (0.6/min),  mean rate E 0.0721 I 0.0719,
%                SFA term E 0.0354 I 0.0355, depression factor 0.5571
%   original at 10 + 50 s: 11 events, mean rate 0.128643, x corr 0.0517183,
%                realized indegree 21.371, event times 5.075 ... 48.595 s:
%                run13 exactly (events_run13.csv)
% One seed pair fixes everything, so a rerun gives these numbers exactly.
%
% SAVED when save_figs is true, in out_dir:
%   run13_chain_mts.fig / .png   the figure
%   traces.mat                   one struct per condition, named
%                                traces_original, traces_chain_mts,
%                                traces_product_mts, each with
%     condition      condition name
%     t              1 x nt time, s (200 Hz, 0..100 s)
%     rate           n x nt_r per-neuron rate r, single, 100 Hz (every
%                    second sample of t; downsampled to keep the file small)
%     t_rate         1 x nt_r time of rate, s
%     E_idx, I_idx   neuron indices of each type (rows of rate)
%     groups         n x 1 module label (1..5)
%     rate_E_mean, rate_I_mean   1 x nt population mean rate (200 Hz)
%     module_rate_E  5 x nt mean E rate of each module (200 Hz)
%     sfa_E_mean, sfa_I_mean     1 x nt mean (c/K) sum_k a_k (200 Hz)
%     depression_EE_mean         1 x nt mean depression factor, E->E route
%                                (= E->I; std_route_factor) (200 Hz)
%     events         table, one row per detected event, ied_event_details'
%                    columns as in docs/IED_model/events_run13.csv: Time (s),
%                    DominantGroup, DetectedGroups (count), Recruitment
%                    (fraction of all neurons, the event's size),
%                    RecruitmentWidth (s, detector width), WaveformWidth (s,
%                    FWHM of the dominant module's mean x, the event's
%                    duration; NaN if censored), WaveformAmplitude (a.u.),
%                    CensoredWidth, PreLocalGrowth / PreExpansion (NaN: no
%                    Lyapunov tracking here)
%     detector_metrics  ied_detect_events' metrics struct
%     within_corr, between_corr  mean x correlation within / between modules
%     params         every parameter value used (run13_condition_model's
%                    info.params) plus connectivity
%
% See also: run13_condition_model, ied_base_config, ied_build_network,
%           ied_detect_events, SRNNCellTypePairs.std_route_factor

%% Settings
duration     = 100;                   % s observed after the 10-s settling interval
seed         = 42;                    % network seed (run13's); the noise seed follows it
fig_visible  = true;                  % false: draw without a window (saved .fig still opens visible)
save_figs    = false;                 % true: write the figure and traces.mat to out_dir
out_dir      = '';                    % '' = figs/explore_B/run13_chain_mts under the repo root

%% Simulate the three conditions on one network
conds  = {'original', 'chain_mts', 'product_mts'};
titles = {'original', 'chain MTS', 'product MTS'};
res = cell(1, numel(conds));
for k = 1:numel(conds)
    started = tic;
    [m, groups, info] = run13_condition_model(conds{k}, 'duration', duration, 'seed', seed);
    assert(all(m.u_ex == 0, 'all'), 'External input must be zero, as in run13.');
    m.run();
    res{k} = summarize_run(m, groups, info);
    res{k}.name = titles{k};
    fprintf(['%-12s %3d events (%5.1f/min), median recruitment %.3f, mean rate E %.4f I %.4f, ' ...
        'SFA term E %.4f I %.4f, depression factor %.4f, x corr %.4f, ladders clipped %d; %.1f s\n'], ...
        titles{k}, res{k}.det.metrics.event_count, res{k}.det.metrics.events_per_min, ...
        res{k}.det.metrics.median_recruitment, res{k}.rE, res{k}.rI, mean(res{k}.sfa, 2), ...
        mean(res{k}.dep), res{k}.det.metrics.mean_corr, sum(info.params.tau_a_clipped), toc(started));
end

%% Regression check: 'original' re-run at run13's own length, 10 + 50 s
% run13 is documented in docs/IED_model/progress_2026_09_30.md and
% events_run13.csv: 11 events, mean rate 0.128643, mean x correlation
% 0.0517183, realized indegree 21.371. (See the header for why this needs
% its own 50-s run.)
[m50, g50, info50] = run13_condition_model('original', 'duration', 50, 'seed', 42);
m50.run();
r13 = summarize_run(m50, g50, info50);
fprintf('original at 10 + 50 s: %d events (run13: 11), mean rate %.6f (0.128643), corr %.7f (0.0517183), indegree %.3f (21.371)\n', ...
    r13.det.metrics.event_count, r13.det.metrics.mean_rate, r13.det.metrics.mean_corr, ...
    info50.connectivity.realized_indegree);
fprintf('  event times: %s\n', mat2str([r13.det.events.time], 5));

%% Figure, and the saved traces
fig = draw_figure(res, fig_visible);
if save_figs
    if isempty(out_dir) %#ok<UNRCH> % save_figs is a setting
        out_dir = fullfile(fileparts(which('setup_paths')), 'figs', 'explore_B', 'run13_chain_mts');
    end
    if ~isfolder(out_dir), mkdir(out_dir); end
    file = fullfile(out_dir, 'run13_chain_mts');
    fig.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';   % the .fig opens visible even if drawn hidden
    savefig(fig, [file '.fig']);
    exportgraphics(fig, [file '.png'], 'Resolution', 100);
    tr = struct();
    for k = 1:numel(conds)
        tr.(['traces_' conds{k}]) = trace_struct(res{k}, conds{k});
    end
    save(fullfile(out_dir, 'traces.mat'), '-struct', 'tr', '-v7');
    fprintf('Saved %s.fig, .png and traces.mat\n', file);
end

%% Local functions
function s = summarize_run(m, groups, info)
p = m.plot_data; params = m.get_params();
s.t = p.t(:)'; s.x = [p.x.E; p.x.I]; s.r = [p.r.E; p.r.I];
s.groups = groups; s.info = info; s.E = m.type_indices{1}; s.I = m.type_indices{2};
s.rE = mean(p.r.E, 'all'); s.rI = mean(p.r.I, 'all');
s.rate = [mean(p.r.E, 1); mean(p.r.I, 1)];
isE = false(m.n, 1); isE(s.E) = true;
G = max(groups); s.module_rate_E = zeros(G, numel(s.t));
for g = 1:G, s.module_rate_E(g, :) = mean(s.r(isE & groups == g, :), 1); end
s.sfa = zeros(2, numel(s.t));
for q = 1:2
    a = p.a.(m.cell_type_names{q});                 % n_q x K x nt
    K = size(a, 2);
    s.sfa(q, :) = info.params.c(q) / K * reshape(sum(a, [1 2]), 1, []) / size(a, 1);
end
D = SRNNCellTypePairs.std_route_factor(p.b.E.E, params, 1, 1);   % n_E x 1 x nt
s.dep = reshape(mean(D, 1), 1, []);
s.det = ied_detect_events(s.t, s.x, s.r, groups, struct());
[s.events, s.spatial] = ied_event_details(struct('t', s.t, 'x', s.x, 'groups', groups), s.det);
end

function tr = trace_struct(s, name) %#ok<DEFNU> % called when save_figs is set
keep = 1:2:numel(s.t);
events = s.events;
params = s.info.params; params.connectivity = s.info.connectivity;
tr = struct('condition', name, 't', s.t, 'rate', single(s.r(:, keep)), 't_rate', s.t(keep), ...
    'E_idx', s.E(:)', 'I_idx', s.I(:)', 'groups', s.groups, ...
    'rate_E_mean', s.rate(1, :), 'rate_I_mean', s.rate(2, :), 'module_rate_E', s.module_rate_E, ...
    'sfa_E_mean', s.sfa(1, :), 'sfa_I_mean', s.sfa(2, :), 'depression_EE_mean', s.dep, ...
    'events', events, 'detector_metrics', s.det.metrics, ...
    'within_corr', s.spatial.within_corr, 'between_corr', s.spatial.between_corr, 'params', params);
end

function fig = draw_figure(res, visible)
nc = numel(res); G = max(res{1}.groups); pal = lines(G);
cE = [.85 .25 .1]; cI = [.1 .4 .7];
fig = figure('Color', 'w', 'Position', [30 30 1700 1100], 'Visible', visible, ...
    'Name', 'run13 network: original, chain MTS, product MTS', 'NumberTitle', 'off');
tl = tiledlayout(fig, 6, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(6, nc);
for k = 1:nc
    s = res{k}; t = s.t; ev = [s.det.events.time];
    [~, order] = sortrows([s.groups (1:numel(s.groups))'], [1 2]);
    ax(1, k) = nexttile(tl, k);
    img = 1:2:numel(t);    % 100 Hz and single, which keeps the .fig near 50 MB
    imagesc(ax(1, k), t(img), 1:numel(order), single(s.x(order, img))); axis(ax(1, k), 'xy');
    clim(ax(1, k), [-.3 .8]); hold(ax(1, k), 'on');
    edges = cumsum(accumarray(s.groups, 1)) + .5;
    yline(ax(1, k), edges(1:end-1), 'w-', 'LineWidth', 0.8);
    title(ax(1, k), sprintf('%s: %d events (%.1f/min), r_E %.3f, r_I %.3f', ...
        s.name, numel(ev), s.det.metrics.events_per_min, s.rE, s.rI));
    if k == 1, ylabel(ax(1, k), {'x, neuron', 'by module'}); end
    if k == nc, colorbar(ax(1, k)); end

    ax(2, k) = nexttile(tl, nc + k); hold(ax(2, k), 'on');
    for g = 1:G, plot(ax(2, k), t, mean(s.x(s.groups == g, :), 1), 'Color', pal(g, :)); end
    if k == 1, ylabel(ax(2, k), 'module mean x'); end

    ax(3, k) = nexttile(tl, 2*nc + k); hold(ax(3, k), 'on');
    plot(ax(3, k), t, s.rate(1, :), 'Color', cE); plot(ax(3, k), t, s.rate(2, :), 'Color', cI);
    if k == 1, ylabel(ax(3, k), 'mean rate'); legend(ax(3, k), {'E', 'I'}, 'Location', 'northwest'); end

    ax(4, k) = nexttile(tl, 3*nc + k); hold(ax(4, k), 'on');
    plot(ax(4, k), t, s.sfa(1, :), 'Color', cE); plot(ax(4, k), t, s.sfa(2, :), 'Color', cI);
    if k == 1, ylabel(ax(4, k), {'mean (c/K) \Sigma_k a_k'}); legend(ax(4, k), {'E', 'I'}, 'Location', 'northwest'); end

    ax(5, k) = nexttile(tl, 4*nc + k);
    plot(ax(5, k), t, s.dep, 'Color', [0 .55 .4]);
    if k == 1, ylabel(ax(5, k), {'mean depression', 'factor, E outputs'}); end

    ax(6, k) = nexttile(tl, 5*nc + k); hold(ax(6, k), 'on');
    for g = 1:G, plot(ax(6, k), t, s.det.group_recruitment(g, :), 'Color', pal(g, :)); end
    for e = 1:numel(ev)
        g = s.events.DominantGroup(e);
        plot(ax(6, k), ev(e), 0.95, 'v', 'MarkerFaceColor', pal(g, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
    end
    ylim(ax(6, k), [0 1]); xlabel(ax(6, k), 'time (s)');
    if k == 1, ylabel(ax(6, k), {'group recruitment', '(\nabla = event)'}); end
    for row = 1:2
        if ~isempty(ev), xline(ax(row, k), ev, ':', 'Color', [.75 .3 .75], 'HandleVisibility', 'off'); end
    end
end
for row = 2:5, linkaxes(ax(row, :), 'y'); end
linkaxes(ax(:), 'x'); xlim(ax(1, 1), [res{1}.t(1) res{1}.t(end)]);
set(ax, 'FontSize', 9, 'Box', 'off');
end
