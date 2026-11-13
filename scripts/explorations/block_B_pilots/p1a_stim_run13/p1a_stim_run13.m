% P1A_STIM_RUN13 Block B pilot 1a: does stimulation suppress the IED-like events of run13 under chained MTS adaptation?
%
% The stimulation experiment of 2026-10-05 (docs/IED_model/stim_2026_10_05.md,
% ied_stim_experiment) showed that a sustained uniform DC input cuts the
% candidate IED-like events of the run13 network (single-timescale adaptation)
% from 18.6 to 2.4 per minute at amplitude 0.20. Pilot 0a found that run13
% with chained multiple-timescale (MTS) adaptation still has localized events
% (7.2 per min against 18). This pilot asks whether the same stimulus still
% suppresses them, whether adaptation is recruited as before, and whether the
% relative reduction is similar.
%
% DESIGN. Conditions 'original' and 'chain_mts' of run13_condition_model
% (scripts/examples/run13_chain_mts_example, the one definition; product MTS
% is left out because it has about one event per 100 s unstimulated).
% Eight network seeds, ied_stim_experiment's (42 73 101 137 211 307 409 503;
% noise seed 224779 + (seed - 42)). For each condition and seed: an OFF arm
% and one ON arm per amplitude, each 10 s settling + 150 s observed, all of
% the SAME length, so they share the network, setpoints, SFA ladders, x0 and
% the Wiener path (the SRA1 noise depends on the run's length). The documented
% experiment's windows are 50 s, under the 100 s the brief requires, so the
% paired-arm design replaces its within-run halves. Stimulus: the documented
% one, dc_staircase_stimulus, uniform over all 350 neurons, no ramp, no input
% noise, from t = 0 to the end in the ON arm (zero during settling); main
% amplitude 0.20 (the documented figure's), plus 0.10 and 0.30. Measures on
% 0..150 s with the IED pilot's fixed detector v1 (its own adaptive threshold
% per arm; the fixed-criterion control and the window 10..150 s, which drops
% the stimulus onset, are secondary). Nothing is tuned.
%
% REGRESSION CHECK. The documented protocol itself (one 10 + 120 s run, input
% 0 then 0.20 at 60 s, detector on 10..60 s and 70..120 s) is re-run for
% 'original' on the eight seeds and compared with
% docs/IED_model/stim_2026_10_05_summary.csv, row by row.
%
% Run setup_paths once first (cwd at the repo root). With a pool (pool_set(4))
% the 64 simulations take about 8 min; serially about 25 min. Writes into
% out_dir: results.mat, per_seed.csv, regression.csv, stats.txt,
% p1a_example_pairs.png/.fig, p1a_summary.png/.fig, example_traces.mat.
%
% WHAT TO EXPECT (main at a33fa6f, pool of 4, about 8.5 min):
%   regression: largest relative difference from the documented csv 4.7e-15
%   (every seed, both windows: events, fixed criterion, recruitment, rate,
%   SFA feedback, resources); x(0) differs between arms by at most 0.0038
%   (the model's linear input interpolant ramps the step over the last
%   settling step)
%   events per minute, mean over 8 seeds, OFF -> ON at 0.10 / 0.20 / 0.30:
%     original   20.3 -> 9.3 / 3.55 / 1.55   (ON lower in 8 of 8 at each)
%     chain MTS  15.3 -> 5.7 / 2.55 / 1.00   (ON lower in 8 of 8 at each)
%   at 0.20: SFA term E 0.0622 -> 0.0941 (original), 0.0431 -> 0.0693
%   (chain); depression factor 0.616 -> 0.527, 0.603 -> 0.485
% One seed list fixes everything, so a rerun gives these numbers exactly
% (the task functions pin the generator type, so pool and client agree).
%
% ONE THING TO KNOW: the model seeds the Wiener path and the setpoints with
% rng(seed), which keeps the session's generator type, 'threefry' on a pool
% worker; p1a_run_pair and p1a_run_documented pin 'twister' (the client's),
% without which a pool run builds different networks from the same seeds.

%% Settings
conditions   = {'original', 'chain_mts'};
seeds        = [42 73 101 137 211 307 409 503];   % ied_stim_experiment's
amplitudes   = [0.10 0.20 0.30];                  % ON-arm uniform DC levels
main_amp     = 0.20;                              % the documented figure's amplitude
duration     = 150;                               % s observed per arm, after 10 s settling
late_start   = 10;                                % secondary window [late_start duration]
example_seed = 42;                                % run13's own seed
run_regression = true;
repo_root    = fileparts(which('setup_paths'));
% Run from a git worktree (.claude/worktrees/<name>), write into the main
% checkout's figure folder, where the other block B pilots are.
wt = strfind(repo_root, [filesep '.claude' filesep 'worktrees' filesep]);
if ~isempty(wt), repo_root = repo_root(1:wt(1) - 1); end
out_dir      = fullfile(repo_root, 'figs', 'explore_B', 'p1a_stim_run13');
if ~isfolder(out_dir), mkdir(out_dir); end
started_all = tic;

%% Simulate: one task per condition and seed (OFF + every ON amplitude)
pool_ensure();
[ci, si] = ndgrid(1:numel(conditions), 1:numel(seeds));
nT = numel(ci); task_out = cell(1, nT);
parfor k = 1:nT
    task_out{k} = p1a_run_pair(conditions{ci(k)}, seeds(si(k)), amplitudes, duration, ...
        late_start, seeds(si(k)) == example_seed);
end
rows = table(); x0_diff = []; traces = struct([]);
for k = 1:nT
    rows = [rows; task_out{k}.rows]; %#ok<AGROW>
    x0_diff = [x0_diff task_out{k}.x0_diff]; %#ok<AGROW>
    traces = [traces task_out{k}.traces]; %#ok<AGROW>
end
fprintf('Simulated %d arms in %.0f s. Max |x(0)| difference between arms: %.2g\n', ...
    height(rows), toc(started_all), max(x0_diff));

%% Regression: the documented protocol, 'original', every seed at 0.20
reg = table();
if run_regression
    reg_out = cell(1, numel(seeds));
    parfor k = 1:numel(seeds)
        reg_out{k} = p1a_run_documented(seeds(k), main_amp);
    end
    for k = 1:numel(seeds), reg = [reg; reg_out{k}]; end %#ok<AGROW>
    doc = readtable(fullfile(fileparts(which('setup_paths')), 'docs', 'IED_model', 'stim_2026_10_05_summary.csv'));
    doc = doc(abs(doc.amplitude - main_amp) < 1e-9, :);
    cols = {'event_count', 'events_per_min', 'events_per_min_fixed', 'median_recruitment', ...
        'mean_rate', 'sfa_feedback', 'resources'};
    worst = 0;
    for k = 1:height(reg)
        d = doc(doc.seed == reg.seed(k) & strcmp(doc.segment, reg.segment(k)), :);
        assert(height(d) == 1);
        for c = 1:numel(cols)
            a = reg.(cols{c})(k); b = d.(cols{c});
            if ~(isnan(a) && isnan(b)), worst = max(worst, abs(a - b) / max(abs(b), 1e-12)); end
            reg.([cols{c} '_doc'])(k) = b;
        end
    end
    fprintf('Regression vs documented csv (amplitude %.2f, %d seeds): largest relative difference %.2g\n', ...
        main_amp, numel(seeds), worst);
end

%% Statistics per condition and amplitude
stats = table();
for c = 1:numel(conditions)
    for A = amplitudes
        off = rows(rows.condition == conditions{c} & rows.amplitude == 0, :);
        on  = rows(rows.condition == conditions{c} & abs(rows.amplitude - A) < 1e-9, :);
        assert(isequal(off.seed, on.seed));
        d = off.events_per_min - on.events_per_min;
        df = off.events_per_min - on.events_per_min_fixed;
        dl = off.events_per_min_late - on.events_per_min_late;
        ratio = on.events_per_min ./ off.events_per_min;
        st = table(string(conditions{c}), A, numel(d), ...
            mean(off.events_per_min), mean(on.events_per_min), mean(d), nnz(d > 0), nnz(d < 0), signflip_p(d), ...
            1 - mean(on.events_per_min) / mean(off.events_per_min), median(ratio, 'omitnan'), ...
            1 - mean(on.events_per_min_fixed) / mean(off.events_per_min), nnz(df > 0), signflip_p(df), ...
            1 - mean(on.events_per_min_late) / mean(off.events_per_min_late), nnz(dl > 0), signflip_p(dl), ...
            mean(off.rate_E), mean(on.rate_E), mean(off.rate_I), mean(on.rate_I), ...
            mean(off.sfa_E), mean(on.sfa_E), mean(off.sfa_I), mean(on.sfa_I), ...
            mean(off.dep), mean(on.dep), nnz(on.sfa_E > off.sfa_E), nnz(on.dep < off.dep), ...
            median(off.median_recruitment, 'omitnan'), median(on.median_recruitment, 'omitnan'), ...
            median(off.median_amplitude, 'omitnan'), median(on.median_amplitude, 'omitnan'), ...
            'VariableNames', {'condition', 'amplitude', 'n_seeds', ...
            'epm_off', 'epm_on', 'epm_diff', 'n_on_lower', 'n_on_higher', 'p_signflip', ...
            'rel_reduction', 'median_ratio', ...
            'rel_reduction_fixed', 'n_on_lower_fixed', 'p_signflip_fixed', ...
            'rel_reduction_late', 'n_on_lower_late', 'p_signflip_late', ...
            'rate_E_off', 'rate_E_on', 'rate_I_off', 'rate_I_on', ...
            'sfa_E_off', 'sfa_E_on', 'sfa_I_off', 'sfa_I_on', 'dep_off', 'dep_on', ...
            'n_sfa_up', 'n_dep_down', 'recruit_off', 'recruit_on', 'amp_off', 'amp_on'});
        stats = [stats; st]; %#ok<AGROW>
    end
end
txt = regexprep(evalc('disp(stats)'), '</?strong>', '');
fprintf('%s', txt);

%% Save
[~, commit] = system(sprintf('git -C "%s" rev-parse --short HEAD', fileparts(which('p1a_stim_run13'))));
settings = struct('conditions', {conditions}, 'seeds', seeds, 'amplitudes', amplitudes, ...
    'main_amp', main_amp, 'duration', duration, 'settling', 10, 'late_start', late_start, ...
    'example_seed', example_seed, 'stimulus', 'dc_staircase_stimulus, uniform, no ramp, no input noise', ...
    'commit', strtrim(commit));
run_seconds = toc(started_all);
save(fullfile(out_dir, 'results.mat'), 'rows', 'stats', 'reg', 'x0_diff', 'settings', 'run_seconds', '-v7');
writetable(rows, fullfile(out_dir, 'per_seed.csv'));
writetable(stats, fullfile(out_dir, 'stats.csv'));
if ~isempty(reg), writetable(reg, fullfile(out_dir, 'regression.csv')); end
ex = traces([traces.amplitude] == 0 | abs([traces.amplitude] - main_amp) < 1e-9);
save(fullfile(out_dir, 'example_traces.mat'), 'ex', '-v7');
fid = fopen(fullfile(out_dir, 'stats.txt'), 'w'); fprintf(fid, '%s', txt); fclose(fid);

%% Figures
fig1 = draw_example_pairs(ex, conditions, example_seed, main_amp);
save_fig(fig1, fullfile(out_dir, 'p1a_example_pairs'));
fig2 = draw_summary(rows, conditions, amplitudes, main_amp);
save_fig(fig2, fullfile(out_dir, 'p1a_summary'));
fprintf('Done in %.0f s; outputs in %s\n', toc(started_all), out_dir);

%% Local functions
function p = signflip_p(d)
% Exact two-sided paired sign-flip test of mean(d) = 0 over all 2^n patterns.
d = d(:); n = numel(d);
S = 1 - 2 * (dec2bin(0:2^n - 1, n) - '0');      % 2^n x n of +-1
null = S * d / n;
p = mean(abs(null) >= abs(mean(d)) - 1e-12);
end

function save_fig(fig, file)
fig.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';
savefig(fig, [file '.fig']);
exportgraphics(fig, [file '.png'], 'Resolution', 96);
close(fig);
end

function fig = draw_example_pairs(ex, conditions, seed, A)
titles = containers.Map({'original', 'chain_mts'}, {'original', 'chain MTS'});
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [20 20 1750 1050], 'Visible', 'off');
tl = tiledlayout(fig, 6, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Example pairs, network seed %d: stimulation OFF and ON (uniform input %.2f from t = 0)', seed, A));
cE = [.85 .25 .1]; cI = [.1 .4 .7]; pal = lines(5);
ax = gobjects(6, 4); col = 0;
for c = 1:numel(conditions)
    for on = [false true]
        col = col + 1;
        sel = strcmp({ex.condition}, conditions{c}) & [ex.seed] == seed;
        amps = [ex.amplitude];
        tr = ex(find(sel & (on & abs(amps - A) < 1e-9 | ~on & amps == 0), 1));
        t = tr.t; ev = tr.events.Time; dg = tr.events.DominantGroup;
        arm = 'OFF'; if on, arm = 'ON'; end
        ax(1, col) = nexttile(tl, col);
        plot(ax(1, col), t, tr.u, 'k', 'LineWidth', 1.2); ylim(ax(1, col), [-0.02 0.35]);
        title(ax(1, col), sprintf('%s, stimulation %s: %d events (%.1f/min)', titles(conditions{c}), arm, ...
            numel(ev), numel(ev) * 60 / (t(end) - t(1))));
        ax(2, col) = nexttile(tl, 4 + col); hold(ax(2, col), 'on');
        for g = 1:5, plot(ax(2, col), t, tr.module_x(g, :), 'Color', pal(g, :)); end
        ax(3, col) = nexttile(tl, 8 + col); hold(ax(3, col), 'on');
        plot(ax(3, col), t, tr.rate_E_mean, 'Color', cE); plot(ax(3, col), t, tr.rate_I_mean, 'Color', cI);
        ax(4, col) = nexttile(tl, 12 + col); hold(ax(4, col), 'on');
        plot(ax(4, col), t, tr.sfa_E_mean, 'Color', cE); plot(ax(4, col), t, tr.sfa_I_mean, 'Color', cI);
        ax(5, col) = nexttile(tl, 16 + col);
        plot(ax(5, col), t, tr.depression_EE_mean, 'Color', [0 .55 .4]);
        ax(6, col) = nexttile(tl, 20 + col); hold(ax(6, col), 'on');
        for g = 1:5, plot(ax(6, col), t, tr.group_recruitment(g, :), 'Color', pal(g, :)); end
        for e = 1:numel(ev)
            plot(ax(6, col), ev(e), 0.95, 'v', 'MarkerFaceColor', pal(dg(e), :), 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
        end
        ylim(ax(6, col), [0 1]); xlabel(ax(6, col), 'time (s)');
        if ~isempty(ev), xline(ax(2, col), ev, ':', 'Color', [.75 .3 .75]); end
        if col == 1
            ylabel(ax(1, col), 'input u'); ylabel(ax(2, col), 'module mean x');
            ylabel(ax(3, col), 'mean rate'); legend(ax(3, col), {'E', 'I'}, 'Location', 'northwest');
            ylabel(ax(4, col), {'mean SFA term', '(c/K) \Sigma_k a_k'}); legend(ax(4, col), {'E', 'I'}, 'Location', 'northwest');
            ylabel(ax(5, col), {'mean depression', 'factor, E outputs'});
            ylabel(ax(6, col), {'module recruitment', '(\nabla = event)'});
        end
    end
end
for row = 2:6, linkaxes(ax(row, :), 'y'); end
linkaxes(ax(:), 'x'); xlim(ax(1, 1), [t(1) t(end)]);
set(ax, 'FontSize', 9, 'Box', 'off');
end

function fig = draw_summary(rows, conditions, amplitudes, A)
titles = {'original', 'chain MTS'}; cc = [.85 .45 .1; .0 .45 .7];
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [20 20 1500 880], 'Visible', 'off');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Pilot 1a: OFF against ON (uniform input %.2f), one line per network seed (8 seeds), 150 s per arm', A));
% Paired event rates, one panel per condition
ymax = max(rows.events_per_min) * 1.08;
for c = 1:numel(conditions)
    ax = nexttile(tl, c); hold(ax, 'on');
    off = rows(rows.condition == conditions{c} & rows.amplitude == 0, :);
    on  = rows(rows.condition == conditions{c} & abs(rows.amplitude - A) < 1e-9, :);
    for k = 1:height(off), plot(ax, [1 2], [off.events_per_min(k) on.events_per_min(k)], '-o', ...
            'Color', [.55 .55 .55], 'MarkerFaceColor', cc(c, :), 'MarkerEdgeColor', 'none'); end
    plot(ax, [1 2], [mean(off.events_per_min) mean(on.events_per_min)], '-', 'Color', cc(c, :), 'LineWidth', 3);
    xlim(ax, [.6 2.4]); xticks(ax, [1 2]); xticklabels(ax, {'OFF', 'ON'}); ylim(ax, [0 ymax]);
    ylabel(ax, 'events per minute');
    title(ax, sprintf('%s: mean %.1f to %.1f /min, ON lower in %d of %d', titles{c}, ...
        mean(off.events_per_min), mean(on.events_per_min), nnz(on.events_per_min < off.events_per_min), height(off)));
end
% Dose: events per minute against amplitude
ax = nexttile(tl, 3); hold(ax, 'on'); lv = [0 amplitudes]; h = gobjects(1, 2);
for c = 1:numel(conditions)
    M = zeros(numel(unique(rows.seed)), numel(lv));
    for j = 1:numel(lv)
        r = rows(rows.condition == conditions{c} & abs(rows.amplitude - lv(j)) < 1e-9, :);
        M(:, j) = r.events_per_min;
    end
    plot(ax, lv, M', '-', 'Color', [cc(c, :) .25]);
    h(c) = plot(ax, lv, mean(M, 1), '-o', 'Color', cc(c, :), 'LineWidth', 2.5, 'MarkerFaceColor', cc(c, :));
end
legend(ax, h, titles, 'Location', 'northeast'); xlabel(ax, 'uniform input amplitude'); ylabel(ax, 'events per minute');
title(ax, 'events against amplitude (thin: seeds; thick: mean)');
% State variables, OFF against ON, both conditions
q = {'sfa_E', 'mean SFA term, E  (c/K) \Sigma_k a_k'; 'dep', 'mean depression factor, E outputs'; 'rate_E', 'mean E rate'};
for j = 1:3
    ax = nexttile(tl, 3 + j); hold(ax, 'on');
    for c = 1:numel(conditions)
        off = rows(rows.condition == conditions{c} & rows.amplitude == 0, :);
        on  = rows(rows.condition == conditions{c} & abs(rows.amplitude - A) < 1e-9, :);
        x0 = 3 * (c - 1);
        for k = 1:height(off), plot(ax, x0 + [1 2], [off.(q{j, 1})(k) on.(q{j, 1})(k)], '-o', ...
                'Color', [.55 .55 .55], 'MarkerFaceColor', cc(c, :), 'MarkerEdgeColor', 'none'); end
        plot(ax, x0 + [1 2], [mean(off.(q{j, 1})) mean(on.(q{j, 1}))], '-', 'Color', cc(c, :), 'LineWidth', 3);
    end
    xlim(ax, [.5 5.5]); xticks(ax, [1 2 4 5]);
    xticklabels(ax, {'orig. OFF', 'orig. ON', 'chain OFF', 'chain ON'}); title(ax, q{j, 2});
end
set(findall(fig, 'Type', 'axes'), 'FontSize', 10, 'Box', 'off');
end
