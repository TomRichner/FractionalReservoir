function out = p3_ablation_pilot(opts)
% P3_ABLATION_PILOT The work behind run_p3_ablation.m: run every ablation, measure, draw, save.
%
%   out = P3_ABLATION_PILOT(name, value, ...)
%
% Read run_p3_ablation.m first: it states the question, the ablations, the
% measures, the regime rule and what is saved. Kept as a function so the
% script's settings reach it as arguments and it can run with the figures
% hidden (a MATLAB with no desktop).
%
% out.T          the results table, one row per run
% out.runs       cell of p3_run_one results (full traces)
% out.figures    struct array: name, h
% out.files      files written (empty unless save_figs)
%
% See also: run_p3_ablation, p3_ablation_args, p3_run_one

arguments
    opts.networks    (1,:) cell    = {'bursting', 'run13'}
    opts.T_show      (1,1) double  = 100
    opts.seed        (1,1) double  = 42
    opts.fig_visible (1,1) logical = true
    opts.save_figs   (1,1) logical = false
    opts.out_dir     (1,:) char    = ''
    opts.scan_fractions (1,:) double = [0.2 0.4 0.6 0.8]
end
networks = opts.networks; T_show = opts.T_show; seed = opts.seed;
fig_visible = opts.fig_visible; out_dir = opts.out_dir;

chain_ablations    = {'sfa_removed', 'std_removed', 'both_removed', 'mean_matched', ...
                      'mean_matched_rw', 'mean_matched_pn', 'slow_removed', 'slow_removed_share'};
original_ablations = {'sfa_removed', 'std_removed', 'mean_matched'};

%% The sign and exactness of the constant input that stands in for SFA
offset_err = p3_check_offset();
fprintf('Constant input -abar against a frozen adaptation term: max difference %.2g (wrong sign: %.2g)\n', offset_err);

%% Run
started = tic;
all_runs = {};
for i = 1:numel(networks)
    net = networks{i};
    for cond = {'chain_mts', 'original'}
        if strcmp(cond{1}, 'chain_mts'), abl = chain_ablations; else, abl = original_ablations; end
        if ~strcmp(net, 'run13'), abl = setdiff(abl, {'mean_matched_pn'}, 'stable'); end
        I = p3_run_one(net, cond{1}, 'intact', [], 'T_show', T_show, 'seed', seed);
        I.meas.regime = 'comparable';
        all_runs{end+1} = I; %#ok<AGROW>
        report_line(I);
        for k = 1:numel(abl)
            R = p3_run_one(net, cond{1}, abl{k}, I.ref, 'T_show', T_show, 'seed', seed);
            R.meas.regime = regime_word(R.meas, I.meas);
            all_runs{end+1} = R; %#ok<AGROW>
            report_line(R);
        end
        % The mean-matched constant between its two definitions: route
        % scale = rate-weighted average + f * (plain - rate-weighted).
        if strcmp(cond{1}, 'chain_mts')
            for f = opts.scan_fractions
                ref = I.ref;
                ref.D_route = I.ref.D_route_rw + f * (I.ref.D_route - I.ref.D_route_rw);
                R = p3_run_one(net, cond{1}, 'mean_matched', ref, 'T_show', T_show, 'seed', seed);
                R.ablation = sprintf('mean_matched_scan_%.2f', f);
                R.spec.label = sprintf('mean-matched, scale %.3f', ref.D_route(1, 1));
                R.spec.how = sprintf(['As mean_matched, with the route scale %.0f %% of the way from the ' ...
                    'rate-weighted to the plain average.'], 100 * f);
                R.meas.regime = regime_word(R.meas, I.meas);
                all_runs{end+1} = R; %#ok<AGROW>
                report_line(R);
            end
        end
    end
end
fprintf('Total %.0f s\n', toc(started));

%% Table
T = results_table(all_runs);
disp(T(:, {'network', 'condition', 'ablation', 'count', 'per_min', 'size_median', 'duration_median', ...
    'mean_rE', 'mean_rI', 'frac_time_E_above_0p5', 'frac_pinned', 'frac_silent', 'regime'}));

%% Figures
figs = struct('name', {}, 'h', {});
pick = @(net, cond, abl) all_runs{cellfun(@(R) strcmp(R.network, net) && strcmp(R.condition, cond) && strcmp(R.ablation, abl), all_runs)};
for i = 1:numel(networks)
    net = networks{i};
    main = {'intact', 'sfa_removed', 'std_removed', 'both_removed', 'mean_matched', 'slow_removed'};
    runs = cellfun(@(a) pick(net, 'chain_mts', a), main, 'UniformOutput', false);
    titles = cellfun(@(R) R.spec.label, runs, 'UniformOutput', false);
    titles{1} = 'chained MTS, intact';
    figs(end+1) = struct('name', [net '_ablations'], 'h', p3_figure(runs, titles, ...
        sprintf('%s network, chained MTS, no stimulation: intact and with adaptation removed (seed %d)', net, seed), fig_visible)); %#ok<AGROW>

    var = {'mean_matched_rw', 'mean_matched_pn', 'slow_removed_share'};
    if ~strcmp(net, 'run13'), var = setdiff(var, {'mean_matched_pn'}, 'stable'); end
    runs = [cellfun(@(a) pick(net, 'chain_mts', a), var, 'UniformOutput', false), ...
            cellfun(@(a) pick(net, 'original', a), [{'intact'}, original_ablations], 'UniformOutput', false)];
    titles = cellfun(@(R) sprintf('%s: %s', strrep(R.condition, '_mts', ' MTS'), R.spec.label), runs, 'UniformOutput', false);
    figs(end+1) = struct('name', [net '_variants'], 'h', p3_figure(runs, titles, ...
        sprintf('%s network, no stimulation: variants on chained MTS, and the original condition (seed %d)', net, seed), fig_visible)); %#ok<AGROW>
end
figs(end+1) = struct('name', 'summary', 'h', p3_summary_figure(T, fig_visible));

%% Save
files = {};
if opts.save_figs
    if isempty(out_dir)
        out_dir = fullfile(fileparts(which('setup_paths')), 'figs', 'explore_B', 'p3_ablation');
    end
    if ~isfolder(out_dir), mkdir(out_dir); end
    for k = 1:numel(figs)
        file = fullfile(out_dir, figs(k).name);
        exportgraphics(figs(k).h, [file '.png'], 'Resolution', 96);
        figs(k).h.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';   % the .fig opens visible even if drawn hidden
        savefig(figs(k).h, [file '.fig']);
        files = [files, {[file '.png'], [file '.fig']}]; %#ok<AGROW>
    end
    writetable(T, fullfile(out_dir, 'results.csv'));
    runs = cellfun(@light_run, all_runs);
    settings = struct('networks', {networks}, 'T_show', T_show, 'T_settle', 10, 'seed', seed, ...
        'offset_check', offset_err);
    save(fullfile(out_dir, 'results.mat'), 'T', 'runs', 'settings', '-v7');
    files = [files, {fullfile(out_dir, 'results.csv'), fullfile(out_dir, 'results.mat')}];
    fprintf('Saved %d figures (.png, .fig), results.csv and results.mat in %s\n', numel(figs), out_dir);
end
out = struct('T', T, 'runs', {all_runs}, 'figures', figs, 'files', {files});
end

%% Local functions
function word = regime_word(M, I)
if M.frac_pinned >= 0.5 || M.frac_time_E_above_0p5 >= 0.5
    word = 'saturation';
elseif M.frac_time_burst >= 0.5
    word = 'continuous high firing';
elseif M.count == 0 && M.mean_rE < 0.5 * I.mean_rE
    word = 'quiet';
elseif M.count == 0
    word = 'steady, no events';
elseif M.count > 2 * I.count
    word = 'more';
elseif M.count < 0.5 * I.count
    word = 'fewer';
else
    word = 'comparable';
end
end

function report_line(R)
M = R.meas;
fprintf('%-8s %-9s %-20s %3d (%5.1f/min) size %.2f dur %.2f s | r_E %.3f r_I %.3f | E>0.5 %.2f pinned %.2f silent %.2f | %s; %.0f s\n', ...
    R.network, R.condition, R.ablation, M.count, M.per_min, M.size_median, M.duration_median, ...
    M.mean_rE, M.mean_rI, M.frac_time_E_above_0p5, M.frac_pinned, M.frac_silent, M.regime, R.run_seconds);
end

function T = results_table(all_runs)
n = numel(all_runs);
rows = cell(n, 1);
for k = 1:n
    R = all_runs{k}; M = R.meas;
    D = R.ref.D_route(1, 1); Drw = R.ref.D_route_rw(1, 1);
    used = NaN;      % the route scale an ablation was given on E->E
    sc = [];
    idx = find(strcmp(R.spec.model_args(1:2:end), 'synapse_config'));
    if ~isempty(idx), sc = R.spec.model_args{2 * idx}; end
    if ~isempty(sc) && isfield(sc.E.E, 'scale'), used = sc.E.E.scale; end
    drive_E = NaN;
    if ~isempty(R.spec.drive), drive_E = mean(R.spec.drive(R.E_idx)); end
    rows{k} = table({R.network}, {R.condition}, {R.ablation}, {R.spec.label}, ...
        M.count, M.per_min, M.size_median, M.duration_median, M.n_bursts, M.frac_time_burst, ...
        M.mean_rE, M.mean_rI, M.peak_mean_rE, M.frac_time_E_above_0p5, M.frac_pinned, M.frac_silent, ...
        {M.regime}, mean(R.sfa_E), mean(R.sfa_I), D, Drw, used, drive_E, {R.spec.how}, ...
        'VariableNames', {'network', 'condition', 'ablation', 'label', ...
        'count', 'per_min', 'size_median', 'duration_median', 'n_bursts', 'frac_time_burst', ...
        'mean_rE', 'mean_rI', 'peak_mean_rE', 'frac_time_E_above_0p5', 'frac_pinned', 'frac_silent', ...
        'regime', 'sfa_term_E', 'sfa_term_I', 'dep_factor_EE', 'dep_factor_EE_rate_weighted', ...
        'scale_given_EE', 'mean_input_given_E', 'how'});
end
T = vertcat(rows{:});
end

function L = light_run(R)
% What is kept of a run in results.mat.
L = struct('network', R.network, 'condition', R.condition, 'ablation', R.ablation, ...
    'label', R.spec.label, 'how', R.spec.how, 't', R.t, ...
    'mean_rE', mean(R.rE, 1), 'mean_rI', mean(R.rI, 1), ...
    'time_avg_rate', [mean(R.rE, 2); mean(R.rI, 2)], ...
    'sfa_E', R.sfa_E, 'sfa_I', R.sfa_I, 'dep', R.dep, ...
    'meas', R.meas, 'model_args', {R.spec.model_args}, 'drive', R.spec.drive, ...
    'col_scale', R.spec.col_scale, 'ref', R.ref);
end
