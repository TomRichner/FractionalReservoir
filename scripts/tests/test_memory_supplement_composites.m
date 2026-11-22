%% TEST_MEMORY_SUPPLEMENT_COMPOSITES The composites of Plan 2 package F build from a small run, and their two new pieces work.
%
% Main Figure 7 (fig_main7_memory) and supplementary S1, S2, S4, S5, S9
% (fig_supp<k>_grouped), the responsive-share row of fig_mc_gain_bias_maps
% and the steps of composite_tools. No simulation.
%
% THE FIXTURE. The stages of the 100-neuron regression run
% (data/regression_fast/reference_post_refactor: three conditions,
% no_adaptation / sfa1_std1 / sfa3_std2, the tau sweep holding sfa3_std2),
% COPIED to a temporary directory, because fig_param_space_allStd writes a
% prep folder into the run directory it is given; beside them two synthetic
% memory grids written the way the mc_gain_sweep stage writes them (the gain
% by self-coupling grid, from make_synthetic_mc_gain_sweep, and a 1-D sweep
% of the input noise). The run is found as in test_sweep_figures_conditions
% (this checkout, else the main checkout); without it the composite checks
% are SKIPPED and the others still run.
%
% Checks
%   1. fig_mc_gain_bias_maps: the default rows are the ones of before; rows
%      {'responsive'} draws one map per condition on [0, 1] whose cells are
%      the median over seeds of 1 - pinned - silent, computed here
%      independently; a spread row alone with one seed is NoRows; an unknown
%      row name is refused.
%   2. composite_tools: copy keeps the pixel positions and the legend with
%      its axes; place maps a block's bounding box onto the rectangle; grid
%      reads rows and columns; letter writes '(A)'; condition_title gives
%      sentence case, on two lines when asked.
%   3. every composite builds with save false (one figure, its panel
%      letters, out.panels naming a source for each letter, no file) and
%      with save true (the PNG whole, i.e. as tall as its canvas allows, the
%      .fig and <tag>_table.md with a section per panel); no source figure
%      is left open; the run directory is left as it was.
%   4. S4's seventh panel (the gain sweep, drawn in the composite) holds
%      the medians computed here from the saved sweep.
%   5. a run directory without the stages is an error, not a blank figure.
% About two minutes. Figures invisible.

n_fail = 0;
T = composite_tools();
three = {'no_adaptation', 'sfa1_std1', 'sfa3_std2'};
figs_before = findall(groot, 'Type', 'figure');

%% 1. the responsive-share row
[Ts, settings] = make_synthetic_mc_gain_sweep(three);
D = read_mc_gain_sweep(struct('T', Ts, 'settings', settings));
quiet = {'save', false, 'visible', false};
o = fig_mc_gain_bias_maps('data', D, quiet{:});
n_fail = n_fail + check(isequal(o.rows, {'mc', 'lle', 'iqr'}) && numel(findall(o.figs, 'Type', 'axes')) == 9, ...
    'maps: the default rows are memory capacity, lambda_1 and the spread');
close(o.figs);
o = fig_mc_gain_bias_maps('data', D, quiet{:}, 'rows', {'responsive'});
ax = findall(o.figs, 'Type', 'axes');
ok = isequal(o.rows, {'responsive'}) && numel(ax) == 3 && all(arrayfun(@(a) isequal(a.CLim, [0 1]), ax));
n_fail = n_fail + check(ok, 'maps: rows {''responsive''} draws one map per condition on [0, 1]');
R = o.responsive_share;
want = nan(numel(D.second_levels), numel(D.first_levels), 3);
for c = 1:3
    Tc = D.T(strcmp(D.T.condition, three{c}), :);
    for b = 1:numel(D.second_levels)
        for g = 1:numel(D.first_levels)
            rows = Tc.bias == D.second_levels(b) & Tc.gain == D.first_levels(g) & ~Tc.failed;
            want(b, g, c) = median(1 - Tc.frac_pinned(rows) - Tc.frac_silent(rows));
        end
    end
end
n_fail = n_fail + check(isequal(R.conditions, three) && max(abs(R.median(:) - want(:))) < 1e-12, ...
    'maps: the responsive share is the median over seeds of 1 - pinned - silent');
% The saturated corner of the synthetic sweep (80 % pinned, 2 % silent) reads 0.18 there.
n_fail = n_fail + check(abs(R.median(1, end, 1) - 0.18) < 1e-12, 'maps: the saturated corner has a responsive share of 0.18');
close(o.figs);
D1 = read_mc_gain_sweep(struct('T', Ts(Ts.seed == 1, :), 'settings', settings));
n_fail = n_fail + check(strcmp(error_of(@() fig_mc_gain_bias_maps('data', D1, quiet{:}, 'rows', {'iqr'})), ...
    'fig_mc_gain_bias_maps:NoRows'), 'maps: a spread row alone with one seed -> NoRows');
n_fail = n_fail + check(~isempty(error_of(@() fig_mc_gain_bias_maps('data', D, quiet{:}, 'rows', {'nonsense'}))), ...
    'maps: an unknown row name is refused');

%% 2. composite_tools
src = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 600 400]);
a = gobjects(2, 2);
for r = 1:2
    for c = 1:2
        a(r, c) = axes(src, 'Units', 'pixels', 'Position', [60 + (c - 1) * 280, 230 - (r - 1) * 180, 200, 120]);
        plot(a(r, c), 1:3, [r c r], 'DisplayName', sprintf('line %d %d', r, c));
    end
end
legend(a(1, 2), 'Location', 'northeast');
G = T.grid(a(:));
n_fail = n_fail + check(isequal(G, a), 'tools: grid reads rows (top first) and columns from the positions');
dst = T.canvas(900, 500, false);
C = T.copy(src, dst);
ok = numel(C.axes) == 4 && isscalar(C.legends) && isempty(C.colorbars);
for k = 1:numel(C.axes)
    twin = a(arrayfun(@(h) strcmp(h.Children(1).DisplayName, C.axes(k).Children(1).DisplayName), a));
    ok = ok && max(abs(C.axes(k).Position - twin.Position)) < 1e-6;
end
n_fail = n_fail + check(ok, 'tools: copy keeps every axes at its pixel position, with the legend');
one = T.copy(a(2, :)', dst, false);
n_fail = n_fail + check(numel(one.axes) == 2 && isempty(one.legends), 'tools: copy of some axes, without their legends');
with = T.copy(a(1, :)', dst);
n_fail = n_fail + check(numel(with.axes) == 2 && isscalar(with.legends), 'tools: copy of some axes brings their legend');
T.place(C.axes, [100 50 700 400]);
P = vertcat(C.axes.Position);
box = [min(P(:, 1)), min(P(:, 2)), max(P(:, 1) + P(:, 3)) - min(P(:, 1)), max(P(:, 2) + P(:, 4)) - min(P(:, 2))];
n_fail = n_fail + check(max(abs(box - [100 50 700 400])) < 1e-6, 'tools: place maps the block''s bounding box onto the rectangle');
h = T.letter(C.axes(1), 'A');
n_fail = n_fail + check(strcmp(h.String, '(A)') && h.FontSize == 14 && strcmp(h.FontWeight, 'normal'), ...
    'tools: the panel letter is ''(A)'', 14 pt, normal weight');
tl = T.condition_title('sfa3_stdc2', true);
n_fail = n_fail + check(strcmp(T.condition_title('no_adaptation'), 'No adaptation') && ...
    isequal(tl, {'Multiple-timescale', 'adaptation'}), 'tools: condition titles in sentence case, two lines when asked');
n_fail = n_fail + check(isequal(T.sentence({'noise-free', 'n_a = 3', '\lambda_1', 'Kept', 'top-20 mismatch'}), ...
    {'Noise-free', 'n_a = 3', '\lambda_1', 'Kept', 'Top-20 mismatch'}) && strcmp(T.sentence('slope 2'), 'Slope 2'), ...
    'tools: sentence capitalizes a first word and leaves a symbol');
m = T.margin(dst, 4);
n_fail = n_fail + check(isscalar(findall(dst, 'Tag', 'composite_crop_margin')) && ~isequal(m.Color, [1 1 1]) && min(m.Color) > 0.98, ...
    'tools: margin puts one almost-white line on the canvas (the crop counts it, the eye does not)');
close([src dst]);

%% The fixture run
src_run = find_regression_run();
if isempty(src_run)
    fprintf('SKIP: no regression run directory (data/regression_fast/reference_post_refactor); the composites were not built\n');
    finish(n_fail);
    return
end
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));
run_dir = fullfile(tmp, 'run'); mkdir(run_dir);
for pat = {'1D_sensitivity_*', 'param_space_*', 'tau_sensitivity_*', 'memory_capacity', 'mc_example', ...
        'single_neuron', 'numerics_verification'}
    d = dir(fullfile(src_run, '*'));
    d = d([d.isdir] & ~cellfun(@isempty, regexp({d.name}, ['^' regexptranslate('wildcard', pat{1}) '$'], 'once')));
    for k = 1:numel(d); copyfile(fullfile(d(k).folder, d(k).name), fullfile(run_dir, d(k).name)); end
end
copyfile(fullfile(src_run, 'run_manifest.mat'), run_dir);
write_grid(run_dir, 'selfcoupling', self_coupling_grid(three));
write_grid(run_dir, 'noise1d_wiener', noise_grid(three));
listing_before = sort({dir(run_dir).name});
reg_preset = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
chain_preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
S = load(fullfile(run_dir, 'single_neuron', 'single_neuron_data.mat'), 'results');
pulse = {S.results.name};
pulse = pulse(1:min(3, numel(pulse)));

%% 3. every composite builds
base = {'run_dir', run_dir, 'visible', false, 'verbose', 'near-none'};
jobs = { ...   function, its own arguments, the panel letters, the file tag
    @fig_main7_memory,  {'preset_name', reg_preset, 'conditions', three, 'run_mode', 'fast'}, 'ABCDEF', 'Fig_Main7_Grouped'; ...
    @fig_supp1_grouped, {'preset_name', chain_preset, 'pulse_conditions', pulse},              'ABC',    'Fig_Supp1_Grouped'; ...
    @fig_supp2_grouped, {'conditions', three},                                                 'ABC',    'Fig_Supp2_Grouped'; ...
    @fig_supp4_grouped, {'preset_name', reg_preset, 'conditions', three},                      'AB',     'Fig_Supp4_Grouped'; ...
    @fig_supp5_grouped, {'preset_name', reg_preset, 'conditions', three},                      'ABC',    'Fig_Supp5_Grouped'; ...
    @fig_supp9_grouped, {'preset_name', reg_preset, 'conditions', three, 'run_mode', 'fast'},  'ABCDE',  'Fig_Supp9_Grouped'};
ppi = get(groot, 'ScreenPixelsPerInch');
for j = 1:size(jobs, 1)
    fn = jobs{j, 1}; args = [base, jobs{j, 2}]; letters = jobs{j, 3}; tag = jobs{j, 4};
    name = func2str(fn);
    open_before = findall(groot, 'Type', 'figure');
    try
        o = fn(args{:}, 'save', false);
    catch err
        n_fail = n_fail + check(false, sprintf('%s builds (save false): %s', name, err.message));
        continue
    end
    found = sort(arrayfun(@(t) t.String(2), findall(o.figs, 'Tag', 'composite_panel_letter'))');
    ok = isscalar(o.figs) && isgraphics(o.figs) && isequal(found, letters) && isempty(o.files) && ...
        isequal(sort(fieldnames(o.panels))', num2cell(letters));
    n_fail = n_fail + check(ok, sprintf('%s builds: one figure, letters %s, a source named for each, no file', name, letters));
    if strcmp(name, 'fig_supp4_grouped'); n_fail = n_fail + gain_panel_check(o.figs, run_dir, three); end
    n_fail = n_fail + polish_checks(name, o.figs, pulse);
    canvas = o.figs.Position(3:4);
    close(o.figs);
    n_fail = n_fail + check(isempty(setdiff(findall(groot, 'Type', 'figure'), open_before)), ...
        sprintf('%s leaves no source figure open', name));

    out_dir = fullfile(tmp, 'figs', name);
    o = fn(args{:}, 'save', true, 'out_dir', out_dir);
    close(o.figs);
    png = fullfile(out_dir, [tag '.png']);
    ok = isfile(png) && isfile(fullfile(out_dir, [tag '.fig'])) && isfile(fullfile(out_dir, [tag '_table.md'])) && ...
        any(strcmp(o.files, [tag '_table.md']));
    n_fail = n_fail + check(ok, sprintf('%s saves its PNG, .fig and _table.md', name));
    if isfile(png)
        info = imfinfo(png);
        res = info.Width / (canvas(1) / ppi);          % dots per inch of canvas, at least (the PNG is cropped to the content)
        % Whole: the content spans most of the canvas height at that
        % resolution. A figure cut at the screen would be far shorter.
        n_fail = n_fail + check(info.Height > 0.8 * canvas(2) / ppi * res * 0.9, ...
            sprintf('%s: the PNG is as tall as its canvas (%d x %d px for %g x %g)', name, info.Width, info.Height, canvas));
    end
    if isfile(fullfile(out_dir, [tag '_table.md']))
        txt = fileread(fullfile(out_dir, [tag '_table.md']));
        n_fail = n_fail + check(all(arrayfun(@(L) contains(txt, ['## Panel ' L]), letters)), ...
            sprintf('%s: the table has a section for every panel', name));
    end
end
n_fail = n_fail + check(isequal(sort({dir(run_dir).name}), listing_before), ...
    'the run directory holds what it held before (no prep folder left)');

%% 5. a run without the stages
empty_run = fullfile(tmp, 'empty'); mkdir(empty_run);
ids = cellfun(@(f) error_of(@() f('run_dir', empty_run, 'save', false, 'conditions', three, 'verbose', 'near-none')), ...
    {@fig_main7_memory, @fig_supp2_grouped, @fig_supp4_grouped, @fig_supp5_grouped, @fig_supp9_grouped, ...
     @(varargin) fig_supp1_grouped(varargin{:}, 'preset_name', chain_preset)}, 'UniformOutput', false);
n_fail = n_fail + check(all(~cellfun(@isempty, ids)), ...
    sprintf('a run directory without the stages is an error for every composite (%s)', strjoin(unique(ids), ', ')));
close(setdiff(findall(groot, 'Type', 'figure'), figs_before));

finish(n_fail);

%% ------------------------------------------------------------------------
function n = polish_checks(name, fig, pulse)
% What the polish of 2026-11-22 (Plan 2, group 3) put into the composites.
n = 0;
st = manuscript_style();
cond = cell2mat(cellfun(@(k) st.condition_color(k), {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}', 'UniformOutput', false));
ax = findall(fig, 'Type', 'axes');
switch name
    case 'fig_supp1_grouped'
        % (A): a single mechanism is in a colour of its own, full strength,
        % and its column is titled by the mechanism.
        for want = {'sfa1_std0', 'SFA alone'; 'sfa0_std1', 'Depression alone'}'
            if ~ismember(want{1}, pulse); continue; end
            a = ax(arrayfun(@(h) strcmp(strjoin(cellstr(h.Title.String), ' '), want{2}), ax));
            ok = isscalar(a);
            if ok
                col = a.Title.Color;
                g = composite_grid_column(ax, a);
                traces = findall(g, 'Type', 'line');
                cols = unique(vertcat(traces.Color), 'rows');
                ok = ~ismember(col, cond, 'rows') && max(col) - min(col) > 0.4 && ...       % not a condition's, not pale
                    all(ismember(cols, [col; 0 0 0], 'rows'));
            end
            n = n + check(ok, sprintf('S1 (A): the %s column is titled ''%s'' and drawn in one colour that is no condition''s', want{1}, want{2}));
        end
    case 'fig_supp2_grouped'
        % (B): the QR lines are neutral (gray), never a condition's hue.
        thin = findall(fig, 'Type', 'line', 'LineWidth', 0.8);
        c = vertcat(thin.Color);
        n = n + check(~isempty(c) && all(max(c, [], 2) - min(c, [], 2) < 1e-12), 'S2 (B): the QR local rates are drawn in grays');
    case 'fig_supp4_grouped'
        n = n + check(isscalar(findall(fig, 'Tag', 'composite_crop_margin')), 'S4: a margin is kept under the bottom row''s x labels');
    case 'fig_supp5_grouped'
        vals = ax(arrayfun(@(h) isstruct(h.UserData) && isfield(h.UserData, 'values'), ax));
        ok = numel(vals) == 3;
        for k = 1:numel(vals)
            h = findall(vals(k), 'Type', 'histogram');
            v = vals(k).UserData.values;
            ok = ok && isscalar(h) && abs(sum(h.Values) - 1) < 1e-9 && ...                  % nothing gathered, nothing dropped
                max(abs(v)) <= 10^max(abs(vals(k).XLim)) * 0.1 && ~isempty(vals(k).XTickLabel);
        end
        n = n + check(ok, 'S5 (A): three symmetric-log histograms that hold every value of the sample');
        lg = findall(fig, 'Type', 'legend');
        ok = isscalar(lg);
        if ok
            lg.Units = 'pixels'; own = lg.Axes; own.Units = 'pixels';
            ok = lg.Position(2) >= own.Position(2) + own.Position(4) - 1 && isequal(lg.String, {'Local', 'Accumulated'});
        end
        n = n + check(ok, 'S5 (C): one key, above its panel, clear of the data');
    case 'fig_supp9_grouped'
        boxed = ax(arrayfun(@(h) contains(strjoin(cellstr(h.XLabel.String), ' '), 'Lyapunov') && strcmp(h.Box, 'on'), ax));
        n = n + check(isempty(boxed), 'S9 (E): the scatter panels are open, not boxed');
        key = findall(fig, 'Tag', 'mc_sweep_scatter_size_legend');
        ok = isscalar(key);
        if ok
            drawnow;
            tx = findall(key, 'Type', 'text');
            tx = tx(arrayfun(@(t) strcmp(t.VerticalAlignment, 'top'), tx));
            set(tx, 'Units', 'pixels');
            e = sortrows(cell2mat(arrayfun(@(t) t.Extent, tx, 'UniformOutput', false)), 1);
            ok = size(e, 1) >= 2 && all(e(2:end, 1) - (e(1:end - 1, 1) + e(1:end - 1, 3)) > 2);
        end
        n = n + check(ok, 'S9 (E): the labels of the marker-size key do not run into each other');
        lg = findall(fig, 'Type', 'legend');
        state = lg(arrayfun(@(l) any(strcmp(l.String, 'SFA')), lg));
        ok = isscalar(state);
        if ok
            state.Units = 'pixels'; own = state.Axes; own.Units = 'pixels';
            ok = state.Position(1) >= own.Position(1) + own.Position(3) && ...
                state.Position(2) >= own.Position(2) && state.Position(2) + state.Position(4) <= own.Position(2) + own.Position(4) + 1;
        end
        n = n + check(ok, 'S9 (B): its key is beside it, level with it');
end
end

function g = composite_grid_column(ax, a)
% The axes of ax in the column of a (the same left edge, within 2 px).
set(ax, 'Units', 'pixels');
g = ax(arrayfun(@(h) abs(h.Position(1) - a.Position(1)) < 2 && h.Position(3) == a.Position(3), ax));
end

function finish(n_fail)
if n_fail > 0
    error('test_memory_supplement_composites:Failed', '%d check(s) failed.', n_fail);
end
fprintf('test_memory_supplement_composites: all checks PASS\n');
end

function n = gain_panel_check(fig, run_dir, names)
% The 'All Weights' panel of block (A): its three lines are the medians of
% lambda_1 over the repetitions at each level of the gain sweep.
ax = findall(fig, 'Type', 'axes');
ax = ax(arrayfun(@(a) strcmp(strjoin(cellstr(a.XLabel.String), ' '), 'All Weights'), ax));
ok = numel(ax) == 2;
if ok
    [~, top] = max(arrayfun(@(a) a.Position(2), ax));
    lines = flipud(findall(ax(top), 'Type', 'line'));
    lines = lines(arrayfun(@(h) ~isempty(h.DisplayName), lines));
    listing = dir(fullfile(run_dir, '1D_sensitivity_*level_of_chaos*'));
    psa = ParamSpaceAnalysis2.from_dir(fullfile(listing(1).folder, listing(1).name));
    x = psa.param_vectors{strcmp(psa.grid_params, 'level_of_chaos')};
    ok = numel(lines) == numel(names);
    for c = 1:numel(names)
        if ~ok; break; end
        want = nan(1, numel(x));
        for k = 1:numel(x)
            v = ParamSpaceAnalysis2.collect_level_values(psa, 'level_of_chaos', k, names{c}, 'LLE');
            if ~isempty(v); want(k) = median(v); end
        end
        ok = isequaln(lines(c).YData(:)', want) && isequal(lines(c).XData(:)', x(:)');
    end
end
n = check(ok, 'fig_supp4_grouped: the gain panel holds the medians of the saved gain sweep, per condition');
end

function D = self_coupling_grid(conditions)
% The synthetic gain by bias table, with its second axis named as the
% self-coupling grid names it.
[T, s] = make_synthetic_mc_gain_sweep(conditions);
T.self_coupling = round(10 * T.bias - 2);          % -4 .. 6, the preset at 0
T.bias = [];
T = movevars(T, 'self_coupling', 'After', 'gain');
s = rmfield(s, 'preset_bias');
s.first_axis = 'level_of_chaos'; s.second_axis = 'self_coupling';
s.preset_self_coupling = 0;
D = struct('T', T, 'settings', s);
end

function D = noise_grid(conditions)
% A 1-D sweep of the input noise at the preset's gain: a zero level, then
% log-spaced levels, three networks.
levels = [0, logspace(-4, -1.5, 6)];
rows = {};
for c = 1:numel(conditions)
    for k = 1:numel(levels)
        for seed = 1:3
            mc = (2 + 3 * c) * exp(-30 * levels(k)) + 0.05 * seed;
            rows(end + 1, :) = {1, levels(k), seed, conditions{c}, mc, 1, -0.1 * c, -0.1 * c, 0.2, 0.2, 0.02, 0.02, false}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'gain', 'sigma_u_noise', 'seed', 'condition', 'mc_total', 'horizon', ...
    'lle_input', 'lle_noinput', 'rate_E', 'rate_I', 'frac_pinned', 'frac_silent', 'failed'});
s = struct('conditions', {conditions}, 'preset_name', '', 'first_axis', 'level_of_chaos', ...
    'second_axis', 'sigma_u_noise', 'preset_gain', 1, 'preset_sigma_u_noise', 0, 'synthetic', true);
D = struct('T', T, 'settings', s);
end

function write_grid(run_dir, grid, D)
% Where the mc_gain_sweep stage writes a grid's result (resolve_mc_gain_sweep).
d = fullfile(run_dir, mc_gain_sweep_folder(grid, ''), 'mc_gain_sweep');
mkdir(d);
T = D.T; settings = D.settings; %#ok<NASGU>
save(fullfile(d, 'mc_gain_sweep_results.mat'), 'T', 'settings');
end

function d = find_regression_run()
% <project root>/data/..., else the main checkout's (a worktree's data/ is empty).
rel = fullfile('data', 'regression_fast', 'reference_post_refactor');
root = fileparts(which('setup_paths'));
d = fullfile(root, rel);
if isfolder(d); return; end
[s, common] = system(sprintf('git -C "%s" rev-parse --path-format=absolute --git-common-dir', root));
d = '';
if s == 0
    main_root = fileparts(strtrim(common));
    if isfolder(fullfile(main_root, rel)); d = fullfile(main_root, rel); end
end
end

function id = error_of(fn)
id = '';
try
    fn();
catch err
    id = err.identifier;
    if isempty(id); id = 'error:without:identifier'; end
end
end

function f = check(ok, what)
if ok
    fprintf('PASS  %s\n', what);
    f = 0;
else
    fprintf('FAIL  %s\n', what);
    f = 1;
end
end
