%% TEST_MEMORY_SUPPLEMENT_COMPOSITES The composites of Plan 2 package F build from a small run, and their two new pieces work.
%
% Main Figure 7 (fig_main7_memory), supplementary S1, S2, S4, S5, S9
% (fig_supp<k>_grouped) and the leading-direction supplement
% (fig_supp_leading_grouped; not yet numbered),
% the responsive-share row of fig_mc_gain_bias_maps
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
    @fig_supp1_grouped, {'preset_name', chain_preset, 'pulse_conditions', pulse},              'AB',     'Fig_Supp1_Grouped'; ...
    @fig_supp2_grouped, {'conditions', three},                                                 'ABC',    'Fig_Supp2_Grouped'; ...
    @fig_supp4_grouped, {'preset_name', reg_preset, 'conditions', three},                      'AB',     'Fig_Supp4_Grouped'; ...
    @fig_supp5_grouped, {'preset_name', reg_preset, 'conditions', three},                      'AB',     'Fig_Supp5_Grouped'; ...
    @fig_supp9_grouped, {'preset_name', reg_preset, 'conditions', three, 'run_mode', 'fast'},  'ABC',    'Fig_Supp9_Grouped'; ...
    @fig_supp_leading_grouped, {'preset_name', reg_preset, 'conditions', three},               'ABC',    'Fig_Supp_Leading_Grouped'};
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
    if strcmp(name, 'fig_main7_memory'); n_fail = n_fail + main7_polish_check(o.figs, three); end
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
    {@fig_main7_memory, @fig_supp2_grouped, @fig_supp4_grouped, @fig_supp5_grouped, @fig_supp9_grouped, @fig_supp_leading_grouped, ...
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
        % Tom's review of 2026-11-22 (Plan2_figure_feedback.md, S1): two
        % panels. (B) is the two widest axes, one above the other.
        set(ax, 'Units', 'pixels');
        P = vertcat(ax.Position);
        is_b = P(:, 3) > 1.5 * min(P(:, 3));
        b = ax(is_b); a = ax(~is_b);
        ok = numel(b) == 2 && abs(b(1).Position(1) - b(2).Position(1)) < 1 && abs(b(1).Position(3) - b(2).Position(3)) < 1 && ...
            abs(b(1).Position(2) - b(2).Position(2)) > b(1).Position(4);
        n = n + check(ok, 'S1 (B): the two step responses, one above the other');
        % (A): every row on one y axis from 0, two ticks, the first at 0.
        [~, ~, row] = unique(round(P(~is_b, 2)));
        ok = numel(a) == 3 * numel(pulse) && max(row) == 3;
        for r = 1:max(row)
            g = a(row == r);
            ok = ok && all(arrayfun(@(h) isequal(h.YLim, g(1).YLim) && h.YLim(1) == 0 && numel(h.YTick) == 2 && h.YTick(1) == 0, g));
        end
        n = n + check(ok, 'S1 (A): the columns of a row share one y axis that starts at 0, with ticks at 0 and one round value');
        lg = findall(fig, 'Type', 'legend');
        ok = isscalar(lg);
        if ok
            s = cellstr(lg.String);
            kids = findall(lg.Axes, 'Type', 'line');
            ok = any(strcmp(s, 'Rate r(t)')) && any(strcmp(s, 'Steady state')) && ~any(contains(s, 'Undepressed')) && ...
                ~any(contains(s, 'n_a')) && any(strcmp(s, 'Reserve b_2 (chain)')) == any(strcmp({kids.Tag}, 'std_step_reserve'));
            lg.Units = 'pixels';
            ok = ok && lg.Position(1) >= b(1).Position(1) + b(1).Position(3);
        end
        n = n + check(ok, 'S1 (B): one key, beside the plots, with the entries they use (no undepressed line, no SFA entries)');
    case 'fig_supp2_grouped'
        % (B): the QR lines are neutral (gray), never a condition's hue.
        thin = findall(fig, 'Type', 'line', 'LineWidth', 0.8);
        c = vertcat(thin.Color);
        n = n + check(~isempty(c) && all(max(c, [], 2) - min(c, [], 2) < 1e-12), 'S2 (B): the QR local rates are drawn in grays');
        % Tom's review of 2026-11-22 (Plan2_figure_feedback.md, S2).
        lg = findall(fig, 'Type', 'legend');
        all_keys = {};
        for k = 1:numel(lg); all_keys = [all_keys, cellstr(lg(k).String)]; end %#ok<AGROW>
        n = n + check(any(strcmp(all_keys, 'Slope 1.5 (SRA1 lower bound)')) && ~any(contains(all_keys, 'strong')), ...
            'S2 (A): the 1.5 reference is the lower bound of SRA1');
        band = findall(fig, 'Tag', 'numerics_slope2_band');
        ok = numel(band) == 3;
        for k = 1:numel(band)
            kids = findall(band(k).Parent, 'Type', 'line');          % top of the stack first
            data = kids(arrayfun(@(h) ~strcmp(h.Marker, 'none'), kids));
            ok = ok && isequal(kids(end), band(k)) && band(k).LineWidth > 2 * max([data.LineWidth]) && min(band(k).Color) > 0.6;
        end
        n = n + check(ok, 'S2 (A): the slope-2 reference is a wide pale band behind the data, in every column');
        b_key = lg(arrayfun(@(l) any(strcmp(l.String, 'QR \lambda_1')), lg));
        n = n + check(isscalar(b_key) && isequal(b_key.String, {'Benettin local \lambda_1', 'Benettin \lambda_1', 'QR local \lambda_1', 'QR \lambda_1'}), ...
            'S2 (B): four traces in the key, Benettin''s two and the QR method''s two');
        ax_b = ax(arrayfun(@(h) strcmp(strjoin(cellstr(h.XLabel.String), ' '), 'Time (s)'), ax));
        ok = numel(ax_b) == 3;
        n_wide = 0;
        for k = 1:numel(ax_b)
            kids = findall(ax_b(k), 'Type', 'line');
            kids = kids(arrayfun(@(h) any(isfinite(h.YData)), kids));   % not the key's samples (NaN proxies)
            is_thin = arrayfun(@(h) h.LineWidth == 0.8, kids);
            ok = ok && nnz(is_thin) == 1 && nnz(~is_thin) == 1 && find(is_thin) < find(~is_thin) && ...   % the QR line over Benettin's
                (isequal(ax_b(k).YLim, [-4 4]) || isequal(ax_b(k).YLim, [-15 5]));
            n_wide = n_wide + isequal(ax_b(k).YLim, [-15 5]);
        end
        n = n + check(ok && n_wide == 1, 'S2 (B): one thin QR line on top of Benettin''s thick one; -15 to 5 without adaptation, plus or minus 4 with');
        % Plan 3 review: the keys' samples in neutral grays (they apply to
        % every column), and at least two labelled decades on each axis of (A).
        ok = ~isempty(lg);
        for k = 1:numel(lg)
            s = lg(k).PlotChildren;
            s = s(isprop(s, 'Color'));
            c = vertcat(s.Color);
            ok = ok && all(max(c, [], 2) - min(c, [], 2) < 0.05);
        end
        n = n + check(ok, 'S2: the keys'' samples are gray, not a condition''s colour');
        ax_a = ax(arrayfun(@(h) strcmp(strjoin(cellstr(h.XLabel.String), ' '), 'Step (s)'), ax));
        n = n + check(numel(ax_a) == 3 && all(arrayfun(@(h) nnz(h.YTick >= h.YLim(1) & h.YTick <= h.YLim(2)) >= 2, ax_a)), ...
            'S2 (A): at least two labelled decades on every y axis');
        ax_c = ax(arrayfun(@(h) strcmp(strjoin(cellstr(h.XLabel.String), ' '), 'Snapshot'), ax));
        ok = numel(ax_c) == 3 && all(arrayfun(@(h) isequal(h.YLim, [1e-12 1e-8]), ax_c)) && ...
            isempty(findall(fig, 'Tag', 'numerics_jacobian_threshold')) && isempty(findall(ax_c, 'Type', 'constantline')) && ...
            ~any(contains(lower(all_keys), 'threshold'));
        n = n + check(ok, 'S2 (C): three panels on 1e-12 to 1e-8, named by snapshot, with no threshold in the panel or the key');
        words = arrayfun(@(h) string(strjoin(cellstr(h.XLabel.String), ' ')), ax);
        n = n + check(~any(contains(lower(words), 'state')), 'S2: no axis calls a sampled time point a state');
    case 'fig_supp4_grouped'
        n = n + check(isscalar(findall(fig, 'Tag', 'composite_crop_margin')), 'S4: a margin is kept under the bottom row''s x labels');
    case 'fig_supp5_grouped'
        % (A) is the firing-rate histogram stacked by E:I weight ratio, with
        % one color bar (Plan 3: the lambda_1 histogram is main Figure 4B).
        rate = ax(arrayfun(@(h) contains(strjoin(cellstr(h.XLabel.String), ' '), 'Firing Rate'), ax));
        cb = findall(fig, 'Type', 'colorbar');
        ok = numel(rate) == 3 && all(arrayfun(@(h) ~isempty(findall(h, 'Type', 'patch')), rate)) && ...
            isscalar(cb) && strcmp(cb.Label.String, 'E:I weight ratio') && isempty(findall(fig, 'Type', 'image')) && ...
            ~any(arrayfun(@(h) contains(strjoin(cellstr(h.XLabel.String), ' '), 'Lyapunov'), ax));
        n = n + check(ok, 'S5 (A): three firing-rate histograms colored by E:I weight ratio, one color bar, no lambda_1 histogram');
        lg = findall(fig, 'Type', 'legend');
        ok = isscalar(lg);
        if ok
            lg.Units = 'pixels'; own = lg.Axes; own.Units = 'pixels';
            ok = lg.Position(2) >= own.Position(2) + own.Position(4) - 1 && isequal(lg.String, {'Local', 'Accumulated'});
        end
        n = n + check(ok, 'S5 (B): one key, above its panel, clear of the data');
    case 'fig_supp9_grouped'
        boxed = ax(arrayfun(@(h) contains(strjoin(cellstr(h.XLabel.String), ' '), 'Lyapunov') && strcmp(h.Box, 'on'), ax));
        n = n + check(isempty(boxed), 'S9 (C): the scatter panels are open, not boxed');
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
        n = n + check(ok, 'S9 (C): the labels of the marker-size key do not run into each other');
        % Tom's review of 2026-11-22 (Plan2_figure_feedback.md, S9): memory
        % only. The leading-direction shares and the cumulative capacity
        % have left; the responsive-share maps mark self-coupling 0.
        words = arrayfun(@(h) string(strjoin(cellstr(h.XLabel.String), ' ')) + " " + ...
            string(strjoin(cellstr(h.YLabel.String), ' ')), ax);
        lgs = findall(fig, 'Type', 'legend');
        n = n + check(~any(contains(words, 'SFA Time')) && ~any(contains(words, 'Cumulative')) && ...
            isscalar(lgs) && isequal(cellstr(lgs.String), {'Input', 'Reconstruction'}), ...
            'S9: memory only (no leading-direction shares, no cumulative capacity); one key, (A)''s input and reconstruction');
        % Plan 3 review: (A)'s traces inside their axes, the scale bar in
        % a.u., and no explanatory note on the figure.
        tr_ax = ax(arrayfun(@(h) startsWith(strjoin(cellstr(h.Title.String), ' '), 'delay'), ax));
        inside = true;
        for k = 1:numel(tr_ax)
            for l = findall(tr_ax(k), 'Type', 'line')'
                y = l.YData(isfinite(l.YData));
                if numel(y) > 10; inside = inside && min(y) > tr_ax(k).YLim(1) && max(y) < tr_ax(k).YLim(2); end
            end
        end
        txt = arrayfun(@(t) strjoin(cellstr(t.String), ' '), findall(fig, 'Type', 'text'), 'UniformOutput', false);
        n = n + check(numel(tr_ax) == 4 && inside && any(strcmp(txt, '1 a.u.')) && ~any(contains(txt, 'AU')) && ...
            ~any(contains(txt, 'drawn at the edge')), 'S9 (A): traces not cut, "1 a.u."; (C): no note on the clipped mark');
        maps = ax(arrayfun(@(h) ~isempty(findall(h, 'Type', 'surface')), ax));     % fig_mc_gain_bias_maps draws its cells as a surface
        n = n + check(numel(maps) == 3 && all(arrayfun(@(h) any(h.YTick == 0), maps)), ...
            'S9 (B): the self-coupling axis of every map is marked at 0');
    case 'fig_supp_leading_grouped'
        % The three panels on the leading direction (Plan 3, question 10).
        lg = findall(fig, 'Type', 'legend');
        ok = isscalar(lg) && any(startsWith(lg.String, 'Spike-frequency')) && ~any(strcmp(lg.String, 'SFA'));
        n = n + check(ok, 'leading (B): one key, spike-frequency adaptation written out');
        tau = ax(arrayfun(@(h) contains(strjoin(cellstr(h.XLabel.String), ' '), 'SFA Time'), ax));
        ok = numel(tau) == 2 && isequal(tau(1).XLim, tau(2).XLim);
        n = n + check(ok, 'leading (A), (B): both against the longest SFA time constant, on one x range');
        c = ax(startsWith(get(ax, {'Tag'}), 'leading_direction_c'));
        ok = numel(c) == 3 && isempty(findall(fig, 'Tag', 'uncoupled_reference')) && ...
            all(arrayfun(@(h) isscalar(findall(h, 'Tag', 'leading_direction_points')), c));
        n = n + check(ok, 'leading (C): one panel per condition, no line at the uncoupled rate');
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

function n = main7_polish_check(fig, names)
% Polish of 2026-11-22. (C): the conditions by name, not '1TS' and 'MTS'.
% (E): gain tick labels under the memory-capacity maps too. (D): the
% labelled noise ticks read their values, and the axis has a break.
% Tom's review of 2026-11-22 (Figure 7 notes 1 to 6): four panels in the
% top row, no lambda_1 against the SFA time constant, self-coupling 0
% marked on every map.
ticks = findall(fig, 'Tag', 'main7_horizon_tick');
words = arrayfun(@(t) string(strjoin(cellstr(t.String), '')), ticks);
full = string(erase(grouped_condition_label(names, 'tick'), ' '));
n = check(numel(ticks) == numel(names) && all(ismember(full, words)) && ...
    isempty(findall(fig, 'Type', 'text', 'String', '1TS')), '(C) names its conditions in full, one text per tick');
top = arrayfun(@(c) findobj(fig, 'Tag', sprintf('main_mc_r1_c%d', c)), 1:3);
bottom = arrayfun(@(c) findobj(fig, 'Tag', sprintf('main_mc_r2_c%d', c)), 1:3);
n = n + check(all(arrayfun(@(a, b) ~isempty(a.XTickLabel) && isequal(a.XTickLabel, b.XTickLabel), top, bottom)), ...
    '(E) carries the gain tick labels of (F)');
n = n + check(all(arrayfun(@(a) any(a.YTick == 0) && any(a.YTick < 0) && any(a.YTick > 0), [top bottom])), ...
    '(E), (F): the self-coupling axis is marked at 0 as well as at its ends');
ax = findall(fig, 'Type', 'axes');
maps = [top bottom];
others = ax(~ismember(ax, maps));
set(others, 'Units', 'pixels');
y = arrayfun(@(a) a.Position(2), others);
n = n + check(numel(others) == 4 && max(y) - min(y) < 1 && min(y) > max(arrayfun(@(a) getpixelposition(a, true) * [0; 1; 0; 1], top)), ...
    'the top row holds four panels, level, above the maps');
words = arrayfun(@(h) string(strjoin(cellstr(h.XLabel.String), ' ')), ax);
n = n + check(~any(contains(words, 'SFA Time')), 'lambda_1 against the longest SFA time constant has left Figure 7');
brk = findall(fig, 'Tag', 'main7_noise_break');
ok = ~isempty(brk);
if ok
    ax = ancestor(brk(1), 'axes');
    levels = getappdata(ax, 'noise_levels');
    lbl = cellstr(ax.XTickLabel); shown = find(~cellfun(@(s) isempty(strtrim(s)), lbl));
    ok = ~isempty(levels) && numel(shown) >= 2 && ...
        all(arrayfun(@(k) strcmp(lbl{k}, sprintf('%g', levels(k))), shown)) && ~any(contains(lbl, 'times'));
end
n = n + check(ok, '(D) labels its noise ticks with their values and marks the break after 0');
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
