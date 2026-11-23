% TEST_FIG_MAIN2_B1 Main Fig 2 row (F) plots b1, the first depression variable (FD-031 item 12).
%
% Quick, no simulation: a synthetic representative_dynamics file with three
% conditions in which b1 and the route factor `std` differ (a product route
% has std = b1 b2). Checks that
%   1. row (F) of every depressed column draws the saved b1 of the displayed
%      neurons (8 E + 8 I of 25 + 25: about 60 percent of the 13 + 13 drawn
%      before, Tom 2026-11-22), not std, with y limits [0 1];
%   2. the no-adaptation column, and a column whose std_model is 'none'
%      (SFA only), leave row (F) blank;
%   3. the provenance notes name b1;
%   4. a saved file without b1 (predating item 12) is refused with
%      fig_main2_dynamics:NoB1;
%   5. the y axes are as Tom asked on 2026-11-22 (0 and 1 on unit axes, the
%      input row 0 to 0.35 with 0 and 0.25, tick labels on the left column
%      only, row (G) on +-5 without adaptation and +-1 with), the type is at
%      its role's size, and the column dividers stand in a gap, clear of the
%      tick labels.
% Asserts; prints one PASS line.
%
% See also: fig_main2_dynamics, run_paper_illustrations

tmp = tempname; mkdir(fullfile(tmp, 'representative_dynamics'));
cleanup = onCleanup(@() rmdir(tmp, 's'));
f = fullfile(tmp, 'representative_dynamics', 'representative_dynamics_data.mat');

nt = 301; t = linspace(0, 20, nt); n = 50; sel = {1:25, 1:25};
names = {'no_adaptation', 'sfa1_std1', 'sfa3_std2'};
models = {'none', 'product', 'product'};
rng(7);
for c = 1:3
    b1 = 0.2 + 0.6 * rand(n, nt); b2 = 0.5 * ones(n, nt);
    if c == 1, b1 = ones(n, nt); b2 = b1; end
    r = struct('name', names{c}, 'title', names{c}, 't', t, 'u', zeros(n, nt), ...
        'x', zeros(n, nt), 'r', zeros(n, nt), 'syn', zeros(n, nt), 'sfa', zeros(n, nt), ...
        'std', b1 .* b2, 'b1', b1, 'std_model', models{c}, 'selected', {sel}, 'n', n, ...
        'cell_type_names', {{'E', 'I'}}, 't_lya', t(:), 'local_lambda', zeros(nt, 1), ...
        'finite_lambda', zeros(nt, 1), 'LLE', 0);
    results(c) = r; %#ok<SAGROW>
end
settings = struct('display_window', [0 20]);
save(f, 'results', 'settings');

%% 1-3. Row (F) is b1
out = fig_main2_dynamics(tmp, tmp, names);
fig = out.figs;
shown = [1:8, 26:33];                  % round(0.3*25) = 8 of each type, as the figure draws
for c = 2:3
    ax = findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r6_c%d', c));
    L = findobj(ax, 'Type', 'line');
    assert(numel(L) == 16, 'row F of column %d draws %d lines, not 16', c, numel(L));
    Y = sortrows(vertcat(L.YData));
    assert(isequal(Y, sortrows(results(c).b1(shown, :))), 'row F of column %d is not b1', c);
    assert(~isequal(Y, sortrows(results(c).std(shown, :))));
    assert(isequal(ax.YLim, [0 1]), 'row F y limits must be [0 1]');
end
ax = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r6_c1');
assert(isempty(findobj(ax, 'Type', 'line')), 'no adaptation: row F must be blank');
assert(any(contains(out.notes, 'b1, the first depression variable')), 'the provenance must name b1');
% Row (G): the accumulated lambda_1 in the column's condition color, the local
% rate a lighter tint of it; neither black nor gray (polish list, 2026-11-22).
st = manuscript_style();
for c = 1:3
    ax = findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r7_c%d', c));
    acc = findobj(ax, 'Tag', 'accumulated_lambda'); loc = findobj(ax, 'Tag', 'local_lambda');
    col = st.condition_color(names{c});
    assert(isscalar(acc) && isscalar(loc), 'row G of column %d: one accumulated and one local trace', c);
    assert(max(abs(acc.Color - col)) < 1e-12, 'row G of column %d: the accumulated trace is not in the condition color', c);
    assert(max(abs(loc.Color - (.6 + .4 * col))) < 1e-12 && all(loc.Color >= acc.Color), ...
        'row G of column %d: the local trace is not a lighter tint of the condition color', c);
end
% The E/I key: one line at the top right of row (B)'s last axes, inside them
% (second round, Tom 2026-11-22: further left), the axes unmoved.
lg = findobj(fig, 'Type', 'legend'); axb = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r2_c3');
assert(isscalar(lg) && isequal(lg.String, {'E', 'I'}) && strcmp(lg.Orientation, 'horizontal'), 'one E/I key, on one line');
lg.Units = 'normalized'; axb.Units = 'normalized';
assert(abs(lg.Position(1) + lg.Position(3) - axb.Position(1) - axb.Position(3)) < 1e-9 && ...
    abs(lg.Position(2) + lg.Position(4) - axb.Position(2) - axb.Position(4)) < 1e-9 && lg.Position(1) > axb.Position(1), ...
    'the E/I key must stand at the top right inside row B''s last axes');
% The time bar: 5 s on a 20 s window, labelled "5 s", under row G's axes (Plan 3 fix X4).
bar = findobj(fig, 'Tag', 'time_bar'); lab = findobj(fig, 'Tag', 'time_bar_label');
assert(isscalar(bar) && abs(diff(bar.XData) - 5) < 1e-12 && strcmp(lab.String, '5 s'), 'a 5 s bar on a 20 s window');
axg = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r7_c1');
assert(all(bar.YData < axg.YLim(1)), 'the time bar must lie under row G''s axes');
axa = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r1_c3');
assert(max(abs(axb.Position([1 3]) - axa.Position([1 3]))) < 1e-9, 'the key must not resize row B''s axes');

% 5. The y axes (Tom's notes on Figure 2, 2026-11-22).
want_lim = {[0 .35], [-6 6], [0 1], [0 1], [0 .6], [0 1]};
want_tick = {[0 .25], [-5 0 5], [0 1], [0 1], [0 .5], [0 1]};
for c = 1:3
    for j = 1:6
        ax = findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r%d_c%d', j, c));
        assert(isequal(ax.YLim, want_lim{j}) && isequal(ax.YTick, want_tick{j}), ...
            'row %d of column %d: y limits %s, ticks %s', j, c, mat2str(ax.YLim), mat2str(ax.YTick));
        assert((c == 1) ~= isempty(ax.YTickLabel), 'row %d: y tick labels belong on the left column only (column %d)', j, c);
    end
    % Row (G): no adaptation on +-5, the adapting conditions on +-1, each
    % with its own tick labels since the axis is not the left column's.
    ax = findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r7_c%d', c));
    lim = 1 + 4 * (c == 1);
    assert(isequal(ax.YLim, lim * [-1 1]) && isequal(ax.YTick, lim * [-1 0 1]), 'row G of column %d: limits %s', c, mat2str(ax.YLim));
    assert(~isempty(ax.YTickLabel), 'row G of column %d must carry tick labels', c);
end
% The type by role, and the sheet as wide as the common canvas.
ax = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r4_c1');
assert(ax.FontSize == st.tick_fs && ax.YLabel.FontSize == st.label_fs && strcmp(ax.YLabel.Interpreter, 'tex'), ...
    'ticks at tick_fs, y labels at label_fs in the figure typeface');
ax = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r1_c2');
assert(ax.Title.FontSize == st.title_fs && iscell(ax.Title.String) && numel(ax.Title.String) == 2, ...
    'a long column title is on two lines at title_fs');
px = getpixelposition(fig);
assert(abs(px(3) / get(groot, 'ScreenPixelsPerInch') - st.canvas_w) < .02 && abs(out.width_scale - 1) < 1e-9, ...
    'three columns: the sheet is canvas_w wide');
% The dividers: placed again without a NoGap warning, each between the
% drawn extent of its left column and the tick labels of its right column.
lastwarn('');
column_dividers(fig);
[~, wid] = lastwarn;
assert(~strcmp(wid, 'column_dividers:NoGap'), 'the columns must leave a gap for the dividers');
dv = sort(arrayfun(@(l) l.XData(1), findall(fig, 'Tag', 'column_divider'))) * px(3);
assert(numel(dv) == 2, 'two dividers for three columns');
for c = 1:2
    lp = getpixelposition(findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r7_c%d', c)));
    ra = findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r7_c%d', c + 1));
    rp = getpixelposition(ra); ra.Units = 'pixels'; ti = ra.TightInset;
    assert(dv(c) > lp(1) + lp(3) + 2 && dv(c) < rp(1) - ti(1) - 2, ...
        'divider %d must lie between column %d and the tick labels of column %d', c, c, c + 1);
end
close(fig);

results(2).std_model = 'none';                          % stands in for an SFA-only column
save(f, 'results', 'settings');
out = fig_main2_dynamics(tmp, tmp, names); fig = out.figs;
ax = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r6_c2');
assert(isempty(findobj(ax, 'Type', 'line')), 'a condition without depression: row F must be blank');
close(fig);

%% 4. Saved data without b1 is refused
results = rmfield(results, 'b1'); results(2).std_model = 'product';
save(f, 'results', 'settings');
before = findall(groot, 'Type', 'figure');
try
    out = fig_main2_dynamics(tmp, tmp, names); close(out.figs);
    error('test_fig_main2_b1:NoError', 'data without b1 must be refused');
catch ME
    close(setdiff(findall(groot, 'Type', 'figure'), before));   % the half-drawn figure
    assert(strcmp(ME.identifier, 'fig_main2_dynamics:NoB1'), 'wrong error: %s', ME.message);
end
fprintf('test_fig_main2_b1: PASS\n');
