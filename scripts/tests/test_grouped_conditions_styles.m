%% TEST_GROUPED_CONDITIONS_STYLES The grouped main figures label by condition NAME (FD-006).
%
% Unit level, no run data, a few seconds. Checks
%   1. grouped_condition_label: every wording the grouped figures paint, for
%      the paper's three conditions and the product MTS, and the inverse
%      (title -> name), which errors on a title no condition carries;
%   2. style_main4_grouped, style_main5_grouped and style_main6_grouped on
%      synthetic figures shaped like the inner figures they style, with three
%      conditions (the paper's set, chain in the MTS slot) and four (the
%      review set, product last): the right count of panels is accepted, every
%      tick label / title / legend entry names the condition drawn there, in
%      order, coloured by name, and nothing cycles (the three-label overwrite
%      the styles used to make repeated labels with four ticks, inventory R5).
% Figures are invisible and closed.

n_fail = 0;
st = manuscript_style();
three = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
four  = [three, {'sfa3_std2'}];

%% 1. labels
L = @(k) grouped_condition_label(four, k);
n_fail = n_fail + check(isequal(L('title'), {'No Adaptation', 'Single-Timescale Adaptation', ...
    'Multiple-Timescale Adaptation', 'Multiple-Timescale Adaptation (product)'}), 'title');
n_fail = n_fail + check(isequal(L('sentence'), {'No adaptation', 'Single-timescale adaptation', ...
    'Multiple-timescale adaptation', 'Multiple-timescale adaptation (product)'}), 'sentence (groups 2, 3)');
n_fail = n_fail + check(isequal(L('tick'), {'No adaptation', 'Single-timescale', ...
    'Multiple-timescale', 'Multiple-timescale (product)'}), 'tick (group 5 margin)');
n_fail = n_fail + check(isequal(L('abbrev'), {'No Adaptation', '1TS', 'MTS', 'prodMTS'}), 'abbrev (group 6 horizon)');
tl = L('two_line');
n_fail = n_fail + check(isequal(tl{1}, 'No Adaptation') && isequal(tl{2}, {'Single-Timescale', 'Adaptation'}) ...
    && isequal(tl{4}, {'Multiple-Timescale', 'Adaptation (product)'}), 'two_line (group 5 gains)');
n_fail = n_fail + check(isequal(grouped_condition_label('sfa3_std0', 'sentence'), {'SFA 3\tau'}), ...
    'sentence keeps acronyms');
n_fail = n_fail + check(isequal(grouped_condition_label(L('title'), 'name_of_title'), four), 'inverse of title');
n_fail = n_fail + check(isequal(grouped_condition_label({{'Single-Timescale', 'Adaptation'}}, 'name_of_title'), ...
    {'sfa1_std1'}), 'inverse of a two-line title');
colored = sprintf('\\color[rgb]{%.6f,%.6f,%.6f}%s', [0 0.7 0], 'Multiple-Timescale Adaptation (product)');
n_fail = n_fail + check(isequal(grouped_condition_label({colored}, 'name_of_title'), {'sfa3_std2'}), ...
    'inverse of a coloured legend label');
n_fail = n_fail + check(strcmp(error_of(@() grouped_condition_label({'Not A Condition'}, 'name_of_title')), ...
    'grouped_condition_label:UnknownTitle'), 'unknown title -> UnknownTitle');
% Every title in the map names one condition (the inverse is well defined).
m = srnn_condition_titles(); v = values(m);
n_fail = n_fail + check(numel(unique(v)) == numel(v), 'srnn_condition_titles titles are unique');

%% 2. the styles, three and four conditions
for set_ = {three, four}
    names = set_{1}; N = numel(names); tag = sprintf('%d conditions', N);

    % ---- group 4: six sensitivity panels (overlaid conditions) + N rate panels
    f = figure('Visible', 'off', 'Position', [40 40 1300 950]);
    for k = 1:6
        ax = axes(f, 'Position', [.05 + mod(k-1, 3)*.32, .40 + (k > 3)*.30, .25, .2]); hold(ax, 'on');
        h = gobjects(1, N);
        for c = 1:N
            h(c) = plot(ax, 1:3, rand(1, 3), 'DisplayName', st.condition_title(names{c}));
        end
        yline(ax, 0, '--', 'HandleVisibility', 'off');
        if k == 3, legend(ax, h); end
    end
    for c = 1:N
        ax = axes(f, 'Position', [.05 + (c-1)*.9/N, .05, .8/N, .2]);
        plot(ax, 1:3, rand(1, 3)); title(ax, st.condition_title(names{c}));
    end
    style_main4_grouped(f);
    rate = findall(f, 'Type', 'axes'); rate = rate(arrayfun(@(a) a.Position(2) < .3, rate));
    [~, o] = sort(arrayfun(@(a) a.Position(1), rate)); rate = rate(o);
    ok = numel(rate) == N;
    for c = 1:N
        ok = ok && strcmp(strjoin(cellstr(rate(c).Title.String), ' '), st.condition_title(names{c})) && ...
            isequal(rate(c).Title.Color, st.condition_color(names{c}));
    end
    n_fail = n_fail + check(ok, ['style 4: rate panels titled and coloured by name, ' tag]);
    lg = findall(f, 'Type', 'legend');
    s = lg.String;
    ok = numel(s) == N + 1 && strcmp(s{end}, 'Zero Growth');
    for c = 1:N
        ok = ok && endsWith(s{c}, st.condition_title(names{c})) && ...
            contains(s{c}, sprintf('%.6f,%.6f,%.6f', st.condition_color(names{c})));
    end
    n_fail = n_fail + check(ok, ['style 4: legend entries by name, in order, plus Zero Growth, ' tag]);
    close(f);

    % ---- group 5: one margin panel + N gain panels
    f = figure('Visible', 'off');
    ax = axes(f, 'Position', [.05 .5 .2 .4]); plot(ax, 1:N, rand(1, N));
    set(ax, 'XTick', 1:N); title(ax, 'non-normal margin per state (5-95%, IQR, median)');
    for c = 1:N
        ax = axes(f, 'Position', [.3 + (c-1)*.6/N, .5, .5/N, .4]);
        h = plot(ax, 1:3, rand(1, 3), 'DisplayName', 'Worst case');
        title(ax, st.condition_title(names{c}));
        if c == 1, legend(ax, h); end
    end
    style_main5_grouped(f, names, names);
    axs = findall(f, 'Type', 'axes');
    margin = axs(arrayfun(@(a) any(contains(string(a.Title.String), 'non-normal')), axs));
    gains = axs(axs ~= margin); [~, o] = sort(arrayfun(@(a) a.Position(1), gains)); gains = gains(o);
    ticks = grouped_condition_label(names, 'tick');
    ok = numel(margin.XTickLabel) == N && numel(gains) == N;
    for c = 1:N
        ok = ok && endsWith(margin.XTickLabel{c}, ticks{c}) && ...
            contains(margin.XTickLabel{c}, sprintf('%.6f,%.6f,%.6f', st.condition_color(names{c}))) && ...
            strcmp(strjoin(cellstr(gains(c).Title.String), ' '), st.condition_title(names{c})) && ...
            isequal(gains(c).Title.Color, st.condition_color(names{c}));
    end
    n_fail = n_fail + check(ok, ['style 5: margin ticks and gain titles by name, no cycling, ' tag]);
    ok = all(arrayfun(@(a) a.Position(1) + a.Position(3) <= 1, axs)) && ...
        all(diff(arrayfun(@(a) a.Position(1), gains)) > 0);
    n_fail = n_fail + check(ok, ['style 5: 1+N panels inside the canvas, left to right, ' tag]);
    % restyling the styled figure (names read back from the gain titles)
    style_main5_grouped(f);
    n_fail = n_fail + check(endsWith(margin.XTickLabel{end}, ticks{end}), ['style 5: restyle reads names back, ' tag]);
    close(f);

    % ---- group 6: tau panels (one LLE + one block per tau condition) + 3 MC axes
    tau = names(end);          % the sweep holds the MTS of the set (chain, or the product last)
    if N == 4, tau = names(3:4); end
    nt = numel(tau);
    f = figure('Visible', 'off');
    for c = 1:nt
        ax = axes(f, 'Position', [.05 + (c-1)*.2, .6, .15, .3]);
        plot(ax, 1:3, rand(1, 3)); yline(ax, 0, '-'); xlabel(ax, 'slowest tau (s)');
        if nt > 1, title(ax, st.condition_title(tau{c})); end
        ax = axes(f, 'Position', [.55 + (c-1)*.2, .6, .15, .3]); hold(ax, 'on');
        for b = {'SFA', 'STD', 'x'}, plot(ax, 1:3, rand(1, 3), 'DisplayName', b{1}); end
        xlabel(ax, 'slowest tau (s)');
        if c == 1, legend(ax); end
    end
    hb = gobjects(1, N);
    for k = 1:3
        ax = axes(f, 'Position', [.05 + (k-1)*.32, .1, .25, .3]); hold(ax, 'on');
        for c = 1:N, hb(c) = plot(ax, 1:3, rand(1, 3)); end
        if k == 3, set(ax, 'XTick', 1:N, 'XTickLabel', mc_display_names(names)); end
        if k == 2, legend(ax, hb, mc_display_names(names)); end
    end
    style_main6_grouped(f, tau);
    axs = findall(f, 'Type', 'axes');
    top = axs(arrayfun(@(a) a.Position(2) > .5, axs));
    titles = arrayfun(@(a) strjoin(cellstr(a.Title.String), ' '), top, 'UniformOutput', false);
    want = cellfun(@(n) st.condition_title(n), tau, 'UniformOutput', false);
    ok = numel(top) == 2*nt && all(ismember(titles, want)) && ...
        ~any(strcmp(titles, 'Multiple-Timescale Adaptation') & ~ismember('sfa3_stdc2', tau));
    n_fail = n_fail + check(ok, ['style 6: tau panels titled by name, ' tag]);
    bottom = axs(arrayfun(@(a) a.Position(2) < .5, axs));
    [~, o] = sort(arrayfun(@(a) a.Position(1), bottom)); bottom = bottom(o);
    n_fail = n_fail + check(isequal(bottom(3).XTickLabel(:)', grouped_condition_label(names, 'abbrev')), ...
        ['style 6: horizon ticks by name, no cycling, ' tag]);
    close(f);
end

if n_fail == 0
    fprintf('\nALL PASS (test_grouped_conditions_styles)\n');
else
    error('test_grouped_conditions_styles:Failed', '%d check(s) FAILED', n_fail);
end

%% ------------------------------------------------------------------------
function n = check(ok, what)
if ok
    fprintf('PASS  %s\n', what); n = 0;
else
    fprintf('FAIL  %s\n', what); n = 1;
end
end

function id = error_of(fn)
id = '';
try
    fn();
catch ME
    id = ME.identifier;
end
end
