% TEST_FIG_MAIN2_B1 Main Fig 2 row (F) plots b1, the first depression variable (FD-031 item 12).
%
% Quick, no simulation: a synthetic representative_dynamics file with three
% conditions in which b1 and the route factor `std` differ (a product route
% has std = b1 b2). Checks that
%   1. row (F) of every depressed column draws the saved b1 of the displayed
%      neurons (13 E + 13 I of 25 + 25), not std, with y limits [0 1];
%   2. the no-adaptation column, and a column whose std_model is 'none'
%      (SFA only), leave row (F) blank;
%   3. the provenance notes name b1;
%   4. a saved file without b1 (predating item 12) is refused with
%      fig_main2_dynamics:NoB1.
% Asserts; prints one PASS line.
%
% See also: fig_main2_dynamics, run_paper_illustrations

tmp = tempname; mkdir(fullfile(tmp, 'representative_dynamics'));
cleanup = onCleanup(@() rmdir(tmp, 's'));
f = fullfile(tmp, 'representative_dynamics', 'representative_dynamics_data.mat');

nt = 301; t = linspace(0, 30, nt); n = 50; sel = {1:25, 1:25};
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
settings = struct('display_window', [0 30]);
save(f, 'results', 'settings');

%% 1-3. Row (F) is b1
out = fig_main2_dynamics(tmp, tmp, names);
fig = out.figs;
shown = [1:13, 26:38];                % round(25/2) = 13 of each type, as the figure draws
for c = 2:3
    ax = findobj(fig, 'Type', 'axes', 'Tag', sprintf('dynamics_r6_c%d', c));
    L = findobj(ax, 'Type', 'line');
    assert(numel(L) == 26, 'row F of column %d draws %d lines, not 26', c, numel(L));
    Y = sortrows(vertcat(L.YData));
    assert(isequal(Y, sortrows(results(c).b1(shown, :))), 'row F of column %d is not b1', c);
    assert(~isequal(Y, sortrows(results(c).std(shown, :))));
    assert(isequal(ax.YLim, [0 1]), 'row F y limits must be [0 1]');
end
ax = findobj(fig, 'Type', 'axes', 'Tag', 'dynamics_r6_c1');
assert(isempty(findobj(ax, 'Type', 'line')), 'no adaptation: row F must be blank');
assert(any(contains(out.notes, 'b1, the first depression variable')), 'the provenance must name b1');
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
