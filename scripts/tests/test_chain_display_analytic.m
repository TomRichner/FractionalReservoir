% TEST_CHAIN_DISPLAY_ANALYTIC The figures and tables that state depression analytically, on chain routes (FD-006).
%
% Quick (about half a minute): no network is simulated. Checks that
%   * fig_STD_steady_state's chain steady state is Plan A section 1's
%     numbers, its single-pool and product curves are the closed forms, and
%     its chain step response agrees with an INDEPENDENT ode45 integration of
%     (C1), (C2) written out here (not through route_pool);
%   * the chain's steady state agrees with a long ode45 integration at
%     constant rate, so (C4) is checked against the equations, not only
%     against itself;
%   * with no selection on a three-condition PRODUCT preset the figure takes
%     the product path and writes the same three file names as before;
%   * fig_main1_concepts' depression endpoint is b1 of the chain (0.2 at
%     r_ref = 0.25) on the chain preset, and the product's 1/9 on the product
%     preset;
%   * the generated tables state (C1)-(C3) and the chain's parameters for a
%     chain condition, print the "product SQUARES the depression" prose only
%     where a product of two pools is covered, and name each condition's
%     depression in the conditions table.
%
% Asserts; silent until it throws, then prints one PASS line.
% See also: fig_STD_steady_state, write_manuscript_tables, fig_main1_concepts,
%           route_pool, test_std_chain_model

chain_preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
prod_preset  = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
three = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
four  = [three, {'sfa3_std2'}];
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));

%% fig_STD_steady_state, four conditions: the numbers
o = fig_STD_steady_state('preset_name', chain_preset, 'conditions', four, ...
    'save', false, 'visible', false);
close(o.figs);
assert(isequal(o.conditions, {'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}));
assert(isequal(o.skipped, {'no_adaptation'}));
S = o.steady_state;
r = S(1).probe.r;
assert(isequal(r, [0.1 0.25 0.5 1]));
by = @(n) S(strcmp({S.name}, n));
% Plan A section 1 (three decimals) and the closed forms.
assert(max(abs(by('sfa3_stdc2').probe.D - [0.420 0.200 0.106 0.054])) < 6e-4, 'chain b1');
assert(max(abs(by('sfa1_std1').probe.D - 1 ./ (1 + 8 * r))) < 1e-14, 'single pool');
assert(max(abs(by('sfa3_std2').probe.D - 1 ./ (1 + 8 * r) .^ 2)) < 1e-14, 'product');
assert(max(abs(by('sfa3_stdc2').probe.theta - r .* by('sfa3_stdc2').probe.D)) < 1e-15);
% The chain's output rises to the ceiling 1/18 and has no peak; the product peaks.
assert(~by('sfa3_stdc2').has_peak && abs(by('sfa3_stdc2').ceiling - 1/18) < 1e-8);
assert(~by('sfa1_std1').has_peak && abs(by('sfa1_std1').ceiling - 1/8) < 1e-8);
assert(by('sfa3_std2').has_peak && isnan(by('sfa3_std2').ceiling));

%% The chain's steady state against the equations, integrated
tau_rec1 = 2; tau_rec2 = 5; tau_rel = 0.25; n2 = 2;
c12 = @(b, rr) [b(2) * (1 - b(1)) / tau_rec1 - b(1) * rr / tau_rel; ...
                (1 - b(2)) / tau_rec2 - b(2) * (1 - b(1)) / (n2 * tau_rec1)];
ode = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);
for k = 1:numel(r)
    [~, Y] = ode45(@(~, b) c12(b, r(k)), [0 400], [1; 1], ode);
    assert(abs(Y(end, 1) - by('sfa3_stdc2').probe.D(k)) < 1e-8, 'b1 at r = %g', r(k));
    assert(abs(Y(end, 2) - by('sfa3_stdc2').probe.b2(k)) < 1e-8, 'b2 at r = %g', r(k));
end

%% The chain's step response against an independent integration of (C1), (C2)
tr = o.traces(strcmp({o.traces.name}, 'sfa3_stdc2'));
t = o.t; rt = o.r_t;
edges = [0, find(diff(rt) ~= 0), numel(rt)];      % segments of constant rate
Y = zeros(numel(t), 2); b = [1; 1]; t_prev = 0;
for s = 1:numel(edges) - 1
    idx = edges(s) + 1:edges(s + 1);
    [~, y] = ode45(@(~, bb) c12(bb, rt(idx(1))), [t_prev, t(idx)], b, ode);
    Y(idx, :) = y(2:end, :);
    b = y(end, :)'; t_prev = t(idx(end));
end
err = max(abs(tr.s - Y), [], 'all');
assert(err < 1e-9, 'chain step response differs from ode45 of (C1),(C2) by %g', err);
assert(max(abs(tr.D - Y(:, 1)')) < 1e-9, 'the depression factor is b1');
assert(max(abs(tr.theta - rt .* Y(:, 1)')) < 1e-9);
% The two share the releasable-pool equation and the reserve starts full: at
% the first onset the chain and the single pool coincide (to O(dt) in time
% after the step, while b2 is still ~1).
t_on = find(rt > 0, 1);
one = o.traces(strcmp({o.traces.name}, 'sfa1_std1'));
assert(abs(tr.D(t_on + 20) - one.D(t_on + 20)) < 1e-3, 'chain and single pool at onset');
assert(tr.D(edges(3)) < one.D(edges(3)) - 0.05, 'chain deeper by the end of the first step');

%% The three-condition set draws the chain and the single pool only
o3 = fig_STD_steady_state('preset_name', chain_preset, 'conditions', three, ...
    'out_dir', fullfile(tmp, 'std3'), 'visible', false);
close(o3.figs);
assert(isequal(o3.conditions, {'sfa1_std1', 'sfa3_stdc2'}));
assert(isequal(sort(o3.files), sort({'Fig_STD_steady_state.fig', 'Fig_STD_steady_state.png', ...
    'Fig_STD_steady_state.svg'})));

%% No selection on a product preset: the product figure, the same file names
op = fig_STD_steady_state('preset_name', prod_preset, 'out_dir', fullfile(tmp, 'stdP'), ...
    'visible', false);
close(op.figs);
assert(isfield(op, 'rho') && isequal(op.rho, [0.125 0.125]), 'product path');
assert(isequal(sort(op.files), sort(o3.files)));

%% fig_main1_concepts: the depression endpoint from the route model
intro = fullfile(tmp, 'intro');
mkdir(fullfile(intro, 'eigenspectra')); mkdir(fullfile(intro, 'statetraces'));
for f = {fullfile(intro, 'eigenspectra', 'panelA_eigenspectrum.fig'), ...
         fullfile(intro, 'statetraces', 'panelA_bottom_traces.fig')}
    h = figure('Visible', 'off');
    for k = 1:3; subplot(1, 3, k, 'Parent', h); end
    savefig(h, f{1}); close(h);
end
m1 = fig_main1_concepts(chain_preset, intro); close(m1.figs);
assert(abs(m1.b_levels(end) - 0.2) < 1e-12, 'chain endpoint b1(0.25) = 0.2');
assert(contains(m1.notes{2}, 'chained vesicle pools') && contains(m1.notes{2}, 'n2=2'));
assert(~contains(m1.notes{3}, 'two-timescale product'));
m1p = fig_main1_concepts(prod_preset, intro); close(m1p.figs);
assert(abs(m1p.b_levels(end) - 1/9) < 1e-15, 'product endpoint 1/9');
assert(contains(m1p.notes{3}, 'two-timescale product'));

%% The generated tables
read = @(d, f) fileread(fullfile(d, f));
d3 = fullfile(tmp, 'tab3');
write_manuscript_tables('preset_name', chain_preset, 'out_dir', d3, 'conditions', three);
eq3 = read(d3, 'equation_table.md');
assert(contains(eq3, '\frac{db_{i2}}{dt} &= \frac{1 - b_{i2}}{\tau_{rec_2}}'), '(C2) stated');
assert(contains(eq3, '\frac{db_{i1}}{dt} &= \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}}'), '(C1) stated');
assert(contains(eq3, 'D_{i} = b_{i1} \quad \text{(chained pools)}'), '(C3) stated');
assert(contains(eq3, '| $\tau_{rec_2}$ (E to E) | Chained STD: reserve refill | 5 | s |'));
assert(contains(eq3, '| $n_2$ (E to E) | Chained STD: reserve capacity / releasable-pool capacity | 2 | -- |'));
assert(contains(eq3, '| $\tau_{rel}$ (E to E) | Chained STD: release (releasable pool only) | 0.25 | s |'));
assert(contains(eq3, '$b_1$ = 0.420, 0.200, 0.106, 0.054'), 'chain steady state');
assert(~contains(eq3, 'SQUARE'), 'no product of two pools is covered: no squaring prose');
cond3 = read(d3, 'adaptation_conditions_table.md');
assert(contains(cond3, 'chained pools: $\tau_{rec}$ = [2, 5] s, $\tau_{rel}$ = 0.25 s, $n_2$ = 2'));
assert(~contains(cond3, 'product of 2 pools'));

d4 = fullfile(tmp, 'tab4');
write_manuscript_tables('preset_name', chain_preset, 'out_dir', d4, 'conditions', four);
eq4 = read(d4, 'equation_table.md');
assert(contains(eq4, 'SQUARE') && contains(eq4, 'independent** pools (not a chain), the depression'));
assert(contains(read(d4, 'adaptation_conditions_table.md'), ...
    'product of 2 pools: $\tau_{rec}$ = [2, 4] s, $\tau_{rel}$ = [0.25, 0.5] s'));

dp = fullfile(tmp, 'tabP');
write_manuscript_tables('preset_name', prod_preset, 'out_dir', dp);
eqp = read(dp, 'equation_table.md');
assert(contains(eqp, '\theta_i &= r_{i} \prod_{m=1}^{M} b_{im}'));
assert(contains(eqp, 'SQUARE') && ~contains(eqp, 'chain') && ~contains(eqp, 'b_{i2}'));

disp('PASS: test_chain_display_analytic -- chain steady state and step response, product path, concept endpoint, tables.');
