% demo_chain_1TS_vs_MTS.m
%
% A candidate definition of the single- versus multiple-timescale conditions
% for short-term depression, built on the chained-pool model:
%   db1/dt = b2*(1 - b1)/tau_rec1 - b1*r/tau_rel               (releasable pool)
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2)*b2*(1 - b1)/tau_rec1   (reserve)
%   theta  = r*b1
% n2 = reserve capacity / releasable-pool capacity.
%   - 1TS: the reserve is held full, b2 = 1 for all time. The b1 equation is
%     then exactly the current single-pool depression equation.
%   - MTS: the same b1 equation, with the reserve b2 dynamic.
% The two conditions are therefore NESTED: they differ by one mechanism, and
% MTS reduces to 1TS as n2 -> infinity.
%
% Two prescribed inputs, chosen by input_type in the parameter block:
%   'steps' - r held 5 s at each of [0 0.25 0 0.5 0 1 0].
%   'sine'  - r(t) = (r_max/2)*(1 - cos(2*pi*f_r*t)), 0 to r_max, four cycles.
% Left: the time series of r, b1, b2 and theta for 1TS and MTS at several
% reserve sizes. Right: the steady-state releasable pool and synaptic output
% as functions of a constant rate (closed form). For the sine input the last
% cycle of each condition is overlaid as a loop, showing how far it departs
% from its own steady-state curve and in which direction it goes round.
% With show_product = true (default) the current, unnormalized implementation
% is drawn as a grey dashed reference: two independent variables (2 s, 0.25 s)
% and (4 s, 0.5 s), D = b1*b2, theta = r*D. It has no reserve, so it is absent
% from the b2 panel; elsewhere "D" means b1 for the chain and b1*b2 for it.
% Each run makes a NEW figure, so running both inputs leaves both open.
%
% Caveat the figure makes visible: bringing in the reserve adds a slower
% timescale AND deepens steady-state depression, by an amount set by n2. The
% two conditions are nested but not matched in strength: the steady-state
% output differs between them, and its high-rate ceiling falls from
% tau_rel/tau_rec1 (1TS) to 1/(tau_rec1/tau_rel + tau_rec2/(n2*tau_rel)).
%
% Integration: ode45, RelTol 1e-10, AbsTol 1e-12, output on the 400 Hz grid.
% Steps: one call per hold, so the rate is constant within every call and the
% discontinuities never sit inside a step. Sine: one call over the whole span
% with r a smooth function of time.
%
% Standalone, base MATLAB only. Leaves its figure open.

%% Parameters
input_type   = 'steps';      % 'steps' | 'sine'
n2_list      = [1 3 10];               % MTS reserve sizes; any number of entries
fs           = 400;                    % output sample rate (Hz)
tau_rec1     = 2;                      % releasable pool refill (s)
tau_rec2     = 4;                      % reserve refill (s)
tau_rel      = 0.25;                   % release (s)
% steps input
r_levels     = [0 0.25 0 0.5 0 1.0 0]; % prescribed rate in each hold
T_hold       = 5;                      % hold duration (s)
% sine input
r_max        = 0.5;                    % peak rate
f_r          = 0.1;                    % rate modulation frequency (Hz)
T_sine       = 40;                     % duration (s): four cycles at 0.1 Hz
% optional reference
show_product = true;                   % also draw the current, unnormalized two-timescale product
tau_rec_slow = 4;                      % (product only) second variable
tau_rel_slow = 0.5;

%% Integrate: y = [b1_1TS; b_slow (product); b1, b2 per n2], all start at 1
n_n2 = numel(n2_list);
n_y  = 2 + 2*n_n2;
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
rhs  = @(y, rr) chain_rhs(y, rr, n2_list, tau_rec1, tau_rec2, tau_rel, tau_rec_slow, tau_rel_slow);

switch input_type
    case 'steps'
        input_label = 'step input';
        n_hold = numel(r_levels);
        n_per  = round(T_hold*fs);
        t = (0:n_hold*n_per).'/fs;
        Y = ones(numel(t), n_y);
        r = zeros(numel(t), 1);
        r(1) = r_levels(1);
        for h = 1:n_hold
            idx = (h-1)*n_per + (1:n_per+1);
            [~, y_seg] = ode45(@(~, y) rhs(y, r_levels(h)), t(idx), Y(idx(1), :).', opts);
            Y(idx(2:end), :) = y_seg(2:end, :);
            r(idx(2:end)) = r_levels(h);
        end
        r_ss = linspace(0, 1, 401);
    case 'sine'
        input_label = sprintf('sinusoidal input, %g Hz, 0 to %g', f_r, r_max);
        r_fun = @(tt) (r_max/2)*(1 - cos(2*pi*f_r*tt));
        t = (0:round(T_sine*fs)).'/fs;
        [~, Y] = ode45(@(tt, y) rhs(y, r_fun(tt)), t, ones(n_y, 1), opts);
        r = r_fun(t);
        r_ss = linspace(0, r_max, 401);
    otherwise
        error('demo_chain_1TS_vs_MTS:input_type', 'input_type must be ''steps'' or ''sine''.');
end
b1_1TS = Y(:, 1);
D_prod = Y(:, 1) .* Y(:, 2);
B1 = Y(:, 3:2:end);
B2 = Y(:, 4:2:end);

%% Steady states (closed form)
ss_1TS = @(rr) 1 ./ (1 + rr*tau_rec1/tau_rel);
ss_MTS = @(rr, n2) chain_b1_ss(rr*tau_rec1/tau_rel, rr*tau_rec2/(n2*tau_rel));
b1ss_1TS = ss_1TS(r_ss);
b1ss_MTS = zeros(numel(r_ss), n_n2);
for v = 1:n_n2
    b1ss_MTS(:, v) = ss_MTS(r_ss(:), n2_list(v));
end
ss_prod = @(rr) 1 ./ ((1 + rr*tau_rec1/tau_rel) .* (1 + rr*tau_rec_slow/tau_rel_slow));
b1ss_prod = ss_prod(r_ss);     % = 1/(1 + 8r)^2 here, both pairs having tau_rec/tau_rel = 8

%% Printed tables
% Rows: 1TS, MTS per n2, and (if show_product) the current product, whose
% "D" is b1*b2 of two independent variables and which has no reserve.
names  = [{'1TS (b2 = 1)'}, compose('MTS n2 = %g', n2_list(:)).'];
D_all  = [b1_1TS, B1];
B2_all = [ones(size(t)), B2];
ss_row = @(rr, v) chain_or_1TS_ss(rr, v, n2_list, ss_1TS, ss_MTS);
if show_product
    names{end+1}  = 'current MTS (product)';
    D_all(:, end+1)  = D_prod;
    B2_all(:, end+1) = NaN;                  % no reserve
    ss_row = @(rr, v) prod_or_chain_ss(rr, v, n_n2, ss_prod, ss_row);
end
n_rows = numel(names);
fprintf('\n=== Chained-pool STD, 1TS vs MTS: %s ===\n', input_label);
is_sine = input_is(input_type, 'sine');

if ~is_sine
    nz = find(r_levels > 0);
    r_test = r_levels(nz);
    fprintf('\nSteady state at the tested rates: D_ss (= b1_ss for the chain), theta_ss, and D_ss/D_ss(1TS)\n');
    fprintf('%-24s', '');
    for rr = r_test, fprintf('|  r = %-4.2g: D     theta  ratio  ', rr); end
    fprintf('\n');
    for v = 1:n_rows
        fprintf('%-24s', names{v});
        for rr = r_test
            b = ss_row(rr, v);
            fprintf('|        %.4f  %.4f  %.3f  ', b, rr*b, b/ss_1TS(rr));
        end
        fprintf('\n');
    end

    r_big = 1e4;
    fprintf('\nHigh-rate ceiling of theta_ss: analytic vs r*D_ss at r = %g\n', r_big);
    for v = 1:n_rows
        if v == 1
            ceil_a = tau_rel/tau_rec1;
        elseif v <= n_n2 + 1
            ceil_a = 1/(tau_rec1/tau_rel + tau_rec2/(n2_list(v-1)*tau_rel));
        else
            ceil_a = 0;                      % product: theta_ss ~ 1/(64 r) -> 0
        end
        val = r_big*ss_row(r_big, v);
        if ceil_a > 0
            fprintf('  %-24s analytic %.5f   at r = %g: %.5f   (rel. diff %.1e)\n', ...
                names{v}, ceil_a, r_big, val, (val - ceil_a)/ceil_a);
        else
            fprintf('  %-24s no ceiling, theta_ss -> 0   at r = %g: %.2e\n', names{v}, r_big, val);
        end
    end

    end_idx = (1:n_hold)*n_per + 1;
    h_last = find(r_levels == max(r_levels), 1, 'last');   % the r = 1 hold
    fprintf('\nRecovery in the %g s of silence after r = %g: (D_end - D_start)/(1 - D_start)\n', ...
        T_hold, r_levels(h_last));
    for v = 1:n_rows
        d0 = D_all(end_idx(h_last), v); d1 = D_all(end_idx(h_last+1), v);
        fprintf('  %-24s %.3f   (D from %.3f to %.3f)\n', names{v}, (d1 - d0)/(1 - d0), d0, d1);
    end
else
    n_cyc = round(fs/f_r);                   % samples per cycle
    last  = numel(t) - n_cyc + (0:n_cyc);    % last full cycle, both ends included
    prev  = last - n_cyc;                    % the cycle before it
    [~, i_rpk] = max(r(last));
    fprintf('\nLast cycle (t = %g to %g s). D = b1 for the chain. Lead = t(r peak) - t(theta peak); positive = theta peaks first.\n', ...
        t(last(1)), t(last(end)));
    fprintf('%-24s| D min    max    mean  | theta mean  peak  | lead (s)  (deg) | b2 min  max   | max|dD| vs prev cycle\n', '');
    for v = 1:n_rows
        dc  = D_all(last, v);
        thc = r(last).*dc;
        [th_pk, i_tpk] = max(thc);
        lead_s = t(last(i_rpk)) - t(last(i_tpk));
        b2c = B2_all(last, v);
        dcyc = max(abs(dc - D_all(prev, v)));
        fprintf('%-24s| %.4f %.4f %.4f | %.4f  %.4f    | %+.3f  %+6.1f  | %5.3f  %5.3f | %.1e\n', ...
            names{v}, min(dc), max(dc), mean(dc), mean(thc), th_pk, ...
            lead_s, 360*f_r*lead_s, min(b2c), max(b2c), dcyc);
    end
end
fprintf('\n');

%% Colours
col_1TS  = [0.72 0.50 0.00];               % darkened Okabe-Ito orange
oi_blue  = [0.00 0.45 0.70];               % Okabe-Ito blue
blue_lo  = [0.60 0.80 0.93];               % light end (smallest n2)
blue_hi  = [0.00 0.22 0.40];               % dark end (largest n2)
w = linspace(0, 1, n_n2).'; if n_n2 == 1, w = 0.5; end
cols = zeros(n_n2, 3);
for v = 1:n_n2
    if w(v) <= 0.5
        a = w(v)/0.5;         cols(v, :) = (1 - a)*blue_lo + a*oi_blue;
    else
        a = (w(v) - 0.5)/0.5; cols(v, :) = (1 - a)*oi_blue + a*blue_hi;
    end
end
col_prod  = [0.55 0.55 0.55];
col_rmark = [0.80 0.80 0.80];
lw_MTS = 1.6; lw_1TS = 2;

%% Figure (always a new one)
fig_chain_ts = figure('Name', ['Chained-pool STD: 1TS vs MTS, ' input_label], 'Color', 'w', ...
    'Position', [100 100 1300 750]);
tl = tiledlayout(fig_chain_ts, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, ['Chained-pool depression, 1TS vs MTS: ' input_label]);
axT = gobjects(4, 1);

% 1. rate
axT(1) = nexttile(tl, 1, [1 2]);
plot(axT(1), t, r, 'k', 'LineWidth', 1.5);
ylim(axT(1), [0 1.1*max(r)]);
ylabel(axT(1), 'r');
title(axT(1), ['Prescribed firing rate r(t), ' input_label]);

% 2-4. releasable pool, reserve, synaptic output
ts_data = {B1, B2, r.*B1};
ts_1TS  = {b1_1TS, ones(size(t)), r.*b1_1TS};
ts_prod = {D_prod, [], r.*D_prod};
ts_ylab = {'D', 'b_2', '\theta = r D'};
ts_ttl  = {'Depression factor D(t)  (= b_1 for the chain)', ...
           'Reserve pool b_2(t)  (1TS: held at 1; the product has no reserve)', ...
           'Synaptic output \theta(t) = r(t) D(t)'};
for k = 1:3
    axT(k+1) = nexttile(tl, 3*k + 1, [1 2]);
    hold(axT(k+1), 'on');
    hP = gobjects(0);
    if show_product && ~isempty(ts_prod{k})  % underneath, as a reference
        hP = plot(axT(k+1), t, ts_prod{k}, '--', 'Color', col_prod, 'LineWidth', 1.25);
    end
    hM = gobjects(n_n2, 1);
    for v = 1:n_n2
        hM(v) = plot(axT(k+1), t, ts_data{k}(:, v), 'Color', cols(v, :), 'LineWidth', lw_MTS);
    end
    h1 = plot(axT(k+1), t, ts_1TS{k}, 'Color', col_1TS, 'LineWidth', lw_1TS);  % on top
    hold(axT(k+1), 'off');
    ylabel(axT(k+1), ts_ylab{k});
    title(axT(k+1), ts_ttl{k});
    if k < 3, ylim(axT(k+1), [0 1]); else, ylim(axT(k+1), [0 inf]); end
    if k == 1   % keep these handles for the one shared legend
        leg_h = [h1; hM; hP(:)];
    end
end
for k = 1:4
    grid(axT(k), 'on'); box(axT(k), 'off');
    if k < 4, axT(k).XTickLabel = []; end
end
xlabel(axT(4), 'Time (s)');
linkaxes(axT, 'x');
xlim(axT(4), [t(1) t(end)]);

% 5-6. steady state vs rate (+ last-cycle loops for the sine input)
axS = gobjects(2, 1);
ss_data = {b1ss_MTS, r_ss(:).*b1ss_MTS};
ss_1    = {b1ss_1TS, r_ss.*b1ss_1TS};
ss_p    = {b1ss_prod, r_ss.*b1ss_prod};
ss_ylab = {'D_{ss}', '\theta_{ss} = r D_{ss}'};
ss_ttl  = {'Steady-state depression factor', 'Steady-state synaptic output'};
if is_sine
    ss_ttl = strcat(ss_ttl, ' + last cycle');
    loop_y = {D_all(last, :), r(last).*D_all(last, :)};
    i_rise = round(n_cyc/4) + 1;             % quarter cycle in: r = r_max/2, rising
    % Loop drawing order: product (if shown) first, then MTS, then 1TS on top.
    % Columns of D_all: 1 = 1TS, 2..n_n2+1 = MTS, n_n2+2 = product.
    cc = [cols; col_1TS];
    order = [2:n_n2+1, 1];
    if show_product
        cc = [col_prod; cc];
        order = [n_n2+2, order];
    end
end
tiles = [3 9];
for k = 1:2
    axS(k) = nexttile(tl, tiles(k), [2 1]);
    hold(axS(k), 'on');
    if ~is_sine
        xline(axS(k), r_test, '-', 'Color', col_rmark, 'LineWidth', 1);
    end
    if show_product                          % underneath, as a reference
        plot(axS(k), r_ss, ss_p{k}, '--', 'Color', col_prod, 'LineWidth', 1.25);
    end
    for v = 1:n_n2
        plot(axS(k), r_ss, ss_data{k}(:, v), 'Color', cols(v, :), 'LineWidth', lw_MTS);
    end
    plot(axS(k), r_ss, ss_1{k}, 'Color', col_1TS, 'LineWidth', lw_1TS);
    if is_sine
        for j = 1:numel(order)
            v = order(j);
            plot(axS(k), r(last), loop_y{k}(:, v), ':', 'Color', cc(j, :), 'LineWidth', 1.1);
            plot(axS(k), r(last(i_rise)), loop_y{k}(i_rise, v), 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', cc(j, :), 'MarkerEdgeColor', 'none');
        end
    end
    hold(axS(k), 'off');
    ylabel(axS(k), ss_ylab{k});
    title(axS(k), ss_ttl{k});
    grid(axS(k), 'on'); box(axS(k), 'off');
end
ylim(axS(1), [0 1]);
ylim(axS(2), [0 inf]);
axS(1).XTickLabel = [];
if is_sine
    xlabel(axS(2), {'Rate r', 'solid: steady state; dotted: last cycle;', 'dot: rising phase at r = r_{max}/2'});
else
    xlabel(axS(2), {'Constant rate r', 'grey lines: tested rates'});
end
linkaxes(axS, 'x');
xlim(axS(2), [r_ss(1) r_ss(end)]);

% One legend for the whole figure
leg_txt = [{'1TS (reserve held full)'}, compose('MTS, n_2 = %g', n2_list(:)).'];
if show_product, leg_txt{end+1} = 'current MTS (product, unnormalized)'; end
lgd = legend(axT(2), leg_h, leg_txt);
lgd.Layout.Tile = 'east';
title(lgd, {'Condition', '(MTS: light = small,', 'dark = large reserve)'});

%% Local functions
function tf = input_is(input_type, name)
% A function rather than an inline strcmp so the Code Analyzer does not fold
% the parameter-block constant and flag the other input's branches as unreachable.
tf = strcmp(input_type, name);
end

function D = chain_or_1TS_ss(rr, v, n2_list, ss_1TS, ss_MTS)
% Steady-state D for table row v: 1 = 1TS, 2.. = MTS at n2_list(v-1).
if v == 1
    D = ss_1TS(rr);
else
    D = ss_MTS(rr, n2_list(v-1));
end
end

function D = prod_or_chain_ss(rr, v, n_n2, ss_prod, ss_chain)
% Adds the product row (v = n_n2 + 2) to the chain/1TS rows.
if v == n_n2 + 2
    D = ss_prod(rr);
else
    D = ss_chain(rr, v);
end
end

function dy = chain_rhs(y, r, n2_list, tau_rec1, tau_rec2, tau_rel, tau_rec_slow, tau_rel_slow)
% y = [b1_1TS; b_slow; b1, b2 per n2]
dy = zeros(size(y));
dy(1) = (1 - y(1))/tau_rec1     - y(1)*r/tau_rel;      % 1TS: b2 = 1
dy(2) = (1 - y(2))/tau_rec_slow - y(2)*r/tau_rel_slow; % product's slow variable
for j = 1:numel(n2_list)
    i1 = 1 + 2*j; i2 = i1 + 1;
    refill = y(i2)*(1 - y(i1))/tau_rec1;
    dy(i1) = refill - y(i1)*r/tau_rel;
    dy(i2) = (1 - y(i2))/tau_rec2 - refill/n2_list(j);
end
end

function b1 = chain_b1_ss(k1, k2)
% Root with b2 >= 0 of k2*b1^2 - (1 + k1 + k2)*b1 + 1 = 0, in the
% cancellation-free form 2/(s + sqrt(s^2 - 4 k2)) (equals 1/(1 + k1) at k2 = 0).
s = 1 + k1 + k2;
b1 = 2 ./ (s + sqrt(s.^2 - 4*k2));
end
