% demo_std_variants_step_rate.m
%
% The short-term-depression (STD) subsystem on its own, driven by a
% PRESCRIBED firing rate r(t) (steps of 5 s, back to zero between steps), to
% compare ways of combining depression timescales into one factor D(t). The
% synaptic output is theta = r*D. Six variants:
%   1. 1TS            one variable (tau_rec 2 s, tau_rel 0.25 s), D = b
%   2. MTS product    two variables (2, 0.25) and (4, 0.5), D = b1*b2   (paper now)
%   3. MTS geo. mean  same two variables, D = sqrt(b1*b2)
%   4. MTS arith mean same two variables, D = (b1 + b2)/2
%   5. same-count     two IDENTICAL variables (2, 0.25), D = b^2
%   6. chained pools  releasable pool b1 refilled from a slower reserve b2
%
% Questions it helps answer:
%   - Which combinations keep the single-timescale steady state? Both pairs
%     have tau_rel/tau_rec = 0.125, so every b_m rests at 1/(1 + 8r). Variants
%     1, 3 and 4 therefore share that steady state; the product (2) squares
%     it, and so does the same-count control (5), which has the product's
%     depth but only one timescale.
%   - How do onset transients and recovery differ between variants with the
%     same steady state? 1, 3 and 4 agree at rest and in steady state but
%     the slow variable makes 3 and 4 approach both more slowly; 2 and 5 agree
%     in steady state, and the difference between them is purely the slow
%     timescale.
%   - What does the chain do differently? Depletion of the reserve slows the
%     refill of the releasable pool, so the chain is deeper than variant 1
%     and its steady state is not of the 1/(1 + k r) form.
%
% The chain (variant 6) is ONE illustrative chain, chosen because it reduces
% exactly to variant 1 as the reserve size n2 -> infinity (b2 stays at 1); it
% is not a fitted vesicle-pool model:
%   db1/dt = b2*(1 - b1)/tau_rec1 - b1*r/tau_rel
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) * b2*(1 - b1)/tau_rec1
% The second term of db2/dt is the flux that refills pool 1, so resources
% are conserved.
%
% Note: variants 1, 2, 3, 4 and 5 all use a variable with constants (2, 0.25)
% obeying the same ODE from the same initial condition, so it is ONE state
% here (b_fast), shared by all five; variants 2-4 add b_slow (4, 0.5).
%
% Integration: ode45 segment by segment (one segment per 5 s hold, so the
% rate is constant within every call and the discontinuity never sits inside
% a step), RelTol 1e-10, AbsTol 1e-12, output on the 400 Hz grid.
%
% Standalone, base MATLAB only. Leaves its figure open.

%% Parameters
r_levels   = [0 0.25 0 0.5 0 1.0 0];   % prescribed rate in each hold
T_hold     = 5;                        % hold duration (s)
fs         = 400;                      % output sample rate (Hz)
tau_rec_MTS = [2 4];                   % two-timescale STD recovery (s)
tau_rel_MTS = [0.25 0.5];              % two-timescale STD release (s)
tau_rec_1TS = tau_rec_MTS(1);          % single timescale = first of the pair
tau_rel_1TS = tau_rel_MTS(1);
ch_tau_rec1 = 2;                       % chain: releasable pool refill (s)
ch_tau_rec2 = 4;                       % chain: reserve refill (s)
ch_tau_rel  = 0.25;                    % chain: release (s)
ch_n2       = 1;                       % chain: reserve size / releasable size

%% Integrate: state y = [b_fast; b_slow; chain b1; chain b2], all start at 1
rhs = @(~, y, r) [ ...
    (1 - y(1))/tau_rec_1TS    - y(1)*r/tau_rel_1TS; ...
    (1 - y(2))/tau_rec_MTS(2) - y(2)*r/tau_rel_MTS(2); ...
    y(4)*(1 - y(3))/ch_tau_rec1 - y(3)*r/ch_tau_rel; ...
    (1 - y(4))/ch_tau_rec2 - (1/ch_n2)*y(4)*(1 - y(3))/ch_tau_rec1];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

n_hold = numel(r_levels);
n_per  = round(T_hold*fs);
t = (0:n_hold*n_per).'/fs;
Y = zeros(numel(t), 4);
Y(1, :) = 1;
r = zeros(numel(t), 1);
r(1) = r_levels(1);
for h = 1:n_hold
    idx = (h-1)*n_per + (1:n_per+1);           % includes both segment ends
    [~, y_seg] = ode45(@(tt, y) rhs(tt, y, r_levels(h)), t(idx), Y(idx(1), :).', opts);
    Y(idx(2:end), :) = y_seg(2:end, :);
    r(idx(2:end)) = r_levels(h);               % r at a boundary = the hold it ends
end
b_fast = Y(:, 1); b_slow = Y(:, 2); ch_b1 = Y(:, 3); ch_b2 = Y(:, 4);

%% The six variants
names = {'1  1TS', '2  MTS product', '3  MTS geometric mean', ...
         '4  MTS arithmetic mean', '5  same-count product', '6  chained pools'};
D = [b_fast, b_fast.*b_slow, sqrt(b_fast.*b_slow), (b_fast + b_slow)/2, ...
     b_fast.^2, ch_b1];
theta = r .* D;
n_var = size(D, 2);

%% Steady states
ss_b = @(r_, rec, rel) 1 ./ (1 + r_*rec/rel);
ss_analytic = @(r_) [ ...
    ss_b(r_, tau_rec_1TS, tau_rel_1TS), ...
    ss_b(r_, tau_rec_MTS(1), tau_rel_MTS(1)) * ss_b(r_, tau_rec_MTS(2), tau_rel_MTS(2)), ...
    sqrt(ss_b(r_, tau_rec_MTS(1), tau_rel_MTS(1)) * ss_b(r_, tau_rec_MTS(2), tau_rel_MTS(2))), ...
    (ss_b(r_, tau_rec_MTS(1), tau_rel_MTS(1)) + ss_b(r_, tau_rec_MTS(2), tau_rel_MTS(2)))/2, ...
    ss_b(r_, tau_rec_1TS, tau_rel_1TS)^2];
% Chain: numerically, with fzero. Setting both derivatives to zero gives
% b2 = 1 - k2*b1 and b2*(1 - b1) = k1*b1, with k1 = tau_rec1*r/tau_rel and
% k2 = tau_rec2*r/(n2*tau_rel); solve the second for b1 in (0, min(1, 1/k2)).
ss_chain = @(r_) fzero(@(b1) (1 - (ch_tau_rec2*r_/(ch_n2*ch_tau_rel))*b1).*(1 - b1) ...
    - (ch_tau_rec1*r_/ch_tau_rel)*b1, ...
    [0, min(1, ch_n2*ch_tau_rel/(ch_tau_rec2*r_))]);

%% Printed table
end_idx = (1:n_hold)*n_per + 1;                % sample at the end of each hold
fprintf('\nD at the end of each %g-s hold vs its steady state  (end / steady, %% off)\n', T_hold);
nz = find(r_levels > 0);
fprintf('%-24s', 'variant');
for h = nz
    fprintf('|     r = %-4.2g               ', r_levels(h));
end
fprintf('\n');
for v = 1:n_var
    fprintf('%-24s', names{v});
    for h = nz
        if v <= 5
            s = ss_analytic(r_levels(h)); s = s(v); tag = ' ';
        else
            s = ss_chain(r_levels(h)); tag = '*';
        end
        d_end = D(end_idx(h), v);
        fprintf('| %.4f / %.4f%s %+6.1f%%  ', d_end, s, tag, 100*(d_end - s)/s);
    end
    fprintf('\n');
end
fprintf('  * chain steady state solved numerically (fzero); others closed form.\n');

fprintf('\nRecovery during each zero-rate hold: fraction of the deficit recovered\n');
fprintf('  (D_end - D_start)/(1 - D_start), D_start = value when the rate returns to 0\n');
zh = find(r_levels == 0 & (1:n_hold) > 1);
fprintf('%-24s', 'variant');
for h = zh
    fprintf('| after r = %-4.2g    ', r_levels(h-1));
end
fprintf('\n');
for v = 1:n_var
    fprintf('%-24s', names{v});
    for h = zh
        d0 = D(end_idx(h-1), v); d1 = D(end_idx(h), v);
        fprintf('| %.3f (D=%.3f)  ', (d1 - d0)/(1 - d0), d1);
    end
    fprintf('\n');
end
fprintf('\n');

%% Figure
cols = [0.00 0.00 0.00;     % 1 1TS
        0.84 0.37 0.00;     % 2 MTS product
        0.00 0.45 0.70;     % 3 geometric mean
        0.00 0.62 0.45;     % 4 arithmetic mean
        0.80 0.47 0.65;     % 5 same-count
        0.90 0.62 0.00];    % 6 chain
styles = {'-', '-', '-', '--', ':', '-.'};

fig_std = figure('Name', 'STD variants under a step rate', 'Color', 'w');
tl = tiledlayout(fig_std, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(4, 1);

ax(1) = nexttile(tl);
plot(ax(1), t, r, 'k', 'LineWidth', 1.5);
ylim(ax(1), [0 1.1*max(r_levels)]);
ylabel(ax(1), 'r');
legend(ax(1), {'prescribed r(t)'}, 'Location', 'eastoutside');
title(ax(1), 'Prescribed firing rate');

ax(2) = nexttile(tl);
hold(ax(2), 'on');
plot(ax(2), t, b_fast, 'Color', cols(1, :), 'LineWidth', 1.5);
plot(ax(2), t, b_slow, 'Color', cols(3, :), 'LineWidth', 1.5);
plot(ax(2), t, ch_b1,  'Color', cols(6, :), 'LineWidth', 1.5);
plot(ax(2), t, ch_b2,  'Color', cols(6, :), 'LineWidth', 1.5, 'LineStyle', ':');
hold(ax(2), 'off');
ylim(ax(2), [0 1]);
ylabel(ax(2), 'b');
legend(ax(2), {sprintf('b (%g s, %g s): 1TS b = MTS b_1', tau_rec_1TS, tau_rel_1TS), ...
               sprintf('MTS b_2 (%g s, %g s)', tau_rec_MTS(2), tau_rel_MTS(2)), ...
               'chain b_1 (releasable)', sprintf('chain b_2 (reserve, n_2 = %g)', ch_n2)}, ...
       'Location', 'eastoutside');
title(ax(2), 'Depression state variables');

ax(3) = nexttile(tl);
hold(ax(3), 'on');
for v = 1:n_var
    plot(ax(3), t, D(:, v), 'Color', cols(v, :), 'LineStyle', styles{v}, 'LineWidth', 1.5);
end
hold(ax(3), 'off');
ylim(ax(3), [0 1]);
ylabel(ax(3), 'D');
legend(ax(3), names, 'Location', 'eastoutside');
title(ax(3), 'Combined depression factor D(t)');

ax(4) = nexttile(tl);
hold(ax(4), 'on');
for v = 1:n_var
    plot(ax(4), t, theta(:, v), 'Color', cols(v, :), 'LineStyle', styles{v}, 'LineWidth', 1.5);
end
hold(ax(4), 'off');
ylim(ax(4), [0 inf]);
ylabel(ax(4), '\theta = r D');
legend(ax(4), names, 'Location', 'eastoutside');
title(ax(4), 'Synaptic output \theta(t) = r(t) D(t)');

for k = 1:4
    grid(ax(k), 'on');
    box(ax(k), 'off');
    if k < 4
        ax(k).XTickLabel = [];
    end
end
xlabel(ax(4), 'Time (s)');
linkaxes(ax, 'x');
xlim(ax(4), [t(1) t(end)]);
