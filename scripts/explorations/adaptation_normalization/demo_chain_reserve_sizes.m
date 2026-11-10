% demo_chain_reserve_sizes.m
%
% The chained-pool depression model (variant 6 of
% demo_std_variants_step_rate.m) at several reserve sizes, driven by a
% PRESCRIBED step rate r(t). A releasable pool b1 is refilled from a slower
% reserve b2, and the synapse transmits from b1 only:
%   db1/dt = b2*(1 - b1)/tau_rec1 - b1*r/tau_rel
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) * b2*(1 - b1)/tau_rec1
%   D = b1,  theta = r*D
% n2 = N2/N1 is the reserve's capacity relative to the releasable pool.
% Transfer between pools is one-for-one in vesicles; b1 and b2 are each a
% fraction of their OWN pool's capacity, so the factor 1/n2 only converts the
% refill flux from pool-1 units into pool-2 units. Resources are conserved.
%
% What it shows:
%   - With a small reserve (n2 = 1) the reserve itself depletes, which slows
%     the refill of the releasable pool, so the chain is much deeper than a
%     single pool: at n2 = 1 its steady-state depth lies between the single
%     pool and the current two-timescale product, and its recovery in
%     silence is as slow as the product's or slower. As the reserve grows, its fractional depletion
%     shrinks and the chain moves toward the single-pool curve, which it
%     reaches exactly as n2 -> infinity (b2 stays at 1).
%   - The reserve recovers on tau_rec2 = 4 s and is drained by every refill,
%     so with a small reserve its depletion carries over from one step to the
%     next: each later step starts from a lower reserve, deepens further, and
%     recovers less in the following 5 s of silence.
% References: the single pool D = b (tau_rec 2 s, tau_rel 0.25 s), the
% n2 -> infinity limit; and the current two-timescale product D = b1*b2 with
% (2 s, 0.25 s) and (4 s, 0.5 s).
%
% This is ONE illustrative chain, chosen because it reduces to the single
% pool as n2 -> infinity; it is not a fitted vesicle-pool model.
%
% Integration: ode45 segment by segment (one segment per hold, so the rate is
% constant within every call), RelTol 1e-10, AbsTol 1e-12, output on the
% 400 Hz grid.
%
% Standalone, base MATLAB only. Leaves its figure open.

%% Parameters
n2_list   = [1 10 100];             % reserve size N2/N1; any number of entries
r_levels  = [0 0.25 0 0.5 0 1.0 0]; % prescribed rate in each hold
T_hold    = 5;                      % hold duration (s)
fs        = 400;                    % output sample rate (Hz)
tau_rec1  = 2;                      % releasable pool refill (s)
tau_rec2  = 4;                      % reserve refill (s)
tau_rel   = 0.25;                   % release (s)
tau_rec_slow = 4;                   % product reference, second variable (s)
tau_rel_slow = 0.5;

%% Integrate: y = [b_single; b_slow; b1(n2_1); b2(n2_1); b1(n2_2); ...]
n_n2 = numel(n2_list);
rhs = @(~, y, r) chain_rhs(y, r, n2_list, tau_rec1, tau_rec2, tau_rel, tau_rec_slow, tau_rel_slow);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

n_hold = numel(r_levels);
n_per  = round(T_hold*fs);
t = (0:n_hold*n_per).'/fs;
n_state = 2 + 2*n_n2;
Y = zeros(numel(t), n_state);
Y(1, :) = 1;
r = zeros(numel(t), 1);
r(1) = r_levels(1);
for h = 1:n_hold
    idx = (h-1)*n_per + (1:n_per+1);
    [~, y_seg] = ode45(@(tt, y) rhs(tt, y, r_levels(h)), t(idx), Y(idx(1), :).', opts);
    Y(idx(2:end), :) = y_seg(2:end, :);
    r(idx(2:end)) = r_levels(h);
end
b_single = Y(:, 1);
D_prod   = Y(:, 1) .* Y(:, 2);
B1 = Y(:, 3:2:end);                 % releasable pool, one column per n2
B2 = Y(:, 4:2:end);                 % reserve, one column per n2

%% Closed-form steady state of the chain
b1_ss = @(r_, n2) chain_b1_ss(r_*tau_rec1/tau_rel, r_*tau_rec2/(n2*tau_rel));

%% Printed tables
end_idx = (1:n_hold)*n_per + 1;
nz = find(r_levels > 0);
row_names = [compose('chain n2 = %g', n2_list(:)).', {'single pool', 'MTS product'}];
D_cols = [B1, b_single, D_prod];

fprintf('\nD at the end of each %g-s hold vs closed-form steady state  (end / steady, %% off)\n', T_hold);
fprintf('%-18s', '');
for h = nz, fprintf('|  r = %-4.2g                  ', r_levels(h)); end
fprintf('\n');
for v = 1:size(D_cols, 2)
    fprintf('%-18s', row_names{v});
    for h = nz
        k1 = r_levels(h)*tau_rec1/tau_rel;
        if v <= n_n2
            s = b1_ss(r_levels(h), n2_list(v));
        elseif v == n_n2 + 1
            s = 1/(1 + k1);
        else
            s = 1/(1 + k1)^2;   % both product variables have tau_rec/tau_rel = 8
        end
        d = D_cols(end_idx(h), v);
        fprintf('| %.4f / %.4f  %+6.1f%%  ', d, s, 100*(d - s)/s);
    end
    fprintf('\n');
end

fprintf('\nReserve b2 at the end of each hold\n');
fprintf('%-18s', '');
for h = 1:n_hold, fprintf('| r=%-4.2g ', r_levels(h)); end
fprintf('\n');
for v = 1:n_n2
    fprintf('%-18s', row_names{v});
    fprintf('| %.4f ', B2(end_idx, v));
    fprintf('\n');
end

fprintf('\nRecovery in each zero-rate hold: (D_end - D_start)/(1 - D_start), D at end in brackets\n');
zh = find(r_levels == 0 & (1:n_hold) > 1);
fprintf('%-18s', '');
for h = zh, fprintf('| after r = %-4.2g   ', r_levels(h-1)); end
fprintf('\n');
for v = 1:size(D_cols, 2)
    fprintf('%-18s', row_names{v});
    for h = zh
        d0 = D_cols(end_idx(h-1), v); d1 = D_cols(end_idx(h), v);
        fprintf('| %.3f (%.3f)    ', (d1 - d0)/(1 - d0), d1);
    end
    fprintf('\n');
end

% Check the closed form against a long integration (n2 = first entry, r = 0.25)
r_chk = 0.25; n2_chk = n2_list(1);
[~, y_long] = ode45(@(~, y) chain_rhs([1; 1; y], r_chk, n2_chk, tau_rec1, tau_rec2, ...
    tau_rel, tau_rec_slow, tau_rel_slow, 'chain_only'), [0 500], [1; 1], opts);
k2_chk = r_chk*tau_rec2/(n2_chk*tau_rel);
b1_cf = b1_ss(r_chk, n2_chk);
fprintf(['\nCheck, n2 = %g, r = %g, 500-s integration: b1 = %.8f, closed form %.8f (diff %.1e);\n' ...
         '        b2 = %.8f, closed form 1 - k2*b1 = %.8f\n\n'], ...
    n2_chk, r_chk, y_long(end, 1), b1_cf, y_long(end, 1) - b1_cf, y_long(end, 2), 1 - k2_chk*b1_cf);

%% Figure
ramp_lo = [0.62 0.79 0.88];         % light blue (smallest n2)
ramp_hi = [0.03 0.19 0.42];         % dark blue (largest n2)
w = linspace(0, 1, n_n2).';
if n_n2 == 1, w = 1; end
cols = (1 - w).*ramp_lo + w.*ramp_hi;
col_single = [0 0 0];
col_prod   = [0.55 0.55 0.55];
lbl_n2 = compose('chain, n_2 = %g', n2_list(:)).';

fig_chain = figure('Name', 'Chained-pool STD: reserve size', 'Color', 'w');
tl = tiledlayout(fig_chain, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(4, 1);

ax(1) = nexttile(tl);
plot(ax(1), t, r, 'k', 'LineWidth', 1.5);
ylim(ax(1), [0 1.1*max(r_levels)]);
ylabel(ax(1), 'r');
legend(ax(1), {'prescribed r(t)'}, 'Location', 'eastoutside');
title(ax(1), 'Prescribed firing rate');

ax(2) = nexttile(tl);
hold(ax(2), 'on');
for v = 1:n_n2
    plot(ax(2), t, B1(:, v), 'Color', cols(v, :), 'LineWidth', 1.75);
end
plot(ax(2), t, b_single, 'Color', col_single, 'LineWidth', 1);
plot(ax(2), t, D_prod, '--', 'Color', col_prod, 'LineWidth', 1.25);
hold(ax(2), 'off');
ylim(ax(2), [0 1]);
ylabel(ax(2), 'b_1 = D');
legend(ax(2), [lbl_n2, {'single pool (n_2 \rightarrow \infty)', 'MTS product b_1 b_2'}], ...
    'Location', 'eastoutside');
title(ax(2), 'Releasable pool b_1(t) = D(t)');

ax(3) = nexttile(tl);
hold(ax(3), 'on');
for v = 1:n_n2
    plot(ax(3), t, B2(:, v), 'Color', cols(v, :), 'LineWidth', 1.75);
end
hold(ax(3), 'off');
ylim(ax(3), [0 1]);
ylabel(ax(3), 'b_2');
legend(ax(3), lbl_n2, 'Location', 'eastoutside');
title(ax(3), 'Reserve pool b_2(t)');

ax(4) = nexttile(tl);
hold(ax(4), 'on');
for v = 1:n_n2
    plot(ax(4), t, r.*B1(:, v), 'Color', cols(v, :), 'LineWidth', 1.75);
end
plot(ax(4), t, r.*b_single, 'Color', col_single, 'LineWidth', 1);
plot(ax(4), t, r.*D_prod, '--', 'Color', col_prod, 'LineWidth', 1.25);
hold(ax(4), 'off');
ylim(ax(4), [0 inf]);
ylabel(ax(4), '\theta = r D');
legend(ax(4), [lbl_n2, {'single pool', 'MTS product'}], 'Location', 'eastoutside');
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

%% Local functions
function dy = chain_rhs(y, r, n2_list, tau_rec1, tau_rec2, tau_rel, tau_rec_slow, tau_rel_slow, mode)
% y = [b_single; b_slow; b1, b2 per n2]. With mode 'chain_only' only the
% chain derivatives are returned (used by the long-time check).
dy = zeros(size(y));
dy(1) = (1 - y(1))/tau_rec1     - y(1)*r/tau_rel;
dy(2) = (1 - y(2))/tau_rec_slow - y(2)*r/tau_rel_slow;
for j = 1:numel(n2_list)
    i1 = 1 + 2*j; i2 = i1 + 1;
    refill = y(i2)*(1 - y(i1))/tau_rec1;
    dy(i1) = refill - y(i1)*r/tau_rel;
    dy(i2) = (1 - y(i2))/tau_rec2 - refill/n2_list(j);
end
if nargin > 8 && strcmp(mode, 'chain_only')
    dy = dy(3:end);
end
end

function b1 = chain_b1_ss(k1, k2)
% Smaller root of k2*b1^2 - (1 + k1 + k2)*b1 + 1 = 0 (the root with b2 >= 0),
% ((s - sqrt(s^2 - 4 k2)) / (2 k2)), written in the algebraically identical
% form 2/(s + sqrt(s^2 - 4 k2)), which has no cancellation at large n2
% (small k2) and gives 1/(1 + k1) at k2 = 0.
s = 1 + k1 + k2;
b1 = 2 / (s + sqrt(s^2 - 4*k2));
end
