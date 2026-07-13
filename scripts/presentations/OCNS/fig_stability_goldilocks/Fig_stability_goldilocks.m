% Fig_stability_goldilocks.m - Three example reservoirs across the edge of chaos
%
% "Goldilocks" stability figure for the OCNS presentation. Three bare
% SRNNModel2 reservoirs with NO adaptation (no SFA, no STD) are shown at three
% points relative to the edge of chaos, dialed in via the connectivity scale
% (level_of_chaos):
%
%   Stable         largest Lyapunov exponent (LLE) ~ -0.25
%   Edge of chaos  LLE ~ +0.05
%   Chaotic        LLE ~ +0.25
%
% The three networks share the same fraction-excitatory (f = 0.5) and the same
% network RNG seed, so only level_of_chaos (a linear multiplier on the spectral
% radius R) differs between them -- isolating the transition to chaos as the one
% knob that is varied. The panels overlay a subset of excitatory firing-rate
% traces r_i(t): the stable network settles, the edge network is marginal, and
% the chaotic network sustains irregular activity.
%
% The (level_of_chaos, f, seed) triples below were found by the search block
% guarded under DO_SEARCH (set it true to re-derive), then hardcoded so the
% figure regenerates deterministically and fast. Achieved LLEs (default
% benettin method): stable -0.225, edge +0.063, chaotic +0.281.

close all; clear; clc;
setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% Shared network settings (no adaptation, fixed structure)
f_shared    = 0.5;        % fraction excitatory
seed_shared = [1 2];      % rng_seeds = [network, stimulus]

%% Hardcoded winners (found via the DO_SEARCH block below)
configs = struct( ...
    'label',          {'Stable',       'Edge of chaos', 'Chaotic'}, ...
    'target_LLE',     {-0.25,           0.05,            0.25}, ...
    'level_of_chaos', { 1.50,           1.55,            1.85});
n_net = numel(configs);

%% -------- Optional tuning search (documents how the winners were found) --------
% Sweeps level_of_chaos + f + seed, runs each no-adaptation SRNNModel2, and
% reports the combo closest to each target LLE. Off by default (slow: one
% benettin run ~12 s). The hardcoded configs above are the result of this search
% at f = 0.5, seed = [1 2].
DO_SEARCH = false;
if DO_SEARCH
    targets    = [-0.25, 0.05, 0.25];
    tol_frac   = 0.15;                 % accept |LLE - target| <= 15% of |target|
    f_grid     = [0.3, 0.5, 0.8];
    loc_grid   = 1.40:0.05:1.90;
    seed_grid  = 1:3;
    best = repmat(struct('err', inf, 'loc', NaN, 'f', NaN, 'seed', NaN, 'LLE', NaN), 1, numel(targets));
    for s = seed_grid
        for ff = f_grid
            for loc = loc_grid
                ms = SRNNModel2('n_a_E',0,'n_a_I',0,'n_b_E',0,'n_b_I',0, ...
                                'level_of_chaos',loc,'f',ff,'rng_seeds',[s 2]);
                ms.build(); ms.run();
                L = ms.lya_results.LLE;
                for ti = 1:numel(targets)
                    e = abs(L - targets(ti));
                    if e < best(ti).err
                        best(ti) = struct('err',e,'loc',loc,'f',ff,'seed',s,'LLE',L);
                    end
                end
            end
        end
    end
    fprintf('\n==== search winners ====\n');
    for ti = 1:numel(targets)
        b = best(ti);
        ok = b.err <= tol_frac*abs(targets(ti));
        fprintf('target %+.2f: loc=%.2f f=%.1f seed=%d  LLE=%+.4f  (err=%.4f%s)\n', ...
            targets(ti), b.loc, b.f, b.seed, b.LLE, b.err, ternary(ok,' OK',' OUT'));
    end
end

%% -------- Build + run each network --------
runs = struct('t', {}, 'rE', {}, 'LLE', {});
for k = 1:n_net
    m = SRNNModel2('n_a_E',0,'n_a_I',0,'n_b_E',0,'n_b_I',0, ...
                   'level_of_chaos', configs(k).level_of_chaos, ...
                   'f', f_shared, 'rng_seeds', seed_shared);
    m.build();
    m.run();
    runs(k).t   = m.plot_data.t;
    runs(k).rE  = m.plot_data.r.E;
    runs(k).LLE = m.lya_results.LLE;
    fprintf('%-14s loc=%.2f  target LLE=%+.2f  achieved LLE=%+.4f\n', ...
        configs(k).label, configs(k).level_of_chaos, configs(k).target_LLE, runs(k).LLE);
end

%% -------- Figure: 1x3 overlaid excitatory firing-rate traces --------
n_traces = 15;             % E neurons drawn per panel
t_disp   = [10, 40];       % display window (skip the initial transient)
cmap     = warm_colormap(n_traces);   % excitatory-style light-orange -> dark-red gradient

fig = figure('Color', 'w', 'Position', [80, 200, 1500, 420]);
tl  = tiledlayout(fig, 1, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:n_net
    ax = nexttile(tl);
    t  = runs(k).t;
    rE = runs(k).rE;
    idx  = round(linspace(1, size(rE, 1), n_traces));
    keep = t >= t_disp(1) & t <= t_disp(2);
    hold(ax, 'on');
    for j = 1:n_traces
        plot(ax, t(keep), rE(idx(j), keep), 'Color', cmap(j, :), 'LineWidth', 1.0);
    end
    hold(ax, 'off');
    xlim(ax, t_disp);
    ylim(ax, [0, 1]);
    yticks(ax, [0, 1]);
    box(ax, 'off');
    xlabel(ax, 'time (s)');
    if k == 1
        ylabel(ax, 'firing rate');
    end
    title(ax, sprintf('%s  (\\lambda_1 = %+.2f)', configs(k).label, runs(k).LLE), ...
        'FontWeight', 'normal', 'FontSize', 14);
end

%% -------- Save --------
save_some_figs_to_folder_2(this_dir, 'Fig_stability_goldilocks', fig.Number, {'fig', 'png', 'svg'});

%% Local helpers
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function cmap = warm_colormap(n)
% WARM_COLORMAP n-row light-orange -> dark-red gradient for excitatory traces.
cmap = [linspace(0.95, 0.45, n)', linspace(0.55, 0.05, n)', linspace(0.25, 0.05, n)'];
end
