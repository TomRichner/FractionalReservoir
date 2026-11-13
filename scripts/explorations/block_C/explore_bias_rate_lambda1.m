function out = explore_bias_rate_lambda1(root)
% EXPLORE_BIAS_RATE_LAMBDA1 How bias, mean population rate and lambda_1 relate on the main gain sweep (block C).
%
%   out = EXPLORE_BIAS_RATE_LAMBDA1()       % the result and figs/ of this checkout
%   out = EXPLORE_BIAS_RATE_LAMBDA1(root)   % ... of the checkout at root (read-only on its data/)
%
% Tom (2026-11-13): "It would be nice to see if bias and mean population
% spike rate are correlated or not, and whether bias and lambda_1 are
% correlated." The reason it is asked: by section 10.7 of
% docs/EquationsParametersDocs/chained_MTS/Equations_chained_MTS.md the bias
% enters the no-adaptation Jacobian only through the slopes phi'(h_i), and
% the recurrent input spreads the operating points h_i more widely than the
% nonlinearity's sloped window, so sliding them by the bias should change
% the rate and the share of neurons on the slope only a little, and
% lambda_1 hardly at all, while the gain moves lambda_1 directly.
%
% Six figures, each fig_mc_sweep_scatter over every run of the main sweep
% (11 gains x 7 biases x 2 seeds), four conditions:
%   1. bias against rate_mean, colored by gain
%   2. bias against lambda_1 (input on), colored by gain
%   3. rate_mean against lambda_1, colored by gain
%   4. gain against frac_on_slope, colored by bias: a first look at the
%      section's test 1. The explanation predicts a share on the slope near
%      1/(g R_W)^2 at the no-adaptation edge, roughly independent of bias at
%      a fixed gain; that curve is drawn on the no-adaptation panel, with
%      R_W the bulk radius of W at gain 1, read off the preset's model.
%   5. rate_mean against lambda_1, colored by gain, marker area by bias
%      (filled for bias >= 0, a ring for bias < 0; Tom, 2026-11-13)
%   6. memory capacity against lambda_1, the same encoding
% In 5 and 6 a clipped lambda_1 is an open square at the axis edge.
% Bias has seven levels, so in figures 1 and 2 the points stand in seven
% columns and Spearman's rho is computed with many tied ranks (averaged).
%
% Printed: each figure's statistics, and for no adaptation, at each bias,
% the two gains that bracket the sign change of the median lambda_1 with
% the mean share on the slope there against 1/(g R_W)^2.
%
% Output: <root>/figs/explore_C/bias_rate_lambda1/Fig_*.{png,svg,fig}
%
% See also: fig_mc_sweep_scatter, mc_sweep_quantity, explore_mc_vs_lambda1_scatter

if nargin < 1 || isempty(root); root = fileparts(which('setup_paths')); end
data_file = fullfile(root, 'data', 'mc_gain_sweep_main', 'mc_gain_sweep', 'mc_gain_sweep_results.mat');
out_dir = fullfile(root, 'figs', 'explore_C', 'bias_rate_lambda1');
conds = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
D = read_mc_gain_sweep(data_file);

% R_W: the bulk radius of W at gain 1 (the model's R scales with the gain).
d = srnn_param_preset(D.preset_name);
args = [struct2namevalue(d), D.settings.model_args, {'level_of_chaos', 1}];
R_W = SRNNCellTypePairs(args{:}).R;

g_curve = logspace(log10(1 / R_W), log10(max(D.gains)), 200);
edge_curve = struct('x', g_curve, 'y', 1 ./ (g_curve * R_W).^2, ...
    'label', sprintf('1/(g R_W)^2, R_W = %.2f', R_W), 'conditions', {{'no_adaptation'}});

common = {'data', D, 'conditions', conds, 'out_dir', out_dir, 'visible', false};
% Each statistics box placed where its panels have fewest points.
specs = {
    {'x', 'bias',      'y', 'rate_mean',     'color_by', 'gain', 'fig_tag', 'Fig_bias_vs_rate_mean', ...
     'stats_position', [0.56 0.98]}
    {'x', 'bias',      'y', 'lle_input',     'color_by', 'gain', 'fig_tag', 'Fig_bias_vs_lambda1', ...
     'stats_position', [0.56 0.27]}
    {'x', 'lle_input', 'y', 'rate_mean',     'color_by', 'gain', 'fig_tag', 'Fig_rate_mean_vs_lambda1', ...
     'stats_position', [0.56 0.99]}
    {'x', 'gain',      'y', 'frac_on_slope', 'color_by', 'bias', 'fig_tag', 'Fig_gain_vs_frac_on_slope', ...
     'overlay', edge_curve, 'stats_position', [0.50 0.22]}
    {'x', 'lle_input', 'y', 'rate_mean',     'color_by', 'gain', 'size_by', 'bias', ...
     'fig_tag', 'Fig_rate_mean_vs_lambda1_size_bias', 'stats_position', [0.56 0.99]}
    {'x', 'lle_input', 'y', 'mc_total',      'color_by', 'gain', 'size_by', 'bias', ...
     'fig_tag', 'Fig_MC_vs_lambda1_size_bias', 'stats_position', [0.56 0.99]}};
out = struct('stats', {cell(numel(specs), 1)}, 'files', {{}}, 'R_W', R_W, 'edge', []);
for k = 1:numel(specs)
    o = fig_mc_sweep_scatter(common{:}, specs{k}{:});
    close(o.figs);
    out.stats{k} = o.stats;
    out.files = [out.files, o.files];
    fprintf('\n%s\n', specs{k}{find(strcmp(specs{k}, 'fig_tag')) + 1});
    disp(o.stats);
end

% The no-adaptation edge: per bias, the gains bracketing the sign change
% of the median lambda_1, and the mean share on the slope at each.
L  = mc_gain_grid(D, 'no_adaptation', 'lle_input');
on = 1 - mc_gain_grid(D, 'no_adaptation', 'frac_pinned', @mean) ...
       - mc_gain_grid(D, 'no_adaptation', 'frac_silent', @mean);
rows = cell(0, 7);
for b = 1:numel(D.biases)
    i = find(L(b, 1:end - 1) < 0 & L(b, 2:end) > 0, 1);
    if isempty(i); continue; end
    g = D.gains([i, i + 1]);
    rows(end + 1, :) = {D.biases(b), g(1), g(2), on(b, i), on(b, i + 1), ...
        1 / (g(1) * R_W)^2, 1 / (g(2) * R_W)^2}; %#ok<AGROW>
end
out.edge = cell2table(rows, 'VariableNames', {'bias', 'gain_below', 'gain_above', ...
    'on_slope_below', 'on_slope_above', 'pred_below', 'pred_above'});
fprintf('\nNo adaptation, R_W = %.3f: share on the slope at the gains bracketing lambda_1 = 0\n', R_W);
disp(out.edge);
out.files = fullfile(out_dir, out.files);
fprintf('%s\n', out.files{endsWith(out.files, '.png')});
end
