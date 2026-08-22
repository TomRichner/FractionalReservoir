function out = fig_introductory_concepts(cfg)
% FIG_INTRODUCTORY_CONCEPTS Figure 1 panel A: chaos onset in a random network.
%
%   out = FIG_INTRODUCTORY_CONCEPTS()
%   out = FIG_INTRODUCTORY_CONCEPTS('gammas', [0.9 1.6 2.5])
%
% Two figures, saved into subfolders of this one:
%   statetraces/  membrane-potential time series at three levels of gain
%   eigenspectra/ the Jacobian eigenvalue disk at the same three gains
%
% Reproduces Sompolinsky, Crisanti & Sommers (1988): three networks sharing ONE
% underlying random weight matrix, scaled by three gains. Only the gain differs,
% which isolates the effect of regulation -- "same brain, different regulation".
% The spectral radius is R = gamma exactly, so chaos onset sits at gamma = 1:
% below it the network relaxes to the trivial fixed point (LLE < 0), above it
% the dynamics are chaotic.
%
% PORTED to SRNNCellTypePairs behind the 'sompolinsky_pairs' preset. The network
% is DELIBERATELY not the paper's operating point -- it is Dale-free, zero-mean,
% tanh, fully connected, with no adaptation and no input -- and that is the
% point of the figure: it is the textbook baseline the rest of the paper departs
% from. Naming it as a preset is what stops those values living as loose
% assignments inside a plotting script.
%
% TWO CELL TYPES, named A and B, both statistically identical. SRNNCellTypePairs
% cannot build a one-type model (see the preset's own comment), and calling two
% indistinguishable zero-mean populations E and I would be a lie. The traces
% concatenate both types, which together are the whole network.
%
% See also: srnn_param_preset, build_from_preset, fig_energy_landscape

arguments
    cfg.preset_name (1,:) char    = 'sompolinsky_pairs'
    cfg.gammas      (1,:) double  = [0.9, 1.6, 2.5]
    cfg.out_dir     (1,:) char    = ''
    cfg.n_traces    (1,1) double  = 15
    cfg.T_range     (1,2) double  = [0, 60]
    cfg.lya_interval (1,2) double = [30, 60]
    cfg.fs          (1,1) double  = 200
    cfg.net_seed    (1,1) double  = 0
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));

gammas   = cfg.gammas;
n_traces = cfg.n_traces;
T_plot   = cfg.T_range;
n_cases  = numel(gammas);
results  = struct('gamma', {}, 't', {}, 'x', {}, 'LLE', {}, 'R', {}, 'W', {});

%% Simulate the three networks
% Same net_seed for every gamma, so W_raw is identical and only the scale
% differs. level_of_chaos multiplies the assembled W, so R = gamma exactly.
for k = 1:n_cases
    gamma = gammas(k);
    fprintf('\n=== Case %d/%d : gamma = %.2f ===\n', k, n_cases, gamma);

    % 'no_adaptation' is the right condition here: this network has none, and
    % the condition is what carries n_a and synapse_config.
    model = build_from_preset(cfg.preset_name, 'no_adaptation', ...
        'level_of_chaos',       gamma, ...
        'T_range',              cfg.T_range, ...
        'fs',                   cfg.fs, ...
        'rng_seeds',            [cfg.net_seed, 1], ...
        'lya_method',           'benettin', ...
        'lya_T_interval',       cfg.lya_interval, ...
        'filter_local_lya',     false, ...
        'store_full_state',     false, ...
        'store_decimated_state', true, ...
        'plot_freq',            cfg.fs, ...
        'u_ex_scale',           0);
    model.run();

    LLE = model.lya_results.LLE;
    fprintf('  R (spectral radius) = %.3f | LLE = %.4f -> %s\n', ...
        model.R, LLE, ternary(LLE > 0, 'CHAOTIC', 'non-chaotic'));

    results(k).gamma = gamma;
    results(k).t     = model.plot_data.t;
    results(k).x     = all_type_states(model.plot_data.x);   % N x nt, both types
    results(k).LLE   = LLE;
    results(k).R     = model.R;
    results(k).W     = model.W;      % gamma already folded in (W = gamma*W_raw)
end
% Shared across cases; the eigenspectrum panel needs both. N comes from the
% model rather than a script-level constant, so it follows the preset.
tau_d = model.tau_d;
N     = model.n;

%% Plot: 1 x 3 linked time-series panels
fig = figure('Color', 'white', 'Position', [100, 300, 1200, 320]);
tl  = tiledlayout(1, n_cases, 'TileSpacing', 'compact', 'Padding', 'compact');

idx  = round(linspace(1, N, n_traces));   % which neurons to show
cols = lines(n_traces);
ax   = gobjects(1, n_cases);

for k = 1:n_cases
    ax(k) = nexttile; hold on;
    t = results(k).t;
    x = results(k).x;
    for j = 1:n_traces
        plot(t, x(idx(j), :), 'LineWidth', 0.6, 'Color', cols(j, :));
    end
    box off;
    set(ax(k), 'Color', 'none');
    ax(k).XAxis.Visible = 'off';   % hide time axis (scale bar added instead)
    if k == 1
        ylabel('x (AU)');
    end
end

linkaxes(ax, 'xy');       % shared x and y scale across the three panels
xlim(ax(1), T_plot);
set(ax, 'YTick', [-3, 0, 3]);

% 10 s time scale bar in the lower-right of the 1st subplot
bar_len = 10;                        % seconds
xr = xlim(ax(1)); yr = ylim(ax(1));
wx = xr(2) - xr(1); wy = yr(2) - yr(1);
x_end   = xr(2) - 0.03*wx;
x_start = x_end - bar_len;
y_bar   = yr(1) + 0.10*wy;
plot(ax(1), [x_start, x_end], [y_bar, y_bar], 'k', 'LineWidth', 4);
text(ax(1), mean([x_start, x_end]), y_bar - 0.03*wy, sprintf('%d s', bar_len), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);
hold(ax(1), 'off');


%% Eigenspectrum of the Jacobian at x = 0:  J = (-I + gamma*W)/tau_d
% (tanh'(0) = 1, so the linearization uses W directly; gamma is already in W.)
% Eigenvalues fill a disk centered at -1/tau_d with radius gamma/tau_d.
% The vertical (Im) axis at Re = 0 is the stability boundary: once the disk
% crosses it (gamma > 1) the fixed point is unstable.
%
% Styling matches ConnectivityAdaptation/RandomMatrixTheory (RMT.plot_spectrum /
% Fig_1_RMT_examples.m): black open-circle eigenvalues, solid black theoretical
% radius circle, axis off with hand-drawn Re/Im axes through the origin, equal
% aspect, common scaling, and (a)/(b)/(c) letters.
fig2  = figure('Color', 'white', 'Position', [100, 300, 977, 380]);
tl2   = tiledlayout(1, n_cases, 'TileSpacing', 'tight', 'Padding', 'compact'); %#ok<NASGU>
axe   = gobjects(1, n_cases);
theta = linspace(0, 2*pi, 200);
xc    = -1 / tau_d;                 % common disk center (shift from -I/tau_d)
mSize = 4;

% Plot eigenvalues + theoretical radius circle (RMT style)
for k = 1:n_cases
    axe(k) = nexttile; hold(axe(k), 'on');
    ev = eig((-eye(N) + results(k).W) / tau_d);
    Rk = results(k).gamma / tau_d;                 % theoretical radius
    unstable = real(ev) > 0;                        % right of the imaginary axis
    plot(axe(k), real(ev(~unstable)), imag(ev(~unstable)), 'ko', ...   % stable modes
        'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
    plot(axe(k), real(ev(unstable)), imag(ev(unstable)), 'o', ...      % unstable modes (Re>0)
        'MarkerSize', mSize, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', [0.7 0 0], 'LineWidth', 0.5);
    plot(axe(k), xc + Rk*cos(theta), Rk*sin(theta), 'k-', 'LineWidth', 2);
end

% Common scaling centered on the disk center (RMT-style equal aspect)
% margin 1.15, then zoomed out a further 20% so the Im axis/label clears the disk
max_R = max([results.gamma]) / tau_d;
common_radius = max_R * 1.15 * 1.20;
for k = 1:n_cases
    xlim(axe(k), [xc - common_radius, xc + common_radius]);
    ylim(axe(k), [-common_radius, common_radius]);
    daspect(axe(k), [1 1 1]);
end

% Format axes: hide box, draw Re/Im axes through the origin (RMT style)
for k = 1:n_cases
    axes(axe(k)); %#ok<LAXES>
    x_lim = xlim; y_lim = ylim;
    y_lim_axis = min(0.75*y_lim, 1.1*max_R);
    axis off;
    hold on;
    h_x = plot(x_lim, [0, 0], 'k', 'LineWidth', 1.5);      % Re axis
    h_y = plot([0, 0], y_lim_axis, 'k', 'LineWidth', 1.5); % Im axis = Re=0 stability boundary
    uistack([h_x, h_y], 'bottom');
    text(x_lim(2), 0, ' Re', 'Interpreter', 'latex', ...
        'VerticalAlignment', 'middle', 'FontSize', 16);
    text(0, y_lim_axis(2), 'Im', 'Interpreter', 'latex', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 16);
    xlim(x_lim); ylim(y_lim);
    hold off;
end

% (a), (b), (c) letters omitted for the eigenspectrum panels
drawnow;


if ~cfg.visible; set([fig, fig2], 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
% Two subfolders, matching the paths the manuscript references.
traces_dir = fullfile(out_dir, 'statetraces');
spec_dir   = fullfile(out_dir, 'eigenspectra');
out = struct('figs', [fig, fig2], 'files', {{}}, ...
             'source', ['preset: ' cfg.preset_name]);
if cfg.save
    save_figure_stable(traces_dir, 'panelA_bottom_traces', fig);
    save_figure_stable(spec_dir,   'panelA_eigenspectrum', fig2);
    out.files = [ ...
        strcat('statetraces/',  existing_outputs(traces_dir, 'panelA_bottom_traces')), ...
        strcat('eigenspectra/', existing_outputs(spec_dir,   'panelA_eigenspectrum'))];
    capture_git_provenance(out_dir, project_root);

    lle_note = sprintf('gamma = %s  ->  LLE = %s', ...
        mat2str(gammas, 3), mat2str([results.LLE], 3));

    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_introductory_concepts', ...
        'title',  'Stability_Manuscript figure 1 panel A: chaos onset in a random network', ...
        'script', 'fig_introductory_concepts.m', ...
        'what',   ['Two figures. statetraces: membrane-potential time series ' ...
                   'at three levels of synaptic gain. eigenspectra: the ' ...
                   'Jacobian eigenvalue disk at the same three gains, with the ' ...
                   'imaginary axis marking the stability boundary. Three ' ...
                   'networks share ONE underlying random weight matrix and ' ...
                   'differ only in gain, so the comparison isolates regulation ' ...
                   'from architecture.'], ...
        'how',    ['Reproduces Sompolinsky, Crisanti & Sommers (1988) on this ' ...
                   'project''s model class. Purely random Dale-free Gaussian ' ...
                   'connectivity, tanh, fully connected, no adaptation and no ' ...
                   'external input, so the spectral radius is R = gamma exactly ' ...
                   'and chaos onset sits at gamma = 1. The Jacobian at x = 0 is ' ...
                   'J = (-I + W)/tau_d, since tanh''(0) = 1 and gamma is ' ...
                   'already folded into W.'], ...
        'source', struct('preset', cfg.preset_name, 'net_seed', cfg.net_seed, ...
                         'gammas', gammas), ...
        'settings', figure_settings(model), ...
        'figures', {out.files}, ...
        'sections', struct( ...
            'heading', {'measured exponents', 'why two cell types'}, ...
            'body',    {lle_note, readme_two_types()})));
end
end

%% ------------------------------------------------------------------------
function X = all_type_states(s)
% Concatenate a per-cell-type state struct into one N x nt matrix.
%
% plot_data.x is keyed by cell type. This network's two types are statistically
% identical -- they exist only because SRNNCellTypePairs cannot build a
% one-type model -- so the whole network is their concatenation, and splitting
% the traces by type would imply a distinction that is not there.
f = fieldnames(s);
X = vertcat(s.(f{1}));
for i = 2:numel(f)
    X = [X; s.(f{i})]; %#ok<AGROW>
end
end

function s = readme_two_types()
s = ['WHY TWO CELL TYPES, NAMED A AND B. SRNNCellTypePairs cannot build a ' ...
     'one-cell-type model: build_W assigns the RMTBlocks generator piecemeal ' ...
     'where set_types is the only supported way to change the number of types, ' ...
     'so a scalar f is expanded back to two and the 1x1 mu_tilde then fails ' ...
     'validation. Two types with IDENTICAL zero-mean blocks is the same ' ...
     'network and builds today. They are named A and B rather than E and I ' ...
     'because they are statistically indistinguishable and the weights take ' ...
     'both signs -- calling them excitatory and inhibitory would be a lie. The ' ...
     'traces concatenate both types, which together are the whole network.'];
end

function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end
