function out = fig_adaptation_methods(cfg)
% FIG_ADAPTATION_METHODS Single-neuron demonstration of SFA, STD and STF.
%
%   out = FIG_ADAPTATION_METHODS()                       % SFA + STD  (figure A)
%   out = FIG_ADAPTATION_METHODS('variant', 'sfa_std_stf')  %          (figure B)
%
% One unconnected neuron driven by a step in external input, under a set of
% mechanism combinations shown as COLUMNS. Rows are the stimulus u, the
% dendritic state x, the firing rate r, the synaptic output, and whichever
% mechanism states are active: adaptation a (SFA), depression b (STD),
% facilitation g (STF).
%
% TWO VARIANTS, one function:
%
%   'sfa_std'      4 columns -- none, SFA, STD, SFA+STD -- on the PAPER'S
%                  preset, so the mechanism cartoon uses the same c, tau_a and
%                  depression timescales as the network figures. This is the
%                  manuscript figure.
%   'sfa_std_stf'  7 columns, adding facilitation, on the 'single_neuron_stf'
%                  preset. Rebuilds the figure produced by the deleted
%                  test_single_neuron_stf.m (last version at commit 60c2992),
%                  which called a class family that no longer exists.
%
% REBUILT ON SRNNCellTypePairs. The original was SRNNModel2 with n = 1 and
% W = 0, and it assembled the figure by COPYING AXES out of model.plot(). Two
% changes:
%
%   * n = 2, not n = 1. SRNNCellTypePairs enforces n >= n_cellTypes and rejects
%     indegree = 0 (it requires 0 < indegree <= n), and it cannot build a
%     one-cell-type model at all (build_W assigns RMTBlocks piecemeal where
%     set_types is required to change D). Two neurons with identically zero
%     weights is the smallest expressible unconnected network; only the E
%     neuron is plotted, and with W == 0 the I neuron cannot influence it.
%   * The traces are plotted from model.plot_data DIRECTLY rather than by
%     harvesting axes from model.plot(). Copying axes tied the figure to that
%     method's row layout, which differs between the two model classes and
%     would have had to be re-derived; reading the data is both stabler and
%     lets a column show only the mechanisms it actually has.
%
% NOISE IS ON in the 'sfa_std' variant, because the preset carries
% sigma_u_noise = 0.025 and is taken literally (TR's decision). On a single
% neuron there is no population averaging, so the jitter is fully visible; that
% is accepted as honest about the model the paper characterises, not a defect
% to smooth away.
%
% See also: build_from_preset, srnn_param_preset, fig_example_timeseries

arguments
    cfg.variant     (1,:) char {mustBeMember(cfg.variant, {'sfa_std','sfa_std_stf'})} = 'sfa_std'
    cfg.preset_name (1,:) char    = ''     % '' -> per variant
    cfg.out_dir     (1,:) char    = ''
    cfg.step_amp    (1,1) double  = 0.5
    cfg.t_on        (1,1) double  = 5
    cfg.t_off       (1,1) double  = 15
    cfg.T_range     (1,2) double  = [0 25]
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''     % unused; accepted for a uniform call
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

%% Variant -> preset, columns, figure tag
switch cfg.variant
    case 'sfa_std'
        preset_name = pick(cfg.preset_name, 'celltype_pairs_Sc0p2_noise0p025_dualStd');
        combos = {'none', 'sfa', 'std', 'sfa+std'};
        fig_tag = 'Fig_single_neuron_SFA_STD';
    case 'sfa_std_stf'
        preset_name = pick(cfg.preset_name, 'single_neuron_stf');
        combos = {'none', 'sfa', 'std', 'stf', 'sfa+std', 'std+stf', 'sfa+std+stf'};
        fig_tag = 'Fig_single_neuron_SFA_STD_STF';
end
titles = cellfun(@combo_title, combos, 'UniformOutput', false);

%% The routes and SFA count this preset switches on
[~, model_class, conditions] = srnn_param_preset(preset_name);
assert(strcmp(model_class, 'SRNNCellTypePairs'), ...
    'fig_adaptation_methods needs a SRNNCellTypePairs preset; ''%s'' is %s.', ...
    preset_name, model_class);
ci   = find(cellfun(@(c) strcmp(c.name, 'sfa_and_std'), conditions), 1);
sc0  = conditions{ci}.synapse_config;
n_a0 = conditions{ci}.n_a;

% The preset's own tau_a, if it carries one. Most presets do NOT -- tau_a is the
% class default, auto-filled per n_a inside build() -- but single_neuron_stf
% carries an explicit single 3 s timescale. That matters here: validate()
% requires numel(tau_a{1}) == n_a(1) EXACTLY, so a column that switches SFA off
% must also shorten tau_a, or it fails with
%   "tau_a{1} must contain n_a(1) positive values"
% Presets without tau_a keep [] and let the auto-fill handle every column.
preset_struct = srnn_param_preset(preset_name);
if isfield(preset_struct, 'tau_a') && ~isempty(preset_struct.tau_a)
    tau_a0 = preset_struct.tau_a;
else
    tau_a0 = [];
end
has_stf = isfield(sc0, 'E') && isfield(sc0.E, 'E') && isfield(sc0.E.E, 'stf');
if strcmp(cfg.variant, 'sfa_std_stf') && ~has_stf
    error('fig_adaptation_methods:NoSTF', ...
        ['Preset ''%s'' declares no STF on E->E, so the STF columns would be ' ...
         'identical to the STD ones. Use the single_neuron_stf preset.'], preset_name);
end

%% Deterministic step stimulus, shared by every column
input_config = struct();
input_config.intrinsic_drive = 0;
input_config.generator = @(params, T, fs, seed, ic) ...
    step_generator(params, T, fs, cfg.t_on, cfg.t_off, cfg.step_amp);

%% Run every column
n_col = numel(combos);
P = cell(1, n_col);
for k = 1:n_col
    [n_a, sc, tau_a] = combo_config(combos{k}, n_a0, sc0, tau_a0);
    extra = {};
    if ~isempty(tau_a); extra = {'tau_a', tau_a}; end
    model = build_from_preset(preset_name, 'sfa_and_std', ...
        'n_a', n_a, 'synapse_config', sc, extra{:}, ...
        'input_config', input_config, ...
        'T_range', cfg.T_range, 'fs', 400, ...
        'lya_method', 'none', 'plot_deci', 1);
    model.run();
    P{k} = model.plot_data;
end

%% Which rows to draw
% Fixed row order; a mechanism row is included only if SOME column uses it, so
% the sheet has no all-blank rows.
rows = {'u', 'x', 'r', 'syn'};
if any(contains(combos, 'sfa')); rows{end+1} = 'a'; end
if any(contains(combos, 'std')); rows{end+1} = 'b'; end
if any(contains(combos, 'stf')); rows{end+1} = 'g'; end
n_row = numel(rows);

row_label = containers.Map( ...
    {'u', 'x', 'r', 'syn', 'a', 'b', 'g'}, ...
    {'u  (input)', 'x', 'r', 'synaptic out', 'a  (SFA)', 'b  (STD)', 'g  (STF)'});

%% Draw
fig = figure('Color', 'w', 'Position', [80 60 260*n_col 130*n_row]);
tl  = tiledlayout(fig, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
ax  = gobjects(n_row, n_col);

for rr = 1:n_row
    for cc = 1:n_col
        ax(rr, cc) = nexttile(tl, (rr-1)*n_col + cc);
        y = row_trace(P{cc}, rows{rr});
        if isempty(y)
            axis(ax(rr, cc), 'off');       % this column has no such mechanism
            continue
        end
        plot(ax(rr, cc), P{cc}.t, y, 'LineWidth', st.line_lw, 'Color', [0 0 0]);
        box(ax(rr, cc), 'off');
        set(ax(rr, cc), 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw);
        xlim(ax(rr, cc), cfg.T_range);
        if cc == 1
            ylabel(ax(rr, cc), row_label(rows{rr}), 'FontSize', st.label_fs, ...
                'Interpreter', 'none');
        end
        if rr == 1
            title(ax(rr, cc), titles{cc}, 'FontWeight', 'normal', ...
                'FontSize', st.title_fs, 'Interpreter', 'none');
        end
        if rr < n_row
            set(ax(rr, cc), 'XTickLabel', []);
        else
            xlabel(ax(rr, cc), 'time (s)', 'FontSize', st.label_fs);
        end
    end
end

% One y-scale per ROW, so a mechanism can be compared across columns. Done per
% row rather than globally because u, x, r and the gates live on different
% scales entirely.
for rr = 1:n_row
    live = ax(rr, arrayfun(@(a) isvalid(a) && strcmp(a.Visible, 'on'), ax(rr, :)));
    if numel(live) > 1; linkaxes(live, 'y'); end
end

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', ['preset: ' preset_name]);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    capture_git_provenance(out_dir, project_root);

    write_figure_readme(out_dir, struct( ...
        'tag',    ['fig_adaptation_methods_' cfg.variant], ...
        'title',  sprintf('Stability_Manuscript figure: single-neuron adaptation mechanisms (%s)', cfg.variant), ...
        'script', 'fig_adaptation_methods.m', ...
        'what',   ['One unconnected neuron driven by a step in external input, ' ...
                   'under each mechanism combination as a column. Rows: the ' ...
                   'stimulus, the dendritic state x, the firing rate r, the ' ...
                   'synaptic output, and the active mechanism states.'], ...
        'how',    ['Built from the named preset with n_a and synapse_config ' ...
                   'overridden per column, so every column shares the preset''s ' ...
                   'timescales and differs only in which mechanisms are on. ' ...
                   'Traces are read from model.plot_data directly rather than ' ...
                   'harvested from model.plot() axes.'], ...
        'source', struct('preset', preset_name, 'variant', cfg.variant, ...
                         'columns', {combos}), ...
        'settings', figure_settings(model), ...
        'figures', {out.files}, ...
        'sections', struct('heading', {'why n = 2', 'stf mapping'}, ...
                           'body', {readme_n2(), readme_stf(cfg.variant, sc0, has_stf)})));
end
end

%% ------------------------------------------------------------------------
function v = pick(a, b)
if isempty(a); v = b; else; v = a; end
end

function t = combo_title(c)
switch c
    case 'none',        t = 'No adaptation';
    case 'sfa',         t = 'SFA only';
    case 'std',         t = 'STD only';
    case 'stf',         t = 'STF only';
    case 'sfa+std',     t = 'SFA + STD';
    case 'std+stf',     t = 'STD + STF';
    case 'sfa+std+stf', t = 'SFA + STD + STF';
    otherwise,          t = c;
end
end

function [n_a, sc, tau_a] = combo_config(combo, n_a0, sc0, tau_a0)
% Adaptation counts, synapse routes and tau_a for one column.
%
% n_a and synapse_config are CONDITION-owned fields, which is exactly why they
% can be overridden per column: they never live in the preset, so setting them
% here does not fight anything.
n_a = zeros(size(n_a0));
if contains(combo, 'sfa')
    n_a(1) = n_a0(1);            % SFA on the first (E) type only
end

% tau_a must be TRUNCATED to match, when the preset carries one at all:
% validate() demands numel(tau_a{1}) == n_a(1) exactly, so an SFA-off column
% with a length-1 tau_a would be rejected. Empty means "let build() auto-fill",
% which is what every preset without an explicit tau_a wants.
tau_a = [];
if ~isempty(tau_a0)
    tau_a = tau_a0;
    tau_a{1} = tau_a0{1}(1:n_a(1));
end

sc = struct();
route = sc0.E.E;
if contains(combo, 'std') && isfield(route, 'std')
    sc.E.E.std = route.std;
end
if contains(combo, 'stf') && isfield(route, 'stf')
    sc.E.E.stf = route.stf;
end
end

function y = row_trace(p, name)
% One row's trace for the E neuron, or [] when this column lacks that mechanism.
y = [];
switch name
    case 'u',   y = first_row(getfield_or(p.u, 'E'));
    case 'x',   y = first_row(getfield_or(p.x, 'E'));
    case 'r',   y = first_row(getfield_or(p.r, 'E'));
    case 'syn'
        if isfield(p.synaptic_output, 'E') && isfield(p.synaptic_output.E, 'E')
            y = first_row(p.synaptic_output.E.E);
        end
    case 'a'
        A = getfield_or(p.a, 'E');
        % a is n x n_a x nt; show the SUM over timescales, which is what the
        % rate actually subtracts (r = phi(x - c*sum_k a_k)).
        if ~isempty(A); y = squeeze(sum(A(1, :, :), 2))'; end
    case 'b'
        if isfield(p.b, 'E') && isfield(p.b.E, 'E') && ~isempty(p.b.E.E)
            % b is n x n_b x nt; show the PRODUCT, which is what the synapse
            % sees when there is more than one depression timescale.
            y = squeeze(prod(p.b.E.E(1, :, :), 2))';
        end
    case 'g'
        if isfield(p.g, 'E') && isfield(p.g.E, 'E') && ~isempty(p.g.E.E)
            y = squeeze(prod(p.g.E.E(1, :, :), 2))';
        end
end
end

function v = getfield_or(s, name)
if isstruct(s) && isfield(s, name); v = s.(name); else; v = []; end
end

function y = first_row(M)
if isempty(M); y = []; else; y = M(1, :); end
end

function s = readme_n2()
s = ['WHY n = 2 AND NOT 1. SRNNCellTypePairs enforces n >= n_cellTypes and ' ...
     'rejects indegree = 0 (it requires 0 < indegree <= n), and it cannot ' ...
     'build a one-cell-type model at all -- build_W assigns the RMTBlocks ' ...
     'generator piecemeal where set_types is the only supported way to change ' ...
     'the number of types, so a scalar f is expanded back to two types and the ' ...
     '1x1 mu_tilde then fails validation. Two neurons with identically zero ' ...
     'weights is the smallest expressible unconnected network. Only the E ' ...
     'neuron is plotted, and with W == 0 the second neuron cannot influence it.'];
end

function s = readme_stf(variant, sc0, has_stf)
if ~strcmp(variant, 'sfa_std_stf') || ~has_stf
    s = ['This variant shows SFA and STD only. The facilitation columns are in ' ...
         'the sfa_std_stf variant, which uses the single_neuron_stf preset.'];
    return
end
f = sc0.E.E.stf;
s = sprintf([ ...
    'STF PARAMETERS map EXACTLY from the deleted test_single_neuron_stf.m. That ' ...
    'model carried a release probability p resting at p0 with gain p/p0 and ' ...
    'dp/dt = (p0-p)/tau_f + kappa*(1-p)*r. Substituting the gain variable ' ...
    'u = p/p0 gives du/dt = (1-u)/tau_f + kappa*(1/p0 - u)*r, which IS this ' ...
    'class''s dg/dt = (1-g)/tau_dec + (G-g)*r/tau_fac with tau_dec = tau_f, ' ...
    'tau_fac = 1/kappa and G = 1/p0. Here tau_dec = %g, tau_fac = %g, G = %.4g, ' ...
    'i.e. the old tau_f = 6, kappa = 0.4, p0 = 0.35. ' ...
    'WHAT DOES NOT MAP: the old STD depleted in proportion to p ' ...
    '(db/dt = (1-b)/tau_rec - (p*b*r)/tau_rel, the Tsodyks-Markram coupling), ' ...
    'where this class''s depletion is independent of g. No tau_rel reproduces ' ...
    'that -- it is a constant factor at rest but ACCELERATES as p rises. ' ...
    'tau_rel is carried literally, so this figure will NOT match the archived ' ...
    'PNG (about 2.9x stronger depression at rest). That is expected. Whether ' ...
    'the decoupling is correct at all is recorded in UserNotes.md.'], ...
    f.tau_dec, f.tau_fac, f.G);
end

%% ------------------------------------------------------------------------
function [u_ex, t_ex] = step_generator(params, T, fs, t_on, t_off, step_amp)
% Deterministic single step, applied to every neuron.
%
% The signature is fixed by SRNNCellTypePairs: the class calls
% config.generator(params, T, fs, rng_seed, config), so the extra arguments are
% bound by the anonymous wrapper at the call site rather than passed through.
dt   = 1 / fs;
t_ex = (0:dt:T)';
u_ex = zeros(params.n, numel(t_ex));
u_ex(:, (t_ex >= t_on) & (t_ex < t_off)) = step_amp;
end
