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
%   'sfa_std'      4 columns -- none, SFA, STD, SFA+STD -- on the
%                  'single_neuron_dualStd' preset, which carries the PAPER'S c,
%                  tau_a and dual depression timescales on one unconnected
%                  neuron. This is the manuscript figure.
%   'sfa_std_stf'  7 columns, adding facilitation, on the 'single_neuron_stf'
%                  preset. Rebuilds the figure produced by the deleted
%                  test_single_neuron_stf.m (last version at commit 60c2992),
%                  which called a class family that no longer exists.
%
% BOTH VARIANTS ARE GENUINELY n = 1, C = 1, W = 0 -- one cell type, one neuron,
% no recurrence. Each preset says so itself; this function overrides no network
% parameter, which is what keeps the figure honest about what it simulated.
%
% THIS WAS NOT ALWAYS TRUE, and the failure is worth recording. Variant A used
% to name the paper's celltype_pairs_Sc0p2_noise0p025_dualStd preset directly
% and override nothing, so it built that preset WHOLE: n = 500, indegree = 100,
% fully recurrent. It then plotted neuron 1 and called it "one unconnected
% neuron". The generated README printed the contradiction outright -- "One
% unconnected neuron" above "n 500" -- because figure_settings reads off the
% BUILT OBJECT rather than the preset struct. The figure demonstrated none of
% the mechanisms: the no-adaptation column was network chaos rather than a step
% response, and the STD columns were silent with b ~ 1. The fix was a preset
% (single_neuron_dualStd) that states the single-neuron network once, rather
% than a figure that reaches into a network preset and forgets to shrink it.
%
% REBUILT ON SRNNCellTypePairs. The original was SRNNModel2 and assembled the
% figure by COPYING AXES out of model.plot(). The traces are now plotted from
% model.plot_data DIRECTLY: copying axes tied the figure to that method's row
% layout, which differs between the two model classes and would have had to be
% re-derived; reading the data is both stabler and lets a column show only the
% mechanisms it actually has.
%
% NOISE IS ON in the 'sfa_std' variant, because single_neuron_dualStd inherits
% sigma_u_noise = 0.025 from the paper's preset and is taken literally (TR's
% decision). On a single neuron there is no population averaging, so the jitter
% is fully visible -- x_noise_std = 0.056 against a 0.5 step, about 11%. That is
% accepted as honest about the model the paper characterises, not a defect to
% smooth away. single_neuron_stf carries no noise, so variant B is clean.
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
st           = manuscript_style();

%% Variant -> preset, columns, figure tag
switch cfg.variant
    case 'sfa_std'
        preset_name = pick(cfg.preset_name, 'single_neuron_dualStd');
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

% The SFA timescales this preset switches on, taken from the CONDITION. Every
% preset's conditions now state tau_a explicitly, so the figure reads the
% timescales rather than a count -- and a column that switches SFA off simply
% empties them.
%
% This used to read tau_a off the PRESET and truncate it per column, because
% validate() demanded numel(tau_a{1}) == n_a(1) and the two lived in different
% places. Presets no longer carry tau_a at all.
tau_a0 = conditions{ci}.tau_a;
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
    [n_a, sc, tau_a] = combo_config(combos{k}, sc0, tau_a0);
    model = build_from_preset(preset_name, 'sfa_and_std', ...
        'tau_a', tau_a, 'synapse_config', sc, ...
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

function [n_a, sc, tau_a] = combo_config(combo, sc0, tau_a0)
% Adaptation timescales and synapse routes for one column.
%
% tau_a and synapse_config are CONDITION-owned, which is exactly why they can be
% overridden per column: they never live in the preset, so setting them here
% does not fight anything.
%
% A column is SFA-on or SFA-off, nothing between, so tau_a is either the
% condition's timescales or empty. n_a follows from it and is passed only
% because it is still a settable property; it is numel(tau_a) by construction.
tau_a = repmat({zeros(1,0)}, 1, numel(tau_a0));
if contains(combo, 'sfa')
    tau_a{1} = tau_a0{1};        % SFA on the first (E) type only
end
n_a = cellfun(@numel, tau_a);
assert(isequal(n_a(1) > 0, contains(combo, 'sfa')), ...
    'combo_config: SFA state disagrees with the column name ''%s''.', combo);

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
% The one neuron. Both presets are n = 1, so this is the whole population
% rather than a choice among neurons; it stays a function because the row
% dimension is still there.
if isempty(M); y = []; else; y = M(1, :); end
end


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
