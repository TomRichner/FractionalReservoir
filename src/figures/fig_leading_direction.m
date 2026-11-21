function out = fig_leading_direction(cfg)
% FIG_LEADING_DIRECTION What decays slowly: the dendritic share of the leading perturbation against lambda_1.
%
%   out = FIG_LEADING_DIRECTION('run_dir', d)
%   out = FIG_LEADING_DIRECTION('run_dir', d, 'layout', 'column', 'reference_file', f)
%
% Panel C of the paper's Figure 4 (figure plan of 2026-11-21; the reason is
% section 2.3 of the manuscript repository's
% docs/Advancing_the_narrative/Report.md). One panel per condition, one
% point per configuration of the joint parameter-space sample
% (<run_dir>/param_space_*): horizontally the run's leading growth rate
% lambda_1, vertically the share of the leading perturbation's squared norm
% that lies in the dendritic coordinates (lead_frac_x, saved per run by the
% sweep; the rest lies in the adaptation and depression variables). A
% dashed line at zero growth, and a dotted line at the growth rate of the
% nearly uncoupled network (the uncoupled_reference stage).
%
% WHAT IT SHOWS. Where an adapting network contracts, its points sit near
% the dotted line with almost no dendritic share: what decays slowly there
% is each neuron's own adaptation, not a mode of the network. Where it
% expands, the leading perturbation is mostly dendritic. Without adaptation
% the state is dendritic only, so every point is at 1.
%
% THE HORIZONTAL AXIS is symmetric-logarithmic (symlog_lle): linear within
% lambda0 of zero, logarithmic beyond, alike on both sides, so the adapting
% conditions' rates near -0.1 1/s and the unadapted network's near +-5 1/s
% are both readable. |lambda_1| beyond 'clip' is drawn at the clip; the
% table says how many.
%
% THE REFERENCE. <run_dir>/uncoupled_reference/uncoupled_reference_data.mat
% (run_uncoupled_reference_stage), or the file given as 'reference_file'
% (the escape hatch for a run made before the stage existed: the caller
% says which reference it is). Each condition's line is the median over the
% stage's repetitions at its first (lowest) gain. A run WITHOUT the stage
% is drawn without the dotted line, with a printed note and a sentence in
% the table: the panel's points do not depend on it. 'require_reference',
% true makes its absence an error.
%
% OPTIONS
%   conditions   names to draw, in order; empty = all the run holds
%                (select_conditions; FD-006)
%   layout       'row' (default; panels side by side) or 'column' (stacked,
%                one shared x axis: the form fig_main4_connectivity places
%                beside panels A and B)
%   lambda0      linear scale of the axis, 1/s (default 0.05)
%   clip         |lambda_1| drawn at most here, 1/s (default 10)
%   reference_file, require_reference   see above
%
% No text on the panels. Numbers: Fig_Leading_Direction_table.md.
%
% out.figs, out.files, out.source ({joint sample dir, reference file or ''}),
% out.table_text, out.stats (per condition), out.axes (panels, in order)
%
% See also: run_uncoupled_reference_stage, uncoupled_reference_stage_plan,
%           fig_main4_connectivity, fig_lle_vs_rate, symlog_lle

arguments
    cfg.verbose        (1,:) char    = 'minimal'
    cfg.data_file      (1,:) char    = ''     % the joint sample's psa_object.mat, when given explicitly
    cfg.out_dir        (1,:) char    = ''
    cfg.save           (1,1) logical = true
    cfg.visible        (1,1) logical = true
    cfg.run_dir        (1,:) char    = ''
    cfg.preset_name    (1,:) char    = ''     % unused; the run records its preset
    cfg.conditions                   = {}
    cfg.layout         (1,:) char {mustBeMember(cfg.layout, {'row', 'column'})} = 'row'
    cfg.lambda0        (1,1) double {mustBePositive} = 0.05
    cfg.clip           (1,1) double {mustBePositive} = 10
    cfg.reference_file (1,:) char    = ''
    cfg.require_reference (1,1) logical = false
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st = manuscript_style();

%% The joint sample
searched = {fullfile(cfg.run_dir, 'param_space')};   % named so that a run without one errors through resolve_data_file
if ~isempty(cfg.run_dir)
    ps = dir(fullfile(cfg.run_dir, 'param_space_*'));
    ps = ps([ps.isdir]);
    if ~isempty(ps)
        [~, order] = sort([ps.datenum], 'descend');
        searched = arrayfun(@(d) fullfile(d.folder, d.name), ps(order), 'UniformOutput', false)';
    end
end
psa_file = resolve_data_file(cfg.data_file, cfg.run_dir, searched, 'psa_object.mat', ...
    'Run the joint parameter-space sample first (the sweeps stage with ''param_space'')');
psa_dir = fileparts(psa_file);
psa = ParamSpaceAnalysis2.from_dir(psa_dir);
held = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
[names, ~] = select_conditions(held, cfg.conditions);
titles = grouped_condition_label(names, 'sentence');
n_cond = numel(names);

S = struct('name', names, 'title', titles, 'lambda', [], 'frac_x', [], 'n', [], 'n_clipped', [], ...
    'reference', NaN, 'ref_range', [NaN NaN], 'ref_n', 0, 'ref_gain', NaN, 'ref_all_gains', [], 'ref_all_medians', []);
for i = 1:n_cond
    if ~isfield(psa.results, names{i})
        error('fig_leading_direction:NoJobs', 'The joint sample in %s holds no results for %s.', psa_dir, names{i});
    end
    R = psa.results.(names{i});
    ok = cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success && isfield(r, 'LLE') && isfinite(r.LLE) && ...
        isfield(r, 'lead_frac_x') && isfinite(r.lead_frac_x), R);
    R = R(ok);
    if isempty(R)
        error('fig_leading_direction:NoJobs', ['No successful run of %s in %s carries lead_frac_x ' ...
            '(the sweep saves it since the top-K Lyapunov metrics).'], names{i}, psa_dir);
    end
    S(i).lambda = cellfun(@(r) r.LLE, R(:));
    S(i).frac_x = cellfun(@(r) r.lead_frac_x, R(:));
    S(i).n = numel(R);
    S(i).n_clipped = nnz(abs(S(i).lambda) > cfg.clip);
end

%% The uncoupled reference
ref_file = cfg.reference_file;
if isempty(ref_file) && ~isempty(cfg.run_dir)
    candidate = fullfile(cfg.run_dir, 'uncoupled_reference', 'uncoupled_reference_data.mat');
    if isfile(candidate); ref_file = candidate; end
elseif ~isempty(ref_file) && ~isfile(ref_file)
    error('fig_leading_direction:NoSuchFile', 'The reference file given explicitly does not exist:\n  %s', ref_file);
end
ref_settings = struct();
if isempty(ref_file)
    if cfg.require_reference
        error('fig_leading_direction:NoReference', ['%s holds no uncoupled reference ' ...
            '(uncoupled_reference/uncoupled_reference_data.mat). Run the uncoupled_reference stage, or give ''reference_file''.'], ...
            cfg.run_dir);
    end
    vprintf(cfg.verbose, 'minimal', ['[fig_leading_direction] no uncoupled reference in %s (stage uncoupled_reference): ' ...
        'the dotted line is not drawn.\n'], cfg.run_dir);
else
    U = load(ref_file, 'ref', 'settings');
    ref_settings = U.settings;
    for i = 1:n_cond
        k = find(strcmp({U.ref.name}, names{i}));
        if isempty(k); continue; end       % a condition the stage did not run: no line, said in the table
        S(i).reference = U.ref(k).reference;
        S(i).ref_range = [U.ref(k).lambda_1_min(1), U.ref(k).lambda_1_max(1)];
        S(i).ref_n = U.ref(k).n_ok(1);
        S(i).ref_gain = U.ref(k).gains(1);
        S(i).ref_all_gains = U.ref(k).gains;
        S(i).ref_all_medians = U.ref(k).lambda_1_median;
    end
end

%% Draw
[~, tick_pos, tick_lab, xlim_y] = symlog_lle(0, cfg.lambda0, cfg.clip);
to_x = @(v) symlog_lle(v, cfg.lambda0, cfg.clip);
is_row = strcmp(cfg.layout, 'row');
if is_row
    fig = figure('Color', 'w', 'Position', [150, 150, 130 + 330 * n_cond, 330]);
    tl = tiledlayout(fig, 1, n_cond, 'TileSpacing', 'compact', 'Padding', 'loose');
else
    fig = figure('Color', 'w', 'Position', [150, 80, 340, 115 + 200 * n_cond]);
    tl = tiledlayout(fig, n_cond, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
end
if ~cfg.visible; set(fig, 'Visible', 'off'); end
ax = gobjects(1, n_cond);
for i = 1:n_cond
    ax(i) = nexttile(tl); hold(ax(i), 'on');
    col = st.condition_color(names{i});
    plot(ax(i), [0 0], [-0.05 1.05], '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw, ...
        'Tag', 'zero_growth_reference');
    if isfinite(S(i).reference)
        plot(ax(i), to_x(S(i).reference) * [1 1], [-0.05 1.05], ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.6, ...
            'Tag', 'uncoupled_reference');
    end
    scatter(ax(i), to_x(S(i).lambda), S(i).frac_x, 22, col, 'filled', 'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeColor', 'none', 'Tag', 'leading_direction_points');
    set(ax(i), 'XLim', xlim_y * 1.03, 'YLim', [-0.05 1.05], 'XTick', tick_pos, 'XTickLabel', tick_lab, ...
        'YTick', [0 0.5 1], 'TickDir', 'out', 'Box', 'off', 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, ...
        'XTickLabelRotation', 0, 'Tag', sprintf('leading_direction_c%d', i));
    title(ax(i), titles{i}, 'FontWeight', 'normal', 'FontSize', st.title_fs, 'Color', col);
    if is_row || i == n_cond
        xlabel(ax(i), 'Leading growth rate \lambda_1 (s^{-1})', 'FontSize', st.tick_fs);
    else
        set(ax(i), 'XTickLabel', {});
    end
end
if is_row
    ylabel(ax(1), {'Dendritic share of the', 'leading perturbation'}, 'FontSize', st.tick_fs);
else
    % On the middle panel, not on the layout: a layout's label is not an axes
    % child and would be lost when the panels are copied into a composite.
    ylabel(ax(ceil(n_cond / 2)), 'Dendritic share of the leading perturbation', 'FontSize', st.tick_fs, 'Clipping', 'off');
end

%% Table
fig_tag = 'Fig_Leading_Direction';
L = {};
L{end + 1} = '# What decays slowly: the dendritic share of the leading perturbation against the leading growth rate';
L{end + 1} = '';
L{end + 1} = sprintf(['Joint parameter-space sample: `%s`. One point per configuration that ran. Vertical: `lead_frac_x`, the ' ...
    'share of the leading tangent vector''s squared norm in the dendritic coordinates (the rest is in the adaptation and ' ...
    'depression variables). Horizontal: the leading growth rate, on a symmetric-logarithmic axis (linear within %g 1/s of ' ...
    'zero; |lambda_1| above %g 1/s drawn at %g). Dashed: zero growth. Dotted: the nearly uncoupled network.'], ...
    psa_dir, cfg.lambda0, cfg.clip, cfg.clip);
L{end + 1} = '';
if isempty(ref_file)
    L{end + 1} = ['**No uncoupled reference in this run** (stage `uncoupled_reference`): the dotted line is not drawn and ' ...
        'the reference columns are empty.'];
else
    L{end + 1} = sprintf(['Uncoupled reference: `%s`. The same network with every recurrent weight scaled to %s of the ' ...
        'preset''s (level_of_chaos), %d repetitions (network seeds) per value, n = %g, T = [%g, %g] s, Lyapunov window ' ...
        '[%g, %g] s, %g Hz, %s. The line is the median at the first value.'], ref_file, mat2str(ref_settings.gains), ...
        ref_settings.n_reps, ref_settings.n, ref_settings.T_range, ref_settings.lya_T_interval, ref_settings.fs, ...
        ref_settings.ode_solver);
end
L{end + 1} = '';
L{end + 1} = ['| Condition | runs | median lambda_1 (1/s) | lambda_1 < 0 | median dendritic share | dendritic share < 0.1 | ' ...
    'same, among lambda_1 < 0 | same, among lambda_1 > 0 | uncoupled lambda_1 (1/s): median [min, max], runs | ' ...
    'at each gain | runs within 25% of it | runs with lambda_1 < 0 within 25% of it | drawn at the clip |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|---|---|---|';
for i = 1:n_cond
    s = S(i); lam = s.lambda; fx = s.frac_x; neg = lam < 0; pos = lam > 0;
    if isfinite(s.reference)
        near = abs(lam - s.reference) < 0.25 * abs(s.reference);
        ref_txt = sprintf('%.3f [%.3f, %.3f], %d', s.reference, s.ref_range, s.ref_n);
        gains_txt = strjoin(arrayfun(@(g, m) sprintf('%g: %.3f', g, m), s.ref_all_gains, s.ref_all_medians, ...
            'UniformOutput', false), '; ');
        near_txt = sprintf('%.2f', mean(near)); near_neg_txt = share(near(neg));
    else
        ref_txt = 'n/a'; gains_txt = 'n/a'; near_txt = 'n/a'; near_neg_txt = 'n/a';
    end
    L{end + 1} = sprintf('| %s | %d | %.3f | %.2f | %.3f | %.2f | %s | %s | %s | %s | %s | %s | %d |', s.title, s.n, ...
        median(lam), mean(neg), median(fx), mean(fx < 0.1), share(fx(neg) < 0.1), share(fx(pos) < 0.1), ...
        ref_txt, gains_txt, near_txt, near_neg_txt, s.n_clipped); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = ['Columns. **dendritic share < 0.1**: share of runs whose leading perturbation is less than 10% dendritic. ' ...
    '**uncoupled lambda_1**: over the reference''s repetitions at its lowest gain. **at each gain**: the median at every ' ...
    'gain the stage ran, to show that the rate has stopped changing. **within 25%**: |lambda_1 - reference| below a ' ...
    'quarter of |reference|.'];
L{end + 1} = '';
L{end + 1} = '## Caption sentence, to adapt';
L{end + 1} = '';
L{end + 1} = ['For each condition, the share of the leading perturbation that lies in the dendritic coordinates against ' ...
    'the leading growth rate, one point per configuration of the joint sample (symmetric-logarithmic axis). Dashed: ' ...
    'zero growth. Dotted: the growth rate of the same network with its recurrent weights scaled to near zero.'];
table_text = strjoin(L, newline);
vprintf(cfg.verbose, 'verbose', '%s\n', table_text);

out = struct('figs', fig, 'files', {{}}, 'source', {{psa_dir, ref_file}}, 'table_text', table_text, ...
    'stats', S, 'axes', ax);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n', table_text);
        fclose(fid);
    end
end
end

function txt = share(tf)
% The mean of a logical vector as text; 'n/a' when it is empty.
if isempty(tf); txt = 'n/a'; else; txt = sprintf('%.2f', mean(tf)); end
end
