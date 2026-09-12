function out = fig_sfa_EOC_allStd(cfg)
% FIG_SFA_EOC_ALLSTD SFA edge of chaos: lambda_1 against the slowest tau_a.
%
%   out = FIG_SFA_EOC_ALLSTD()
%   out = FIG_SFA_EOC_ALLSTD('run_dir', d)
%
% Replots the tau-sensitivity LLE panel -- how the largest Lyapunov exponent
% approaches 0 (the edge of chaos) as the maximum SFA adaptation timescale grows.
% No simulation is re-run: it reloads the saved tau_a_E_max PSA object and
% restyles the single-condition (SFA+STD) panel.
%
% THE RUN IS RESOLVED, NOT HARDCODED. This script used to carry
%   data_root = .../run_all_aug_14_26_17_25
% and three sibling figures carried the same line while a fourth pointed at a
% different, older run on a different model class -- so "regenerate the figures"
% silently built the manuscript from two runs. resolve_run_dir picks the newest
% run whose manifest names the requested preset, and errors if there is none.
%
% See also: resolve_run_dir, ParamSpaceAnalysis2/plot_sensitivity, paper_config

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% --- Locate the tau sensitivity subfolder ----------------------------------
% By GLOB rather than by timestamped name, so re-pointing at another run needs
% no other edit.
tau_listing = dir(fullfile(run_dir, 'tau_sensitivity_*'));
tau_listing = tau_listing([tau_listing.isdir]);
if isempty(tau_listing)
    error('fig_sfa_EOC_allStd:NoTauDir', ...
        'No tau_sensitivity_* subfolder found in:\n  %s', run_dir);
end
if numel(tau_listing) > 1
    warning('fig_sfa_EOC_allStd:MultipleTauDirs', ...
        'Found %d tau_sensitivity_* subfolders; using the newest.', numel(tau_listing));
    [~, newest] = max([tau_listing.datenum]);
    tau_listing = tau_listing(newest);
end
tau_dir = fullfile(tau_listing.folder, tau_listing.name);

% --- Presentation constants ------------------------------------------------
% Kept identical to the original figure by decision, so the two read the same.
% Worth knowing when reading the output: on this run the tau LLEs span
% -0.26 .. +0.29 with a median of +0.008, so slightly over half the
% distribution is POSITIVE and 32% of it sits above the y_view ceiling of 0.05
% -- in the +inf overflow band, which is itself off-screen. The panel therefore
% reads as almost entirely sub-zero when it is not. Widen lle_range / y_view /
% y_ticks together to show it.
lle_range = [-0.3, 0.1];    % histogram range passed to plot_sensitivity
y_view    = [-0.25, 0.05];  % visible y-range (crops the overflow bands)
y_ticks   = -0.2:0.1:0.1;
fig_position = [457 637 264 252];

tick_fs   = st.tick_fs;
label_fs  = st.label_fs;
clim_frac = 0.8;   % darken imagesc: cap CLim at total_reps*clim_frac
% Colormap ramps white (0 counts) -> 90% black (max), not pure black, so the
% blue median line stays visible over the darkest cells.
dark_cmap    = repmat(linspace(1, 0.1, 256)', 1, 3);
median_alpha = 0.35;   % blue median line transparency
median_lw    = st.median_lw;
zeroline_lw  = st.zeroline_lw;

% Start clean: plot_sensitivity operates on the current session's figures.
close all force;

% --- Reload the tau PSA and regenerate the single LLE panel ----------------
psa = ParamSpaceAnalysis2.from_dir(tau_dir);

psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_range);
cf = gcf;

% --- Presentation restyle (same knobs as Fig_sensitivity_analysis_allStd) ---
ax = findobj(cf, 'Type', 'axes');
set(ax, 'FontSize', tick_fs);
box(ax, 'off');
colormap(ax, dark_cmap);

% Darken the histogram density (shared black->white scale from plot_sensitivity).
cl = get(ax, 'CLim');
set(ax, 'CLim', [0, cl(2) * clim_frac]);

% Blue median line: more transparent + thinner. (imagesc is Type 'image', the
% zero line is 'constantline', so 'line' is the median.)
ml = findobj(ax, 'Type', 'line');
for m = 1:numel(ml)
    mc = get(ml(m), 'Color');
    if numel(mc) < 4; mc(4) = 1; end
    mc(4) = median_alpha;
    set(ml(m), 'Color', mc, 'LineWidth', median_lw);
end
% Green dashed zero line: thinner.
set(findobj(ax, 'Type', 'constantline'), 'LineWidth', zeroline_lw);

% Labels: keep lambda_1; x -> "max \tau_a (s)"; drop the title. latex here (not
% tex) so the symbol matches the latex xlabel below.
ylabel(ax, '$\lambda_1$', 'Interpreter', 'latex', 'FontSize', label_fs);
xlabel(ax, 'max $\tau_a$ (s)', 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax, '');

% Fewer y-ticks + a tighter view around the data.
ylim(ax, y_view);
yticks(ax, y_ticks);

set(cf, 'Position', fig_position);

if ~cfg.visible; set(cf, 'Visible', 'off'); end

%% --- Second figure: WHERE the leading direction lives, vs max tau_a ----------
% The mechanistic form of the claim above. Per level of the sweep, the median
% over reps of the leading Lyapunov vector's norm^2 fractions in the x, SFA and
% STD blocks (stored per job by the top-K sweeps; see
% SRNNCellTypePairs.lya_summary). If the slowest SFA timescale sets lambda_1,
% the SFA fraction dominates until the STD mode takes over at the knee.
bf = gobjects(0);
param = 'tau_a_E';
has_blocks = isfield(psa.vector_param_lookup, param) && ...
    ~isempty(psa.results) && any(cellfun(@(c) has_field_in(psa.results, c, 'lead_frac_sfa'), ...
    cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false)));
if has_blocks
    lookup = psa.vector_param_lookup.(param);
    x_tau  = cellfun(@(v) v(end), lookup);
    cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
    bf = figure('Color', 'w', 'Position', [457 300 300 * numel(cond_names), 260]);
    blocks = {'lead_frac_sfa', 'lead_frac_std', 'lead_frac_x'};
    block_lab = {'SFA', 'STD', 'x'};
    block_col = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.30 0.30 0.30];
    for ci = 1:numel(cond_names)
        ax_b = subplot(1, numel(cond_names), ci, 'Parent', bf); hold(ax_b, 'on');
        Y = nan(numel(x_tau), numel(blocks));
        for li = 1:numel(x_tau)
            for bi = 1:numel(blocks)
                v = ParamSpaceAnalysis2.collect_level_values(psa, param, li, cond_names{ci}, blocks{bi});
                if ~isempty(v); Y(li, bi) = median(v); end
            end
        end
        for bi = 1:numel(blocks)
            plot(ax_b, x_tau, Y(:, bi), '-o', 'Color', block_col(bi, :), 'LineWidth', 2, ...
                'MarkerFaceColor', block_col(bi, :), 'MarkerSize', 4, 'DisplayName', block_lab{bi});
        end
        hold(ax_b, 'off');
        set(ax_b, 'XScale', 'log', 'FontSize', tick_fs); box(ax_b, 'off');
        ylim(ax_b, [0 1]);
        xlabel(ax_b, 'max $\tau_a$ (s)', 'Interpreter', 'latex', 'FontSize', label_fs);
        if ci == 1; ylabel(ax_b, 'leading-vector fraction', 'FontSize', label_fs); end
        if numel(cond_names) > 1
            title(ax_b, st.condition_title(cond_names{ci}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
        end
        if ci == 1; legend(ax_b, 'Location', 'east', 'FontSize', 10); end
    end
    if ~cfg.visible; set(bf, 'Visible', 'off'); end
else
    warning('fig_sfa_EOC_allStd:NoBlockFractions', ...
        ['%s has no lead_frac_* fields (a Benettin-era run?); the leading-vector ' ...
         'block-fraction figure is skipped.'], tau_dir);
end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_SFA_EOC_allStd';
out = struct('figs', [cf, bf], 'files', {{}}, 'source', tau_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, cf);
    out.files = existing_outputs(out_dir, fig_tag);
    if isgraphics(bf)
        save_figure_stable(out_dir, 'Fig_SFA_EOC_blocks', bf);
        out.files = [out.files, existing_outputs(out_dir, 'Fig_SFA_EOC_blocks')];
    end
end
end

function tf = has_field_in(results, cond, field)
tf = false;
if ~isfield(results, cond); return; end
R = results.(cond); R = R(~cellfun(@isempty, R));
tf = ~isempty(R) && isfield(R{1}, field);
end



