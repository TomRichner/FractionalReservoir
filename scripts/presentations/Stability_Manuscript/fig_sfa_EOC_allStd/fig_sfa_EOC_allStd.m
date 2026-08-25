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
project_root = fileparts(which('setup_paths'));
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

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_SFA_EOC_allStd';
out = struct('figs', cf, 'files', {{}}, 'source', tau_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, cf);
    out.files = existing_outputs(out_dir, fig_tag);
    capture_git_provenance(out_dir, project_root);

    [~, tau_name] = fileparts(tau_dir);
    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_sfa_EOC_allStd', ...
        'title',  'Stability_Manuscript figure: SFA edge of chaos (tau_a sensitivity)', ...
        'script', 'fig_sfa_EOC_allStd.m', ...
        'what',   ['Single panel: the largest Lyapunov exponent lambda_1 against ' ...
                   'the maximum SFA adaptation timescale tau_a (the last and ' ...
                   'largest of the log-spaced tau_a_E vector), under SFA + STD.'], ...
        'how',    ['Presentation replot -- no simulation is re-run. The saved ' ...
                   'tau_a_E_max PSA object is reloaded via ' ...
                   'ParamSpaceAnalysis2.from_dir and plot_sensitivity rebuilds ' ...
                   'the panel, which is then restyled: imagesc CLim capped at ' ...
                   'total_reps*clim_frac; colormap white to 90 percent black so ' ...
                   'the blue median stays visible; axis box removed; the ' ...
                   'condition title dropped; x relabelled to max tau_a (s).'], ...
        'source', struct('run_dir', run_dir, 'tau_subfolder', tau_name, ...
                         'preset', cfg.preset_name), ...
        'settings', struct('lle_range', lle_range, 'y_view', y_view, ...
                           'y_ticks', y_ticks, 'clim_frac', clim_frac, ...
                           'median_alpha', median_alpha), ...
        'figures', {out.files}, ...
        'notes',  ['THE HISTOGRAM RANGE AND Y-VIEW CROP REAL DATA, and are kept ' ...
                   'identical to the original figure by choice. On the aug_14 ' ...
                   'preset the tau LLEs spanned about -0.26 to +0.29 with a ' ...
                   'MEDIAN of +0.008, so slightly over half the distribution ' ...
                   'was positive and about a third of it sat above the y_view ' ...
                   'ceiling -- in the +inf overflow band, which is itself ' ...
                   'off-screen. The panel therefore reads as almost entirely ' ...
                   'sub-zero when it is not. Widen lle_range, y_view and ' ...
                   'y_ticks together to show the full spread. Re-check those ' ...
                   'numbers against the run actually plotted; they are ' ...
                   'preset-dependent.']));
end
end
