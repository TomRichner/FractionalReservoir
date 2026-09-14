function out = fig_eig_heatmap_imbalance(cfg)
% FIG_EIG_HEATMAP_IMBALANCE Jacobian-eigenvalue occupancy across E:I imbalance.
%
%   out = FIG_EIG_HEATMAP_IMBALANCE('run_dir', d)
%   out = FIG_EIG_HEATMAP_IMBALANCE('run_dir', d, 'density_scale', 'log')
%
% Rows = the examples run_eig_heatmap sampled (by default mu_EE_relative at
% 0.5x, 1x and 1.5x the preset's value: inhibition-dominant, reference,
% excitation-dominant), columns = adaptation regimes. Every example x regime
% shares one structural seed, so the picture is a DEFINED comparison rather
% than a chosen realization. Each panel is annotated with the numbers that
% make it interpretable: its matched finite-time lambda_1 (top-K, over the
% sampled window), the mean firing rate over the same window, and the
% realised E:I weight balance B_E = |sum W_E| / (|sum W_E| + |sum W_I|).
%
% Panels share square, density-trimmed axis limits (keep_frac of ALL pooled
% eigenvalues, the fig_eig_heatmap approach) and one colour scale, so a shift
% of the cloud between rows is a shift in the network, not in the axes. The
% zero-real-part line is drawn on every panel; it marks LOCAL instability of
% the instantaneous Jacobian, and the panel's lambda_1 -- not the cloud -- is
% the stability statement.
%
% density_scale as in fig_eig_heatmap: 'log' -> log10(1 + D), 'loglog' ->
% log10(1 + log10(1 + D)) (default here: the dense core saturates otherwise).
%
% A markdown table beside the figure lists, per panel, the override value,
% lambda_1, the mean rate, B_E and the median spectral and numerical
% abscissae of the dendritic block.
%
% Errors fig_eig_heatmap_imbalance:NoExamples on a .mat from before the
% examples existed (run run_eig_heatmap again).
%
% See also: run_eig_heatmap, fig_eig_heatmap, fig_EI_weights_param_space

arguments
    cfg.verbose       (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file     (1,:) char    = ''
    cfg.out_dir       (1,:) char    = ''
    cfg.save          (1,1) logical = true
    cfg.visible       (1,1) logical = true
    cfg.run_dir       (1,:) char    = ''
    cfg.preset_name   (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.density_scale (1,:) char {mustBeMember(cfg.density_scale, {'log', 'loglog'})} = 'loglog'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'eig_heatmap')}, ...
    'eig_heatmap_data.mat', ...
    'Run run_eig_heatmap first');
D = load(data_file);
if ~isfield(D, 'examples') || isempty(D.examples)
    error('fig_eig_heatmap_imbalance:NoExamples', ...
        ['%s has no `examples` (it predates the imbalance examples, 2026-09-14). ' ...
         'Run run_eig_heatmap again.'], data_file);
end
ex         = D.examples;
titles     = D.condition_titles;
cond_names = D.cond_names;
n_ex       = numel(ex);
n_cond     = numel(titles);
lle_window = D.lle_window;

%% ---- Heatmap parameters (as fig_eig_heatmap) --------------------------------
grid_res   = 250;
sigma_bins = 1.25;
keep_frac  = 0.999;
switch cfg.density_scale
    case 'log'
        scale_fn = @(M) log10(1 + M);
        cb_label = 'log_{10}(1 + eigenvalue density)';
    case 'loglog'
        scale_fn = @(M) log10(1 + log10(1 + M));
        cb_label = 'log_{10}(1 + log_{10}(1 + eigenvalue density))';
end

%% ---- Shared square limits from EVERY panel's eigenvalues ------------------
all_evals = [];
for e = 1:n_ex
    all_evals = [all_evals; vertcat(ex(e).evals_by_cond{:})]; %#ok<AGROW>
end
re = real(all_evals); im = imag(all_evals);
tail   = (1 - keep_frac) / 2;
re_lo  = quantile(re, tail); re_hi = quantile(re, 1 - tail);
re_ctr = (re_lo + re_hi) / 2; re_span = re_hi - re_lo;
im_span = 2 * quantile(abs(im), keep_frac);
span   = max(re_span, im_span);
re_lim = re_ctr + [-0.5, 0.5] * span;
im_lim = [-0.5, 0.5] * span;
re_edges = linspace(re_lim(1), re_lim(2), grid_res + 1);
im_edges = linspace(im_lim(1), im_lim(2), grid_res + 1);

Dm = cell(n_ex, n_cond);
cmax = 0;
for e = 1:n_ex
    for i = 1:n_cond
        Dm{e, i} = scale_fn(compute_eigenvalue_density(ex(e).evals_by_cond{i}, re_edges, im_edges, sigma_bins));
        cmax = max(cmax, max(Dm{e, i}(:)));
    end
end
clim = [0, max(cmax, 0.5)];   % the colorbar is ticked at 0, 0.25, 0.5

%% ---- Figure ---------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [120, 80, 400 * n_cond + 60, 360 * n_ex]);
tl  = tiledlayout(fig, n_ex, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
ax  = gobjects(n_ex, n_cond);
rows = cell(1, n_ex * n_cond);
k = 0;
for e = 1:n_ex
    for i = 1:n_cond
        ax(e, i) = nexttile(tl);
        plot_eigenvalue_heatmap_helper(ax(e, i), Dm{e, i}, re_edges, im_edges, clim, false);
        if e == 1
            title(ax(e, i), titles{i}, 'FontWeight', 'normal', 'FontSize', st.title_fs);
        end
        if i == 1
            ylabel(ax(e, i), {row_label(ex(e)); 'Im \lambda'}, 'FontSize', st.label_fs);
        else
            ylabel(ax(e, i), '');
        end
        if e == n_ex
            xlabel(ax(e, i), 'Re \lambda', 'FontSize', st.label_fs);
        else
            xlabel(ax(e, i), '');
        end
        text(ax(e, i), 0.03, 0.96, sprintf('\\lambda_1 = %+.3f\n\\langle r\\rangle = %.2f\nB_E = %.2f', ...
            ex(e).lle_by_cond(i), ex(e).mean_rate_by_cond(i), ex(e).B_E_by_cond(i)), ...
            'Units', 'normalized', 'Color', [0.25 0.25 0.25], 'FontSize', 10, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
        k = k + 1;
        rows{k} = sprintf('| %s | %s | %s | %+.4f | %.3f | %.3f | %+.3f | %+.3f | %d |', ...
            ex(e).label, override_text(ex(e)), cond_names{i}, ex(e).lle_by_cond(i), ...
            ex(e).mean_rate_by_cond(i), ex(e).B_E_by_cond(i), ...
            median(ex(e).spec_abscissa_by_cond{i}), median(ex(e).num_abscissa_by_cond{i}), ...
            numel(ex(e).evals_by_cond{i}));
    end
end
cb = colorbar(ax(end, end));
cb.Layout.Tile = 'east';
cb.Label.String = cb_label;
cb.Ticks = [0 0.25 0.5];
cb.Box = 'off';
title(tl, {'Jacobian eigenvalue occupancy across E:I imbalance and adaptation regime', ...
    sprintf('one structural seed; \\lambda_1 = finite-time top-K exponent over the last %g s; \\langle r\\rangle mean rate; B_E = excitatory share of summed weight', lle_window)}, ...
    'FontWeight', 'normal', 'FontSize', 11);

hdr = '| Example | Override | Condition | lambda_1 | mean rate | B_E | median alpha(J_xx) | median omega(J_xx) | eigenvalues |';
sep = '|---|---|---|---|---|---|---|---|---|';
vprintf(cfg.verbose, 'verbose', '%s\n%s\n%s\n', hdr, sep, strjoin(rows, newline));

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% ---- Save -----------------------------------------------------------------
fig_tag = 'Fig_Eig_Heatmap_Imbalance';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '# Jacobian occupancy examples\n\nSource: `%s`\n\nn (reference) = %d; non-reference examples at n = %d; %s.\n\n', ...
            data_file, D.settings.n, ex_n_nonref(ex), density_note(cfg.density_scale));
        fprintf(fid, '%s\n%s\n%s\n', hdr, sep, strjoin(rows, newline));
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function s = row_label(e)
% Row label from the override: '\mu_{EE} x 0.5' style, or 'reference'.
if isempty(fieldnames(e.overrides))
    s = 'reference';
    return;
end
f = fieldnames(e.overrides);
parts = cell(1, numel(f));
for k = 1:numel(f)
    v = e.overrides.(f{k});
    parts{k} = sprintf('%s = %s', tex_name(f{k}), mat2str(v, 3));
end
s = sprintf('%s (%s)', strrep(e.label, '-', ' '), strjoin(parts, ', '));
end

function s = override_text(e)
if isempty(fieldnames(e.overrides))
    s = 'none';
    return;
end
f = fieldnames(e.overrides);
parts = cell(1, numel(f));
for k = 1:numel(f)
    parts{k} = sprintf('%s = %s', f{k}, mat2str(e.overrides.(f{k}), 4));
end
s = strjoin(parts, '; ');
end

function s = tex_name(name)
switch name
    case 'mu_EE_relative', s = '\mu_{EE}';
    case 'mu_EI_relative', s = '\mu_{EI}';
    case 'mu_IE_relative', s = '\mu_{IE}';
    case 'mu_II_relative', s = '\mu_{II}';
    otherwise,             s = strrep(name, '_', '\_');
end
end

function n = ex_n_nonref(ex)
nn = [ex(~strcmp({ex.label}, 'reference')).n];
if isempty(nn); n = NaN; else; n = nn(1); end
end

function s = density_note(scale)
switch scale
    case 'log',    s = 'colour = log10(1 + density)';
    case 'loglog', s = 'colour = log10(1 + log10(1 + density))';
end
end
