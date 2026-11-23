function out = fig_burst_mode(cfg)
%FIG_BURST_MODE Figure 8 G and H on their own: why the bursting network bursts, and why sustained input stops it.
%
%   out = FIG_BURST_MODE('run_dir', d)
%
% Draws what the opt-in stage burst_mode saved (its first seed: the
% bursting figures' network); simulates nothing. The panels are
% bm_fig8_panels' (it documents what each shows and how the five states of
% C are chosen); the combined Figure 8 (ied_stim_fig_combined) draws the
% same function into its top third. Here they stand alone, 7 in wide, at
% the combined figure's font sizes.
%
% Writes Fig_Burst_Mode.{png,svg,fig} and Fig_Burst_Mode_table.md: the
% network, the burst and states drawn, and the per-network numbers of the
% stage (bursts, crossings before bursts and their lead, the share of time
% beyond the boundary before and during the step, gains frozen and in the
% full Jacobian).
%
% See also: bm_fig8_panels, burst_mode_stage_plan, fig_supp10_grouped

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''      % the stage's burst_mode_index.mat, instead of run_dir
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''      % unused; accepted for a uniform call
    cfg.conditions                = {}      % unused: one bursting network under its most-adapted condition
end
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
B = bm_load_runs(cfg.run_dir, cfg.data_file, 'first');
R = B.R{1};
if isempty(R.intact.ev)
    error('fig_burst_mode:NoEigenvalues', ['The first network of %s kept no eigenvalues (the stage''s ' ...
        'keep_ev was 0); panel C needs them.'], B.index_file);
end
vis = 'on'; if ~cfg.visible; vis = 'off'; end
W = 7.0; Ht = 4.9;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 W Ht], 'Visible', vis);
P = bm_fig8_panels(fig, [0 0.05 W Ht - 0.3], R, 'left', 0.75);   % room for the input's name, written across
fs_letter = 12;
% Figure 8's letters since 2026-11-22: G the time course with its zoom, H the planes.
panel_letter(fig, P.ax_A(1), 'G', 0.02, fs_letter);
panel_letter(fig, P.ax_C(1), 'H', 0.02, fs_letter, 0.2);

fig_tag = 'Fig_Burst_Mode';
out = struct('figs', fig, 'files', {{}}, 'source', B.index_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    bm_write_table(fullfile(out_dir, [fig_tag '_table.md']), B, P.info, 'Figure 8 A to C');
    out.files{end + 1} = fullfile(out_dir, [fig_tag '_table.md']);
end
end

function panel_letter(fig, ax, ch, x_in, fs, raise)
if nargin < 6; raise = 0.04; end
u = ax.Units; ax.Units = 'inches'; p = ax.Position; ax.Units = u;
annotation(fig, 'textbox', 'Units', 'inches', 'Position', [x_in, p(2) + p(4) + raise, 0.4, 0.2], ...
    'String', ['(' ch ')'], 'FontSize', fs, 'FontName', 'Helvetica', 'LineStyle', 'none', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Margin', 0, 'Tag', 'letter');
end
