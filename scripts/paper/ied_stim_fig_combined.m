function out = ied_stim_fig_combined(opts)
% IED_STIM_FIG_COMBINED Fig 8 (fig_stim_combined, mirror layout) as a figure of the ied_stim bundle (FD-031 item 5).
%
%   out = ied_stim_fig_combined('data_dir', run_dir, 'out_dir', d, ...
%           'model_fig', f, 'clinical_dir', c, ...)
%
% Registered in ied_stim_config; ied_stim_figures calls it after
% fig_ied_stim_blocks, whose .fig (in the run's own figure root) is
% model_fig. clinical_dir is the bundle's one clinical setting, a results
% folder of cscs_stim_stability's run_manuscript_stim_awake; the panels are
% read from its panels/ subfolder (clinical_files). clinical_psd_fig, when
% not empty, takes the PSD panel from that full path instead (testing only).
%
% The output holds clinical recordings: it lands under figs/ (gitignored);
% never force-add it to git. Replaces scripts/paper/make_fig_stim_combined*.m
% and their 'stage' switch (retired 2026-11-17).
%
% THE NEW ORDER (Tom, 2026-11-21; Plan 2, package E). With burst_mode_dir
% set (a run directory holding the burst_mode stage, burst_mode/), the
% figure is Figure 8 of figures_main_supp_appen.md: A to C the bursting
% network (bm_fig8_panels: input, rates and the network mode; one burst;
% the complex plane as rates), D to F the clinical panels (the mirror
% layout's A to C), G to I the modular model (its D to F). The mirror
% figure is drawn by fig_stim_combined as before, unsaved, and moved whole
% into the lower part of a taller canvas (stim_fig8_compose), so the
% clinical and modular panels are exactly the present figure's. Written as
% <tag>.{png,fig,svg} (tag Fig_Main8_Stim by default), <tag>_provenance.md
% and <tag>_table.md (the burst_mode numbers). With burst_mode_dir empty
% the present figure (Fig_Stim_Combined) is drawn, unchanged.
arguments
    opts.data_dir (1,:) char = ''       % unused; the standard field of a bundle entry
    opts.out_dir (1,:) char {mustBeNonempty}
    opts.visible (1,1) logical = false
    opts.model_fig (1,:) char {mustBeNonempty}
    opts.clinical_dir (1,:) char {mustBeNonempty}
    opts.clinical_files struct = struct('traces','Fig_Awake_Main_PanelD.fig', ...
        'rate','Fig_Awake_Main_PanelA_Normalized.fig','psd','Fig_Awake_Main_PanelB_PSD.fig')
    opts.clinical_psd_fig (1,:) char = ''
    opts.run_name (1,:) char = ''
    opts.sim struct = struct()
    opts.burst_mode_dir (1,:) char = ''  % a run directory with burst_mode/: the new order
    opts.tag (1,:) char = ''             % '' = Fig_Stim_Combined, or Fig_Main8_Stim in the new order
end
panels = fullfile(opts.clinical_dir, 'panels');
assert(isfolder(panels), 'ied_stim_fig_combined:NoPanels', ...
    'No clinical panels folder %s (cfg.clinical_dir in ied_stim_config).', panels);

[~, clin_run] = fileparts(strip_sep(opts.clinical_dir));
clinical_note = sprintf('cscs_stim_stability results/manuscript_stim/%s/panels%s', clin_run, repo_commit(opts.clinical_dir));
if ~isempty(opts.clinical_psd_fig)
    clinical_note = sprintf('%s; panel C PSD from %s', clinical_note, opts.clinical_psd_fig);
end
s = opts.sim;
if all(isfield(s, {'amplitude','seeds','n_no','n_stim','block_len','fs'}))
    model_note = sprintf(['IED model block run %s (ied_stim bundle), uniform input %.2f, %d networks, ' ...
        '%d + %d blocks of %g s each, fs %g Hz (fig_ied_stim_blocks_stats.md beside the model figure)'], ...
        opts.run_name, s.amplitude, numel(s.seeds), s.n_no, s.n_stim, s.block_len, s.fs);
else
    model_note = sprintf('IED model block run %s', opts.run_name);
end
common = {'layout', 'mirror', 'clinical_files', opts.clinical_files, 'clinical_psd_fig', opts.clinical_psd_fig, ...
    'clinical_note', clinical_note, 'model_note', model_note};
if isempty(opts.burst_mode_dir)
    tag = opts.tag; if isempty(tag); tag = 'Fig_Stim_Combined'; end
    r = fig_stim_combined(panels, opts.model_fig, 'out_dir', opts.out_dir, 'tag', tag, common{:}, 'visible', opts.visible);
    out = struct('files', {fullfile(opts.out_dir, r.files)}, 'out_dir', opts.out_dir);
    return
end
tag = opts.tag; if isempty(tag); tag = 'Fig_Main8_Stim'; end
B = bm_load_runs(opts.burst_mode_dir);          % errors (resolve_data_file) when the stage has not run there
r = fig_stim_combined(panels, opts.model_fig, common{:}, 'save', false, 'visible', false);
files = stim_fig8_compose(r, B, opts.out_dir, tag, clinical_note, model_note, opts.visible);
out = struct('files', {files}, 'out_dir', opts.out_dir);
end

function d = strip_sep(d)
while ~isempty(d) && any(d(end) == '/\'), d(end) = []; end
end
function s = repo_commit(d)
% ", commit <hash>" of the git repo holding d, or '' when there is none.
s = '';
[st, h] = system(sprintf('git -C "%s" log -1 --format=%%h', d));
h = strtrim(h);
if st == 0 && ~isempty(regexp(h, '^[0-9a-f]{7,}$', 'once')), s = [' (repo at commit ' h ')']; end
end
