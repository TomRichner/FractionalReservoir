function out = fig_doc_tables(cfg)
% FIG_DOC_TABLES The manuscript's generated tables, as a figure-registry entry.
%
%   out = FIG_DOC_TABLES()
%   out = FIG_DOC_TABLES('out_dir', d, 'preset_name', p)
%
% A thin adapter over write_manuscript_tables, giving it the uniform signature
% every registry entry has -- cfg in, out.figs/.files/.source back -- so
% make_all_paper_figures can run it in the same loop as everything else.
%
% WHY IT EXISTS. The tables used to be called separately at the end of
% make_all_paper_figures, inside their own try/catch, OUTSIDE the loop. That put
% their failures outside the headline count: when the n_a refactor broke
% write_manuscript_tables, the run reported 17/17 while the tables were throwing,
% and nothing in the summary said otherwise. As an ordinary entry a broken table
% reduces the total, which is the only place anyone looks.
%
% run_dir, save and visible are accepted and ignored. That is HONEST here, unlike
% the same comment on fig_eig_heatmap and fig_memory_capacity_example, which
% turned out to be a bug: those two had a run directory to read and did not read
% it. These tables are built from a live model constructed FROM THE PRESET --
% there is no run directory involved, and nothing to fix by wiring one in.
%
% See also: write_manuscript_tables, make_all_paper_figures, paper_config

arguments
    cfg.run_dir     (1,:) char    = ''    % unused; the tables describe the PRESET
    cfg.out_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.save        (1,1) logical = true  % unused; writing IS the work
    cfg.visible     (1,1) logical = true  % unused; no figures are drawn
end

setup_paths();

paths = write_manuscript_tables( ...
    'preset_name', cfg.preset_name, ...
    'out_dir',     cfg.out_dir, ...
    'verbose',     false);

% File NAMES, matching what the figures return via existing_outputs, so the
% master's "did anything land on disk" check and the manifest read alike.
files = cell(1, numel(paths));
for k = 1:numel(paths)
    [~, base, ext] = fileparts(paths{k});
    files{k} = [base ext];
end

% gobjects(0), not [] -- the master calls close() on out.figs, and an empty
% graphics array is what it expects from an entry that draws nothing.
out = struct('figs', gobjects(0), 'files', {files}, 'source', cfg.preset_name);
end
