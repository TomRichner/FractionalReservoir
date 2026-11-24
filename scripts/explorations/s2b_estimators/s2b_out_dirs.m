function [ddir, fdir, root] = s2b_out_dirs()
%S2B_OUT_DIRS data\explorations\s2b_estimators and figs\explorations\s2b_estimators of this worktree.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
ddir = fullfile(root, 'data', 'explorations', 's2b_estimators');
fdir = fullfile(root, 'figs', 'explorations', 's2b_estimators');
if ~isfolder(ddir); mkdir(ddir); end
if ~isfolder(fdir); mkdir(fdir); end
end
