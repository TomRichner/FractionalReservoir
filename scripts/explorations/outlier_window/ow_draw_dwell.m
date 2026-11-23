% OW_DRAW_DWELL fig_outlier_dwell, unchanged, from the copied 200 s fastmedium run (10 seeds) and from the copied 20 s fast run (2 seeds).
%
% Output: figs/explorations/outlier_window/dwell_200s_fastmedium/ and
% figs/explorations/outlier_window/dwell_20s_fast/ (Fig_Outlier_Dwell.png and
% its _table.md). The run directories are copies under
% data/explorations/outlier_window, so nothing is written near the run.

setup_paths();
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dd = fullfile(root, 'data', 'explorations', 'outlier_window');
fd = fullfile(root, 'figs', 'explorations', 'outlier_window');
src = {fullfile(dd, 'fastmedium_copy'), fullfile(dd, 'fast_copy')};
dst = {fullfile(fd, 'dwell_200s_fastmedium'), fullfile(fd, 'dwell_20s_fast')};
for k = 1:2
    out = fig_outlier_dwell('run_dir', src{k}, 'out_dir', dst{k}, 'visible', false, 'verbose', 'verbose');
    close(out.figs);
    fprintf('%s\n', strjoin(fullfile(dst{k}, out.files), newline));
end
close all force
