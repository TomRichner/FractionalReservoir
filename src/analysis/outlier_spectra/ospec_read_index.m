function [runs, settings, folder] = ospec_read_index(run_dir)
%OSPEC_READ_INDEX The index of a run directory's outlier_spectra stage: one row per run it holds.
%
%   [runs, settings, folder] = ospec_read_index(run_dir)
%
%   run_dir   a bundle's run directory (the stage's folder is
%             <run_dir>/outlier_spectra), or that folder itself
%   runs      table, one row per run, in the order the stage's plan lists
%             them (seeds outermost):
%               network          cellstr, as ospec_network names it
%               condition        cellstr, the preset's condition name
%               seed1, seed2     the run's rng_seeds
%               file             the run's .mat file, relative to folder
%               has_full_series  the full-Jacobian series was computed
%               has_dense        the dense full-Jacobian spectra were computed
%               window_s         the measured window after settling (s)
%               dt_sample        the sampling step of the series (s)
%               ok               the run finished and its file is on disk
%               seconds          its wall-clock time
%               message          the error text of a run that failed ('' otherwise)
%   settings  the plan the stage ran (outlier_spectra_stage_plan), without
%             its list of runs
%   folder    the stage's folder, absolute
%
% Errors when there is no index (the stage did not run there) and when a
% row marked ok names a file that is not on disk.
%
% See also: ospec_load_run, run_outlier_spectra_stage, outlier_spectra_stage_plan

arguments
    run_dir (1,:) char
end
folder = run_dir;
if ~isfile(fullfile(folder, 'ospec_index.mat'))
    folder = fullfile(run_dir, 'outlier_spectra');
end
index_file = fullfile(folder, 'ospec_index.mat');
if ~isfile(index_file)
    error('ospec_read_index:NoIndex', ...
        'No ospec_index.mat in %s or its outlier_spectra folder: the outlier_spectra stage has not run there.', run_dir);
end
L = load(index_file, 'runs', 'settings');
runs = L.runs; settings = L.settings;
missing = runs.ok & ~cellfun(@(f) isfile(fullfile(folder, f)), runs.file);
if any(missing)
    error('ospec_read_index:MissingFile', 'The index of %s lists %d run file(s) that are not there, the first %s.', ...
        folder, nnz(missing), runs.file{find(missing, 1)});
end
end
