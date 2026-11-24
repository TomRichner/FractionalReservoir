% OW_TIME_PASSES Time the two per-sample computations of outlier_spectra on saved states, one thread (as a pool worker).
%
% For the reference network under each condition (seed [7 8] of the copied
% fastmedium run): rebuild the model (ospec_rebuild_model, no integration),
% then at each of the run's 12 saved calibration states time
%   frozen   ospec_snapshot(S, params, ws, 3, 30)     (the fast pass, per sample)
%   full     ospec_full_spectrum(S, params, ws, 8, 0.1)  (the full pass, per sample)
% Writes data/explorations/outlier_window/ow_time_passes.md.

setup_paths();
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dd = fullfile(root, 'data', 'explorations', 'outlier_window');
src = fullfile(dd, 'fastmedium_copy', 'outlier_spectra');
conds = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
n_thr0 = maxNumCompThreads(1);
restore = onCleanup(@() maxNumCompThreads(n_thr0));
fid = fopen(fullfile(dd, 'ow_time_passes.md'), 'w');
fprintf(fid, ['# Per-sample cost of the two passes, measured here (one thread, reference network, seed [7 8], ' ...
    'the 12 calibration states)\n\n| condition | state variables | frozen-gain spectrum (s), median [min, max] | ' ...
    'full Jacobian rightmost modes (s), median [min, max] | ratio full / frozen |\n|---|---|---|---|---|\n']);
for c = 1:numel(conds)
    R = ospec_load_run(src, struct('network', 'reference', 'condition', conds{c}, 'seed', [7 8]));
    [~, params, ws] = ospec_rebuild_model(R);
    Sc = R.cal.S; N = R.N;
    clear R
    n_cal = size(Sc, 1);
    tf = zeros(1, n_cal); tj = zeros(1, n_cal);
    for k = 1:n_cal
        s = Sc(k, :)';
        t0 = tic; ospec_snapshot(s, params, ws, 3, 30); tf(k) = toc(t0);
        t0 = tic; ospec_full_spectrum(s, params, ws, 8, 0.1); tj(k) = toc(t0);
    end
    fprintf(fid, '| %s | %d | %.3f [%.3f, %.3f] | %.3f [%.3f, %.3f] | %.1f |\n', conds{c}, N, median(tf), min(tf), max(tf), ...
        median(tj), min(tj), max(tj), median(tj) / median(tf));
    fprintf('%s: N %d, frozen %.3f s, full %.3f s\n', conds{c}, N, median(tf), median(tj));
end
fclose(fid);
