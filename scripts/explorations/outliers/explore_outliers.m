% EXPLORE_OUTLIERS Does adaptation rein in outlier eigenvalues? The runs behind the measures and the figures.
%
% THE PROBLEM. Dale's law gives the weight matrix a low-rank mean, and when
% excitation and inhibition are not balanced that mean puts an eigenvalue
% far outside the bulk (E-to-E times 1.5: about +26 in units of W against
% a bulk radius of 5.6). Such a network should go unstable through that one
% population-wide mode. With adaptation its leading Lyapunov exponent is
% near zero all the same.
%
% THE QUESTION (Tom, 2026-11-19). Does adaptation rein the outliers in?
% Which eigenvalues leave the bulk, how far do they go, how long do they
% stay out, and what of those near the abscissa?
%
% THE HYPOTHESIS. Adaptation scales each neuron's outgoing weights by its
% gain (slope times available vesicles). The outlier of the scaled matrix
% follows the DIFFERENCE of the mean gains by type (the eigenvalues of
% K_eff), the bulk follows their rms. A uniform loss of gain, which is all
% saturation gives the network without adaptation, moves both together and
% leaves the outlier where the structure put it relative to the bulk.
% Adaptation driven by the very excess that makes the outlier can lower
% the gain of the over-driven population alone, and pull the outlier back
% to the edge of the bulk, with the bulk edge itself near the abscissa.
% Slower variables should then hold it there instead of letting it cycle
% across. (docs/future_directions/adaptation_reins_in_outlier_eigenvalues
% in the manuscript repository; nothing here is tuned towards it.)
%
% WHAT IS RUN. Each of `networks` under each of `conditions` on each of
% `seeds`, by outlier_run: 60 s to settle, then 20 s measured, the
% frozen-gain spectrum every 10 ms, null draws for the bulk edge at 12
% states and the whole full-Jacobian spectrum at the same states, the
% rightmost eigenvalues of the full Jacobian every 100 ms for as long as
% its time budget lasts, and the top 3 Lyapunov exponents over
% the same 20 s. One file per run under data/explorations/outliers/.
%
% HOW TO RUN IT. A run takes minutes, so the script does as many missing
% runs as fit in `call_budget_s` and stops; call it again until it says
% every run is done. Runs already on disk are skipped (delete a file to
% redo it). With `smoke = true` the same code runs on a 100-neuron network
% for 7 s into data/explorations/outliers_smoke/: a check that everything
% executes, not a result.
%
% POOL. The script opens nothing. With a pool up (pool_set in this
% session) the fast pass and the null use it; without one they run in the
% client.
%
% THE DOSE RUNS (added after the first results; not looped here). E-to-E
% times 1.1, 1.25 and 1.75 on seed [7 8], lighter sampling, into
% data/explorations/outliers_dose/, which fig_outliers reads when present:
%   outlier_run('EEx1p25', cond, [7 8], 'dt_sample', 0.05, 'n_cal', 8, ...
%       'dense', false, 'dt_full', 20, 'out_dir', dose_dir)
%
% THEN: outlier_table (the measures, one row per run) and fig_outliers
% (the figures and their sidecars).
%
% See also: outlier_run, outlier_measures, outlier_table, fig_outliers

%% Settings
if ~exist('smoke', 'var');         smoke = false;          end
if ~exist('call_budget_s', 'var'); call_budget_s = 420;    end
if ~exist('networks', 'var');      networks = {'reference', 'EEx1p5'}; end
if ~exist('conditions', 'var');    conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}; end
if ~exist('seeds', 'var');         seeds = {[7 8], [1 2]}; end

setup_paths();
root = fileparts(which('setup_paths'));
if smoke
    out_dir = fullfile(root, 'data', 'explorations', 'outliers_smoke');
    run_args = {'n', 100, 'T_settle', 5, 'T_measure', 2, 'n_cal', 2, 'n_null', 5, ...
        'n_shuffle', 5, 'full_budget_s', 20, 'n_kin', 1, 'out_dir', out_dir};
else
    out_dir = fullfile(root, 'data', 'explorations', 'outliers');
    run_args = {'out_dir', out_dir};
end

%% The runs, seeds outermost so the first seed is complete before the second starts
t_call = tic;
n_missing = 0; n_ran = 0;
for s = 1:numel(seeds)
    for w = 1:numel(networks)
        for c = 1:numel(conditions)
            file = fullfile(out_dir, outlier_run_name(networks{w}, conditions{c}, seeds{s}, ''));
            if isfile(file); continue; end
            if n_ran > 0 && toc(t_call) > call_budget_s
                n_missing = n_missing + 1;
                continue;
            end
            outlier_run(networks{w}, conditions{c}, seeds{s}, run_args{:});
            n_ran = n_ran + 1;
        end
    end
end
if n_missing > 0
    fprintf('[explore_outliers] %d run(s) made, %d still missing: call the script again.\n', n_ran, n_missing);
else
    fprintf('[explore_outliers] every run is on disk (%d made in this call).\n', n_ran);
    outlier_table(out_dir);
end
