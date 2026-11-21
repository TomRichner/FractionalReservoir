% EXPLORE_OUTLIER_STRUCTURES Which kinds of outlier eigenvalue does adaptation rein in? The runs.
%
% THE QUESTION (Tom, 2026-11-20). explore_outliers took one example, the
% E-to-E mean times 1.5: a real eigenvalue far outside the bulk, carried by
% a whole population, and adaptation brings it back to the bulk edge. Is
% that special to a block-mean mismatch? What of the other routes, of
% low-rank structure that is not a population, of complex outliers, and of
% the outliers a balanced sparse Dale matrix has without any mean
% imbalance (Harris et al. 2023, Sec. III D 2)?
%
% THE HYPOTHESES (the lead's; outlier_mode_reference and
% outlier_mode_snapshot hold the algebra)
%   H1  Gains scale columns, so an eigenvalue moves by the participation-
%       weighted gain of the neurons that carry it. A uniform loss of gain
%       moves the bulk as much. A mode is reined in relative to the bulk
%       only if its carriers lose more gain than the average neuron.
%   H2  Adaptation and depression follow the rate, which is slow and
%       one-signed, so a real or slowly rotating mode is reined in and a
%       fast-rotating one (beyond about 1/0.25 s) is not.
%   H3  The local outliers of the balanced sparse Dale matrix are the H1
%       question with carriers that are not a population.
%
% THE NETWORKS, in priority order (outlier_network has the name grammar).
% Every one was checked on its structural spectrum first, and its strength
% set so that the intended eigenvalue is 1.3 to 2 structural bulk radii
% from the center on the right; the numbers are what
% outlier_mode_reference reports for the seed used.
%   Rotation sweep (H2): structure ADDED to the reference W, on random
%   patterns of neurons, dense, ignoring Dale's law (OutlierNet).
%     lrA1p5W0p0    rank one, real:      1.40 R
%     lrA1p4W0p15   rank two: (1.36 + 0.11i) R
%     lrA1p4W0p35             (1.37 + 0.36i) R
%     lrA1p4W0p6              (1.39 + 0.62i) R
%     lrA1p4W1p0              (1.41 + 1.02i) R
%     (the rank-one case needs a = 1.5 to land where the pairs do: the
%     balanced mean pulls a real eigenvalue in further)
%   Local outliers (H3), no change to the block means. Of seeds 1 to 150
%   of the reference network 26 have an eigenvalue beyond 1.2 R in
%   modulus and 9 have one on the right (Re > 0.7 R); at indegree 200
%   it is 113 and 38 (denser, not sparser, makes them more common: the
%   block-mean matrix grows against the bulk as sqrt(n alpha / (1 - alpha))).
%     reference, seed [27 28]   real, 1.33 R; the realized block-mean
%                               matrix has its eigenvalues at 0.37 R, so
%                               this is not a block-mean eigenvalue
%     ind200, seed [25 26]      complex, (1.16 + 0.64i) R, the rightmost
%     (reference seed 123 has the largest, 1.45 R, but its realized block
%     means are unbalanced enough to put an eigenvalue at 1.31 R, so it was
%     not used.) These modes have 40% of their eigenvector in the
%     population-uniform patterns: they are fluctuation-driven, and they
%     are still partly made of the balanced mean.
%   An E-I loop (complex, from the block means). Strengthening E-to-I and
%   I-to-E alone gives a pair on the imaginary axis of W (real part at the
%   CENTER of the bulk); E-to-E has to rise with it for a right-side pair.
%     EEx1p6-EIx1p33-IEx1p33    (1.26 + 1.15i) R
%     EEx1p6-EIx1p31-IEx1p31    (1.31 + 0.66i) R   (the slower loop)
%   Other single routes.
%     EIx0p9      I-to-E times 0.9:  +1.78 R (and -1.75 R)
%     IIx1p25     I-to-I times 1.25: +2.00 R (and -3.22 R): more
%                 inhibition of inhibition is a RIGHT-side outlier
%   Added after the first results (the 1.4 R modes come back to the edge
%   through the spread of the gains alone, so they could not tell whether
%   a non-population mode is reined in when that is not enough):
%     lrA2p5W0p0    rank one, real:      2.45 R
%     lrA2p5W1p0    rank two: (2.48 + 0.99i) R
%   Yesterday's two, rerun so that they carry the mode measures.
%     EEx1p5, reference (seed [7 8])
%
% WHAT IS RUN. outlier_run with the sampling of yesterday's dose runs
% (every 50 ms over 20 s after 60 s, 8 calibration states) plus the dense
% full-Jacobian pass at the calibration states, the three conditions, one
% seed. About 1, 1.7 and 3 minutes a run (no adaptation, one timescale,
% multiple) on 3 workers.
%
% HOW TO RUN IT. As explore_outliers: the script does as many missing runs
% as fit in call_budget_s and stops; call it again until it says every run
% is on disk. `networks` (an n x 2 cell of name and seed) can be set
% beforehand to run a subset.
%
% THEN: fig_outlier_structures (the figures, their sidecars, the table).
%
% See also: outlier_run, outlier_network, outlier_mode_measures, fig_outlier_structures

%% Settings
if ~exist('call_budget_s', 'var'); call_budget_s = 840; end
if ~exist('networks', 'var')
    networks = { ...
        'lrA1p5W0p0',  [7 8]; 'lrA1p4W0p15', [7 8]; 'lrA1p4W0p35', [7 8]; ...
        'lrA1p4W0p6',  [7 8]; 'lrA1p4W1p0',  [7 8]; ...
        'reference',   [27 28]; 'ind200',    [25 26]; ...
        'EEx1p6-EIx1p33-IEx1p33', [7 8]; ...
        'EEx1p5',      [7 8]; 'reference',   [7 8]; ...
        'EIx0p9',      [7 8]; 'IIx1p25',     [7 8]; ...
        'EEx1p6-EIx1p31-IEx1p31', [7 8]; ...
        'lrA2p5W0p0',  [7 8]; 'lrA2p5W1p0',  [7 8]};
end
if ~exist('conditions', 'var'); conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}; end

setup_paths();
out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'explorations', 'outlier_structures');
run_args = {'dt_sample', 0.05, 'n_cal', 8, 'dt_full', 20, 'n_kin', 1, 'out_dir', out_dir};

%% The runs, networks outermost so that the priority order is the order of completion
t_call = tic;
n_missing = 0; n_ran = 0;
for w = 1:size(networks, 1)
    for c = 1:numel(conditions)
        file = fullfile(out_dir, outlier_run_name(networks{w, 1}, conditions{c}, networks{w, 2}, ''));
        if isfile(file); continue; end
        if n_ran > 0 && toc(t_call) > call_budget_s
            n_missing = n_missing + 1;
            continue;
        end
        outlier_run(networks{w, 1}, conditions{c}, networks{w, 2}, run_args{:}, 'verbose', false);
        n_ran = n_ran + 1;
        fprintf('[explore_outlier_structures] %s | %s done (%.0f s into the call)\n', ...
            networks{w, 1}, conditions{c}, toc(t_call));
    end
end
if n_missing > 0
    fprintf('[explore_outlier_structures] %d run(s) made, %d still missing: call the script again.\n', n_ran, n_missing);
else
    fprintf('[explore_outlier_structures] every run is on disk (%d made in this call).\n', n_ran);
end
