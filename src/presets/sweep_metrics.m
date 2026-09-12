function M = sweep_metrics(name)
% SWEEP_METRICS The one registry of per-job measures the sweeps store and plot.
%
%   M = SWEEP_METRICS()          every entry, a struct array
%   m = SWEEP_METRICS('hks')     one entry by key ...
%   m = SWEEP_METRICS('h_KS_bits')  ... or by result field name
%
% Each entry describes one scalar that ParamSpaceAnalysis2.run_single_job
% stores per job (see SRNNCellTypePairs.lya_summary for the Lyapunov ones)
% and how the plotters present it:
%
%   key          short name used in 'Metrics' option lists ('lle', 'r', ...)
%   field        the result struct field ('LLE', 'mean_rate', 'h_KS_bits', ...)
%   label        axis label (TeX)
%   stem         the token used in figure tags, e.g. Fig_Sensitivity_<stem>
%   dist_range   bin range for ParamSpaceAnalysis2.plot (the param-space
%                distribution) and the unit histograms
%   sens_range   bin range for plot_sensitivity's imagesc sheets
%   median_ylim  y window for the sensitivity-medians sheet
%   yticks       fixed y ticks on the sheets ([] = automatic)
%   zero_line    draw the zero reference (a sign that matters)
%   inf_both     overflow bins on both sides (else above only)
%   nan_means    what a NaN result means for this measure (for captions)
%   in_sheets    gets its own sensitivity / param-space / E:I figure
%
% WHY ONE TABLE. Until 2026-09-12 the labels, ranges and the LLE-only zero
% line were restated in five places (ParamSpaceAnalysis2.plot,
% .plot_sensitivity, .plot_unit_histograms, load_and_make_unit_histograms,
% and the metric_specs of two figures), two of them drifted copies of each
% other, so adding a measure meant editing all of them. Now a measure is a
% row here and the figures loop over the rows with in_sheets set.
%
% Ranges: h_KS in bit/s and D_KY are extensive (they grow with n), so their
% upper limits are overflow bins, not hard walls -- values past the range
% land in the last (inf) bin as LLE's already do.
%
% See also: SRNNCellTypePairs.lya_summary, ParamSpaceAnalysis2,
%           load_and_make_unit_histograms, fig_sensitivity_medians

E = @(key, field, label, stem, dist, sens, mylim, yt, zl, ib, nanm, sheet) struct( ...
    'key', key, 'field', field, 'label', label, 'stem', stem, ...
    'dist_range', dist, 'sens_range', sens, 'median_ylim', mylim, 'yticks', yt, ...
    'zero_line', zl, 'inf_both', ib, 'nan_means', nanm, 'in_sheets', sheet);

M = [ ...
    E('lle',  'LLE',            '\lambda_1 (1/s)',           'LLE',       [-1.5, 1.5],  [-2, 2],   [-1.75, 1.75], [],     true,  true,  'no estimate', true), ...
    E('r',    'mean_rate',      'Mean Firing Rate',          'mean_rate', [0, 1],       [0, 1],    [0, 1],        [0, 1], false, false, 'no run',      true), ...
    E('br',   'mean_synaptic_output', 'Mean Synaptic Output', 'synaptic_output', [0, 1], [0, 1],  [0, 1],        [0, 1], false, false, 'no run',      false), ...
    E('hks',  'h_KS_bits',      'h_{KS} (bit/s)',            'hKS',       [0, 20],      [0, 20],   [0, 20],       [],     false, false, 'no spectrum (Benettin run)', true), ...
    E('dky',  'D_KY',           'D_{KY}',                    'DKY',       [0, 30],      [0, 30],   [0, 30],       [],     false, false, 'unresolved within K', true), ...
    E('npos', 'n_positive',     'n_+ (positive exponents)',  'npos',      [0, 30],      [0, 30],   [0, 30],       [],     false, false, 'no spectrum', false), ...
    E('gap',  'lambda_gap',     '\lambda_1 - \lambda_2 (1/s)', 'gap',     [0, 2],       [0, 2],    [0, 2],        [],     false, false, 'no spectrum', false), ...
    E('fpos', 'frac_local_positive', 'fraction of time \lambda_1^{local} > 0', 'fpos', [0, 1], [0, 1], [0, 1],   [0, 1], false, false, 'no local series', false), ...
    E('p95',  'p95_finite_0p2s', 'p95 of 0.2 s finite-time \lambda_1 (1/s)', 'p95', [-2, 6], [-2, 6], [-2, 6],  [],     true,  true,  'no local series', false), ...
    E('exc',  'mean_positive_excursion_s', 'mean positive excursion (s)', 'exc', [0, 2], [0, 2],  [0, 2],        [],     false, false, 'never positive', false)];

if nargin == 0
    return;
end
hit = strcmpi({M.key}, name) | strcmp({M.field}, name);
if ~any(hit)
    error('sweep_metrics:UnknownMetric', ...
        'Unknown metric ''%s''. Keys: %s. Fields: %s.', name, ...
        strjoin({M.key}, ', '), strjoin({M.field}, ', '));
end
M = M(find(hit, 1));
end
