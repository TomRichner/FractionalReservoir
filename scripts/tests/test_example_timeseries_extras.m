% TEST_EXAMPLE_TIMESERIES_EXTRAS The opt-in extra figures of fig_example_timeseries.
%
% show_stimulus and show_W_spectrum (2026-11-12) each add one figure after
% the main one. Checks, on the 100-neuron regression network with a short
% window: the default call makes exactly one figure; both options on make
% three valid figures; the stimulus figure's count of stimulated neurons
% (out.stimulus and the title) equals a count computed here from a model
% built the same way; the W spectrum's numbers are those of eig(params.W);
% and with the preset's own input (none) the stimulus figure still builds
% and says so. Nothing is saved.

preset = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
seeds  = [3 4];
T      = [0 4];
ic = struct('n_steps', 3, 'step_density', struct('E', 0.2, 'I', 0.2), 'amp', 0.5, ...
    'no_stim_pattern', [true false true], 'intrinsic_drive', 0, 'positive_only', false);
common = {'preset_name', preset, 'rng_seeds', seeds, 'T_range', T, 'conditions', {'no_adaptation'}, ...
    'save', false, 'visible', false, 'verbose', 'near-none'};

%% Default: one figure, no extra fields
out0 = fig_example_timeseries(common{:});
assert(isscalar(out0.figs) &&all(isgraphics(out0.figs)), 'default call must make one figure');
assert(~isfield(out0, 'stimulus') && ~isfield(out0, 'W_spectrum'), 'default call must add no fields');
close(out0.figs);

%% Both extras, with a stimulus
out = fig_example_timeseries(common{:}, 'model_args', {'input_config', ic}, ...
    'show_stimulus', true, 'show_W_spectrum', true);
assert(numel(out.figs) == 3 && all(isgraphics(out.figs, 'figure')), 'both options must make three figures');

m = build_from_preset(preset, 'no_adaptation', 'rng_seeds', seeds, 'T_range', T, 'fs', 400, ...
    'input_config', ic, 'verbose', 'near-none');
stim = any(m.u_ex ~= 0, 2);
n_E = nnz(stim(m.type_indices{1}));
n_I = nnz(stim(m.type_indices{2}));
assert(n_E + n_I > 0, 'test stimulus must stimulate someone');
assert(out.stimulus.n_stimulated == nnz(stim), 'stimulated count %d, expected %d', ...
    out.stimulus.n_stimulated, nnz(stim));
assert(isequal(out.stimulus.n_stimulated_per_type, [n_E n_I]), 'per-type counts differ');
ttl = strjoin(cellfun(@(c) strjoin(cellstr(c), ' '), ...
    get(findall(out.figs(2), 'Type', 'Text'), 'String'), 'UniformOutput', false), ' | ');
assert(contains(ttl, sprintf('E %d of %d', n_E, numel(m.type_indices{1}))) && ...
    contains(ttl, sprintf('I %d of %d', n_I, numel(m.type_indices{2}))), ...
    'stimulus figure title does not state the counts: %s', ttl);
iv = out.stimulus.on_intervals;
assert(isrow(iv) && abs(iv(1) - T(2) / 3) < 0.01 && abs(iv(2) - 2 * T(2) / 3) < 0.01, ...
    'on interval should be the middle third');

W = m.get_params().W;
lam = eig(full(W));
assert(abs(out.W_spectrum.abscissa - max(real(lam))) < 1e-10, 'abscissa');
assert(abs(out.W_spectrum.radius - max(abs(lam))) < 1e-10, 'radius');
assert(~out.W_spectrum.differs_across_conditions, 'one condition cannot differ');
close(out.figs);

%% No input (the preset's own): the stimulus figure builds and says so
out2 = fig_example_timeseries(common{:}, 'show_stimulus', true);
assert(numel(out2.figs) == 2 && all(isgraphics(out2.figs, 'figure')), 'no-input stimulus figure must build');
assert(out2.stimulus.n_stimulated == 0, 'no-input preset reports stimulated neurons');
txt = get(findall(out2.figs(2), 'Type', 'Text'), 'String');
if ~iscell(txt), txt = {txt}; end
assert(any(cellfun(@(s) contains(strjoin(cellstr(s), ' '), 'No external input'), txt)), ...
    'no-input figure must say so');
close(out2.figs);

fprintf('test_example_timeseries_extras: PASS\n');
