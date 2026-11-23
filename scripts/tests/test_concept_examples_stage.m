%% test_concept_examples_stage.m
% The opt-in stage 'concept_examples' (Plan 2, package B) and its two
% figures, Figure 1 (fig_main1_concept_networks) and S3
% (fig_supp3_connectivity_gallery), on tiny networks.
%
% Checks:
%   1. The plan from mu7chain_AC_config: the four default examples and the
%      eight gallery structures with their levels and seeds, the knobs'
%      defaults, the size at 500; bad settings refused (unknown knob, odd
%      n_image, lya_from after T, bad network name, a preset without
%      'no_adaptation' is not tested: every paper preset has it).
%   2. END TO END at 80 neurons, 2 s: three examples and three gallery
%      structures. The file's fields; the spectrum is (eig(W) - 1)/tau_d of
%      the rebuilt network; the weight piece is the submatrix of W on its
%      neurons, half E and half I; the kept traces are the model's states
%      (rates through the class's nonlinearity in [0, 1]); the E-to-E x 1.5
%      network's rightmost eigenvalue is real and right of its bulk; the
%      low-rank network's added mode is where it was built; lambda_1 is
%      finite.
%   3. Both figures draw from the run (PNG and _table.md written, the
%      far eigenvalue beyond the break counted), and error on an empty run
%      directory.
%   4. fig_main1_concepts still draws its bottom row through
%      draw_adaptation_concepts (test_chain_display_analytic checks its
%      numbers).
%
% Asserts; prints one PASS line. About 20 s. Assumes setup_paths has run.
%
% See also: concept_examples_stage_plan, run_concept_examples_stage, cx_load

%% 1. The plan
cfg = mu7chain_AC_config('fast');
cfg.verbose = 'near-none';
pl = concept_examples_stage_plan(cfg);
assert(strcmp(pl.folder, 'concept_examples') && strcmp(pl.condition, 'no_adaptation'), 'folder, condition');
assert(isequal({pl.examples.key}, {'decay', 'irregular', 'runaway', 'oscillation'}), 'the four examples');
assert(isequal({pl.examples.name}, {'reference', 'reference', 'EEx1p5', 'lrA2p5W1p0'}), 'their networks');
assert(isequal([pl.examples.level_of_chaos], [0.15 1 1 1]), 'their levels');
assert(numel(pl.gallery) == 8 && isequal(pl.gallery(6).seed, [27 28]) && isequal(pl.gallery(1).seed, [7 8]), 'the gallery');
% The gallery shows structures the outlier_structures stage tests (S7 panel
% B summarizes them): the same names and weight seeds, and no 'E-to-E down',
% which that stage does not run.
assert(isequal({pl.gallery.name}, {'reference', 'EEx1p5', 'EIx0p9', 'IIx1p25', 'EEx1p6-EIx1p33-IEx1p33', 'reference', ...
    'lrA1p5W0p0', 'lrA1p4W1p0'}), 'the gallery''s networks');
[g_titles, g_keys] = cx_gallery_names(pl.gallery, pl.seed);
po = outlier_structures_stage_plan(cfg);
assert(all(ismember(g_keys, {po.networks.key})) && strcmp(g_keys{6}, 'reference_w27'), 'every gallery structure is a network of outlier_structures');
assert(isequal(g_titles, {'Reference', ['E-to-E ' char(215) ' 1.5'], 'I-to-E down', 'I-to-I up', 'E-I loop', 'Local outlier', 'Real low-rank mode', ...
    'Rotating low-rank mode'}), 'the tiles'' titles');
assert(pl.T == 12 && pl.dt_sample == 0.01 && pl.n_trace == 5 && pl.n_image == 30 && pl.lya && pl.n_override == 0, 'knobs');
bad = {{'concept_examples_bogus', 1}, {'concept_examples_n_image', 31}, {'concept_examples_lya_from', 20}, ...
    {'concept_examples_examples', {'a', 'notAName', 1, []}}, {'concept_examples_dt_sample', 0.0031}};
for k = 1:numel(bad)
    c = cfg; c.(bad{k}{1}) = bad{k}{2};
    try
        concept_examples_stage_plan(c);
        error('test:NotRefused', 'Not refused: %s', bad{k}{1});
    catch err
        assert(strcmp(err.identifier, 'concept_examples_stage_plan:BadConfig'), 'refused %s: %s', bad{k}{1}, err.message);
    end
end
c = cfg; c.concept_examples_gallery = {};
assert(isempty(concept_examples_stage_plan(c).gallery), 'an empty gallery is a choice');

%% 2. End to end, tiny
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
c = cfg;
c.concept_examples_n_override = 80;
c.concept_examples_T = 2;
c.concept_examples_lya_from = 1;
c.concept_examples_n_trace = 3;
c.concept_examples_n_image = 20;
c.concept_examples_examples = {'decay', 'reference', 0.15, []; 'runaway', 'EEx1p5', 1, []; 'oscillation', 'lrA2p5W1p0', 1, []};
c.concept_examples_gallery = {'balanced', 'reference', 1, []; 'EE_down', 'EEx0p9', 1, []; 'lowrank_rot', 'lrA1p4W1p0', 1, [3 4]};
pl = concept_examples_stage_plan(c);
run_concept_examples_stage(c, tmp, pl);
C = cx_load(tmp);
assert(numel(C.examples) == 3 && numel(C.gallery) == 3, 'both lists saved');
nt = round(2 / 0.01) + 1;
for k = 1:3
    x = C.examples(k);
    assert(x.n == 80 && x.n_E == 40 && numel(x.ev) == 80, 'size');
    assert(numel(x.t) == nt && isequal(size(x.x_trace), [6 nt]) && isequal(x.trace_type', [1 1 1 2 2 2]), 'traces');
    assert(all(x.rate_trace(:) >= 0 & x.rate_trace(:) <= 1), 'rates in [0, 1]');
    assert(isfinite(x.LLE) && isequal(x.lya_interval, [1 2]), 'lambda_1');
    assert(isequal(size(x.sample.W), [20 20]) && isequal(x.sample.type', [ones(1, 10), 2 * ones(1, 10)]), 'the piece');
end
% the spectrum and the piece against the rebuilt network
e = pl.examples(2);
[model, params, ws] = cx_network(pl, e);
ev = (eig(ws.Wf) - 1) / model.tau_d;
x = C.examples(2);
assert(max(abs(sort(ev, 'ComparisonMethod', 'real') - sort(x.ev, 'ComparisonMethod', 'real'))) < 1e-8, 'spectrum');
assert(isequal(x.sample.W, ws.Wf(x.sample.idx, x.sample.idx)), 'the piece is the submatrix');
assert(all(ws.type_of(x.sample.idx) == x.sample.type), 'its types');
assert(abs(imag(x.rightmost)) < 1e-9 && real(x.rightmost) > x.center + x.R_rate, 'E-to-E x 1.5: a real outlier');
assert(abs(x.mean_rate_E - 1) < 0.05, 'E-to-E x 1.5: E runs to the ceiling');
assert(x.center == -1 / params.tau_d, 'center');
xo = C.examples(3);
assert(abs(min(abs(xo.ev - xo.lowrank_rate))) < 0.15 * abs(xo.lowrank_rate - xo.center), 'the added mode is where it was built');
assert(isnan(C.examples(1).lowrank_rate) && isfield(C.gallery, 'sample') && ~isfield(C.gallery, 'LLE'), 'gallery: structure only');
assert(isequal(C.gallery(3).seed, [3 4]), 'a gallery row with its own weight seed');

%% 3. The figures
fd = fullfile(tmp, 'figs');
o1 = fig_main1_concept_networks('run_dir', tmp, 'out_dir', fullfile(fd, 'f1'), 'visible', false, 't_window', [0 2], ...
    'n_E', 3, 'n_I', 3);
assert(any(strcmp(o1.files, 'Fig_Main1_Concept_Networks.png')) && any(strcmp(o1.files, 'Fig_Main1_Concept_Networks_table.md')), 'Figure 1 written');
assert(o1.n_far(2) >= 1, 'the E-to-E outlier beyond the break');
tab = fileread(fullfile(fd, 'f1', 'Fig_Main1_Concept_Networks_table.md'));
assert(contains(tab, 'EEx1p5') && contains(tab, 'lambda_1'), 'its table');
% Tom, 2026-11-22: eigenvalues left of the boundary gray, right of it
% black, an outlier far to the right larger and green; what is not drawn is
% said in the table; the diagrams' axis names and the added mode's title in
% sentence case; letters by row (A to D); the arrows Increasing SFA and
% Increasing STD in row D. Tom, 2026-11-23: row C back to the potential x_i
% (the runaway column with a y range of its own), at most four y ticks;
% row A as S3 (balanced + added where something was added); arrows on the
% discs too.
n_dots = @(hs) sum(arrayfun(@(h) numel(h.XData), hs));
for k = 1:numel(o1.axes.spectrum)
    ax = o1.axes.spectrum(k);
    hc = findobj(ax, 'Tag', 'eig_crossed'); hl = findobj(ax, 'Tag', 'eig_left'); ho = findobj(ax, 'Tag', 'eig_outlier');
    ev_k = C.examples(k).ev;
    assert(o1.n_crossed(k) == nnz(real(ev_k) > 0), 'column %d: the count right of the boundary', k);
    assert(n_dots(hc) + n_dots(hl) + n_dots(ho) == numel(ev_k) - o1.n_left(k) && n_dots(ho) == o1.n_out(k), ...
        'column %d: every eigenvalue the window holds drawn once, the outliers on their own', k);
    assert(n_dots(hc) <= nnz(real(ev_k) > 0), 'column %d: black only right of the boundary', k);
    assert(all(arrayfun(@(h) isequal(h.CData(1, :), [0 0 0]), hc)) && ...
        all(arrayfun(@(h) h.CData(1, 1) > 0 && h.CData(1, 1) < 1 && all(h.CData(1, :) == h.CData(1, 1)), hl)) && ...
        all(arrayfun(@(h) h.CData(1, 2) > h.CData(1, 1) && h.CData(1, 2) > h.CData(1, 3), ho)), ...
        'column %d: black right of the boundary, gray left of it, the outliers green', k);
end
assert(strcmp(o1.axes.activity(1).YLabel.String, 'x_i') && all(arrayfun(@(ax) numel(ax.YTick) <= 4, o1.axes.activity)), ...
    'row C: the potential x_i, at most four ticks');
assert(isequal(isgraphics(o1.axes.added), o1.add_clim > 0) && ...
    numel(findall(o1.figs, 'Tag', 'plus')) == nnz(o1.add_clim > 0), 'row A: an added image and a plus sign where something was added');
assert(~isempty(findall(o1.figs, 'Tag', 'arrow_sfa_discs')) && ~isempty(findall(o1.figs, 'Tag', 'text_std_discs')), 'row D: the arrows on the discs');
assert(all(arrayfun(@(L) ~isempty(findall(o1.figs, 'Tag', ['letter_' L])), 'ABCD')), 'letters by row, A to D');
assert(~isempty(findall(o1.figs, 'Tag', 'text_sfa')) && ~isempty(findall(o1.figs, 'Tag', 'arrow_std')), 'row D: the arrows');
assert(contains(tab, 'Eigenvalues right of the boundary') && ...
    (contains(tab, 'Not drawn (left of the window') == any(o1.n_left > 0)) && ...
    (contains(tab, 'Every eigenvalue is drawn') == ~any(o1.n_left > 0)), 'the table says what is not drawn');
assert(strcmp(o1.axes.concepts(1).XLabel.String, 'Dendritic potential') && ...
    strcmp(o1.axes.concepts(3).YLabel.String, 'Synaptic output'), 'panel B: axis names in sentence case');
assert(any(strcmp(o1.titles, 'Added rotating mode')) == any(strcmp({C.examples.key}, 'oscillation')), 'the added mode''s title');
close(o1.figs);
o3 = fig_supp3_connectivity_gallery('run_dir', tmp, 'out_dir', fullfile(fd, 's3'), 'visible', false);
assert(any(strcmp(o3.files, 'Fig_Supp3_Connectivity_Gallery.png')) && numel(o3.axes.spectrum) == 3, 'S3 written');
% S3 (Tom, 2026-11-22): balanced plus added, the added matrix on its own
% color range; a tile with nothing added shows one image.
has_add = o3.add_clim > 0;
assert(isequal(isgraphics(o3.axes.added), has_add) && ~has_add(1) && any(has_add), ...
    'S3: an added image exactly where something was added, none on the balanced tile');
assert(numel(findall(o3.figs, 'Tag', 'plus')) == nnz(has_add), 'S3: one plus sign per added image');
close(o3.figs);
empty = fullfile(tmp, 'empty'); mkdir(empty);
for f = {@fig_main1_concept_networks, @fig_supp3_connectivity_gallery}
    try
        f{1}('run_dir', empty, 'save', false, 'visible', false);
        error('test:NoError', '%s drew from an empty run directory', func2str(f{1}));
    catch err
        assert(~strcmp(err.identifier, 'test:NoError'), err.message);
    end
end

%% 4. fig_main1_concepts through the shared helper
intro = fullfile(tmp, 'intro');
mkdir(fullfile(intro, 'eigenspectra')); mkdir(fullfile(intro, 'statetraces'));
for f = {fullfile(intro, 'eigenspectra', 'panelA_eigenspectrum.fig'), fullfile(intro, 'statetraces', 'panelA_bottom_traces.fig')}
    h = figure('Visible', 'off');
    for k = 1:3; subplot(1, 3, k, 'Parent', h); end
    savefig(h, f{1}); close(h);
end
m1 = fig_main1_concepts(cfg.preset_name, intro);
tags = get(findall(m1.figs, 'Type', 'axes'), 'Tag');
assert(all(ismember({'sfa_sigmoid', 'sfa_discs', 'std_sigmoid', 'std_discs'}, tags)), 'the bottom row');
close(m1.figs);

fprintf('PASS test_concept_examples_stage\n');
