%% test_composite_spectral_figures.m
% The composite figures of the figure plan of 2026-11-21 that are built
% from the spectral stages (Plan 2, package D): fig_main5_mechanism,
% fig_main6_amplification, fig_supp6_grouped, fig_supp7_grouped,
% fig_supp8_grouped, their helpers (composite_figure, composite_letter,
% composite_part, composite_save), the options added to the source figures
% for them, and the export of a figure taller than the screen.
%
% Checks:
%   1. A FIGURE LARGER THAN THE SCREEN IS SAVED WHOLE, twice in a row
%      (save_figure_stable): the PNG has the figure's own proportions, the
%      mark drawn at its very top and the one at its very bottom are both
%      in it, the second save equals the first in size, and the figure
%      keeps the size it was laid out for.
%   2. ONE tiny run directory with the real stages: outlier_spectra (60
%      neurons, 2 seeds, reference, E-to-E x 1.25 and x 1.5, three
%      conditions), eigenvalue_density and frequency_response on its
%      output, outlier_structures (80 neurons, five structures), and, for
%      panels A and B of Figure 6, MADE-UP eig_heatmap and transient_gain
%      files with the fields the figure reads (those stages have their own
%      tests; here only the drawing is exercised).
%   3. Each composite draws from it: one figure, invisible, its PNG and its
%      _table.md written (the caption, the data files, one section per
%      source figure), the panel letters the plan names, and nothing left
%      open but the figure returned. With save false nothing is written.
%   4. What the plan asks of the panels: Figure 5 has the E and I marks on
%      every inset and names the drawn seed (the stage's first) in its
%      sidecar; Figure 6's frequency panels are the reference network only,
%      carry their own frequency-axis labels, and leave out the middle
%      panel; S6 draws one cumulative panel per network, titled.
%   5. The source figures' new options: fig_freqresp_network networks and
%      panels (the table still has every network's row it was asked for),
%      fig_freqresp_plane networks and sentence-case titles, a network the
%      stage lacks is an error; fig_eigdens_radial panels; `images` false
%      writes the sidecar and no image; fig_ostruct_summary keeps negative
%      detachments on its axis.
%   6. Each composite errors on an empty run directory (resolve_data_file).
%
% Asserts; prints one PASS line. A few minutes (the stages). Assumes
% setup_paths has run and a pool allocation is set (pool_set(0) will do).
%
% NOT TESTED HERE: how the pictures look at 500 neurons (looked at by eye
% on data/narrative_essential_spectral, 2026-11-21), the composites through
% make_all_paper_figures (wired by package A), and panels A and B of
% Figure 6 on real stage output (drawn by eye from data/mu7chain_AC_fastmedium).
%
% See also: composite_figure, fig_main5_mechanism, fig_main6_amplification,
%           test_outlier_figures, test_frequency_response_stage

%% 1. A figure larger than the screen is saved whole
scr = get(groot, 'ScreenSize'); ppi = get(groot, 'ScreenPixelsPerInch');
tall = [8, max(16, 1.4 * scr(4) / ppi)];                   % inches; well above the screen's height
big_dir = tempname; cleanup_big = onCleanup(@() rmdir_if(big_dir));
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 tall], 'Visible', 'off');
ax = axes(fig, 'Units', 'inches', 'Position', [0 0 tall], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 tall(2)]);
hold(ax, 'on');
plot(ax, [0.1 0.9], tall(2) - [0.1 0.1], 'k-', 'LineWidth', 4);    % a tenth of an inch from the top
plot(ax, [0.1 0.9], [0.1 0.1], 'k-', 'LineWidth', 4);              % and from the bottom
sizes = zeros(2, 2);
for pass = 1:2
    save_figure_stable(big_dir, 'Tall', fig);
    A = imread(fullfile(big_dir, 'Tall.png'));
    sizes(pass, :) = [size(A, 1), size(A, 2)];
    ink = find(any(any(A < 128, 3), 2));
    per_in = size(A, 1) / tall(2);
    assert(abs(size(A, 1) / size(A, 2) - tall(2) / tall(1)) < 0.01, 'pass %d: the PNG has the figure''s proportions', pass);
    assert(~isempty(ink) && ink(1) < 0.25 * per_in && ink(end) > size(A, 1) - 0.25 * per_in, ...
        'pass %d: the marks at the top and at the bottom are both in the PNG (ink rows %d to %d of %d)', pass, ink(1), ink(end), size(A, 1));
    assert(max(abs(fig.Position(3:4) - tall)) < 1e-6, 'pass %d: the figure keeps its size', pass);
    svg = fileread(fullfile(big_dir, 'Tall.svg'));
    hgt = regexp(svg, 'height="([\d.]+)mm"', 'tokens', 'once');
    assert(~isempty(hgt) && abs(str2double(hgt{1}) / 25.4 - tall(2)) < 0.05, 'pass %d: the SVG is the figure''s height', pass);
    assert(numel(strfind(svg, '<polyline')) == 2, 'pass %d: the SVG has the top mark and the bottom mark', pass);
    assert(~contains(svg, '#010203') && ~contains(svg, '#030201') && isempty(findall(fig, 'Tag', 'print_whole_check')) && ...
        isscalar(findall(fig, 'Type', 'axes')), 'pass %d: the check squares are gone from the file and from the figure', pass);
    corner = A(1:20, 1:20, :);
    assert(all(corner(:) == 255), 'pass %d: the check square is painted out of the PNG', pass);
end
assert(isequal(sizes(1, :), sizes(2, :)), 'the second save is the size of the first');
close(fig);

%% 2. One tiny run directory
tmp = tempname; empty_run = tempname; fig_dir = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {tmp, empty_run, fig_dir}));
mkdir(empty_run);
conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
t0 = tic;
c = mu7chain_config('essential');
c.verbose = 'near-none';
c.outlier_spectra_n_override = 60;
c.outlier_spectra_networks = {'reference', 'EEx1p25', 'EEx1p5'};
c.outlier_spectra_full_networks = {'reference', 'EEx1p5'};
c.outlier_spectra_conditions = conditions;
c.outlier_spectra_seeds = [1 2; 2 3];
c.outlier_spectra_n_seeds = 2;
c.outlier_spectra_T_settle = 2;
c.outlier_spectra_window_s = 2;
c.outlier_spectra_full_window_s = 2;
c.outlier_spectra_lya_window_s = 1.5;
c.outlier_spectra_n_cal = 3;
c.outlier_spectra_n_null = 4;
c.outlier_spectra_n_shuffle = 3;
c.outlier_spectra_chunk_s = 1;
run_outlier_spectra_stage(c, tmp, outlier_spectra_stage_plan(c));

ce = struct('run_mode', 'essential', 'stages', {{'outlier_spectra', 'eigenvalue_density'}}, 'verbose', 'near-none');
run_eigenvalue_density_stage(ce, tmp, eigenvalue_density_stage_plan(ce));

cf = struct('run_mode', 'essential', 'verbose', 'near-none');
cf.frequency_response = struct('f', logspace(-3, 2, 21), 'f_single', logspace(-3, 2, 41), ...
    'sigma', -12:0.5:8, 'w', 0:1:30, 'n_sample', 60);
run_frequency_response_stage(cf, tmp, frequency_response_stage_plan(cf));

cs = mu7chain_config('essential');
cs.verbose = 'near-none';
cs.outlier_structures_n_override = 80;
cs.outlier_structures_conditions = conditions;
cs.outlier_structures_T_settle = 2;
cs.outlier_structures_window_s = 2;
cs.outlier_structures_lya_window_s = 1.5;
cs.outlier_structures_n_cal = 3;
cs.outlier_structures_n_null = 4;
cs.outlier_structures_n_shuffle = 3;
cs.outlier_structures_chunk_s = 0.7;
cs.outlier_structures_networks = {'lrA1p5W0p0', [], []; 'lrA1p4W0p6', [], []; 'EEx2p0', [], []; 'IIx1p25', [], []};
evalc('run_outlier_structures_stage(cs, tmp, outlier_structures_stage_plan(cs));');

made_up_transient(tmp, [conditions, {'sfa3_std2'}]);     % a fourth condition the selection must leave out
t_stage = toc(t0);

common = {'run_dir', tmp, 'verbose', 'near-none', 'conditions', conditions};
open_before = findall(groot, 'Type', 'figure');

%% 3, 4. The composites
% Figure 5
o = draw(@fig_main5_mechanism, common, fullfile(fig_dir, 'main5'), 'Fig_Main5_Grouped', 'ABCD', ...
    {'Panel A (fig_outlier_plane)', 'Panels B to D (fig_outlier_dose)'});
insets = findall(o.figs, 'Tag', 'fig_outlier_plane_inset');
marks = findall(o.figs, 'Tag', 'fig_outlier_plane_inset_mark');
assert(numel(insets) == 6 && numel(marks) == 12, 'six insets, each with an E and an I');
assert(isequal(sort(get(marks, 'String')), sort(repmat({'E'; 'I'}, 6, 1))));
ei = [excitatory_colormap(1); inhibitory_colormap(1)];
for k = 1:numel(marks)
    assert(isequal(marks(k).Color, ei(1 + strcmp(marks(k).String, 'I'), :)), 'the marks are in the cell-type colors');
end
runs = ospec_read_index(tmp);
first_seed = [runs.seed1(1), runs.seed2(1)];
txt = fileread(fullfile(fig_dir, 'main5', 'Fig_Main5_Grouped_table.md'));
assert(isequal(o.seed, first_seed) && contains(txt, sprintf('network seed [%d %d], the first of the stage', first_seed)), ...
    'the seed drawn in A is the stage''s first and the sidecar says so');
titles = arrayfun(@(a) string(a.Title.String), findall(o.figs, 'Type', 'axes'));
assert(all(ismember(["No adaptation", "Single-timescale adaptation", "Multiple-timescale adaptation"], titles)), 'sentence-case titles');
close(o.figs);
o = fig_main5_mechanism(common{:}, 'seed', [2 3], 'networks', {'EEx1p5'}, 'save', false);
assert(isequal(o.seed, [2 3]) && isempty(o.files) && numel(findall(o.figs, 'Tag', 'fig_outlier_plane_inset')) == 3, 'seed and networks are passed on');
close(o.figs);

% Figure 6
o = draw(@fig_main6_amplification, common, fullfile(fig_dir, 'main6'), 'Fig_Main6_Grouped', 'ABCDE', ...
    {'Panels A and B', 'Panels C and D', 'Panel E'});
net_ax = findall(o.figs, '-regexp', 'Tag', '^freqresp_network_');
assert(isequal(sort(get(net_ax, 'Tag')), {'freqresp_network_reference_c1'; 'freqresp_network_reference_c3'}), ...
    'C and D: the reference network, bulk radius and loop gain; no middle panel');
for k = 1:2
    assert(strcmp(net_ax(k).XLabel.String, 'Frequency (Hz)') && ~isempty(net_ax(k).XTickLabel), 'C and D carry their own frequency axis');
end
plane_ax = findall(o.figs, '-regexp', 'Tag', '^freqresp_plane_');
assert(numel(plane_ax) == 3 && all(startsWith(get(plane_ax, 'Tag'), 'freqresp_plane_reference_')), 'E: the reference network, one panel per condition');
gain_ax = findall(o.figs, '-regexp', 'Tag', '^main6_gain_');
assert(numel(gain_ax) == 3 && isscalar(findall(o.figs, 'Tag', 'main6_margin')) && isscalar(findall(o.figs, 'Type', 'legend')), ...
    'A and B: one margin panel, a gain panel per selected condition, one legend');
txt = fileread(fullfile(fig_dir, 'main6', 'Fig_Main6_Grouped_table.md'));
assert(contains(txt, '| no_adaptation | Worst case |') && ~contains(txt, 'sfa3_std2'), 'the numbers of A and B, the selected conditions only');
close(o.figs);
assert(throws(@() fig_main6_amplification(common{:}, 'network', 'EEx9p9', 'save', false), 'fresp_pick_networks:NoSuchNetwork'));

% S6
o = draw(@fig_supp6_grouped, common, fullfile(fig_dir, 'supp6'), 'Fig_Supp6_Grouped', 'ABC', ...
    {'Panel A (fig_eigdens_gains, reference)', 'Panel B, panel 1', 'Panel B, panel 2', 'Panel C (fig_eigdens_measures, reference)'});
titles = arrayfun(@(a) string(strjoin(cellstr(a.Title.String), ' ')), findall(o.figs, 'Type', 'axes'));
assert(nnz(titles == "Reference") == 1 && nnz(contains(titles, "E-to-E")) == 1, 'B: one cumulative panel per network, titled by it');
assert(numel(findall(o.figs, 'Type', 'axes')) == 2 + 2 + 4, 'A two panels, B two, C four');
close(o.figs);

% S7
o = draw(@fig_supp7_grouped, common, fullfile(fig_dir, 'supp7'), 'Fig_Supp7_Grouped', 'AB', ...
    {'Panel A (fig_outlier_slow_variables)', 'Panel B (fig_ostruct_summary)'});
close(o.figs);

% S8
o = draw(@fig_supp8_grouped, common, fullfile(fig_dir, 'supp8'), 'Fig_Supp8_Grouped', 'AB', ...
    {'Panel A (fig_freqresp_single_neuron)', 'Panel B (fig_freqresp_phase)'});
assert(numel(findall(o.figs, 'Type', 'axes')) == 7, 'A six panels, B one');
close(o.figs);

%% 5. The source figures' new options
o = fig_freqresp_network('run_dir', tmp, 'visible', false, 'out_dir', fullfile(fig_dir, 'net'), 'networks', {'EEx1p5'}, 'panels', 3);
assert(isscalar(findall(o.figs, 'Type', 'axes')) && strcmp(o.figs.Children.Children.Tag, 'freqresp_network_EEx1p5_c3'));
txt = fileread(fullfile(fig_dir, 'net', 'Fig_Freqresp_Network_table.md'));
assert(contains(txt, '| EEx1p5 | sfa3_stdc2 |') && ~contains(txt, '| reference | sfa3_stdc2 |'), 'the table has the networks drawn');
close(o.figs);
o = fig_freqresp_network('run_dir', tmp, 'visible', false, 'save', false);
assert(numel(findall(o.figs, 'Type', 'axes')) == 6, 'the default is unchanged: every network, three panels');
close(o.figs);
assert(throws(@() fig_freqresp_network('run_dir', tmp, 'save', false, 'visible', false, 'networks', {'nope'}), 'fresp_pick_networks:NoSuchNetwork'));
o = fig_freqresp_plane('run_dir', tmp, 'visible', false, 'save', false, 'networks', {'reference'});
ax = findall(o.figs, 'Type', 'axes');
assert(numel(ax) == 3 && any(arrayfun(@(a) strcmp(a.Title.String, 'Single-timescale adaptation'), ax)), 'one row, sentence-case titles');
close(o.figs);
o = fig_eigdens_radial('run_dir', tmp, 'visible', false, 'save', false, 'panels', {'density'});
ax = findall(o.figs, 'Type', 'axes');
assert(isscalar(ax) && strcmp(ax.YScale, 'log'), 'the density panel alone');
close(o.figs);
assert(throws(@() fig_eigdens_radial('run_dir', tmp, 'visible', false, 'save', false, 'panels', {'nope'}), 'fig_eigdens_radial:Panels'));
o = fig_outlier_dose('run_dir', tmp, 'visible', false, 'images', false, 'out_dir', fullfile(fig_dir, 'table_only'));
assert(isequal(o.files, {'Fig_Outlier_Dose_table.md'}) && isempty(dir(fullfile(fig_dir, 'table_only', '*.png'))), ...
    'images false: the sidecar alone');
close(o.figs);
% The summary's detachment panel holds every value, negative ones included.
o = fig_ostruct_summary('run_dir', tmp, 'visible', false, 'save', false);
F = ostruct_figure_data(struct('data_file', '', 'run_dir', tmp, 'conditions', {{}}), 'test');
ax = findall(o.figs, 'Type', 'axes');
top = ax(arrayfun(@(a) contains(a.YLabel.String, 'bulk radius'), ax));
lo = min([F.lo.mode_rho(:); -1], [], 'omitnan'); hi = max([F.hi.mode_rho_uniform(:); F.hi.mode_rho(:); 1], [], 'omitnan');
assert(isscalar(top) && strcmp(top.YScale, 'linear') && top.YLim(1) < lo && top.YLim(2) > hi, ...
    'every detachment is inside the axis (%.3g to %.3g in [%.3g, %.3g])', lo, hi, top.YLim);
close(o.figs);

%% 6. Errors, never a fallback; nothing left open
fns = {@fig_main5_mechanism, @fig_main6_amplification, @fig_supp6_grouped, @fig_supp7_grouped, @fig_supp8_grouped};
for k = 1:numel(fns)
    assert(throws(@() fns{k}('run_dir', empty_run, 'save', false), 'resolve_data_file:NotFound'), 'composite %d on an empty run', k);
    assert(throws(@() fns{k}('save', false), 'resolve_data_file:NoRunDir'), 'composite %d without a run', k);
end
left = setdiff(findall(groot, 'Type', 'figure'), open_before);
assert(isempty(left), '%d figure(s) left open', numel(left));

fprintf('PASS test_composite_spectral_figures (stages %.0f s)\n', t_stage);

%% ------------------------------------------------------------------------
function o = draw(fn, args, out_dir, tag, letters, sections)
% Draw and save one composite; assert its files, its letters and its sidecar's sections.
before = findall(groot, 'Type', 'figure');
o = fn(args{:}, 'out_dir', out_dir);
assert(isscalar(o.figs) && isgraphics(o.figs, 'figure') && strcmp(o.figs.Visible, 'off'), '%s: the figure', tag);
assert(isequal(setdiff(findall(groot, 'Type', 'figure'), before), o.figs), '%s: the source figures are closed', tag);
assert(any(strcmp(o.files, [tag '.png'])) && any(strcmp(o.files, [tag '_table.md'])), '%s: the files', tag);
info = imfinfo(fullfile(out_dir, [tag '.png']));
pos = o.figs.Position;                                   % inches
assert(abs(info.Height / info.Width - pos(4) / pos(3)) < 0.25 * pos(4) / pos(3), '%s: the PNG is about the figure''s shape', tag);
for k = 1:numel(letters)
    h = findall(o.figs, 'Tag', ['panel_letter_' letters(k)]);
    assert(isscalar(h) && strcmp(h.String, sprintf('(%s)', letters(k))) && strcmp(h.FontWeight, 'normal'), '%s: letter %s', tag, letters(k));
end
assert(isempty(findall(o.figs, 'Tag', ['panel_letter_' char(letters(end) + 1)])), '%s: no further letter', tag);
txt = fileread(fullfile(out_dir, [tag '_table.md']));
assert(startsWith(txt, '# ') && contains(txt, 'Data read:') && contains(txt, '|---|'), '%s: the sidecar', tag);
for k = 1:numel(sections)
    assert(contains(txt, ['## ' sections{k}]), '%s: the sidecar''s section "%s"', tag, sections{k});
end
end

function made_up_transient(run_dir, names)
% Files with the fields panels A and B of Figure 6 read, filled with
% made-up numbers: the eig_heatmap stage's abscissas and the transient_gain
% stage's gains. Not a stand-in for those stages' tests.
rs = RandStream('mt19937ar', 'Seed', 3);
n = numel(names);
cond_names = names;
num_abscissa_by_cond = cell(1, n); spec_abscissa_by_cond = cell(1, n);
for i = 1:n
    spec_abscissa_by_cond{i} = -1 + rs.rand(1, 20);
    num_abscissa_by_cond{i} = spec_abscissa_by_cond{i} + 40 + 10 * i + 5 * rs.rand(1, 20);
end
mkdir(fullfile(run_dir, 'eig_heatmap'));
save(fullfile(run_dir, 'eig_heatmap', 'eig_heatmap_data.mat'), 'cond_names', 'num_abscissa_by_cond', 'spec_abscissa_by_cond');
t = linspace(0, 3, 31);
settings = struct('variants', {{'frozen_x', 'frozen_full', 'active'}}, 'horizon_s', 3);
keys = {'G_worst', 'G_noise', 'G_ei_diff', 'G_ei_sum', 'G_lyap'};
results = struct('name', names, 'samples', []);
for i = 1:n
    smp = struct('kind', repmat({'regular'}, 1, 4), 't', t);
    for s = 1:4
        for k = 1:5
            smp(s).(keys{k}) = exp((6 - k) * 0.3 * t / i) .* (1 + 0.1 * rs.rand(3, numel(t)));
        end
    end
    smp(end + 1) = smp(1); smp(end).kind = 'other'; %#ok<AGROW>
    results(i).samples = smp;
end
mkdir(fullfile(run_dir, 'transient_gain'));
save(fullfile(run_dir, 'transient_gain', 'transient_gain_data.mat'), 'results', 'settings');
end

function tf = throws(fn, id)
try
    fn();
    tf = false;
catch err
    tf = strcmp(err.identifier, id);
    if ~tf; fprintf(2, 'expected %s, got %s: %s\n', id, err.identifier, err.message); end
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end
