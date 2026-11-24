function out = fig_freqresp_single_neuron(cfg)
% FIG_FREQRESP_SINGLE_NEURON One neuron's frequency response, magnitude and phase, at three operating points taken from the runs.
%
%   out = FIG_FREQRESP_SINGLE_NEURON('run_dir', d)
%   out = FIG_FREQRESP_SINGLE_NEURON('data_file', f, 'out_dir', o)
%
% Draws what the opt-in stage frequency_response saved (FD-041;
% run_frequency_response_stage). It reads one run directory,
% <run_dir>/frequency_response/frequency_response_data.mat, or the file
% given, errors when neither holds it, and computes nothing beyond medians.
% Supplementary.
%
% THE RESPONSE. The gain from a neuron's dendritic input to its synaptic
% output on the route to the excitatory type, with its slow variables free,
% divided by the frozen gain (slope times available vesicle fraction, the
% limit of high frequency). Magnitude below 1 is attenuation by adaptation
% and depression; 1 where they cannot follow.
%
% TWO ROWS (magnitude; phase lead), THREE COLUMNS, each a neuron chosen
% among those with nonzero slope at the stage's shown calibration state:
%   left    the excitatory neuron of median rate, first network
%   middle  the inhibitory neuron of median rate, first network
%   right   the excitatory neuron at the 90th percentile of rate, LAST
%           network (E-to-E x 1.5 in the bundle: a hard-driven neuron)
% Each seed and condition has its own neuron and operating point; drawn is
% the MEDIAN over seeds of the magnitude and of the phase.
%   color   the condition (the conditions with slow variables)
%   line    solid: spike-frequency adaptation and depression free, with the
%           quartiles over seeds as a band; dashed: spike-frequency
%           adaptation alone; dotted: depression alone (the key says
%           "spike-frequency adaptation", not "adaptation": Tom's notes of
%           2026-11-22, S8 and rule 7)
%
% ON THE PANELS: the drawing, a title over each column saying which neuron
% it is, and one key at the right (the line styles and the conditions'
% colors). Until 2026-11-22 there was neither, and the columns and
% the solid, dashed and dotted lines could only be told from the caption.
%
% ONE NEURON (cfg.columns, a single index): its gain and its phase side by
% side, each with its own frequency axis, no title (the caption names the
% neuron). S8 since 2026-11-22 draws the E neuron of median rate this way
% (Tom: one gain panel and one phase panel); the table keeps all three
% neurons whatever is drawn.
% Fig_Freqresp_Single_Neuron_table.md
% holds the caption and, per column, condition and mechanism, the medians
% over seeds of the operating point, the zero-frequency gain, the corner
% frequencies, the largest phase lead, and the agreement of the closed form
% with the Jacobian's blocks.
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions, in the order given; empty
% draws every condition with slow variables the stage saved.
%
% See also: run_frequency_response_stage, fresp_run_one, fresp_summary,
%           fig_freqresp_phase, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the networks are recorded in the .mat
    cfg.conditions                = {}    % names to draw, in order; empty = all with slow variables (select_conditions)
    cfg.fig_size    (1,2) double  = [NaN NaN]   % inches [w h]; NaN = 1380 x 620 px (fig_supp8_grouped lays it out on the common canvas)
    cfg.key         (1,:) char {mustBeMember(cfg.key, {'east', 'south'})} = 'east'   % where the key goes: right of the panels, or under them in two columns
    cfg.columns     (1,:) double {mustBeMember(cfg.columns, 1:3)} = 1:3   % the neurons drawn (1 E median, 2 I median, 3 E driven hard); the table has all three
end
if numel(unique(cfg.columns)) ~= numel(cfg.columns)
    error('fig_freqresp_single_neuron:Columns', 'Each neuron at most once.');
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
[D, data_file] = fresp_fig_load(cfg.data_file, cfg.run_dir);
slow_runs = D.runs([D.runs.has_slow]);
if isempty(slow_runs)
    error('fig_freqresp_single_neuron:NoSlow', ...
        'None of the stage''s conditions has slow variables: there is no single-neuron response to draw.');
end
available = unique({slow_runs.condition}, 'stable');
if isempty(cfg.conditions)
    conditions = available;
else
    select_conditions(D.settings.conditions, cfg.conditions);
    conditions = select_conditions(available, intersect(cellstr(cfg.conditions), available, 'stable'));
end
networks = D.settings.networks;
columns = {networks{1}, 'E median'; networks{1}, 'I median'; networks{end}, 'E driven'};
modes = {'both', 'sfa', 'std'};
% Rule 7 of Tom's notes (2026-11-22): where a name means spike-frequency
% adaptation as against depression, it says so ("adaptation alone" read as
% any adaptation). The three read alike.
mode_names = {'spike-frequency adaptation and depression', 'spike-frequency adaptation alone', 'depression alone'};
styles = {'-', '--', ':'};
f = D.settings.f_single;
drawn = cfg.columns;
n_d = numel(drawn);

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Position', [80 80 1380 620], 'Visible', vis);
if n_d == 1; fig.Position(3:4) = [1000 480]; end
if all(isfinite(cfg.fig_size)); fig.Units = 'inches'; fig.Position(3:4) = cfg.fig_size; end
% One neuron (S8 since 2026-11-22): its gain and its phase side by side.
% Several: gain on top, phase below, a column per neuron.
if n_d == 1
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
else
    tl = tiledlayout(fig, 2, n_d, 'TileSpacing', 'compact', 'Padding', 'compact');
end
rows = {};
ax1 = gobjects(1, n_d); ax2 = gobjects(1, n_d);
for j = 1:n_d
    if n_d == 1
        ax1(j) = nexttile(tl, 1); ax2(j) = nexttile(tl, 2);
    else
        ax1(j) = nexttile(tl, j); ax2(j) = nexttile(tl, n_d + j);
    end
    hold(ax1(j), 'on'); hold(ax2(j), 'on');
end
y_hi = [1.05 80];
% Every column's numbers go to the table, drawn or not (the neurons S8 no
% longer draws are kept there).
for p = 1:3
    j = find(drawn == p, 1);
    for c = 1:numel(conditions)
        runs = slow_runs(strcmp({slow_runs.network}, columns{p, 1}) & strcmp({slow_runs.condition}, conditions{c}));
        own = arrayfun(@(r) r.single(strcmp({r.single.label}, columns{p, 2})), runs);
        own = own(isfinite([own.neuron]));              % a seed with no neuron of nonzero slope has none
        if isempty(own); continue; end
        col = st.condition_color(conditions{c});
        for m = 1:3
            G = arrayfun(@(o) o.G.(modes{m}), own, 'UniformOutput', false);
            G = cat(1, G{:});                           % seeds x nf
            mag = abs(G); ph = angle(G) * 180 / pi;
            if ~isempty(j)
                if m == 1
                    band(ax1(j), f, mag, col); band(ax2(j), f, ph, col);
                end
                plot(ax1(j), f, median(mag, 1), styles{m}, 'Color', col, 'LineWidth', st.line_lw);
                plot(ax2(j), f, median(ph, 1), styles{m}, 'Color', col, 'LineWidth', st.line_lw);
                y_hi = max(y_hi, [1.05 * max(median(mag, 1)), 10 * ceil(max(median(ph, 1)) / 10)]);
            end
            s = arrayfun(@(o) o.summary.(modes{m}), own);
            rows(end + 1, :) = {columns{p, 1}, columns{p, 2}, conditions{c}, mode_names{m}, numel(own), ...
                median([own.slope]), median([own.rate]), median([own.F]), ...
                median([s.G0]), median([s.dB0]), median([s.f_half], 'omitnan'), median([s.f_3dB], 'omitnan'), ...
                median([s.f_90], 'omitnan'), median([s.phase_max]), median([s.f_phase]), median([s.peak]), ...
                max([own.max_diff])}; %#ok<AGROW>
        end
    end
end
% The columns' titles (several neurons only: with one, the caption names it)
% and the key.
column_titles = {'E neuron, median rate', 'I neuron, median rate', 'E neuron, driven hard'};
if n_d > 1
    for j = 1:n_d
        title(ax1(j), column_titles{drawn(j)}, 'FontSize', st.title_fs, 'FontWeight', 'normal', 'Color', 'k');
    end
end
key = gobjects(1, 0); key_names = {};
style_names = {'Spike-frequency adaptation and depression', 'Spike-frequency adaptation alone', 'Depression alone'};
for m = 1:3
    key(end + 1) = plot(ax1(end), NaN, NaN, styles{m}, 'Color', [0.35 0.35 0.35], 'LineWidth', st.line_lw); %#ok<AGROW>
    key_names{end + 1} = style_names{m}; %#ok<AGROW>
end
cond_titles = grouped_condition_label(conditions, 'sentence');
for c = 1:numel(conditions)
    key(end + 1) = plot(ax1(end), NaN, NaN, '-', 'Color', st.condition_color(conditions{c}), 'LineWidth', st.line_lw); %#ok<AGROW>
    key_names{end + 1} = cond_titles{c}; %#ok<AGROW>
end
lg = legend(ax1(end), key, key_names, 'Box', 'off', 'FontSize', st.tick_fs, 'AutoUpdate', 'off');
lg.Layout.Tile = cfg.key;
% Under the panels: the line styles down the first column, the conditions'
% colors down the second (the entries fill the columns top to bottom).
if strcmp(cfg.key, 'south'); lg.NumColumns = 2; end
ticks = 10 .^ (ceil(log10(f(1)) - 1e-9):floor(log10(f(end)) + 1e-9));
for j = 1:n_d
    set([ax1(j), ax2(j)], 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', ...
        'Layer', 'top', 'XScale', 'log', 'XLim', f([1 end]), 'XTick', ticks, 'XMinorTick', 'off');
    set(ax1(j), 'YLim', [0 y_hi(1)]);
    set(ax2(j), 'YLim', [0 y_hi(2)], 'YTick', 0:40:y_hi(2));   % three ticks on 0 to 80 (Tom, 2026-11-23, G5)
    xlabel(ax2(j), 'Frequency (Hz)', 'FontSize', st.label_fs);
    if n_d == 1
        xlabel(ax1(j), 'Frequency (Hz)', 'FontSize', st.label_fs);   % side by side: each its own axis
    else
        set(ax1(j), 'XTickLabel', []);
    end
    if j == 1 || n_d == 1
        ylabel(ax1(j), 'Gain / frozen gain', 'FontSize', st.label_fs);
        ylabel(ax2(j), 'Phase lead (deg)', 'FontSize', st.label_fs);
    else
        set([ax1(j), ax2(j)], 'YTickLabel', []);
    end
end
for j = 1:n_d
    ax1(j).Tag = sprintf('freqresp_single_gain_%d', drawn(j));
    ax2(j).Tag = sprintf('freqresp_single_phase_%d', drawn(j));
end

r1 = slow_runs(1);
head = {'network', 'neuron picked', 'condition', 'free', 'seeds', 'slope', 'rate', 'output factor F', ...
    'G(0)', 'G(0) dB', 'f_half Hz', 'f_3dB Hz', 'f_90 Hz', 'max phase lead deg', 'at Hz', 'peak', ...
    'closed vs numeric (largest)'};
if n_d == 1
    layout_text = sprintf('LEFT: magnitude; RIGHT: phase lead, of neuron (%d) below. ', drawn);
else
    layout_text = sprintf('TOP: magnitude; BOTTOM: phase lead; COLUMNS, left to right, neurons (%s) below. ', strjoin(string(drawn), ', '));
end
caption = [ ...
    'One neuron''s frequency response: the gain from its dendritic input to its synaptic output (route to the excitatory type) with its slow variables free, ' ...
    'divided by the frozen gain (slope times available vesicle fraction, the limit of high frequency). ' ...
    layout_text ...
    sprintf('The neurons (the table has all three): (1) the excitatory neuron of median rate in the network `%s`; (2) the inhibitory neuron of median rate there; ', networks{1}) ...
    sprintf('(3) a hard-driven excitatory neuron (90th percentile of rate) in the network `%s`; each among the neurons with nonzero slope, ', networks{end}) ...
    sprintf('at calibration state %d (t = %.4g s on the first seed) of every seed, so each seed and condition has its own neuron and operating point. ', r1.shown_state, r1.t_cal(r1.shown_state)) ...
    sprintf('The curves are the MEDIAN over the %d seed(s) of the magnitude and of the phase. ', size(D.settings.seeds, 1)) ...
    sprintf('COLOR: the condition (%s), in the manuscript''s condition colors. ', strjoin(conditions, ', ')) ...
    'LINE: solid spike-frequency adaptation and depression both free, with the quartiles over seeds as a band; ' ...
    'dashed spike-frequency adaptation alone (depression held); dotted depression alone (spike-frequency adaptation held). A value of 1 means the slow variables do nothing at that frequency. ' ...
    'The curves are the closed form (`fresp_gain`); the same response computed from the model''s own Jacobian blocks for that neuron (`fresp_neuron_schur`) agrees to the last column of the table (relative, largest over seeds). ' ...
    'Table: medians over seeds. G(0): gain at zero frequency over the frozen gain. f_half: where half the attenuation in dB is gone. f_3dB, f_90: above these the gain stays within 3 dB, within 10%, of the frozen gain. ' ...
    'peak: the largest magnitude (1 = no resonance of the neuron itself).'];
out = fresp_fig_save(fig, cfg, out_dir, 'Fig_Freqresp_Single_Neuron', data_file, ...
    'Single-neuron frequency response', caption, head, rows);
end

%% ------------------------------------------------------------------------
function band(ax, x, M, col)
% Quartiles down the rows (the seeds); nothing for a single seed.
if size(M, 1) < 2; return; end
lo = prctile(M, 25, 1); hi = prctile(M, 75, 1);
patch(ax, [x, fliplr(x)], [lo, fliplr(hi)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
