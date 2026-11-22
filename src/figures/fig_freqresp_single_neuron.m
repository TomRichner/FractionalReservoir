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
%   line    solid: adaptation and depression free, with the quartiles over
%           seeds as a band; dashed: adaptation alone; dotted: depression
%           alone
%
% ON THE PANELS: the drawing, a title over each column saying which neuron
% it is, and one key at the right (the conditions' colors and the three
% line styles). Until 2026-11-22 there was neither, and the columns and
% the solid, dashed and dotted lines could only be told from the caption.
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
    cfg.key         (1,:) char {mustBeMember(cfg.key, {'east', 'south'})} = 'east'   % where the key goes: right of the panels, or under them in three columns
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
mode_names = {'adaptation and depression', 'adaptation alone', 'depression alone'};
styles = {'-', '--', ':'};
f = D.settings.f_single;

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Position', [80 80 1380 620], 'Visible', vis);
if all(isfinite(cfg.fig_size)); fig.Units = 'inches'; fig.Position(3:4) = cfg.fig_size; end
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
ax1 = gobjects(1, 3); ax2 = gobjects(1, 3);
y_hi = [1.05 80];
for p = 1:3
    ax1(p) = nexttile(tl, p); hold(ax1(p), 'on');
    ax2(p) = nexttile(tl, 3 + p); hold(ax2(p), 'on');
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
            if m == 1
                band(ax1(p), f, mag, col); band(ax2(p), f, ph, col);
            end
            plot(ax1(p), f, median(mag, 1), styles{m}, 'Color', col, 'LineWidth', st.line_lw);
            plot(ax2(p), f, median(ph, 1), styles{m}, 'Color', col, 'LineWidth', st.line_lw);
            y_hi = max(y_hi, [1.05 * max(median(mag, 1)), 10 * ceil(max(median(ph, 1)) / 10)]);
            s = arrayfun(@(o) o.summary.(modes{m}), own);
            rows(end + 1, :) = {columns{p, 1}, columns{p, 2}, conditions{c}, mode_names{m}, numel(own), ...
                median([own.slope]), median([own.rate]), median([own.F]), ...
                median([s.G0]), median([s.dB0]), median([s.f_half], 'omitnan'), median([s.f_3dB], 'omitnan'), ...
                median([s.f_90], 'omitnan'), median([s.phase_max]), median([s.f_phase]), median([s.peak]), ...
                max([own.max_diff])}; %#ok<AGROW>
        end
    end
end
% The columns' titles and the key.
column_titles = {'E neuron, median rate', 'I neuron, median rate', 'E neuron, driven hard'};
for p = 1:3
    title(ax1(p), column_titles{p}, 'FontSize', st.title_fs, 'FontWeight', 'normal', 'Color', 'k');
end
key = gobjects(1, 0); key_names = {};
cond_titles = grouped_condition_label(conditions, 'sentence');
for c = 1:numel(conditions)
    key(end + 1) = plot(ax1(3), NaN, NaN, '-', 'Color', st.condition_color(conditions{c}), 'LineWidth', st.line_lw); %#ok<AGROW>
    key_names{end + 1} = cond_titles{c}; %#ok<AGROW>
end
style_names = {'Adaptation and depression', 'Adaptation alone', 'Depression alone'};
for m = 1:3
    key(end + 1) = plot(ax1(3), NaN, NaN, styles{m}, 'Color', [0.35 0.35 0.35], 'LineWidth', st.line_lw); %#ok<AGROW>
    key_names{end + 1} = style_names{m}; %#ok<AGROW>
end
lg = legend(ax1(3), key, key_names, 'Box', 'off', 'FontSize', st.tick_fs, 'AutoUpdate', 'off');
lg.Layout.Tile = cfg.key;
if strcmp(cfg.key, 'south'); lg.NumColumns = 3; end
ticks = 10 .^ (ceil(log10(f(1)) - 1e-9):floor(log10(f(end)) + 1e-9));
for p = 1:3
    set([ax1(p), ax2(p)], 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', ...
        'Layer', 'top', 'XScale', 'log', 'XLim', f([1 end]), 'XTick', ticks, 'XMinorTick', 'off');
    set(ax1(p), 'YLim', [0 y_hi(1)], 'XTickLabel', []);
    set(ax2(p), 'YLim', [0 y_hi(2)], 'YTick', 0:20:y_hi(2));
    xlabel(ax2(p), 'Frequency (Hz)', 'FontSize', st.label_fs);
    if p == 1
        ylabel(ax1(p), 'Gain / frozen gain', 'FontSize', st.label_fs);
        ylabel(ax2(p), 'Phase lead (deg)', 'FontSize', st.label_fs);
    else
        set([ax1(p), ax2(p)], 'YTickLabel', []);
    end
end

r1 = slow_runs(1);
head = {'network', 'neuron picked', 'condition', 'free', 'seeds', 'slope', 'rate', 'output factor F', ...
    'G(0)', 'G(0) dB', 'f_half Hz', 'f_3dB Hz', 'f_90 Hz', 'max phase lead deg', 'at Hz', 'peak', ...
    'closed vs numeric (largest)'};
caption = [ ...
    'One neuron''s frequency response: the gain from its dendritic input to its synaptic output (route to the excitatory type) with its slow variables free, ' ...
    'divided by the frozen gain (slope times available vesicle fraction, the limit of high frequency). TOP: magnitude; BOTTOM: phase lead. ' ...
    sprintf('COLUMNS, left to right: the excitatory neuron of median rate in the network `%s`; the inhibitory neuron of median rate there; ', networks{1}) ...
    sprintf('a hard-driven excitatory neuron (90th percentile of rate) in the network `%s`; each among the neurons with nonzero slope, ', networks{end}) ...
    sprintf('at calibration state %d (t = %.4g s on the first seed) of every seed, so each seed and condition has its own neuron and operating point. ', r1.shown_state, r1.t_cal(r1.shown_state)) ...
    sprintf('The curves are the MEDIAN over the %d seed(s) of the magnitude and of the phase. ', size(D.settings.seeds, 1)) ...
    sprintf('COLOR: the condition (%s), in the manuscript''s condition colors. ', strjoin(conditions, ', ')) ...
    'LINE: solid adaptation and depression both free, with the quartiles over seeds as a band; ' ...
    'dashed adaptation alone (depression held); dotted depression alone (adaptation held). A value of 1 means the slow variables do nothing at that frequency. ' ...
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
