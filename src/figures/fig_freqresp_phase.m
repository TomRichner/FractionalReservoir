function out = fig_freqresp_phase(cfg)
% FIG_FREQRESP_PHASE The phase lead of a neuron's gain against frequency, one against multiple timescales, at a matched operating point.
%
%   out = FIG_FREQRESP_PHASE('run_dir', d)
%   out = FIG_FREQRESP_PHASE('data_file', f, 'out_dir', o)
%
% Draws what the opt-in stage frequency_response saved (FD-041;
% run_frequency_response_stage). It reads one run directory,
% <run_dir>/frequency_response/frequency_response_data.mat, or the file
% given, errors when neither holds it, and computes nothing beyond medians.
% Supplementary. (The phase-lead panel of the exploration's
% freqresp_one_vs_multiple, as a figure of its own; its two ratio panels
% are not drawn: Tom, 2026-11-20.)
%
% ONE PANEL. Every excitatory neuron of the stage's first network is put at
% one rate (settings.r0) with its slow variables at their steady state for
% that rate (fresp_matched), so the conditions are compared at one
% operating point and what differs is the timescale structure. Drawn: the
% phase of the gain on the E-to-E route over the frozen gain, in degrees,
% the MEDIAN over every excitatory neuron of every seed.
%   color   the condition (the conditions with slow variables)
%   line    solid: adaptation and depression free; dashed: adaptation alone;
%           dotted: depression alone
%   gray    the phase lag of the dendrite, atan(omega tau_d), which a
%           mode's gain must cancel
%
% NOTHING BUT THE DRAWING ON THE PANEL. Fig_Freqresp_Phase_table.md holds
% the caption and, per condition and mechanism, the largest phase lead,
% where it is, its quartiles over neurons, and the lead at chosen
% frequencies.
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions, in the order given; empty
% draws every condition with slow variables the stage saved.
%
% See also: run_frequency_response_stage, fresp_matched, fig_freqresp_network,
%           fig_freqresp_single_neuron, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the networks are recorded in the .mat
    cfg.conditions                = {}    % names to draw, in order; empty = all with slow variables (select_conditions)
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
[D, data_file] = fresp_fig_load(cfg.data_file, cfg.run_dir);
if isempty(D.matched)
    error('fig_freqresp_phase:NoMatched', ...
        'The stage saved no matched operating point: none of its conditions has slow variables.');
end
available = unique({D.matched.condition}, 'stable');
if isempty(cfg.conditions)
    conditions = available;
else
    % A condition without slow variables has no response to draw: asking
    % for the run's whole set is not an error, a name the run lacks is.
    select_conditions(D.settings.conditions, cfg.conditions);
    conditions = select_conditions(available, intersect(cellstr(cfg.conditions), available, 'stable'));
end
modes = {'both', 'sfa', 'std'};
mode_names = {'adaptation and depression', 'adaptation alone', 'depression alone'};
styles = {'-', '--', ':'};
f = D.matched(1).f;
tau_d = D.runs(1).tau_d;
f_tab = [0.01 0.1 0.3 1 3];

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Position', [80 80 560 440], 'Visible', vis);
ax = axes(fig); hold(ax, 'on');
plot(ax, f, atan(2 * pi * f * tau_d) * 180 / pi, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
rows = {};
y_hi = 80;
for c = 1:numel(conditions)
    M = D.matched(strcmp({D.matched.condition}, conditions{c}));
    col = st.condition_color(conditions{c});
    for m = 1:3
        G = arrayfun(@(mt) double(mt.G.(modes{m})), M, 'UniformOutput', false);
        ph = angle(cat(1, G{:})) * 180 / pi;              % every neuron of every seed x nf
        med = median(ph, 1);
        plot(ax, f, med, styles{m}, 'Color', col, 'LineWidth', st.line_lw);
        y_hi = max(y_hi, 10 * ceil(max(med) / 10));
        [pk, ipk] = max(med);
        q = prctile(ph(:, ipk), [25 75]);
        rows(end + 1, :) = [{conditions{c}, mode_names{m}, numel(M), size(ph, 1), pk, f(ipk), ...
            sprintf('%.3g to %.3g', q(1), q(2))}, num2cell(interp1(log(f), med, log(f_tab)))]; %#ok<AGROW>
    end
end
set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top', ...
    'XScale', 'log', 'XLim', f([1 end]), 'XTick', 10 .^ (ceil(log10(f(1)) - 1e-9):floor(log10(f(end)) + 1e-9)), ...
    'XMinorTick', 'off', 'YLim', [0 y_hi], 'YTick', 0:20:y_hi);
xlabel(ax, 'Frequency (Hz)', 'FontSize', st.label_fs);
ylabel(ax, 'Phase lead (deg)', 'FontSize', st.label_fs);

head = [{'condition', 'free', 'seeds', 'neurons pooled', 'largest median phase lead (deg)', 'at Hz', ...
    'quartiles over neurons there (deg)'}, arrayfun(@(v) sprintf('lead at %g Hz (deg)', v), f_tab, 'UniformOutput', false)];
caption = [ ...
    'The phase lead of one neuron''s gain against frequency, at a matched operating point. ' ...
    sprintf('Every excitatory neuron of the network `%s` is put at rate %.3g (the linear part of the nonlinearity, slope 1) with its slow variables at their steady state for that rate (`fresp_matched_state`), ', ...
        D.matched(1).network, D.matched(1).r0) ...
    'so the conditions are compared at one operating point; on one seed a neuron has the same fast adaptation constant and releasable-pool times in every condition. ' ...
    'Drawn: the phase of the gain from dendritic input to synaptic output on the E-to-E route, relative to the frozen gain (the limit of high frequency), ' ...
    sprintf('the MEDIAN over every excitatory neuron of every seed (%d seed(s)). ', size(D.settings.seeds, 1)) ...
    sprintf('COLOR: the condition (%s), in the manuscript''s condition colors. ', strjoin(conditions, ', ')) ...
    'LINE: solid adaptation and depression both free; dashed adaptation alone (depression held); dotted depression alone (adaptation held). ' ...
    sprintf('GRAY: the phase lag of the dendrite, atan(omega tau_d) with tau_d = %.3g s, which a mode''s gain must cancel for the mode to be marginal at that frequency. ', tau_d) ...
    'A positive value is a lead: the slow variables subtract a delayed copy of the rate, which advances the output. ' ...
    'It is not a state the network visits (there the rates differ by neuron and the slow variables lag): it isolates what the timescale structure does at a given slope and rate.'];
out = fresp_fig_save(fig, cfg, out_dir, 'Fig_Freqresp_Phase', data_file, ...
    'Phase lead at a matched operating point', caption, head, rows);
end
