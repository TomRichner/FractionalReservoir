% test_single_neuron_adaptation.m - Single-neuron demo of SFA vs STD
%
% Quick first pass at a methods figure showing what spike-frequency
% adaptation (SFA) and short-term synaptic depression (STD) each do to a
% single neuron's output. A lone excitatory neuron (no recurrent
% connectivity, W = 0) is driven with a step in external input u under 4
% conditions: no adaptation, SFA only, STD only, and SFA+STD. model.plot()
% shows the stimulus, dendritic state x, firing rate r, synaptic output
% b.*r, and the adaptation/depression variable itself for each condition;
% the 4 conditions are then assembled into a single 4-column figure -- based
% on scripts/tests/test_SRNN2_defaults.m, but n=1 and a clean deterministic
% step replaces the default random sparse-step stimulus.

close all; clear; clc;

setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% Step stimulus (shared by both conditions)
step_amp = 0.8;       % external drive amplitude during the step (mid-range, avoids saturating r)
t_on     = 5;         % step onset (s)
t_off    = 15;        % step offset (s); long enough for SFA to fully develop (tau_a_E=3s below)
T_range  = [-10, 20];  % simulate a bit before/after the step

input_config = struct();
input_config.intrinsic_drive = [];  % filled in by generate_stimulus()
input_config.generator = @(params, T, fs, rng_seed, input_config) ...
    single_neuron_step_generator(params, T, fs, t_on, t_off, step_amp);

%% Common single-neuron model settings
base_args = { ...
    'n', 1, ...
    'f', 1, ...              % the lone neuron is excitatory
    'indegree', 0, ...       % no recurrent connectivity (W = 0)
    'mu_E_tilde_relative', 0, 'mu_I_tilde_relative', 0, ...
    'sigma_E_tilde_relative', 0, 'sigma_I_tilde_relative', 0, ...
    'T_range', T_range, ...
    'input_config', input_config, ...
    'lya_method', 'none', ...
    'x0_std', 0, ...         % deterministic x(0) = 0
    'activation', 'piecewise', ...
    'S_c', 0.35, ...         % nonlinearity center
    'S_a', 1.0, ...          % hard sigmoid (piecewise linear)
    'tau_b_E_rec', 2};

%% Conditions: no adaptation, SFA only, STD only, SFA+STD
% Single adaptation timescale with a larger-than-production c_E so the
% firing-rate decay is clearly visible in this illustrative demo (the
% production default splits a smaller budget over 3 timescales).
condition_titles = {'No adaptation', 'SFA only', 'STD only', 'SFA + STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...
    {'n_a_E', 1, 'tau_a_E', 3, 'c_E', 1.0, 'n_b_E', 0}, ...
    {'n_a_E', 0, 'n_b_E', 1}, ...
    {'n_a_E', 1, 'tau_a_E', 3, 'c_E', 1.0, 'n_b_E', 1}};
n_cond = numel(condition_titles);

fig_list = gobjects(n_cond, 1);
ax_cols = cell(n_cond, 1);
for i = 1:n_cond
    model = SRNNModel2(base_args{:}, condition_args{i}{:});
    model.build();
    model.run();
    [fig_list(i), ax_cols{i}] = model.plot('T_plot', [0, 20]);
    strip_scale_bar(ax_cols{i});
    set_traces_black(ax_cols{i});
end

%% Assemble all conditions into a single figure (one column per condition)
fig_combined = combine_columns(ax_cols, condition_titles);
close(fig_list);

%% Save the combined figure
save_some_figs_to_folder_2(this_dir, 'sfa_std_single_neuron_example', fig_combined.Number, {'fig', 'png', 'svg'});

%% Local functions
function strip_scale_bar(ax_handles)
% STRIP_SCALE_BAR Remove the per-panel scale bar drawn by plot_SRNN_tseries
% (a thick line + "X seconds" text). Must run BEFORE set_traces_black,
% since that recolors the scale-bar line and would make it unidentifiable
% by LineWidth.
for i = 1:numel(ax_handles)
    delete(findobj(ax_handles(i), 'Type', 'line', 'LineWidth', 4));
    txt = findobj(ax_handles(i), 'Type', 'text');
    for k = 1:numel(txt)
        if contains(txt(k).String, 'seconds')
            delete(txt(k));
        end
    end
end
end

function set_traces_black(ax_handles)
% SET_TRACES_BLACK Recolor all line traces on the given axes to black.
for i = 1:numel(ax_handles)
    lines = findobj(ax_handles(i), 'Type', 'line');
    set(lines, 'Color', 'k', 'LineWidth', 1.5);

end
end

function fig_combined = combine_columns(ax_cols, titles)
% COMBINE_COLUMNS Copy N column vectors of axes into a single figure,
% arranged as an N-column tiledlayout (one column per condition). Panels
% are assigned to a fixed row by their ylabel (stim, dendrite, firing
% rate, synaptic output, adaptation, depression), so e.g. row 5 is always
% the SFA/adaptation panel and row 6 is always the STD/depression panel,
% regardless of which panels a given condition happens to have.
row_order = {'stim', 'dendrite', 'firing rate', 'synaptic output', 'adaptation', 'depression'};
n_cols = numel(ax_cols);
n_rows = numel(row_order);
fig_combined = figure();
tl = tiledlayout(fig_combined, n_rows, n_cols);

no_adapt_col = find(strcmpi(titles, 'No adaptation'), 1);
if isempty(no_adapt_col)
    no_adapt_col = 1;
end

first_ax = gobjects(n_cols, 1);
for c = 1:n_cols
    ax_new = gobjects(numel(ax_cols{c}), 1);
    for i = 1:numel(ax_cols{c})
        row = find(strcmp(get(get(ax_cols{c}(i), 'YLabel'), 'String'), row_order));
        ax_new(i) = copyobj(ax_cols{c}(i), tl);
        ax_new(i).Layout.Tile = (row - 1) * n_cols + c;
        ax_new(i).YAxis.LineWidth = 1.0;
        ylim(ax_new(i), [-0.05, 1]);
        yticks(ax_new(i), [0, 1]);
        ax_new(i).YLabel.FontSize = 12;
    end
    linkaxes(ax_new, 'x');
    title(ax_new(1), titles{c}, 'FontSize', 16, 'FontWeight', 'normal');
    if isempty(first_ax(c)) || ~isgraphics(first_ax(c))
        first_ax(c) = ax_new(1);
    end
end

% Single shared scale bar: 5 seconds, in row 6 ("depression" row) under
% "No adaptation" -- an empty, axis-less tile since that condition has no
% panel there. Its x-limits are matched to the rest of the column.
scale_row = find(strcmp('adaptation', row_order));
scale_tile = (scale_row - 1) * n_cols + no_adapt_col;
ax_scale = nexttile(tl, scale_tile);
axis(ax_scale, 'off');
xlim(ax_scale, xlim(first_ax(no_adapt_col)));
add_scale_bar(ax_scale, 10);
end

function add_scale_bar(ax, length_sec)
% ADD_SCALE_BAR Draw a horizontal time scale bar in the lower-left of ax.
hold(ax, 'on');
xlims = xlim(ax);
ylim(ax, [0, 1]);
x_start = xlims(1) + 0.05 * (xlims(2) - xlims(1));
x_end = x_start + length_sec;
y_pos = 0.6;
plot(ax, [x_start, x_end], [y_pos, y_pos], 'k-', 'LineWidth', 4);
text(ax, (x_start + x_end) / 2, y_pos - 0.08, ...
    sprintf('%g seconds', length_sec), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold(ax, 'off');
end

function [u_ex, t_ex] = single_neuron_step_generator(params, T, fs, t_on, t_off, step_amp)
% SINGLE_NEURON_STEP_GENERATOR Deterministic single step input for 1 neuron.
dt = 1 / fs;
t_ex = (0:dt:T)';
u_ex = zeros(params.n, length(t_ex));
on_mask = (t_ex >= t_on) & (t_ex < t_off);
u_ex(:, on_mask) = step_amp;
end
