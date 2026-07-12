% test_single_neuron_stf.m - Single-neuron demo of SFA vs STD vs STF
%
% Companion to test_single_neuron_adaptation.m, extended to short-term
% synaptic facilitation (STF). A lone neuron (no recurrent connectivity,
% W = 0) is driven with a step in external input u under 7 conditions:
% no adaptation, SFA only, STD only, STF only, SFA+STD, STD+STF, and
% SFA+STD+STF. Each column shows the stimulus, dendritic state x, firing
% rate r, synaptic output eff.*r, and the mechanism variables that are
% active: adaptation a (SFA), depression b (STD), facilitation p (STF).
% Facilitation raises the gain eff=(p/p0).*b above 1, so the synaptic-output
% row uses its own (across-column uniform) y-axis; all other rows are [0,1].
%
% Unlike the SFA/STD reference (which uses SRNNModel2), STF only exists in
% SRNNModelCellTypes, so this script integrates the committed static kernel
% SRNNModelCellTypes.dynamics_fast_ct directly with a hand-built n=1, K=1
% params struct (W = 0), and reconstructs traces with unpack_states_ct.
% This keeps the figure faithful to the committed equations, including the
% new p-coupled STD depletion  db/dt = (1-b)/tau_rec - (p.*b.*r)/tau_rel
% and gain  eff = (p/p0).*b.

close all; clear; clc;

setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% Shared step stimulus and model settings (mirror the SFA/STD reference)
step_amp = 0.8;        % external drive during the step (mid-range, avoids saturating r)
t_on     = 5;          % step onset (s)
t_off    = 15;         % step offset (s); long enough for slow SFA/STF to develop
T_range  = [-10, 20];  % simulate before/after the step
fs       = 400;        % sample rate (Hz)
dt       = 1 / fs;
t_ex     = (T_range(1):dt:T_range(2))';

% Deterministic single-neuron step, as an interpolant the kernel can query.
u = single_neuron_step(t_ex, t_on, t_off, step_amp);   % 1 x nt
u_interpolant = griddedInterpolant(t_ex, u', 'linear', 'none');

% Shared (condition-independent) parameters.
base = struct();
base.n     = 1;
base.K     = 1;
base.W     = 0;                 % no recurrence
base.tau_d = 0.1;               % dendritic time constant (base default)
base.type_of = 1;
base.activation_function            = @(x) SRNNModelBase.piecewiseSigmoid(x, 1.0, 0.35);
base.activation_function_derivative = @(x) SRNNModelBase.piecewiseSigmoidDerivative(x, 1.0, 0.35);
base.u_interpolant = u_interpolant;

% Illustrative mechanism values (exaggerated for a clear, readable demo).
vals = struct();
vals.tau_a        = 3;      % SFA timescale (s)
vals.c            = 1.0;    % SFA strength
vals.tau_b_rec    = 2;      % STD recovery (s)
vals.tau_b_rel    = 0.3;    % STD depletion timescale (s)
vals.p0           = 0.35;   % baseline release prob (STF rest / floor)
vals.tau_f        = 6;      % STF facilitation decay (s)
vals.kappa        = 0.4;    % STF facilitation rate

%% Conditions: [n_a n_b n_u] per column
condition_titles = {'No adaptation', 'SFA only', 'STD only', 'STF only', ...
                    'SFA + STD', 'STD + STF', 'SFA + STD + STF'};
condition_flags = [ ...
    0 0 0; ...   % No adaptation
    1 0 0; ...   % SFA only
    0 1 0; ...   % STD only
    0 0 1; ...   % STF only
    1 1 0; ...   % SFA + STD
    0 1 1; ...   % STD + STF
    1 1 1];      % SFA + STD + STF
n_cond = numel(condition_titles);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', dt);

%% Integrate every condition and collect display traces
row_order = {'stim', 'dendrite', 'firing rate', 'synaptic output', ...
             'adaptation', 'depression', 'facilitation'};
n_rows = numel(row_order);
traces = cell(n_cond, n_rows);   % traces{c, r} = [] if that mechanism is inactive

keep = t_ex >= 0;                % display window (trim negative settling)
t_disp = t_ex(keep);

for c = 1:n_cond
    n_a = condition_flags(c, 1);
    n_b = condition_flags(c, 2);
    n_u = condition_flags(c, 3);
    params = build_params(base, n_a, n_b, n_u, vals);

    % Initial state: a=0, b=1, p=p0, x=0  (S = [a(:); b(:); p(:); x])
    S0 = [zeros(params.n * n_a, 1); ...
          ones(params.n * params.K * n_b, 1); ...
          repmat(params.p0_mat(:), n_u, 1); ...
          0];

    [~, S] = ode45(@(t, S) SRNNModelCellTypes.dynamics_fast_ct(t, S, params), t_ex, S0, opts);
    st = SRNNModelCellTypes.unpack_states_ct(S, params);   % x,a,b,p,r,br as 1 x (K x) nt

    traces{c, 1} = u(keep);                          % stim
    traces{c, 2} = st.x(1, keep);                    % dendrite
    traces{c, 3} = st.r(1, keep);                    % firing rate
    traces{c, 4} = st.br(1, keep);                   % synaptic output (eff.*r)
    if n_a > 0, traces{c, 5} = reshape(st.a(1, 1, keep), 1, []); end   % adaptation a
    if n_b > 0, traces{c, 6} = reshape(st.b(1, 1, keep), 1, []); end   % depression b
    if n_u > 0, traces{c, 7} = reshape(st.p(1, 1, keep), 1, []); end   % facilitation p
end

%% Assemble the combined figure (rows x conditions)
fig = figure('Color', 'w', 'Position', [50, 50, 1600, 950]);
tl = tiledlayout(fig, n_rows, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');

no_adapt_col = find(strcmpi(condition_titles, 'No adaptation'), 1);
scale_row = find(strcmp('adaptation', row_order), 1);

% Leftmost column that has a trace in each row (label only there, so the
% repeated row labels don't collide across the narrow multi-column layout).
label_col = zeros(1, n_rows);
for r = 1:n_rows
    populated = find(~cellfun(@isempty, traces(:, r)), 1);
    if ~isempty(populated), label_col(r) = populated; end
end

% Per-row y-limits. All rows are normalized to [0,1] EXCEPT synaptic output:
% facilitation drives the gain (p/p0) above 1, so that row gets its own
% uniform (across-column) axis to show enhancement honestly instead of
% clipping it. STD-only output stays <1; STF-only output rises >1.
syn_row = find(strcmp('synaptic output', row_order), 1);
syn_max = max(cellfun(@(y) max([y, 0]), traces(:, syn_row)));
syn_top = max(1, ceil(syn_max / 0.5) * 0.5);      % round up to a clean 0.5
row_ylim   = repmat([-0.05, 1], n_rows, 1);
row_yticks = repmat({[0, 1]}, n_rows, 1);
row_ylim(syn_row, :) = [-0.05 * syn_top, syn_top];
row_yticks{syn_row}  = unique([0, 1, syn_top]);

for c = 1:n_cond
    for r = 1:n_rows
        y = traces{c, r};
        if isempty(y), continue; end
        ax = nexttile(tl, (r - 1) * n_cond + c);
        plot(ax, t_disp, y, 'k-', 'LineWidth', 1.5);
        xlim(ax, [t_disp(1), t_disp(end)]);
        ylim(ax, row_ylim(r, :));
        yticks(ax, row_yticks{r});
        ax.YAxis.LineWidth = 1.0;
        ax.XAxis.Visible = 'off';
        box(ax, 'off');
        if c == label_col(r)
            ylabel(ax, row_order{r}, 'FontSize', 12);
        end
        if r == 1
            title(ax, condition_titles{c}, 'FontSize', 16, 'FontWeight', 'normal');
        end
    end
end

% Single shared scale bar in an empty tile (adaptation row, No-adaptation col).
ax_scale = nexttile(tl, (scale_row - 1) * n_cond + no_adapt_col);
axis(ax_scale, 'off');
xlim(ax_scale, [t_disp(1), t_disp(end)]);
add_scale_bar(ax_scale, 10);

%% Save
save_some_figs_to_folder_2(this_dir, 'sfa_std_stf_single_neuron_example', fig.Number, {'fig', 'png', 'svg'});

%% Local functions
function u = single_neuron_step(t_ex, t_on, t_off, step_amp)
% SINGLE_NEURON_STEP Deterministic single step input (1 x nt) over t_ex.
u = zeros(1, numel(t_ex));
u((t_ex >= t_on) & (t_ex < t_off)) = step_amp;
end

function params = build_params(base, n_a, n_b, n_u, vals)
% BUILD_PARAMS Assemble the dynamics_fast_ct params struct for one condition
% (n=1, K=1). Only the active-mechanism fields matter, but all are filled.
params = base;
params.n_a = n_a; params.n_b = n_b; params.n_u = n_u;
params.N_sys_eqs = params.n * n_a + params.n * params.K * n_b + ...
                   params.n * params.K * n_u + params.n;
params.tau_a = vals.tau_a;                 % 1 x n_a (scalar broadcasts)
params.c_vec = vals.c * ones(params.n, 1); % n x 1
params.tau_b_rec_mat = vals.tau_b_rec * ones(params.n, params.K);
params.tau_b_rel_mat = vals.tau_b_rel * ones(params.n, params.K);
params.p0_mat        = vals.p0     * ones(params.n, params.K);
params.tau_f_mat     = vals.tau_f  * ones(params.n, params.K);
params.kappa_mat     = vals.kappa  * ones(params.n, params.K);
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
