% EXAMPLE_BURSTING_CHAIN_MTS The bursting network over 100 s: original adaptation, chain MTS, product MTS.
%
% The exploratory look Tom asked for at the start of planning block B
% (2026-11-12): the small, strongly coupled bursting network of the
% manuscript's stimulation figure (preset bursting_pairs, the port of the
% deleted bursting_SRNN_example.m), run with no input and no noise for 100 s
% after 10 s of settling, under three conditions on the SAME network, seed
% and (zero) input. One figure, one column per condition:
%
%   row 1  rates, neurons by time, E (1-35) above I (36-50), white = 1
%   row 2  population mean rate, E and I; bursts shaded, threshold dotted
%   row 3  mean adaptation term (c/K) sum_k a_k over the E neurons
%   row 4  mean depression factor over the depressed routes (E->E, E->I),
%          from SRNNCellTypePairs.std_route_factor
% Rows share their y axis across columns.
%
% THE CONDITIONS are built by bursting_condition_model (in this folder; its
% header lists every value): 'original' is the preset's sfa3_std1; the two
% MTS conditions keep the original's SFA (it already has three rungs,
% 0.25 / 1.58 / 10 s on E, c = 0.5 total) and its depression routes and
% constants (tau_rec 1 s, tau_rel 0.25 s), and add only what the MTS forms
% need, at the mu7chain bundle's current defaults:
%   chain MTS    tau_rec [1 5] s, tau_rel 0.25 s, n2 = 2
%   product MTS  tau_rec [1 4] s, tau_rel [0.25 0.5] s
%
% A BURST is the population mean E rate above 0.2 for at least 0.1 s, dips
% shorter than 0.2 s not splitting it. A neuron is ACTIVE when r_i > 0.1;
% burst size is the fraction of all 50 neurons active at the burst's peak.
%
% Run setup_paths once in the session first (cwd at the repository root);
% like every script in scripts/examples this file has no path code. It takes
% about 20 s. Nothing is saved unless save_figs is set.
%
% WHAT TO EXPECT with the settings below, on the code as of 2026-11-12
% (main at eddd8ee):
%   original (sfa3_std1)  33 bursts (0.33 /s), period about 3.0 s, median
%                         duration 0.46 s, 82 % of neurons active at the
%                         peak; mean r_E 0.0807, r_I 0.0866
%   chain MTS             34 bursts (0.34 /s), period about 2.95 s, median
%                         duration 0.38 s, 70 % active at the peak;
%                         mean r_E 0.0613, r_I 0.0628
%   product MTS           0 bursts: a small regular oscillation (mean E
%                         rate peaks near 0.09, period about 2 s);
%                         mean r_E 0.0370, r_I 0.0392
% The network is deterministic here (sigma_u_noise = 0, rk4), so a rerun
% gives these numbers exactly.
%
% SAVED with save_figs, in out_dir: bursting_chain_mts.fig (reopens
% visible), bursting_chain_mts.png, and traces.mat with one struct per
% condition, named original, chain_mts and product_mts, each holding
%   label             column title
%   t                 1 x nt, s, 100 Hz, [0, 100]
%   r_E, r_I          neurons x nt, per-neuron rate (E 35, I 15)
%   mean_r_E, mean_r_I         1 x nt population means
%   mean_adaptation   1 x nt, mean (c/K) sum_k a_k over E
%   mean_depression   1 x nt, mean depression factor over the routes in
%   depression_routes          (cell of 'E->E', 'E->I')
%   bursts            struct: n_bursts, burst_rate, t_on, t_off, peak_t,
%                     duration, frac_active_peak, medians, mean rates
%   values            struct of every value the condition was built with
%                     (bursting_condition_model's second output)
% plus burst_definition (the four thresholds above).
%
% TO TRY SOMETHING ELSE, edit the settings block: input_config is handed to
% the model as is ([] = no input), and model_args are further name-value
% overrides applied last in every condition, for example
%   model_args = {'sigma_u_noise', 0.025, 'ode_solver', 'sra1'};
%
% See also: bursting_condition_model, bursting_chain_mts_figure,
%           fig_stim_engages_adaptation, srnn_param_preset

%% Settings
rng_seeds    = [42 42];     % network seed pair (the stimulation figure's)
T_settle     = 10;          % s simulated before t = 0, not shown
T_show       = 100;         % s shown and analysed
input_config = [];          % [] = no external input
model_args   = {};          % further name-value overrides, applied last
save_figs    = false;       % true: write the .fig, .png and traces.mat to out_dir
out_dir      = '';          % e.g. fullfile(fileparts(which('setup_paths')), 'figs', 'explore_B', 'bursting_chain_mts')

%% Simulate the three conditions and draw the figure
out = bursting_chain_mts_figure( ...
    'rng_seeds',    rng_seeds, ...
    'T_settle',     T_settle, ...
    'T_show',       T_show, ...
    'input_config', input_config, ...
    'model_args',   model_args, ...
    'visible',      true, ...
    'save_figs',    save_figs, ...
    'out_dir',      out_dir);

%% Say what was found
for s = out.stats
    fprintf(['%-22s %2d bursts (%.2f /s), median duration %.2f s, median %2.0f %% active at peak, ' ...
        'mean r_E %.4f, r_I %.4f\n'], s.name, s.n_bursts, s.burst_rate, s.median_dur, ...
        100 * s.median_frac_peak, s.mean_rE, s.mean_rI);
end
