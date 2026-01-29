close all; clear all; clc;
setup_paths();

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextFontSize', 14);
set(groot, 'DefaultLineLineWidth', 1.25);
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTitleFontWeight', 'normal');

save_figs = false;
save_workspace = false;

level_of_chaos = 1.0; % gamma from somplinsky

u_ex_scale = 1; % can change scale of stimulus.

rng_seeds = [42 42]; % M.O.L., seed 1 is for connection matrix, seed 2 is for sitmulus. Setting these ensures comparison between SFA vs SFA+STD use the same connection matrix and stimulus

time_config.J_periods = [false true true];  % there are three periods: no-stim, stim, no-stim.  The first no-stim period is ignored for stats and plots
time_config.T_range = [-15, 45]; % seconds. simulation time starts at -25 seconds to allow for IC transient to settle and local lyapunov vector to alignt to dynamcis
time_config.T_plot = [7.5, 45];  % seconds. time period to be included in plot. Trim off half of the first no-stim period since it is not used for stats

combined_runs = {}; % cell array to hold multiple simulations

%% Run 1: spike frequency adaptation (SFA) only
close all;
note = 'SFA_only';

n_a_E = 3; % three timeconstants of SFA
n_b_E = 0; % no short-term synaptic depression

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review_revised', note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent();
[~, ~, params_1, lya_1, plot_data_1] = full_SRNN_run_SRNNModel(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config);

run1.plot_data = plot_data_1;
run1.params = params_1;
run1.lya_results = lya_1;
run1.Lya_method = 'benettin'; % reshooting method to get the largest lyapunov exponent
combined_runs{1} = run1;

%% Run 4
% close all;
note = 'STD_and_SFA';
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review_revised', note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent();
[~, ~, params_4, lya_4, plot_data_4] = full_SRNN_run_SRNNModel(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config);
run4.plot_data = plot_data_4;
run4.params = params_4;
run4.lya_results = lya_4;
run4.Lya_method = 'benettin';
combined_runs{2} = run4;

%% Plot Combined
[fig_handle, ~] = plot_SRNN_combined_tseries(combined_runs, 3, {'u_ex', 'x', 'br', 'a', 'b', 'lya'});

% Add letters to subplots
AddLetters2Plots(fig_handle, {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'}, 'FontSize', 16, 'FontWeight', 'normal', 'HShift', -0.06, 'VShift', -0.02);

ylim([-1.9 1.9]) % y limits of the local lyapunov exponent

if save_figs
    save_dir_combined = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review_revised');
    save_name_base = 'combined_comparison';

    % Use the existing helper function
    save_some_figs_to_folder_2(save_dir_combined, save_name_base, [fig_handle.Number], {'fig', 'svg', 'png'});
    fprintf('Combined plot saved to %s\n', save_dir_combined);
end
