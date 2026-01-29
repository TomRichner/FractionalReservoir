close all; clear all; clc;
setup_paths();

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 1.25);
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTitleFontWeight', 'normal');

save_figs = false;
save_workspace = false;

level_of_chaos = 1.0;

% u_ex_scale = 1.6;
u_ex_scale = 0.75;


% rng_seeds = [105 25];
% rng_seeds = [10 10];
% rng_seeds = [100 100];
rng_seeds = [1 1];

time_config.T_range = [-25, 45];
time_config.T_plot = [7, 45];  % Cutoff halfway through second no-stim period
time_config.J_periods = [false true true];  % Only plot J_eff/eigenspectrum for first two periods

combined_runs = {};

% Run 1: no adapt no depression
close all;
note = 'Review_no_adapt_ex';

n_a_E = 3;
n_b_E = 0;

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review', note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent();
[~, ~, params_1, lya_1, plot_data_1] = full_SRNN_run_SRNNModel(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config);

run1.plot_data = plot_data_1;
run1.params = params_1;
run1.lya_results = lya_1;
run1.Lya_method = 'benettin';
combined_runs{1} = run1;

%% Run 4
% close all;
note = 'Review_STD_and_3TS_SFA_ex';
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review', note);
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
AddLetters2Plots(fig_handle, {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'}, 'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.065, 'VShift', -0.025);

ylim([-1.8 1.8])

if save_figs
    save_dir_combined = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'revised_results_review');
    save_name_base = 'combined_comparison_v3';

    % Use the existing helper function
    save_some_figs_to_folder_2(save_dir_combined, save_name_base, [fig_handle.Number], {'fig', 'svg', 'png'});
    fprintf('Combined plot saved to %s\n', save_dir_combined);
end
