close all; clear all; clc;
setup_paths();

save_figs = false;
save_workspace = false;

level_of_chaos = 1.8;

u_ex_scale = 1.7;
rng_seeds = [8 16 3 4 5];

time_config.T_range = [-15, 30];
time_config.T_plot = [0, 30];

combined_runs = {};

%% Run 1: no adapt no depression
close all;
note = 'review_no_adapt_ex';

n_a_E = 0;
n_b_E = 0;

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review', note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); 
[~, ~, params_1, lya_1, plot_data_1] = full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config);

run1.plot_data = plot_data_1;
run1.params = params_1;
run1.lya_results = lya_1;
run1.Lya_method = 'benettin';
combined_runs{1} = run1;

%% Run 4
close all;
note = 'review_STD_and_3TS_SFA_ex';
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review', note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); 
[~, ~, params_4, lya_4, plot_data_4] = full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config);

run4.plot_data = plot_data_4;
run4.params = params_4;
run4.lya_results = lya_4;
run4.Lya_method = 'benettin';
combined_runs{2} = run4;

%% Plot Combined
plot_SRNN_combined_tseries(combined_runs, 5, {'u_ex', 'x', 'br', 'a', 'b', 'lya'});
