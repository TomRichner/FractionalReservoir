close all; clear all; clc;
setup_paths();

save_figs = false;
save_workspace = false;

level_of_chaos = 1.3;

u_ex_scale = 4;
rng_seeds = [8 8 3 4 5];

%% Run 0: no stim
% note = 'no_stim_no_adapt_no_depress';
% u_ex_scale = 0;
% n_a_E = 0;
% n_b_E = 0;
% 
% save_dir = fullfile('/Users/tom/Desktop/local_code/FractionalReservoir', 'results', note);
% 
% fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note);



% %% Run 1: no adapt no depression
% close all;
% note = 'review_no_adapt_ex';
% 
% n_a_E = 0;
% n_b_E = 0;
% 
% save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review', note);
% 
% fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note);

% %% Run 2
% close all;
% note = 'SFA_1_timescale';
% u_ex_scale = 1;
% n_a_E = 0;
% n_b_E = 1;
% 
% save_dir = fullfile(pwd, 'results', note);
% 
% fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note);

% %% Run 3
% close all;
% note = 'SFA_3TS';
% u_ex_scale = 1;
% n_a_E = 3;
% n_b_E = 0;
% 
% save_dir = fullfile(pwd, 'results', note);
% 
% fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note);

%% Run 4
close all;
note = 'review_STD_and_3TS_SFA_ex';
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile('/Users/richner.thomas/Desktop/local_code/FractionalResevoir/figs', 'results_review', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note);