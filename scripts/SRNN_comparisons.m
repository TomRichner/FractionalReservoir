close all; clear all; clc;
setup_paths();

save_figs = true;
save_workspace = false;

level_of_chaos = 2;

%% Run 0: no stim
% note = 'no_stim_no_adapt_no_depress';
% u_ex_scale = 0;
% n_a_E = 0;
% n_b_E = 0;
% 
% save_dir = fullfile('/Users/tom/Desktop/local_code/FractionalReservoir', 'results', note);
% 
% fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);



%% Run 1: no adapt no depression
close all;
note = 'no_adapt_no_depress';
u_ex_scale = 1;
n_a_E = 0;
n_b_E = 0;

save_dir = fullfile('/Users/tom/Desktop/local_code/FractionalReservoir', 'results', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

% %% Run 2
% close all;
% note = 'SFA_1_timescale';
% u_ex_scale = 1;
% n_a_E = 1;
% n_b_E = 0;
% 
% save_dir = fullfile(pwd, 'results', note);
% 
% fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 3
close all;
note = 'SFA_3TS';
u_ex_scale = 1;
n_a_E = 3;
n_b_E = 0;

save_dir = fullfile(pwd, 'results', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 4
close all;
note = 'STD_and_3TS_SFA';
u_ex_scale = 1;
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile(pwd, 'results', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);