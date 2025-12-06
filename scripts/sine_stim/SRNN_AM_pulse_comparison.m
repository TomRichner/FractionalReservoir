close all; clear all; clc;
setup_paths();

save_figs = true;
save_workspace = true;

level_of_chaos = 2.5;
u_ex_scale = 20;

% %% Run 1: no adapt no depression
% close all;
% note = 'AM_no_adapt_no_depress';
% n_a_E = 0;
% n_b_E = 0;
% 
% save_dir = fullfile(pwd, 'results', note);
% 
% fprintf('Running SRNN AM pulse with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); 
% full_SRNN_AM_pulse(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);
% 
% %% Run 2: SFA 3TS
% close all;
% note = 'AM_SFA_3TS';
% n_a_E = 3;
% n_b_E = 0;
% 
% save_dir = fullfile(pwd, 'results', note);
% 
% fprintf('Running SRNN AM pulse with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); 
% full_SRNN_AM_pulse(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 3: STD and 3TS SFA
close all;
note = 'AM_STD_and_3TS_SFA';
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile(pwd, 'results', note);

fprintf('Running SRNN AM pulse with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); 
full_SRNN_AM_pulse(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

