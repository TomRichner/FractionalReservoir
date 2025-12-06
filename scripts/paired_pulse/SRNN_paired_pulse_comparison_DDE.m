close all; clear all; clc;
setup_paths();

save_figs = false;
save_workspace = false;

level_of_chaos = 1.1; % Adjusted for DDE stability usually lower

u_ex_scale = 20;

%% Run 1: no adapt no depression (DDE)
close all;
note = 'PP_DDE_no_adapt_no_depress';
n_a_E = 0;
n_b_E = 0;

save_dir = fullfile(pwd, 'results_DDE', note);

fprintf('Running DDE SRNN paired pulse with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); 
full_SRNN_paired_pulse_DDE(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 2: SFA 3TS (DDE)
% close all;
% note = 'PP_DDE_SFA_3TS';
% n_a_E = 3;
% n_b_E = 0;
% 
% save_dir = fullfile(pwd, 'results_DDE', note);
% 
% fprintf('Running DDE SRNN paired pulse with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); 
% full_SRNN_paired_pulse_DDE(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 3: STD and 3TS SFA (DDE)
% close all;
% note = 'PP_DDE_STD_and_3TS_SFA';
% n_a_E = 3;
% n_b_E = 1;
% 
% save_dir = fullfile(pwd, 'results_DDE', note);
% 
% fprintf('Running DDE SRNN paired pulse with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); 
% full_SRNN_paired_pulse_DDE(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

