close all; clear all; clc;
setup_paths();

save_figs = true;
save_workspace = false;

level_of_chaos = 1;

%% Run 1: No Adaptation, No Depression (DDE Baseline)
close all;
note = 'DDE_paired_no_adapt';
u_ex_scale = 1;
n_a_E = 0;
n_b_E = 0;

save_dir = fullfile(pwd, 'results_DDE', note);

fprintf('Running DDE SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); 
full_SRNN_run_DDE(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

% %% Run 2: With Adaptation (3 timescales)
% close all;
% note = 'DDE_paired_SFA_3TS';
% u_ex_scale = 1;
% n_a_E = 3;
% n_b_E = 0;
% 
% save_dir = fullfile(pwd, 'results_DDE', note);
% 
% fprintf('Running DDE SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
% clear_SRNN_persistent(); 
% full_SRNN_run_DDE(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

