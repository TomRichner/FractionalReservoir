close all; clear all; clc;
setup_paths();

save_figs = true;
save_workspace = true;



u_ex_scale = 0;

%% Run 1: no adapt no depression
close all;
note = 'chaotic';
n_a_E = 0;
n_b_E = 0;
level_of_chaos = 7;

save_dir = fullfile('/Users/tom/Desktop/local_code/FractionalReservoir', 'results', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); minimal_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 3
close all;
note = 'edge';
n_a_E = 1;
n_b_E = 1;
level_of_chaos = 2.3

save_dir = fullfile('/Users/tom/Desktop/local_code/FractionalReservoir', 'results', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); minimal_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

%% Run 4
close all;
note = 'overly_stable';
n_a_E = 1;
n_b_E = 1;
level_of_chaos = 2.1;

save_dir = fullfile('/Users/tom/Desktop/local_code/FractionalReservoir', 'results', note);

fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);
clear_SRNN_persistent(); minimal_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note);

