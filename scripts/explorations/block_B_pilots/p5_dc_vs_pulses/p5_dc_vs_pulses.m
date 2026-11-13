% P5_DC_VS_PULSES Block B pilot 5: DC input versus pulse trains of the same mean, on the bursting network.
%
% QUESTION. At the same time-averaged input, do a DC step and a pulse
% train have the same effect on the spontaneous population bursts and on
% the recruited adaptation, under chained multiple-timescale (MTS)
% adaptation?
%
% WHAT IS RUN. The bursting network (preset bursting_pairs, 50 neurons),
% condition 'chain_mts' from bursting_condition_model (the one definition of
% the conditions), seed [42 42], no noise, rk4. 10 s of settling, then
% 100 s: no input from 0 to 40 s (the baseline, inside the same run), input
% ON from 40 to 100 s, the same on every neuron and in phase. Two mean input
% levels (0.1, 0.2) and four forms at each, all with the same mean over the
% on period:
%   DC        the mean, constant
%   2 Hz      rectangular pulses, 20 % duty (100 ms wide), amplitude 5 x mean
%   10 Hz     20 % duty (20 ms wide)
%   50 Hz     20 % duty (4 ms wide)
% plus the 'original' condition (sfa3_std1) under DC and 2 Hz at both levels,
% as a reference: 12 simulations in all. The input is made by
% pulse_train_input (this folder), a generator for the model's
% input_config.generator.
%
% THE SAMPLING RATE. The model integrates at fs = 400 Hz by default, where
% a 50 Hz pulse at 20 % duty would be 1.6 samples wide. Every case here runs
% at fs = 2000 Hz instead (model_args {'fs', 2000}), so the pulses are 200,
% 40 and 8 samples wide and the periods 1000, 200 and 40 samples. That is
% the one change from pilot 0b's run (fs 400), and its baseline differs from
% 0b's at the level of the step size.
%
% MEASURES, before (10-40 s) and during (50-100 s, leaving 10 s for the
% transient): bursts per second (pilot 0b's definition: population mean E
% rate above 0.2 for at least 0.1 s, dips under 0.2 s merged; restated in
% p5_dc_vs_pulses_compute because 0b's is a local function), mean E and I
% rate, the largest population mean E rate, the mean adaptation term
% (c/K) sum_k a_k over E and the mean depression factor over the depressed
% routes (SRNNCellTypePairs.std_factors, which calls std_route_factor), and
% locking to the pulses: the share of the mean E rate's power at the pulse
% frequency (and at it plus its harmonics), the peak-to-trough swing of the
% pulse-triggered average, and for the bursts the vector strength of their
% onset phase within the pulse cycle and the number of bursts per pulse.
%
% Run setup_paths once in the session first (cwd at the repository root).
% Takes about 12.5 minutes with no pool (twelve 110-s runs at 2000 Hz).
% The figures are drawn invisible and shown at the end; nothing is saved
% unless save_figs is set (then out_dir receives the two chain MTS figures,
% the original's and the summary as .png and .fig, traces.mat and
% measures.csv).
%
% WHAT TO EXPECT with the settings below, on the code as of 2026-11-12
% (main at b051af1). Deterministic, so a rerun gives these exactly.
% Before (10-40 s, identical in every run of a condition): chain MTS
% 0.333 bursts/s, mean r_E 0.0601; original 0.333 bursts/s, r_E 0.0815.
% During (50-100 s), bursts/s | mean r_E | mean depression | mean adaptation:
%   chain MTS  mean 0.1  DC     0 | 0.0969 | 0.631 | 0.0482
%                        2 Hz   2 | 0.1179 | 0.565 | 0.0585  (one event per pulse)
%                        10 Hz  0 | 0.0981 | 0.627 | 0.0488
%                        50 Hz  0 | 0.0969 | 0.631 | 0.0482
%   chain MTS  mean 0.2  DC     0 | 0.1333 | 0.560 | 0.0661
%                        2 Hz   2 | 0.1663 | 0.459 | 0.0823  (one event per pulse)
%                        10 Hz  0 | 0.1364 | 0.551 | 0.0676
%                        50 Hz  0 | 0.1334 | 0.560 | 0.0662
%   original   mean 0.1  DC     0 | 0.1704 | 0.673 | 0.0844
%                        2 Hz   2 | 0.1816 | 0.652 | 0.0899
%   original   mean 0.2  DC     0 | 0.2166 | 0.629 | 0.1071
%                        2 Hz   2 | 0.2460 | 0.572 | 0.1214
% Every generated input has the stated mean over 40-100 s to 1e-14.
%
% See also: pulse_train_input, p5_dc_vs_pulses_compute,
%           bursting_condition_model, bursting_chain_mts_figure

%% Settings
S = struct();
S.conditions        = {'chain_mts', 'original'};
S.original_forms_Hz = [0 2];          % the original is run under these forms only
S.mean_levels       = [0.1 0.2];      % time-averaged input while on
S.forms_Hz          = [0 2 10 50];    % 0 = DC
S.duty              = 0.2;            % pulse duty cycle
S.rng_seeds         = [42 42];
S.T_settle          = 10;             % s before t = 0
S.T_show            = 100;            % s; input on from t_on to T_show
S.t_on              = 40;             % s
S.fs                = 2000;           % Hz, integration rate (see header)
S.plot_deci         = 4;              % 2000 Hz -> 500 Hz plot_data
S.win_before        = [10 40];        % s
S.win_during        = [50 100];       % s
S.zoom              = [60 62];        % s, the input row of the figures
S.burst_thresh      = 0.2;            % pilot 0b's burst definition
S.burst_min_dur     = 0.1;
S.burst_merge       = 0.2;
S.save_figs         = false;
S.out_dir           = '';             % e.g. fullfile(fileparts(which('setup_paths')), 'figs', 'explore_B', 'p5_dc_vs_pulses')

%% Run
out = p5_dc_vs_pulses_compute(S);
disp(out.table(:, {'condition', 'mean_input', 'pulse_Hz', 'bursts_per_s_before', ...
    'bursts_per_s_during', 'mean_rE_during', 'mean_dep_during', 'mean_adapt_during'}));
for f = out.figs, set(f, 'Visible', 'on'); end
