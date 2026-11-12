% RUN_P3_ABLATION Block B pilot 3: the unstimulated networks with adaptation removed.
%
% QUESTION. Block B's causal experiment is a 2 x 2, stimulation off / on
% crossed with adaptation intact / ablated. Removing adaptation also changes
% the baseline. What does each network do, unstimulated, with adaptation
% removed, and is there a way of removing it that leaves a baseline with
% comparable aberrant activity?
%
% WHAT IT DOES. For the bursting network (bursting_condition_model, 50
% neurons, no noise) and the IED-like run13 network (run13_condition_model,
% 350 neurons, noise 0.025), on the same seed, input (none) and noise path
% as the intact run, 10 s of settling and 100 s analysed:
%   chained MTS   intact; SFA removed; depression removed; both removed;
%                 mean-matched (both removed, each adapting variable
%                 replaced by a constant equal to its time average in the
%                 intact run): with the plain average of the depression
%                 factor, with its rate-weighted average, with four
%                 constants between the two, and (run13) per neuron; the
%                 slow parts removed (fastest SFA rung and first depression
%                 pool only), with the total SFA strength kept or with only
%                 the fast rung's share (c/K)
%   original      intact; SFA removed; depression removed; mean-matched
% The ablations are defined in p3_ablation_args (read its header) and reach
% the builders only through their 'model_args'. No model code is touched.
%
% MEASURES (p3_measures): bursts (bursting network: mean E rate above 0.2
% for at least 0.1 s) or detector-v1 events (run13) per minute, with size
% and duration; mean E and I rate; the fraction of time the mean E rate is
% above 0.5; the fraction of neurons with time-averaged rate above 0.9
% (pinned) and below 0.01 (silent). The REGIME word is assigned by rule:
% saturation if at least half the neurons are pinned or the mean E rate is
% above 0.5 at least half the time; continuous high firing if the mean E
% rate is above 0.2 at least half the time; with no events, quiet if the
% mean E rate is below half the intact run's and "steady, no events"
% otherwise; else comparable (0.5 to 2 times the count of the same
% condition's intact run), more, or fewer. A saturated bursting network
% counts as ONE burst that never ends.
%
% Run setup_paths once in the session first (cwd at the repository root).
% About 8 minutes, serial, no pool (33 simulations of 110 s). Nothing is
% saved unless save_figs is set.
%
% WHAT TO EXPECT with the settings below (main at a33fa6f, 2026-11-12);
% count per 100 s, mean E rate, regime:
%   bursting, chained MTS          run13, chained MTS
%   intact            34  0.061    intact             12  0.091
%   SFA removed        0  0.098    SFA removed        11  0.104  comparable
%     steady, no events
%   depression removed    0.829    depression removed  0  0.997  saturation
%     saturation
%   both removed          0.923    both removed        0  0.999  saturation
%     saturation
%   mean-matched          0.769    mean-matched        0  0.965  saturation
%     saturation (scale 0.752)       (scale 0.585)
%   mean-matched, rate-weighted    mean-matched, rate-weighted
%                      0  0.003                        0  0.138  steady
%     quiet (scale 0.527)            (scale 0.332)
%   scales between: quiet up to    scales between: 3 events and 5 % of
%     0.707 (bistable network)       neurons pinned at 0.382, then
%                                    module after module saturates
%                                  mean-matched, per neuron
%                                                      1  0.550  some
%                                    modules saturated
%   slow parts removed 70  0.104   slow parts removed  30  0.125  (this IS
%     more                           the original condition, exactly)
%   slow removed, c/K   0  0.179   slow removed, c/K   23  0.134
%   bursting, original             run13, original
%   intact            33  0.081    intact             30  0.125
%   SFA removed           0.238    SFA removed        17  0.142  comparable
%     continuous high firing
%   depression removed    0.829    depression removed  0  0.997  saturation
%   mean-matched          0.794    mean-matched        0  0.974  saturation
% One seed fixes everything, so a rerun gives these numbers exactly.
%
% SAVED with save_figs, in out_dir:
%   <network>_ablations.png / .fig   chained MTS: intact and the ablations
%   <network>_variants.png / .fig    the mean-matched variants, c/K, and
%                                    the original condition's ablations
%   summary.png / .fig               bars of events per minute and mean E rate
%   results.csv                      the table, one row per run (with the
%                                    constants each mean match was given)
%   results.mat                      T (the table), runs (per run: mean
%                                    rates over time, time-averaged rate of
%                                    every neuron, mean adaptation term and
%                                    depression factor, bursts or events,
%                                    the model_args, the constant input and
%                                    the intact run's averages), settings
%
% See also: p3_ablation_pilot, p3_ablation_args, p3_run_one, p3_measures,
%           p3_figure, p3_summary_figure, p3_check_offset,
%           bursting_condition_model, run13_condition_model

%% Settings
networks    = {'bursting', 'run13'};
T_show      = 100;      % s analysed, after 10 s of settling
seed        = 42;       % both networks' seed (the unstimulated pilots')
fig_visible = true;     % false: draw without windows
save_figs   = false;    % true: write the figures and the results to out_dir
out_dir     = '';       % '' = figs/explore_B/p3_ablation under the repo root

%% Run every ablation, measure, draw, save
out = p3_ablation_pilot('networks', networks, 'T_show', T_show, 'seed', seed, ...
    'fig_visible', fig_visible, 'save_figs', save_figs, 'out_dir', out_dir);
