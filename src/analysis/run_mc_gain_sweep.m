function [out_dir, results_file] = run_mc_gain_sweep(varargin)
% RUN_MC_GAIN_SWEEP Memory capacity and lambda_1 at every point of a sweep of gain against a second axis (the bias by default).
%
%   out_dir = RUN_MC_GAIN_SWEEP()                          % 'essential': a first look
%   out_dir = RUN_MC_GAIN_SWEEP('run_mode', 'fast', 'output_dir', d)
%   [out_dir, results_file] = RUN_MC_GAIN_SWEEP('gain_levels', g, 'bias_levels', b, ...
%       'seeds', s, 'conditions', {'no_adaptation', 'sfa3_stdc2'})
%   out_dir = RUN_MC_GAIN_SWEEP('gain_levels', g, 'second_axis', 'self_feedback', ...
%       'second_levels', k, 'seeds', s, 'output_dir', d)
%
% Block C (FD-010): is the memory-capacity result a property of adaptation
% or of one well-tuned network? For every adaptation condition of the
% preset, on several network seeds, this measures total memory capacity and
% the largest Lyapunov exponent over a grid of
%   gain          level_of_chaos, the multiplier on every connection weight
%   a second axis ONE other scalar model property, named by the option
%                 second_axis. The default is the bias, S_c: the center of
%                 the per-neuron setpoint draw (the per-neuron spread
%                 sigma_S_c stays as the preset has it). Any other settable
%                 property that takes a scalar can be named instead, for
%                 example 'self_coupling' (the diagonal of W, in multiples
%                 of the bulk radius at gain 1; a scalar there means every
%                 cell type). The stage holds no code for a particular
%                 property: mc_sweep_second_axis turns the property into the
%                 names used in the job files and the table, and
%                 mc_sweep_quantity labels it in the figures.
% run_memory_capacity measures the same task at the preset network only.
% Below, "bias" stands for the second axis wherever the default is described.
%
% THIS IS ONE RUN: one network preset, one second axis. It is a set of one
% for run_mc_gain_sweep_set, which holds the code and runs several such
% runs (a bundle's grids and variant networks) as one pool of jobs.
%
% ONE JOB is one (gain, second-axis value, seed). It builds every condition on the same
% seeds -- the same W, input weights and input sequence, only the adaptation
% differs -- so the conditions are PAIRED exactly as run_memory_capacity
% pairs its trials (verify_shared_build asserts it once, on the first job).
% The stimulus seed follows the network seed as it does there
% (seed + seed_stim_base - seed_net_base). The input weights use the
% reservoir's fixed rng_seed_input, so they are the same in every job.
% Per condition:
%   memory capacity   by SRNN_ESN_reservoir.run_memory_capacity_lean: the
%                     task, protocol (mc_protocol: 10 s washout, 600 s
%                     training, 150 s test, delays to 15 s, fs 200, 0.3 s
%                     holds, synaptic readout) and scoring of
%                     run_memory_capacity, integrated segment by segment so
%                     no state trajectory and no n x T input table is held.
%                     Bit-identical to the existing path
%                     (test_mc_lean_path), at tens of MB instead of GB.
%   lambda_1          by a separate short run() of the same reservoir with
%                     the top-K estimator (lya_K exponents, no retry;
%                     lya_T_sec seconds, accumulated over the last
%                     lya_acc_sec): once with the memory-task input on
%                     (lle_input, the value plotted: it is the first
%                     lya_T_sec of the same input sequence) and once with no
%                     input (lle_noinput).
%   rates             each neuron's time-mean rate over the training and
%                     test windows of the memory task: the mean over E and
%                     over I neurons, and the fractions of neurons pinned
%                     (above 0.9) and silent (below 0.01), so a low memory
%                     capacity can be read as saturated, silent or chaotic.
%
% SIZES, by run_mode: mc_gain_sweep_grid holds the one size table (levels
% per axis and networks, with job counts and times) and the default ranges:
% gains log-spaced 0.1 to 4, biases linear -0.6 to 0.6, seeds 3001, 3002,
% ...; 'essential' is 3 x 3 levels on 1 network, 'fast' 7 x 7 on 3. A run
% on another axis states its own levels (second_levels) and takes gains and
% seeds from the table unless given. Explicit gain_levels, bias_levels /
% second_levels and seeds override the table. The protocol is the same in
% every mode. (Until 2026-11-13 the table was a list of hand-picked levels
% per mode, fast = 11 gains x 7 biases x 2 seeds; a result made with it
% still reads, and resumes when its levels are given explicitly.)
%
% RESUMABLE. Each finished job writes <out_dir>/jobs/job_*.mat; a rerun into
% the same output directory skips the jobs that are there (and reruns the
% ones recorded as failed, unless retry_failed is false). A job file made
% with a different preset, protocol, integrator or model override is
% refused, not mixed in. A condition that throws is recorded as failed with
% its message, and the other conditions and jobs carry on.
%
% OUTPUT, in <output_dir>/mc_gain_sweep/ (default output_dir:
% data/mc_gain_sweep_<timestamp>/, so the default lands in data/ and the
% figures' run_dir is that dated folder):
%   mc_gain_sweep_results.mat   T (table) and settings (struct), below
%   mc_gain_sweep_results.csv   T without the R^2 curve
%   jobs/job_*.mat              one per job
%   git_provenance.txt          (+ working_changes.patch when dirty)
% T has one row per (gain, second-axis value, seed, condition). The second
% axis's column is named 'bias' for S_c and by the property otherwise
% (mc_sweep_second_axis):
%   gain, bias, seed, seed_stim, condition, mc_total, horizon (s),
%   R2_delay (1 x n_delay, R^2 against delay), lle_input, lle_noinput (1/s),
%   rate_E, rate_I, frac_pinned, frac_silent, seconds (wall time of the
%   condition), failed (logical), message (the error text of a failed row).
% settings: conditions, delay_s, preset_name, second_axis (the property),
%   preset_gain, preset_bias, run_mode, gain_levels, bias_levels (for
%   another axis preset_<column> and <column>_levels), seeds, protocol, ode_solver, the
%   Lyapunov settings, model_args, n, n_jobs, n_failed, minutes (wall time
%   of the call; the runs of a set share it), git.
%
% Options (mc_gain_sweep_options names them and holds the defaults):
%   preset_name   a SRNNCellTypePairs preset with noise off (default the
%                 mu7chain bundle's noise-free memory-capacity preset)
%   run_mode      the size (default 'essential'); see SIZES
%   gain_levels, seeds   explicit levels of the gain and network seeds
%   second_axis   the model property of the second axis (default 'S_c')
%   second_levels its levels. Required when second_axis is not the bias.
%   bias_levels   the same option under the name it has for the default
%                 axis; refused with another axis
%   conditions    names of the preset's conditions to run (default all)
%   output_dir    where the run's mc_gain_sweep folder goes (default a
%                 dated folder under data/)
%   ode_solver    default 'sra1', as the bundle's memory capacity uses
%   verbose       'minimal' (default) | 'verbose' | 'near-none'
%   model_args    name-value overrides of the preset (cell), applied after
%                 it and before the job's gain and second-axis value
%   protocol      a struct overriding fields of mc_protocol (T_wash_sec,
%                 T_train_sec, ...). FOR TESTS: a shorter training window
%                 makes the readout under-determined at the paper's n.
%   lya_K, lya_T_sec, lya_acc_sec   the Lyapunov run: exponents kept, its
%                 length, and the final stretch accumulated (5, 30 s, 20 s)
%   segment_len   samples per integration segment of the lean path
%   retry_failed  rerun the jobs an earlier run recorded as failed (true)
%   use_parallel  false runs the jobs serially in the client (true)
%
% The pool is the session's allocation (pool_set); with none the jobs run
% serially.
%
% See also: run_mc_gain_sweep_set, mc_gain_sweep_options, mc_gain_sweep_grid,
%           SRNN_ESN_reservoir.run_memory_capacity_lean, mc_protocol,
%           mc_sweep_second_axis, mc_gain_sweep_job_rows, read_mc_gain_sweep,
%           read_mc_gain_sweep_jobs, run_memory_capacity, verify_shared_build

[out_dirs, results_files] = run_mc_gain_sweep_set({varargin});
out_dir = out_dirs{1};
results_file = results_files{1};
end
