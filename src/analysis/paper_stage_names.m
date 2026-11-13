function [stages, sweeps, opt_in] = paper_stage_names()
% PAPER_STAGE_NAMES The stages of run_all_paper_analyses and the sweeps of its first stage, in run order.
%
%   [stages, sweeps] = PAPER_STAGE_NAMES()
%   [stages, sweeps, opt_in] = PAPER_STAGE_NAMES()
%
% THE ONE LIST of each (FD-023). run_all_paper_analyses runs its stages in
% this order and prints these names in its summary; run_all_analyses runs the
% sweeps in this order; paper_stage_plan validates a bundle's cfg.stages and
% cfg.sweeps against them. A stage added to the entry point is added here, or
% the entry point's own check fails before any compute.
%
%   stages  'sweeps' is run_all_analyses (it also creates the run directory
%           and run_manifest.mat, so those exist whether or not any sweep is
%           selected). 'paper_illustrations' exists only for a bundle with
%           cfg.illustrations true.
%   sweeps  'sensitivity' is the seven 1-D sweeps (run_sensitivity_analysis),
%           'tau' the SFA-timescale sweep (run_tau_sensitivity_analysis),
%           'param_space' the joint sample (run_param_space_analysis).
%   opt_in  OPT-IN STAGES (FD-010): valid in a bundle's cfg.stages, never part
%           of "every stage". A bundle that names no stages, or only stages
%           of the first list, runs exactly what it ran before this list
%           existed, and its plan does not mention them, not even as
%           skipped. A bundle that names one runs it after the stages of the
%           first list, in this list's order.
%           'mc_gain_sweep' is block C's gain-by-bias memory-capacity sweep
%           (run_mc_gain_sweep), with the E:I gain-only runs a bundle names.
%           A stage that only some bundles want goes here, not into the
%           first list, so no existing bundle starts running it.
%
% See also: paper_stage_plan, run_all_paper_analyses, run_all_analyses

stages = {'sweeps', 'memory_capacity', 'mc_example', 'eig_heatmap', 'numerics', ...
          'lyapunov_spectrum', 'transient_gain', 'local_lyapunov', 'dc_lle', ...
          'paper_illustrations'};
sweeps = {'sensitivity', 'tau', 'param_space'};
opt_in = {'mc_gain_sweep'};
end
