function registry = paper_opt_in_stages()
% PAPER_OPT_IN_STAGES The stages run_all_paper_analyses runs only for a bundle that names them.
%
%   registry = PAPER_OPT_IN_STAGES()
%
% A struct array, one element per opt-in stage, in run order:
%   name   the stage's name, as a bundle writes it in cfg.stages
%   plan   @(cfg) -> plan. Validates the bundle's settings for this stage
%          and returns whatever its run needs. No compute and no file
%          access: the entry point calls it in its pre-flight, before the
%          run directory exists, so a bad setting fails in milliseconds.
%   run    @(cfg, run_dir, plan) -> char. Runs the stage into run_dir and
%          returns a one-line description of what it wrote. It throws on
%          failure; the entry point wraps it like every stage.
%
% WHY OPT-IN. A stage that only some bundles want (block C's sweep is hours
% of compute) must not start running in every bundle that names no stages.
% So "every stage" means paper_stage_names' first list, and the stages here
% run only when cfg.stages names them, after the others.
%
% ADDING ONE is a row here and its two functions. paper_stage_names,
% paper_stage_plan and run_all_paper_analyses read this list and name no
% opt-in stage themselves.
%
% See also: paper_stage_names, paper_stage_plan, run_all_paper_analyses

% lyapunov_spectrum_n80  the top-240 Lyapunov spectrum of the 80-neuron
%                network per condition (FD-032; Tom, 2026-11-19). First,
%                so in a bundle that also runs block C it follows block A's
%                stages rather than hours of gain sweeps.
% outlier_spectra     frozen-gain and full-Jacobian spectra along settled
%                runs of the bundle's network and of copies with a stronger
%                E-to-E mean, per condition and seed (FD-040). The input of
%                the two stages after it, so it runs before them.
% frequency_response  the network frequency response at outlier_spectra's
%                calibration states (FD-041); simulates nothing
% eigenvalue_density  the radial eigenvalue density of the frozen-gain
%                connectivity at those states (FD-042); simulates nothing.
%                The plans of these two refuse a bundle whose cfg.stages
%                names them without outlier_spectra. The three follow the
%                80-neuron spectrum and come before block C's sweep.
% mc_gain_sweep  block C: memory capacity over gain and a second axis
% noise_switch   copies of one network that share the noise, then stop
%                sharing it (run_noise_switch_rollouts; Tom, 2026-11-14)
registry = struct( ...
    'name', {'lyapunov_spectrum_n80', 'outlier_spectra', 'frequency_response', 'eigenvalue_density', ...
             'mc_gain_sweep', 'noise_switch'}, ...
    'plan', {@lyapunov_spectrum_n80_stage_plan, @outlier_spectra_stage_plan, @frequency_response_stage_plan, ...
             @eigenvalue_density_stage_plan, @mc_gain_sweep_stage_plan, @noise_switch_stage_plan}, ...
    'run',  {@run_lyapunov_spectrum_n80_stage, @run_outlier_spectra_stage, @run_frequency_response_stage, ...
             @run_eigenvalue_density_stage, @run_mc_gain_sweep_stage, @run_noise_switch_stage});
end
