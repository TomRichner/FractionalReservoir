function C = cx_load(run_dir, data_file, who)
%CX_LOAD What the concept_examples stage saved, from a run directory or the file itself.
%
%   C = cx_load(run_dir)
%   C = cx_load('', data_file)
%   C = cx_load(run_dir, '', who)     who names the caller in the error
%
%   C.examples   1 x m struct, one per column of Figure 1A
%   C.gallery    1 x k struct, the structures of S3 (may be empty)
%   C.settings   the stage's plan without its lists
%   C.file       the file read
%
% FIELDS OF EVERY NETWORK (both lists)
%   key, name, label, target, level_of_chaos, seed, n, n_E, n_I, alpha,
%   tau_d
%   ev           n x 1 complex: (eig(W) - 1) / tau_d, the Jacobian's
%                eigenvalues with every gain 1, in 1/s
%   center       -1 / tau_d; R_rate = R_struct / tau_d, the bulk circle;
%                R_theory_rate, the class's ensemble radius over tau_d
%   rightmost    the eigenvalue with the largest real part (1/s)
%   mu, s2       realized block mean and variance, zeros included
%                (post, pre)
%   density, nz_mean, nz_sd   per block: fraction of nonzero entries, and
%                mean and SD of the nonzero entries
%   K_struct, K_eig_rate      the structural mean matrix and its
%                eigenvalues as rates
%   lowrank_rate where the added low-rank structure was built to put its
%                eigenvalue (1/s; NaN when nothing was added)
%   sample       weight_image_sample's struct (W, idx, type, density, ...)
% AND OF EVERY EXAMPLE
%   t (1 x nt), trace_idx, trace_type (1 E, 2 I), x_trace, rate_trace,
%   S_c_trace, S_a, activation, rate_E, rate_I (population means over
%   time), mean_rate_E, mean_rate_I, frac_floor, frac_ceiling, x_sd (over
%   the second half), LLE and lya_interval
%
% Errors through resolve_data_file when the stage has not run.
%
% See also: run_concept_examples_stage, concept_examples_stage_plan

arguments
    run_dir (1,:) char
    data_file (1,:) char = ''
    who (1,:) char = 'cx_load'
end
if isempty(data_file) && isempty(run_dir)
    error([who ':NoData'], 'Give the run directory that holds the concept_examples stage, or its file.');
end
f = resolve_data_file(data_file, run_dir, {fullfile(run_dir, 'concept_examples')}, 'concept_examples.mat', ...
    'the opt-in stage ''concept_examples'' (concept_examples_stage_plan, run_concept_examples_stage)');
L = load(f, 'examples', 'gallery', 'settings');
C = struct('examples', L.examples, 'gallery', L.gallery, 'settings', L.settings, 'file', f);
end
