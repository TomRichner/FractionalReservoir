function T = regression_time_model(opts)
% REGRESSION_TIME_MODEL Time the model's hot functions, alone or against a second copy of the class.
%
%   T = regression_time_model()
%   T = regression_time_model('reference_class', 'SRNNCellTypePairsBefore')
%
% WHY. A rewrite of the model's arithmetic is gated on values by the
% tolerance tier (regression_tol_check) and on speed by this (FD-021;
% Plan A, step 0d: a slowdown above 10 % goes to Tom before the merge).
% It times, on the multiple-timescale condition sfa3_std2 at the regression
% network's size (n = 100) and at the manuscript's (n = 500):
%
%   dynamics_fast            one right-hand-side evaluation
%   compute_Jacobian_fast    one assembled sparse Jacobian
%   jacobian_times (K = 15)  one matrix-free product, the K the sweeps use
%   build + get_params       constructing and building the model from its
%                            preset, then packing params
%
% at a state from the end of a short noise-free run. It uses only the
% public surface the refactor keeps (the constructor, build, run,
% get_params and the three statics).
%
% COMPARE IN ONE CALL, NOT BETWEEN TWO. These calls take tens of
% microseconds to milliseconds, and on this machine two runs of this
% function on the SAME code a minute apart differed by 20 to 30 % at
% n = 100 (measured 2026-11-11, before this option existed). A before /
% after difference of 10 % cannot be read off two such runs. So the
% comparison is made inside one call: REFERENCE_CLASS names a second
% class with the same static interface (a copy of the model class as it
% stood before the edit, with the class name replaced throughout, in a
% folder on the path that is NOT in the repository), both are built from
% the same preset with the same seeds, and the two are timed ALTERNATELY,
% round by round. The reported ratio is the median over rounds of
% (this class / reference class) within a round, which cancels whatever
% the machine was doing during that round.
%
% One "round" is the mean time of REPS consecutive calls. `median_us` is
% the median over rounds, `spread` the (max - min) / median of the rounds.
%
% Options:
%   reps             calls per round (default 300; a tenth of that for the
%                    assembled Jacobian)
%   rounds           rounds per function (default 15)
%   build_reps       builds per round (default 3)
%   K                columns of Y for jacobian_times (default 15)
%   reference_class  see above; '' (default) times the model class alone
%
% Returns a table: case, n, N, what, median_us, spread, and with a
% reference class also ref_us (its median), ratio (see above) and
% ratio_lo / ratio_hi (the quartiles of the per-round ratios).
%
% See also: regression_tol_check, test_jacobian_times
arguments
    opts.reps (1,1) double {mustBePositive, mustBeInteger} = 300
    opts.rounds (1,1) double {mustBePositive, mustBeInteger} = 15
    opts.build_reps (1,1) double {mustBePositive, mustBeInteger} = 3
    opts.K (1,1) double {mustBePositive, mustBeInteger} = 15
    opts.reference_class (1,:) char = ''
end
cases = struct( ...
    'label', {'n100', 'n500'}, ...
    'preset', {'regression_pairs_noStim_noise0p025_dualStd_3cond', ...
               'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium'}, ...
    'condition', {'sfa3_std2', 'sfa3_std2'});
run_args = {'sigma_u_noise', 0, 'ode_solver', 'rk4', 'rng_seeds', [1 2], 'T_range', [0 1], ...
    'fs', 200, 'lya_method', 'none', 'store_full_state', true, 'verbose', 'near-none'};
classes = {'SRNNCellTypePairs'};
if ~isempty(opts.reference_class)
    classes{2} = opts.reference_class;
end
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));

rows = {};
for c = 1:numel(cases)
    [preset, ~, conditions] = srnn_param_preset(cases(c).preset);
    names = cellfun(@(s) s.name, conditions, 'UniformOutput', false);
    cond = rmfield(conditions{strcmp(names, cases(c).condition)}, 'name');
    args = [struct2namevalue(preset), struct2namevalue(cond), run_args];

    fns = cell(4, numel(classes));
    for j = 1:numel(classes)
        cls = classes{j};
        m = feval(cls, args{:});
        m.build();
        m.run();
        p = m.get_params();
        p.u_interpolant = m.u_interpolant;
        S = m.S_out(end, :)';
        t = m.t_out(end);
        N = numel(S);
        n = m.n;
        rng(7, 'twister');
        Y = randn(N, opts.K);
        rhs = str2func([cls '.dynamics_fast']);
        jac = str2func([cls '.compute_Jacobian_fast']);
        jt  = str2func([cls '.jacobian_times']);
        fns(:, j) = {@() rhs(t, S, p); @() jac(S, p); @() jt(S, Y, p); @() built_params(cls, args)};
    end
    what = {'dynamics_fast', 'compute_Jacobian_fast', sprintf('jacobian_times (K = %d)', opts.K), ...
        'build + get_params'};
    reps = [opts.reps, max(1, round(opts.reps / 10)), opts.reps, opts.build_reps];
    for k = 1:numel(what)
        per_call = time_rounds(fns(k, :), reps(k), opts.rounds);
        med = median(per_call, 1);
        row = {cases(c).label, n, N, what{k}, 1e6 * med(1), (max(per_call(:, 1)) - min(per_call(:, 1))) / med(1)};
        if numel(classes) == 2
            ratio = per_call(:, 1) ./ per_call(:, 2);
            row = [row, {1e6 * med(2), median(ratio), prctile(ratio, 25), prctile(ratio, 75)}]; %#ok<AGROW>
        end
        rows(end + 1, :) = row; %#ok<AGROW>
    end
end
vars = {'case', 'n', 'N', 'what', 'median_us', 'spread'};
if numel(classes) == 2
    vars = [vars, {'ref_us', 'ratio', 'ratio_lo', 'ratio_hi'}];
end
T = cell2table(rows, 'VariableNames', vars);
end

function p = built_params(cls, args)
m = feval(cls, args{:});
m.build();
p = m.get_params();
end

function per_call = time_rounds(fns, reps, rounds)
% per_call(r, j): the mean time of REPS calls of fns{j} in round r. The
% functions are timed one after another within a round, in an order that
% alternates between rounds.
for j = 1:numel(fns)
    for k = 1:min(reps, 20); fns{j}(); end          % warm up: JIT, caches
end
per_call = zeros(rounds, numel(fns));
for r = 1:rounds
    order = 1:numel(fns);
    if mod(r, 2) == 0; order = fliplr(order); end
    for j = order
        fn = fns{j};
        t0 = tic;
        for k = 1:reps; fn(); end
        per_call(r, j) = toc(t0) / reps;
    end
end
end
