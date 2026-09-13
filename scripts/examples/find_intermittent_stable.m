function T = find_intermittent_stable(run_dir, opts)
% FIND_INTERMITTENT_STABLE Sweep jobs that are asymptotically stable but intermittent.
%
%   T = find_intermittent_stable()                      % data/topk_med
%   T = find_intermittent_stable(run_dir, 'min_frac', 0.25, 'min_excursion_s', 0.1, 'top', 10)
%
% A query over the per-job Lyapunov scalars the sweeps store (no simulation):
% every successful job in the run directory's 1D_sensitivity_*, tau_sensitivity_*
% and param_space_* sweeps with
%     lambda_1 < 0                              (asymptotically stable)
%     frac_local_positive   >= min_frac         (the leading local rate is positive
%                                                at least that share of the time)
%     mean_positive_excursion_s >= min_excursion_s   (excursions long enough to act on)
% sorted by p95_finite_0p2s (how strongly the 0.2-s finite-time exponent
% peaks). Prints the counts per condition and the top jobs per condition with
% their sweep, grid values (effective_param, so vector parameters decode to
% values), lambda_1, the three transient scalars and D_KY, and returns the
% whole filtered set as a table.
%
% WHY: the stimulation question in the handoff note (Next 3) needs a network
% that is stable but has positive-LLE transients -- bursting_pairs is a limit
% cycle, not that. The medium sweep is the place to look first.
%
% See also: ParamSpaceAnalysis2.from_dir, SRNNCellTypePairs.lya_summary,
%           docs/notes/Session_2026-09-12_topK_and_nonnormal_handoff.md

arguments
    run_dir (1,:) char = ''   % '' -> data/topk_med
    opts.min_frac        (1,1) double = 0.25
    opts.min_excursion_s (1,1) double = 0.1
    opts.top             (1,1) double = 10
    opts.verbose         (1,1) logical = true
end

setup_paths();
if isempty(run_dir)
    run_dir = fullfile(fileparts(which('setup_paths')), 'data', 'topk_med');
end
d = [dir(fullfile(run_dir, '1D_sensitivity_*')); dir(fullfile(run_dir, 'tau_sensitivity_*')); ...
     dir(fullfile(run_dir, 'param_space_*'))];
d = d([d.isdir]);
if isempty(d)
    error('find_intermittent_stable:NoSweeps', 'No sweep directories under %s', run_dir);
end

rows = {};
n_seen = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(d)
    sweep_dir = fullfile(d(k).folder, d(k).name);
    try
        psa = ParamSpaceAnalysis2.from_dir(sweep_dir);
    catch ME
        warning('find_intermittent_stable:LoadFailed', '%s: %s', d(k).name, ME.message);
        continue;
    end
    axes_ = setdiff(psa.grid_params, {'reps'}, 'stable');
    sweep = regexprep(d(k).name, '_nLevs.*$', '');
    for c = fieldnames(psa.results)'
        cname = c{1};
        res = psa.results.(cname);
        res = res(~cellfun(@isempty, res));
        for j = 1:numel(res)
            r = res{j};
            if ~isfield(r, 'success') || ~r.success || ~isfield(r, 'frac_local_positive'); continue; end
            if n_seen.isKey(cname); n_seen(cname) = n_seen(cname) + 1; else; n_seen(cname) = 1; end
            if ~(r.LLE < 0 && r.frac_local_positive >= opts.min_frac && ...
                 r.mean_positive_excursion_s >= opts.min_excursion_s)
                continue;
            end
            vals = cell(1, numel(axes_));
            for a = 1:numel(axes_)
                v = psa.effective_param(r, axes_{a});
                if isnumeric(v) && isscalar(v); vals{a} = sprintf('%s=%.3g', axes_{a}, v);
                else; vals{a} = sprintf('%s=%s', axes_{a}, mat2str(v, 3)); end
            end
            rep = NaN; if isfield(r.config, 'reps'); rep = r.config.reps; end
            rows(end + 1, :) = {cname, sweep, strjoin(vals, ', '), rep, r.LLE, r.frac_local_positive, ...
                r.p95_finite_0p2s, r.mean_positive_excursion_s, r.D_KY, r.K_used}; %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', {'condition', 'sweep', 'grid', 'rep', 'lambda_1', ...
    'frac_local_positive', 'p95_finite_0p2s', 'mean_excursion_s', 'D_KY', 'K_used'});
T = sortrows(T, 'p95_finite_0p2s', 'descend');

if opts.verbose
    fprintf('Intermittent-but-stable jobs in %s\n', run_dir);
    fprintf('  criteria: lambda_1 < 0, frac_local_positive >= %.2f, mean excursion >= %.2f s\n', ...
        opts.min_frac, opts.min_excursion_s);
    for c = n_seen.keys
        n_hit = nnz(strcmp(T.condition, c{1}));
        fprintf('  %-14s %4d of %4d successful jobs\n', c{1}, n_hit, n_seen(c{1}));
    end
    for c = n_seen.keys
        Tc = T(strcmp(T.condition, c{1}), :);
        if isempty(Tc); continue; end
        fprintf('\n  -- %s: top %d by p95 of the 0.2-s finite-time exponent --\n', c{1}, min(opts.top, height(Tc)));
        fprintf('  %-38s %-36s rep  lambda_1  frac+  p95_0.2s  exc(s)  D_KY  K\n', 'sweep', 'grid');
        for i = 1:min(opts.top, height(Tc))
            fprintf('  %-38s %-36s %3d  %+7.3f  %.2f  %7.2f  %6.2f  %4.1f  %d\n', Tc.sweep{i}, Tc.grid{i}, ...
                Tc.rep(i), Tc.lambda_1(i), Tc.frac_local_positive(i), Tc.p95_finite_0p2s(i), ...
                Tc.mean_excursion_s(i), Tc.D_KY(i), Tc.K_used(i));
        end
    end
end
end
