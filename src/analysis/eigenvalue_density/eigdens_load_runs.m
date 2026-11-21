function out = eigdens_load_runs(run_dir, index, k)
%EIGDENS_LOAD_RUNS The ONE place the eigenvalue_density stage reads the outlier_spectra stage's output.
%
%   index = eigdens_load_runs(run_dir)            the runs the index lists
%   X     = eigdens_load_runs(run_dir, index, k)  run k of that index, as this stage needs it
%
% WHY ONE FUNCTION. The input is what the opt-in stage outlier_spectra
% saves (FD-040), and every read of it by this stage is here and nowhere
% else, through that stage's own accessors: ospec_read_index (the index),
% ospec_load_run (a run; ITS HEADER IS THE DATA CONTRACT) and
% ospec_rebuild_model (the run's model, weights and block statistics, built
% and not integrated). Everything downstream (eigdens_run_measures,
% run_eigenvalue_density_stage, the figures) sees the normalized struct X
% below and never the stage's own layout.
%
% ONE RUN DIRECTORY, OR AN ERROR. The input is <run_dir>/outlier_spectra/;
% there is no fallback to the exploration's data (data/explorations/...)
% or to any other folder. A missing run directory or index is an error
% (eigdens_load_runs:NoInput), and so is a run that does not hold what the
% contract says (eigdens_load_runs:Layout) or whose rebuilt network is not
% the one it ran (eigdens_load_runs:WrongNetwork).
%
% THE INDEX is ospec_read_index(run_dir) without the runs it marks as
% failed (ok false: they have no file); run k of it is read with
% ospec_load_run(run_dir, index(k, :)).
%
% WHAT IS READ FROM A RUN R (ospec_load_run's header):
%   R.network, R.condition, R.rng_seeds   checked against the index
%   R.n, R.N, R.tau_d                     checked against the rebuilt model
%   R.cal.t          1 x K calibration times
%   R.cal.gain       n x K PER-NEURON gain d_j (slope of the nonlinearity
%                    times the available vesicle fraction, averaged over the
%                    neuron's routes): the column scaling of M = W D. NaN for
%                    a neuron with no synapse at all (its column of W is
%                    zero, so its gain does not enter M): read as 0.
%   R.cal.cell_type  n x 1 type index of each neuron
%   R.cal.ev_M       n x K, the whole spectrum of M at each state (units of W)
%   R.cal.ev_uniform n x K, the actual W under one uniform gain (the rms gain)
%   R.cal.R_eff      1 x K, the predicted bulk radius at each state (kept for a check)
%   R.null(k).ev_null   n x q, the first q = 3 null draws at state k (the
%                    same gains on W with its block means removed and its
%                    entries permuted within each block)
%   R.dense(k).lam, .xfrac   the whole full-Jacobian spectrum and the
%                    fraction of each eigenvector's norm^2 in the dendritic
%                    rows; R.dense is [] for a run of a network that is not
%                    among the outlier_spectra stage's full_networks
%   R.lya.LLE        the leading Lyapunov exponent (an empty struct when
%                    that stage ran with lya_K 0)
%   R.structural.R_realized   the structural bulk radius of the drawn W
%
% THE WEIGHTS. The block theory needs the realized (post, pre) block
% statistics of W. ospec_rebuild_model(R) rebuilds the run's network and
% returns them (ospec_weight_stats: mu and s2 over all entries of a block,
% zeros included; the type of each neuron; R_struct). The rebuild is
% checked against R.structural.R_realized and the cell types against
% R.cal.cell_type, so a wrong rebuild is an error
% (eigdens_load_runs:WrongNetwork), not a silently different s2.
%
% THE NORMALIZED RUN X (doubles; K states):
%   network, condition, seed (1 x 2), file
%   n, N, tau_d, C, type_of (n x 1), n_per_type (1 x C)
%   s2, mu        C x C (post, pre) realized variance and mean of W's
%                 blocks, zeros included (ospec_weight_stats)
%   R_struct      the structural bulk radius of W (all gains 1)
%   harris        the inputs of eigdens_harris_density for W alone: alpha,
%                 f_e, mu_e, mu_i, sigma_e, sigma_i (the class's normalized
%                 values times sqrt(n) times level_of_chaos), R_class
%   t             1 x K calibration times
%   gain          n x K per-neuron gains (|d| below 1e-10, and NaN, set to exactly 0)
%   ev_M          n x K;  ev_null  n x q x K;  ev_uniform  n x K
%   R_eff_saved   1 x K
%   dense         1 x K struct array with lam and xfrac, or [] (no dense pass)
%   lambda_1      the run's leading Lyapunov exponent, NaN when absent
%
% See also: run_eigenvalue_density_stage, eigdens_run_measures,
%           eigenvalue_density_stage_plan

arguments
    run_dir (1,:) char
    index = []
    k (1,1) double = 0
end

if nargin < 2
    out = read_index(run_dir);
    return;
end
row = index(k, :);
R = ospec_load_run(run_dir, row);
out = normalize(R, row, fullfile(ospec_folder(run_dir), row.file{1}));
end

%% ------------------------------------------------------------------------
function T = read_index(run_dir)
if isempty(run_dir) || ~isfolder(run_dir)
    error('eigdens_load_runs:NoInput', 'The run directory ''%s'' does not exist.', run_dir);
end
try
    T = ospec_read_index(run_dir);
catch err
    if ~strcmp(err.identifier, 'ospec_read_index:NoIndex'); rethrow(err); end
    error('eigdens_load_runs:NoInput', ['%s\nThe eigenvalue_density stage reads the output of the ' ...
        'outlier_spectra stage of the SAME run directory. Run that stage first.'], err.message);
end
T = T(T.ok, :);                          % a failed run has no file; that stage reports it
if height(T) == 0
    error('eigdens_load_runs:NoInput', 'The outlier_spectra index of %s lists no finished run.', run_dir);
end
end

function folder = ospec_folder(run_dir)
[~, ~, folder] = ospec_read_index(run_dir);
end

function X = normalize(R, row, file)
X = struct();
X.network = row.network{1};
X.condition = row.condition{1};
X.seed = [row.seed1, row.seed2];
X.file = file;
if ~strcmp(R.network, X.network) || ~strcmp(R.condition, X.condition) || ~isequal(double(R.rng_seeds(:)'), X.seed)
    error('eigdens_load_runs:Layout', 'The index says %s / %s / seed [%d %d] and the run says %s / %s / seed [%d %d] (%s).', ...
        X.network, X.condition, X.seed, R.network, R.condition, R.rng_seeds, X.file);
end
need = {'gain', 'cell_type', 't', 'ev_M', 'ev_uniform', 'R_eff'};
if ~(isfield(R, 'cal') && isstruct(R.cal) && isscalar(R.cal) && all(isfield(R.cal, need)))
    error('eigdens_load_runs:Layout', 'Run %s: R.cal must be a scalar struct with fields %s (ospec_load_run).', ...
        X.file, strjoin(need, ', '));
end

%% The states: gains, times, types, spectra
X.gain = double(R.cal.gain);
X.gain(isnan(X.gain) | abs(X.gain) < 1e-10) = 0;
[n, K] = size(X.gain);
X.t = double(R.cal.t(:)');
type_of = double(R.cal.cell_type(:));
X.ev_M = double(R.cal.ev_M);
X.ev_uniform = double(R.cal.ev_uniform);
X.R_eff_saved = double(R.cal.R_eff(:)');
if ~(isfield(R, 'null') && numel(R.null) == K && isfield(R.null, 'ev_null'))
    error('eigdens_load_runs:Layout', 'Run %s: R.null must hold the null draws (ev_null) at its %d calibration states.', X.file, K);
end
X.ev_null = double(cat(3, R.null.ev_null));
if ~isequal(size(X.ev_M), [n K]) || ~isequal(size(X.ev_uniform), [n K]) || size(X.ev_null, 1) ~= n || ...
        size(X.ev_null, 3) ~= K || numel(type_of) ~= n || numel(X.t) ~= K || numel(X.R_eff_saved) ~= K
    error('eigdens_load_runs:Layout', 'Run %s: the gains, types and spectra do not have n = %d rows and K = %d states.', ...
        X.file, n, K);
end

%% The full Jacobian at the same states
X.dense = [];
if ~isempty(R.dense)
    if numel(R.dense) ~= K
        error('eigdens_load_runs:Layout', 'Run %s has %d dense spectra for %d states.', X.file, numel(R.dense), K);
    end
    X.dense = struct('lam', arrayfun(@(d) double(d.lam(:)), R.dense, 'UniformOutput', false), ...
        'xfrac', arrayfun(@(d) double(d.xfrac(:)), R.dense, 'UniformOutput', false));
end
X.lambda_1 = NaN;
if isfield(R.lya, 'LLE') && ~isempty(R.lya.LLE); X.lambda_1 = R.lya.LLE; end

%% The network's weights: the run's model, rebuilt
[model, p, ws] = ospec_rebuild_model(R);      % checks n and N against the run
if ~isequal(ws.type_of, type_of)
    error('eigdens_load_runs:WrongNetwork', 'Run %s: the rebuilt network''s cell types differ from R.cal.cell_type.', X.file);
end
if abs(ws.R_struct - R.structural.R_realized) > 1e-6 * ws.R_struct
    error('eigdens_load_runs:WrongNetwork', 'Run %s: the rebuilt W has realized radius %.6g, the run had %.6g.', ...
        X.file, ws.R_struct, R.structural.R_realized);
end
if abs(R.tau_d - p.tau_d) > 1e-12
    error('eigdens_load_runs:WrongNetwork', 'Run %s: tau_d %g saved, %g rebuilt.', X.file, R.tau_d, p.tau_d);
end
X.n = n; X.N = R.N; X.C = p.n_cellTypes; X.type_of = type_of;
X.n_per_type = cellfun(@numel, ws.ti); X.n_per_type = X.n_per_type(:)';
X.tau_d = p.tau_d;
X.s2 = ws.s2; X.mu = ws.mu; X.R_struct = ws.R_struct;
sc = sqrt(p.n) * p.level_of_chaos;
X.harris = struct('alpha', p.alpha, 'f_e', p.n_per_type(1) / p.n, ...
    'mu_e', sc * p.mu_tilde(1, 1), 'mu_i', sc * p.mu_tilde(1, 2), ...
    'sigma_e', sc * p.sigma_tilde(1, 1), 'sigma_i', sc * p.sigma_tilde(1, 2), 'R_class', model.R);
end
