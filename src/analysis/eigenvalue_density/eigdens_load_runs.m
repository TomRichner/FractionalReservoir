function out = eigdens_load_runs(run_dir, index, k)
%EIGDENS_LOAD_RUNS The ONE place the eigenvalue_density stage reads the outlier_spectra stage's output.
%
%   index = eigdens_load_runs(run_dir)            the runs the index lists
%   X     = eigdens_load_runs(run_dir, index, k)  run k of that index, as this stage needs it
%
% WHY ONE FUNCTION. The input is what the opt-in stage outlier_spectra
% saves (FD-040; docs/Plan_spectral_supplement.md, section 5), which was
% being built in parallel with this stage, so its field names may move at
% the merge. Every read of that data is here and nowhere else: a renamed
% field is a change to this file only. Everything downstream
% (eigdens_run_measures, run_eigenvalue_density_stage, the figures) sees
% the normalized struct X below and never the stage's own layout.
%
% ONE RUN DIRECTORY, OR AN ERROR. The input is <run_dir>/outlier_spectra/;
% there is no fallback to the exploration's data (data/explorations/...)
% or to any other folder. A missing folder, index or run file is an error
% (eigdens_load_runs:NoInput, :NoRunFile), and so is a run whose layout
% this function does not know (eigdens_load_runs:Layout).
%
% THE INDEX is read with ospec_read_index(run_dir) when that accessor is on
% the path (FD-040 provides it), else from
% <run_dir>/outlier_spectra/ospec_index.mat, variable runs: a table with
% network, condition, seed1, seed2, file (relative to the folder or
% absolute) and, unused here, has_full_series and window_s. A run is read
% with ospec_load_run(run_dir, index(k, :)) when that is on the path, else
% loaded from its file (variable R).
%
% WHAT IS READ FROM A RUN R (the contract of section 5):
%   R.network, R.condition          strings; also taken from the index
%   R.rng_seeds                     [seed1 seed2] (else the index's)
%   R.preset_name, R.net_args       to rebuild the network's W (below)
%   R.tau_d, R.n, R.N               (n and N fall back to the spectra's sizes)
%   R.cal.t                         1 x K calibration times (else R.cal_t)
%   R.cal.gain                      n x K PER-NEURON gain d_j (slope of the
%                                   nonlinearity times the available vesicle
%                                   fraction): the column scaling of M = W D
%   R.cal.cell_type                 n x 1 (or n x K) type index of each neuron
%   the spectra at the K states, from the outlier_null fields:
%     ev_M         n x 1 per state, the whole spectrum of M (units of W)
%     ev_null      n x q per state, the null draws (same gains on a shuffled,
%                  mean-free W)
%     ev_uniform   n x 1 per state, the actual W under one uniform gain
%     R_eff        the predicted bulk radius at the state (kept for a check)
%   R.dense(k).lam, .xfrac          the whole full-Jacobian spectrum and the
%                                   fraction of each eigenvector's norm^2 in
%                                   the dendritic rows (outlier_dense_spectrum)
%   R.lya.LLE                       the leading Lyapunov exponent, if any
% Two layouts of the per-state fields are accepted, because the contract
% names both R.cal.gain as an n x K matrix and "everything outlier_run
% saves", where R.cal is a 1 x K struct array of outlier_null results:
%   (a) R.cal a SCALAR struct whose fields hold the states as columns
%       (gain n x K, ev_M n x K, ev_null n x q x K, ev_uniform n x K,
%       R_eff 1 x K); or
%   (b) R.cal a 1 x K STRUCT ARRAY, one element per state (gain n x 1,
%       ev_M n x 1, ev_null n x q, ev_uniform n x 1, R_eff scalar).
% In either, the spectra may instead sit in R.null (same two layouts) when
% R.cal holds only the states and gains.
%
% THE WEIGHTS. The block theory needs the realized (post, pre) block
% variance s2 of W, which the contract does not list. W is redrawn from
% the run's own description, without integrating: ospec_rebuild_model(R)
% when that is on the path (FD-040), else build_from_preset(R.preset_name,
% R.condition, R.net_args{:}, 'rng_seeds', R.rng_seeds). The redrawn
% network is checked against R.structural.R_realized when the run carries
% it, and its cell types against R.cal.cell_type, so a wrong rebuild is an
% error (eigdens_load_runs:WrongNetwork), not a silently different s2.
%
% THE NORMALIZED RUN X (doubles; K states):
%   network, condition, seed (1 x 2), file
%   n, N, tau_d, C, type_of (n x 1), n_per_type (1 x C)
%   s2, mu        C x C (post, pre) realized variance and mean of W's
%                 blocks, zeros included (outlier_weight_stats's definitions)
%   R_struct      the structural bulk radius of W (all gains 1)
%   harris        the inputs of eigdens_harris_density for W alone: alpha,
%                 f_e, mu_e, mu_i, sigma_e, sigma_i (the class's normalized
%                 values times sqrt(n) times level_of_chaos), R_class
%   t             1 x K calibration times
%   gain          n x K per-neuron gains (|d| below 1e-10 set to exactly 0)
%   ev_M          n x K;  ev_null  n x q x K;  ev_uniform  n x K
%   R_eff_saved   1 x K (NaN when absent)
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

folder = fullfile(run_dir, 'outlier_spectra');
if nargin < 2
    out = read_index(run_dir, folder);
    return;
end
row = index(k, :);
R = read_run(run_dir, folder, row);
out = normalize(R, row, folder);
end

%% ------------------------------------------------------------------------
function T = read_index(run_dir, folder)
if isempty(run_dir) || ~isfolder(run_dir)
    error('eigdens_load_runs:NoInput', 'The run directory ''%s'' does not exist.', run_dir);
end
if exist('ospec_read_index', 'file') == 2
    T = ospec_read_index(run_dir);
else
    f = fullfile(folder, 'ospec_index.mat');
    if ~isfile(f)
        error('eigdens_load_runs:NoInput', ['No %s: the eigenvalue_density stage reads the output of the ' ...
            'outlier_spectra stage of the SAME run directory. Run that stage first.'], f);
    end
    L = load(f, 'runs');
    if ~isfield(L, 'runs')
        error('eigdens_load_runs:Layout', '%s holds no variable ''runs''.', f);
    end
    T = L.runs;
end
if ~istable(T)
    error('eigdens_load_runs:Layout', 'The outlier_spectra index is a %s, not a table.', class(T));
end
need = {'network', 'condition', 'seed1', 'seed2', 'file'};
missing = setdiff(need, T.Properties.VariableNames);
if ~isempty(missing)
    error('eigdens_load_runs:Layout', 'The outlier_spectra index lacks the column(s) %s.', strjoin(missing, ', '));
end
if height(T) == 0
    error('eigdens_load_runs:NoInput', 'The outlier_spectra index of %s lists no runs.', run_dir);
end
T.network = cellstr(string(T.network));
T.condition = cellstr(string(T.condition));
T.file = cellstr(string(T.file));
end

function R = read_run(run_dir, folder, row)
if exist('ospec_load_run', 'file') == 2
    R = ospec_load_run(run_dir, row);
    return;
end
f = row.file{1};
if ~isfile(f); f = fullfile(folder, f); end
if ~isfile(f)
    error('eigdens_load_runs:NoRunFile', 'The index names %s, which does not exist.', f);
end
L = load(f, 'R');
if ~isfield(L, 'R')
    error('eigdens_load_runs:Layout', '%s holds no variable ''R''.', f);
end
R = L.R;
end

function X = normalize(R, row, folder)
X = struct();
X.network = row.network{1};
X.condition = row.condition{1};
X.seed = [row.seed1, row.seed2];
if isfield(R, 'rng_seeds') && ~isempty(R.rng_seeds); X.seed = double(R.rng_seeds(1:2)); end
X.file = fullfile(folder, row.file{1});
for f = {'network', 'condition'}
    if isfield(R, f{1}) && ~strcmp(char(R.(f{1})), X.(f{1}))
        error('eigdens_load_runs:Layout', 'The index says %s ''%s'' and the run says ''%s'' (%s).', ...
            f{1}, X.(f{1}), char(R.(f{1})), X.file);
    end
end
if ~isfield(R, 'cal') || isempty(R.cal)
    error('eigdens_load_runs:Layout', 'Run %s has no calibration states (R.cal).', X.file);
end

%% The states: gains, times, types, spectra
K = n_states(R.cal);
gain = per_state(R.cal, 'gain', K, X.file);
X.gain = double(gain);
X.gain(abs(X.gain) < 1e-10) = 0;
n = size(X.gain, 1);
if isfield(R.cal, 't')
    X.t = reshape(double(per_state(R.cal, 't', K, X.file)), 1, []);
elseif isfield(R, 'cal_t')
    X.t = double(R.cal_t(:)');
else
    X.t = nan(1, K);
end
ct = double(per_state(R.cal, 'cell_type', K, X.file));
type_of = ct(:, 1);
spec = R.cal;
if ~isfield(spec, 'ev_M')
    if isfield(R, 'null') && isfield(R.null, 'ev_M')
        spec = R.null;
    else
        error('eigdens_load_runs:Layout', 'Run %s holds no frozen-gain spectra (ev_M in R.cal or R.null).', X.file);
    end
end
X.ev_M = double(per_state(spec, 'ev_M', K, X.file));
X.ev_uniform = double(per_state(spec, 'ev_uniform', K, X.file));
if numel(spec) > 1
    X.ev_null = double(cat(3, spec.ev_null));
else
    X.ev_null = double(spec.ev_null);
    if K == 1; X.ev_null = reshape(X.ev_null, size(X.ev_null, 1), [], 1); end
end
X.R_eff_saved = nan(1, K);
if isfield(spec, 'R_eff')
    X.R_eff_saved = reshape(double(per_state(spec, 'R_eff', K, X.file)), 1, []);
end
if size(X.ev_M, 1) ~= n || size(X.ev_M, 2) ~= K || size(X.ev_null, 3) ~= K || numel(type_of) ~= n
    error('eigdens_load_runs:Layout', 'Run %s: the gains, types and spectra do not have n = %d rows and K = %d states.', ...
        X.file, n, K);
end

%% The full Jacobian at the same states
X.dense = [];
if isfield(R, 'dense') && ~isempty(R.dense)
    if numel(R.dense) ~= K
        error('eigdens_load_runs:Layout', 'Run %s has %d dense spectra for %d states.', X.file, numel(R.dense), K);
    end
    X.dense = struct('lam', arrayfun(@(d) double(d.lam(:)), R.dense, 'UniformOutput', false), ...
        'xfrac', arrayfun(@(d) double(d.xfrac(:)), R.dense, 'UniformOutput', false));
end
X.lambda_1 = NaN;
if isfield(R, 'lya') && isfield(R.lya, 'LLE') && ~isempty(R.lya.LLE); X.lambda_1 = R.lya.LLE; end

%% The network's weights, redrawn
model = rebuild_model(R, X);
p = model.get_params();
W = full(p.W);
C = p.n_cellTypes; ti = p.type_indices;
if p.n ~= n
    error('eigdens_load_runs:WrongNetwork', 'Run %s: the rebuilt network has n = %d, the saved gains %d.', X.file, p.n, n);
end
type_model = zeros(n, 1);
mu = zeros(C); s2 = zeros(C);
for a = 1:C
    type_model(ti{a}) = a;
    for b = 1:C
        blk = W(ti{a}, ti{b});
        mu(a, b) = mean(blk(:));
        s2(a, b) = var(blk(:), 1);
    end
end
if ~isequal(type_model, type_of)
    error('eigdens_load_runs:WrongNetwork', 'Run %s: the rebuilt network''s cell types differ from R.cal.cell_type.', X.file);
end
n_per = cellfun(@numel, ti);
R_struct = sqrt(max(0, max(real(eig(s2 .* n_per(:)')))));
if isfield(R, 'structural') && isfield(R.structural, 'R_realized') && ...
        abs(R_struct - R.structural.R_realized) > 1e-6 * R_struct
    error('eigdens_load_runs:WrongNetwork', 'Run %s: the rebuilt W has realized radius %.6g, the run had %.6g.', ...
        X.file, R_struct, R.structural.R_realized);
end
X.n = n; X.C = C; X.type_of = type_of; X.n_per_type = n_per(:)';
X.N = NaN;
if isfield(R, 'N'); X.N = R.N; elseif ~isempty(X.dense); X.N = numel(X.dense(1).lam); end
X.tau_d = p.tau_d;
if isfield(R, 'tau_d') && abs(R.tau_d - p.tau_d) > 1e-12
    error('eigdens_load_runs:WrongNetwork', 'Run %s: tau_d %g saved, %g rebuilt.', X.file, R.tau_d, p.tau_d);
end
X.s2 = s2; X.mu = mu; X.R_struct = R_struct;
sc = sqrt(p.n) * p.level_of_chaos;
X.harris = struct('alpha', p.alpha, 'f_e', p.n_per_type(1) / p.n, ...
    'mu_e', sc * p.mu_tilde(1, 1), 'mu_i', sc * p.mu_tilde(1, 2), ...
    'sigma_e', sc * p.sigma_tilde(1, 1), 'sigma_i', sc * p.sigma_tilde(1, 2), 'R_class', model.R);
end

function K = n_states(cal)
if numel(cal) > 1
    K = numel(cal);
else
    K = size(cal.gain, 2);
end
end

function v = per_state(s, name, K, file)
% The field NAME at the K states, as columns, from either layout.
if ~isfield(s, name)
    error('eigdens_load_runs:Layout', 'Run %s: no field ''%s'' at the calibration states.', file, name);
end
if numel(s) > 1
    v = cat(2, s.(name));
else
    v = s.(name);
    if isvector(v) && numel(v) == K && K > 1 && ~any(strcmp(name, {'gain', 'cell_type', 'ev_M', 'ev_uniform'}))
        v = v(:)';
    end
end
end

function model = rebuild_model(R, X)
if exist('ospec_rebuild_model', 'file') == 2
    model = ospec_rebuild_model(R);
    return;
end
need = {'preset_name', 'net_args'};
missing = need(~isfield(R, need));
if ~isempty(missing)
    error('eigdens_load_runs:Layout', 'Run %s lacks %s, needed to redraw its weights.', X.file, strjoin(missing, ', '));
end
args = R.net_args;
if isstruct(args); args = struct2namevalue(args); end
model = build_from_preset(char(R.preset_name), X.condition, args{:}, 'rng_seeds', X.seed, ...
    'T_range', [0 1], 'verbose', 'near-none');
end
