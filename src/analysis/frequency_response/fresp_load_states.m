function out = fresp_load_states(run_dir, row)
%FRESP_LOAD_STATES The ONE place the frequency-response stage reads what the outlier_spectra stage saved.
%
%   T  = fresp_load_states(run_dir)          the index of the runs
%   st = fresp_load_states(run_dir, row)     one run's calibration states
%
% The frequency response is computed at the calibration states of the
% stage outlier_spectra (FD-040), which owns the files under
% <run_dir>/outlier_spectra/ and their layout (docs/Plan_spectral_supplement.md,
% section 5). EVERY read of that layout by the stage frequency_response
% (FD-041) is in this function, so a renamed field is a change here and
% nowhere else. It reads one run directory or errors: there is no fallback
% to any other data.
%
% WHAT IT READS
%   the index   ospec_index.mat, table `runs` with network, condition,
%               seed1, seed2, file
%   per run     struct R with
%                 R.network, R.condition, R.rng_seeds (1 x 2)
%                 R.preset_name, R.net_args   what builds the network
%                 R.cal.t   1 x n_states: the calibration times (s)
%                 R.cal.S   n_states x N: the state vectors there
%                 R.dense   struct array over the calibration states, with
%                           .cont: the leftover mode, the eigenvalue of the
%                           full Jacobian the frozen-gain outlier continues
%                           into (1/s)
%   the model   of a run, built and not integrated, for its parameters
%
% THROUGH FD-040'S OWN FUNCTIONS WHEN THEY ARE ON THE PATH: ospec_read_index(run_dir),
% ospec_load_run(run_dir, row) and ospec_rebuild_model(R). Until FD-040 is
% merged they are not, and the same files are read directly as the contract
% lays them out (the index and the run file by load; the model by
% build_from_preset(R.preset_name, R.condition, R.net_args{:}, 'rng_seeds',
% R.rng_seeds)). Both routes read the same files of the same run directory.
% AT THE MERGE: check the three call signatures below against FD-040's, and
% the direct route can go.
%
% Inputs
%   run_dir   the bundle's run directory
%   row       one row of the index: a one-row table or a scalar struct with
%             the index's variables (table2struct of a row)
%
% Outputs
%   T    table: network, condition (cellstr), seed1, seed2, file, and
%        whatever else the index holds; rows the index marks as failed
%        (a variable `ok` that is false) are left out
%   st   struct:
%          network, condition, seed (1 x 2), preset_name, net_args
%          t      1 x n_states
%          S      N x n_states (one state per COLUMN)
%          cont   1 x n_states complex, Im >= 0
%          params model.get_params() of the run's network
%
% Errors: fresp_load_states:NoIndex (the stage outlier_spectra has not run
% in this run directory), fresp_load_states:Contract (a file does not hold
% what the contract says; the message names what is missing).
%
% See also: run_frequency_response_stage, fresp_run_one

arguments
    run_dir (1,:) char
    row = []
end

folder = fullfile(run_dir, 'outlier_spectra');
if isempty(row)
    out = read_index(run_dir, folder);
    return;
end
if istable(row); row = table2struct(row); end
if ~(isstruct(row) && isscalar(row))
    error('fresp_load_states:Contract', 'row must be one row of the index (a one-row table or a scalar struct).');
end

R = load_run(run_dir, folder, row);
need(R, {'network', 'condition', 'rng_seeds', 'preset_name', 'net_args', 'cal', 'dense'}, 'R');
if ~(isstruct(R.cal) && isscalar(R.cal))
    error('fresp_load_states:Contract', 'R.cal must be a scalar struct with fields t and S (it is %s of size %s).', ...
        class(R.cal), mat2str(size(R.cal)));
end
need(R.cal, {'t', 'S'}, 'R.cal');
if ~(isstruct(R.dense) && isfield(R.dense, 'cont'))
    error('fresp_load_states:Contract', ...
        'R.dense must be a struct array over the calibration states with field cont (the leftover mode).');
end

model = rebuild_model(R);
params = model.get_params();
N = params.state_layout.N;

t = R.cal.t(:)';
S = R.cal.S;
if size(S, 2) ~= N || size(S, 1) ~= numel(t)
    error('fresp_load_states:Contract', ...
        'R.cal.S is %d x %d; expected n_states x N = %d x %d (network %s, condition %s).', ...
        size(S, 1), size(S, 2), numel(t), N, char(R.network), char(R.condition));
end
cont = [R.dense.cont];
if numel(cont) ~= numel(t)
    error('fresp_load_states:Contract', 'R.dense has %d states and R.cal.t has %d.', numel(cont), numel(t));
end
cont = complex(real(cont(:)'), abs(imag(cont(:)')));

out = struct('network', char(R.network), 'condition', char(R.condition), ...
    'seed', double(R.rng_seeds(:)'), 'preset_name', char(R.preset_name), 'net_args', {R.net_args}, ...
    't', t, 'S', double(S)', 'cont', cont, 'params', params);
end

%% ------------------------------------------------------------------------
function T = read_index(run_dir, folder)
if exist('ospec_read_index', 'file') == 2
    T = ospec_read_index(run_dir);
else
    index_file = fullfile(folder, 'ospec_index.mat');
    if ~isfile(index_file)
        error('fresp_load_states:NoIndex', ...
            ['No %s.\nThe stage frequency_response reads the calibration states the stage ' ...
             'outlier_spectra saves in the SAME run directory; run that stage first.'], index_file);
    end
    L = load(index_file);
    if ~isfield(L, 'runs') || ~istable(L.runs)
        error('fresp_load_states:Contract', '%s must hold a table named runs.', index_file);
    end
    T = L.runs;
end
missing = setdiff({'network', 'condition', 'seed1', 'seed2', 'file'}, T.Properties.VariableNames);
if ~isempty(missing)
    error('fresp_load_states:Contract', 'The index of the outlier_spectra stage lacks the variable(s) %s.', ...
        strjoin(missing, ', '));
end
if ismember('ok', T.Properties.VariableNames)
    T = T(logical(T.ok), :);            % a run that failed has no file; the outlier_spectra stage reports it
end
T.network = cellstr(T.network);
T.condition = cellstr(T.condition);
T.file = cellstr(T.file);
end

function R = load_run(run_dir, folder, row)
if exist('ospec_load_run', 'file') == 2
    R = ospec_load_run(run_dir, row);
    return;
end
name = char(row.file);
file = fullfile(folder, name);
if ~isfile(file) && isfile(name); file = name; end        % an index that stores full paths
if ~isfile(file)
    error('fresp_load_states:Contract', 'The index names %s, which is not in %s.', name, folder);
end
L = load(file, 'R');
if ~isfield(L, 'R')
    error('fresp_load_states:Contract', '%s must hold a struct named R.', file);
end
R = L.R;
end

function model = rebuild_model(R)
if exist('ospec_rebuild_model', 'file') == 2
    model = ospec_rebuild_model(R);
    return;
end
% A network seed must mean one network whatever generator the session (or
% a pool worker, which starts on threefry) is on.
prev = rng;
restore = onCleanup(@() rng(prev));
rng(0, 'twister');
model = build_from_preset(char(R.preset_name), char(R.condition), R.net_args{:}, ...
    'rng_seeds', double(R.rng_seeds(:)'), 'verbose', 'near-none', 'lya_method', 'none');
end

function need(s, names, what)
missing = setdiff(names, fieldnames(s));
if ~isempty(missing)
    error('fresp_load_states:Contract', '%s lacks the field(s) %s (docs/Plan_spectral_supplement.md, section 5).', ...
        what, strjoin(missing, ', '));
end
end
