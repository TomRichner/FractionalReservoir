function out = fresp_load_states(run_dir, row)
%FRESP_LOAD_STATES The ONE place the frequency-response stage reads what the outlier_spectra stage saved.
%
%   T  = fresp_load_states(run_dir)          the index of the runs
%   st = fresp_load_states(run_dir, row)     one run's calibration states
%
% The frequency response is computed at the calibration states of the
% stage outlier_spectra (FD-040), which owns the files under
% <run_dir>/outlier_spectra/ and their layout. EVERY read of that output by
% the stage frequency_response (FD-041) is in this function, and every one
% goes through that stage's own accessors: ospec_read_index (the index),
% ospec_load_run (a run) and ospec_rebuild_model (its model, built and not
% integrated). The header of ospec_load_run is the data contract. It reads
% one run directory or errors: there is no fallback to any other data.
%
% WHAT IT TAKES FROM A RUN R (ospec_load_run's header)
%   R.network, R.condition, R.rng_seeds (1 x 2), R.preset_name, R.net_args
%   R.cal.t   1 x n_cal: the calibration times (s)
%   R.cal.S   n_cal x N: the state vectors there (one state per ROW)
%   R.dense   1 x n_cal struct array, ospec_dense_spectrum at the
%             calibration states, of which .cont: the continuation of the
%             rightmost frozen-gain mode into the full Jacobian (1/s), the
%             leftover mode. Only a run of one of the outlier_spectra
%             stage's full_networks has it (R.has_dense).
%   the model of the run, from ospec_rebuild_model(R), for its parameters
%
% Inputs
%   run_dir   the bundle's run directory
%   row       one row of the index: a one-row table, or a scalar struct
%             with the index's variables (table2struct of a row)
%
% Outputs
%   T    table: the index (ospec_read_index) without the runs it marks as
%        failed (ok false): network, condition (cellstr), seed1, seed2,
%        file, has_full_series, has_dense, window_s, dt_sample, ok,
%        seconds, message
%   st   struct:
%          network, condition, seed (1 x 2), preset_name, net_args
%          t      1 x n_cal
%          S      N x n_cal (one state per COLUMN)
%          cont   1 x n_cal complex, Im >= 0
%          params model.get_params() of the run's network
%
% Errors: fresp_load_states:NoIndex (the stage outlier_spectra has not run
% in this run directory), fresp_load_states:NoDense (the run has no dense
% spectra: its network is not among the outlier_spectra stage's
% full_networks), fresp_load_states:Contract (a run file does not hold what
% the contract says; the message names what is missing).
%
% See also: run_frequency_response_stage, fresp_run_one, ospec_read_index,
%           ospec_load_run, ospec_rebuild_model

arguments
    run_dir (1,:) char
    row = []
end

if isempty(row)
    out = read_index(run_dir);
    return;
end
if istable(row); row = table2struct(row); end
if ~(isstruct(row) && isscalar(row) && all(isfield(row, {'network', 'condition', 'seed1', 'seed2'})))
    error('fresp_load_states:Contract', ...
        'row must be one row of the index (a one-row table or a scalar struct with network, condition, seed1, seed2).');
end

R = ospec_load_run(run_dir, struct('network', char(row.network), 'condition', char(row.condition), ...
    'seed', [row.seed1, row.seed2]));
need(R, {'network', 'condition', 'rng_seeds', 'preset_name', 'net_args', 'cal', 'dense'}, 'R');
if ~(isstruct(R.cal) && isscalar(R.cal))
    error('fresp_load_states:Contract', 'R.cal must be a scalar struct with fields t and S (it is %s of size %s).', ...
        class(R.cal), mat2str(size(R.cal)));
end
need(R.cal, {'t', 'S'}, 'R.cal');
if isempty(R.dense)
    error('fresp_load_states:NoDense', ...
        ['Run %s / %s / seed [%d %d] has no dense full-Jacobian spectra, so no leftover mode: the ' ...
         'frequency_response stage reads the networks the outlier_spectra stage names in its ' ...
         'full_networks (cfg.outlier_spectra_full_networks).'], ...
        char(R.network), char(R.condition), R.rng_seeds(1), R.rng_seeds(2));
end
if ~(isstruct(R.dense) && isfield(R.dense, 'cont'))
    error('fresp_load_states:Contract', ...
        'R.dense must be a struct array over the calibration states with field cont (the leftover mode).');
end

[~, params] = ospec_rebuild_model(R);
N = params.state_layout.N;

t = R.cal.t(:)';
S = R.cal.S;
if size(S, 2) ~= N || size(S, 1) ~= numel(t)
    error('fresp_load_states:Contract', ...
        'R.cal.S is %d x %d; expected n_cal x N = %d x %d (network %s, condition %s).', ...
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
function T = read_index(run_dir)
try
    T = ospec_read_index(run_dir);
catch err
    if ~strcmp(err.identifier, 'ospec_read_index:NoIndex'); rethrow(err); end
    error('fresp_load_states:NoIndex', ...
        ['%s\nThe stage frequency_response reads the calibration states the stage ' ...
         'outlier_spectra saves in the SAME run directory; run that stage first.'], err.message);
end
T = T(T.ok, :);            % a run that failed has no file; the outlier_spectra stage reports it
end

function need(s, names, what)
missing = setdiff(names, fieldnames(s));
if ~isempty(missing)
    error('fresp_load_states:Contract', '%s lacks the field(s) %s (the header of ospec_load_run).', ...
        what, strjoin(missing, ', '));
end
end
