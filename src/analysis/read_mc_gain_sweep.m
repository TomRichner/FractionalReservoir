function D = read_mc_gain_sweep(source)
% READ_MC_GAIN_SWEEP The one reader of the gain-by-bias memory-capacity table (block C).
%
%   D = READ_MC_GAIN_SWEEP(mat_file)
%   D = READ_MC_GAIN_SWEEP(T)            % a table already in memory
%   D = READ_MC_GAIN_SWEEP(struct('T', T, 'settings', s))
%
% Every block C figure and mc_gain_bias_statistics read the stage's
% consolidated result (run_mc_gain_sweep, Plan C section 2) through THIS
% function and nowhere else, so a column the stage names differently is a
% one-line fix in the ALIASES table below.
%
% THE SHAPE (Plan C, section 2): one row per (gain, bias, seed, condition).
% Canonical columns of D.T, after renaming:
%   gain         level_of_chaos of the job
%   bias         S_c, the centre of the per-neuron setpoint draw
%   seed         network seed (the pairing key across conditions)
%   condition    condition name (cellstr)
%   mc_total     total memory capacity
%   horizon      memory horizon (s)
%   R2_delay     R^2 against delay, one row per job (n_delay columns)
%   lle_input    lambda_1 with the memory-task input on (1/s), the one plotted
%   lle_noinput  lambda_1 with no input (1/s)
%   rate_E, rate_I            mean rate of E and I
%   frac_pinned, frac_silent  fractions of neurons above 0.9 / below 0.01
%   failed       logical; a failed row has every value set to NaN here
% Required: gain, bias, seed, condition, mc_total. Any other column that is
% absent is filled with NaN (failed with false) and listed in D.missing.
%
% D also carries: conditions (in the settings' order if it names them, else
% the order of first appearance), gains, biases, seeds (sorted unique),
% delay_s (settings.delay_s, else []), preset_gain / preset_bias (settings
% fields of those names, else 1 and the preset's S_c, else 0.2), preset_name,
% settings, source.
%
% THE FILE: a .mat holding one table (variable T, results, job_table or
% table, else the only table in the file) and optionally a struct settings.
%
% See also: mc_gain_grid, mc_gain_bias_statistics, fig_mc_gain_bias_maps

% ---- canonical name -> names the stage might use (first is canonical).
% THIS IS THE PLACE TO EDIT WHEN THE STAGE'S COLUMN NAMES ARE KNOWN.
ALIASES = { ...
    'gain',        {'level_of_chaos', 'gain_level'}; ...
    'bias',        {'S_c', 'bias_level'}; ...
    'seed',        {'network_seed', 'rng_seed', 'trial'}; ...
    'condition',   {'cond', 'condition_name'}; ...
    'mc_total',    {'MC', 'MC_total', 'total_mc', 'mc'}; ...
    'horizon',     {'H', 'memory_horizon', 'horizon_s'}; ...
    'R2_delay',    {'R2', 'R2_curve', 'r2_delay'}; ...
    'lle_input',   {'lambda1_input', 'lle', 'lambda1', 'lle_task'}; ...
    'lle_noinput', {'lambda1_noinput', 'lle_no_input', 'lle_free'}; ...
    'rate_E',      {'mean_rate_E', 'r_E'}; ...
    'rate_I',      {'mean_rate_I', 'r_I'}; ...
    'frac_pinned', {'frac_sat', 'frac_saturated', 'frac_above'}; ...
    'frac_silent', {'frac_below', 'frac_quiet'}; ...
    'failed',      {'is_failed', 'fail'}};
REQUIRED = {'gain', 'bias', 'seed', 'condition', 'mc_total'};

settings = struct();
src_txt  = '(table in memory)';
if ischar(source) || isstring(source)
    src_txt = char(source);
    S = load(src_txt);
    [T, settings] = table_and_settings(S, src_txt);
elseif istable(source)
    T = source;
elseif isstruct(source) && isfield(source, 'T')
    T = source.T;
    if isfield(source, 'settings'); settings = source.settings; end
else
    error('read_mc_gain_sweep:BadInput', 'Expected a .mat path, a table or a struct with field T.');
end

%% Rename to canonical
names = T.Properties.VariableNames;
missing = {};
for k = 1:size(ALIASES, 1)
    canon = ALIASES{k, 1};
    if ismember(canon, names); continue; end
    hit = find(ismember(ALIASES{k, 2}, names), 1);
    if ~isempty(hit)
        T = renamevars(T, ALIASES{k, 2}{hit}, canon);
        names = T.Properties.VariableNames;
    elseif ismember(canon, REQUIRED)
        error('read_mc_gain_sweep:MissingColumn', ...
            'The table has no column for ''%s'' (accepted: %s). Columns: %s.', ...
            canon, strjoin([{canon}, ALIASES{k, 2}], ', '), strjoin(names, ', '));
    else
        missing{end+1} = canon; %#ok<AGROW>
        if strcmp(canon, 'failed')
            T.failed = false(height(T), 1);
        elseif strcmp(canon, 'R2_delay')
            T.R2_delay = nan(height(T), 0);
        else
            T.(canon) = nan(height(T), 1);
        end
    end
end

%% Types
T.condition = cellstr(string(T.condition));
T.failed    = logical(T.failed) | ~isfinite(T.mc_total);
num_cols = {'gain', 'bias', 'seed', 'mc_total', 'horizon', 'lle_input', 'lle_noinput', ...
            'rate_E', 'rate_I', 'frac_pinned', 'frac_silent'};
for k = 1:numel(num_cols)
    T.(num_cols{k}) = double(T.(num_cols{k}));
end
if iscell(T.R2_delay)    % a cell of rows -> a matrix (rows padded with NaN)
    n_d = max(cellfun(@numel, T.R2_delay));
    M = nan(height(T), n_d);
    for i = 1:height(T); M(i, 1:numel(T.R2_delay{i})) = T.R2_delay{i}; end
    T.R2_delay = M;
end
% A failed job carries no values.
vals = setdiff(num_cols, {'gain', 'bias', 'seed'});
for k = 1:numel(vals)
    T.(vals{k})(T.failed) = NaN;
end
T.R2_delay(T.failed, :) = NaN;

%% What the run holds
in_table = unique(T.condition, 'stable')';
conditions = in_table;
if isfield(settings, 'conditions') && ~isempty(settings.conditions)
    sc = cellstr(string(settings.conditions));
    conditions = [sc(ismember(sc, in_table)), in_table(~ismember(in_table, sc))];
    conditions = reshape(conditions, 1, []);
end

D = struct();
D.T          = T;
D.conditions = conditions;
D.gains      = unique(T.gain)';
D.biases     = unique(T.bias)';
D.seeds      = unique(T.seed)';
D.delay_s    = field_or(settings, 'delay_s', []);
D.preset_name = char(field_or(settings, 'preset_name', ''));
D.preset_gain = field_or(settings, 'preset_gain', 1);
D.preset_bias = field_or(settings, 'preset_bias', NaN);
if isnan(D.preset_bias)
    D.preset_bias = 0.2;
    if ~isempty(D.preset_name)
        try
            d = srnn_param_preset(D.preset_name);
            D.preset_bias = d.S_c;
        catch
        end
    end
end
D.missing  = missing;
D.settings = settings;
D.source   = src_txt;
end

%% ------------------------------------------------------------------------
function [T, settings] = table_and_settings(S, file)
settings = struct();
if isfield(S, 'settings') && isstruct(S.settings); settings = S.settings; end
for nm = {'T', 'results', 'job_table', 'table'}
    if isfield(S, nm{1}) && istable(S.(nm{1})); T = S.(nm{1}); return; end
end
fn = fieldnames(S);
is_t = cellfun(@(f) istable(S.(f)), fn);
if nnz(is_t) ~= 1
    error('read_mc_gain_sweep:NoTable', ...
        '%s holds %d tables; expected one (or one named T / results / job_table).', file, nnz(is_t));
end
T = S.(fn{is_t});
end

function v = field_or(s, f, default)
if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = default; end
end
