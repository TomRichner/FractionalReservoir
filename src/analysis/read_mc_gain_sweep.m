function D = read_mc_gain_sweep(source)
% READ_MC_GAIN_SWEEP The one reader of the gain-by-bias memory-capacity table (block C).
%
%   D = READ_MC_GAIN_SWEEP(mat_file)     % mc_gain_sweep_results.mat
%   D = READ_MC_GAIN_SWEEP(T)            % a table already in memory
%   D = READ_MC_GAIN_SWEEP(struct('T', T, 'settings', s))
%   D = READ_MC_GAIN_SWEEP(D)            % its own output, returned as it is
%
% Every block C figure and mc_gain_bias_statistics take the result of
% run_mc_gain_sweep through THIS function, so what a failed row holds, the
% order of the conditions and where the preset sits are decided once. A
% caller that accepts "a file, a table or a result already read" passes
% whatever it was given straight in.
%
% THE TABLE: one row per (gain, bias, seed, condition), with the stage's
% column names:
%   gain         level_of_chaos of the job
%   bias         S_c, the centre of the per-neuron setpoint draw
%   seed         network seed (the pairing key across conditions)
%   condition    condition name (cellstr)
%   mc_total     total memory capacity
%   horizon      memory horizon (s)
%   R2_delay     R^2 against delay (n_delay columns)
%   lle_input    lambda_1 with the memory-task input on (1/s), the one plotted
%   lle_noinput  lambda_1 with no input (1/s)
%   rate_E, rate_I            mean rate of E and I
%   frac_pinned, frac_silent  fractions of neurons pinned / silent
%   failed       logical; a failed row has every value set to NaN here, and
%                a row with no finite mc_total counts as failed
% Required: gain, bias, seed, condition, mc_total. Any other column that is
% absent is filled with NaN (failed with false, R2_delay with no columns)
% and listed in D.missing, so a hand-made table of the five required
% columns can be drawn and summarised.
%
% D: T (as above), conditions (the settings' order when it names them, else
% the order of first appearance), gains, biases, seeds (sorted unique),
% delay_s (settings.delay_s, else []), preset_name, preset_gain and
% preset_bias (the settings' values; a table that carries none is taken to
% sit at gain 1, bias 0.2, the paper presets' point), missing, settings,
% source.
%
% See also: run_mc_gain_sweep, read_mc_gain_sweep_jobs, mc_gain_grid,
%           mc_gain_bias_statistics

REQUIRED = {'gain', 'bias', 'seed', 'condition', 'mc_total'};
VALUES   = {'mc_total', 'horizon', 'lle_input', 'lle_noinput', ...
            'rate_E', 'rate_I', 'frac_pinned', 'frac_silent'};
% Where a table with no settings is taken to sit: the class default gain
% and the S_c of the paper's presets.
DEFAULT_PRESET_GAIN = 1;
DEFAULT_PRESET_BIAS = 0.2;

settings = struct();
source_text = '(table in memory)';
if isstruct(source) && all(isfield(source, {'T', 'gains', 'biases', 'conditions'}))
    D = source;
    return
elseif ischar(source) || isstring(source)
    source_text = char(source);
    loaded = load(source_text);
    if ~(isfield(loaded, 'T') && istable(loaded.T))
        error('read_mc_gain_sweep:NoTable', '%s holds no table T.', source_text);
    end
    T = loaded.T;
    if isfield(loaded, 'settings') && isstruct(loaded.settings); settings = loaded.settings; end
elseif istable(source)
    T = source;
elseif isstruct(source) && isfield(source, 'T')
    T = source.T;
    if isfield(source, 'settings'); settings = source.settings; end
else
    error('read_mc_gain_sweep:BadInput', ...
        'Expected a .mat path, a table, a struct with field T, or this function''s output.');
end

%% Columns: the required ones present, the others filled
names = T.Properties.VariableNames;
absent = setdiff(REQUIRED, names, 'stable');
if ~isempty(absent)
    error('read_mc_gain_sweep:MissingColumn', ...
        'The table has no column %s. Columns: %s.', strjoin(absent, ', '), strjoin(names, ', '));
end
missing = setdiff([VALUES, {'R2_delay', 'failed'}], names, 'stable');
for k = 1:numel(missing)
    switch missing{k}
        case 'failed',   T.failed = false(height(T), 1);
        case 'R2_delay', T.R2_delay = nan(height(T), 0);
        otherwise,       T.(missing{k}) = nan(height(T), 1);
    end
end

%% Types, and a failed row carries no values
T.condition = cellstr(string(T.condition));
T.failed    = logical(T.failed) | ~isfinite(T.mc_total);
for name = [{'gain', 'bias', 'seed'}, VALUES]
    T.(name{1}) = double(T.(name{1}));
end
for name = VALUES
    T.(name{1})(T.failed) = NaN;
end
T.R2_delay(T.failed, :) = NaN;

%% What the run holds
in_table = unique(T.condition, 'stable')';
conditions = in_table;
if isfield(settings, 'conditions') && ~isempty(settings.conditions)
    named = reshape(cellstr(string(settings.conditions)), 1, []);
    conditions = [named(ismember(named, in_table)), in_table(~ismember(in_table, named))];
end

D = struct();
D.T           = T;
D.conditions  = conditions;
D.gains       = unique(T.gain)';
D.biases      = unique(T.bias)';
D.seeds       = unique(T.seed)';
D.delay_s     = field_or(settings, 'delay_s', []);
D.preset_name = char(field_or(settings, 'preset_name', ''));
D.preset_gain = field_or(settings, 'preset_gain', DEFAULT_PRESET_GAIN);
D.preset_bias = field_or(settings, 'preset_bias', DEFAULT_PRESET_BIAS);
D.missing     = missing;
D.settings    = settings;
D.source      = source_text;
end

%% ------------------------------------------------------------------------
function v = field_or(s, f, default)
if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = default; end
end
