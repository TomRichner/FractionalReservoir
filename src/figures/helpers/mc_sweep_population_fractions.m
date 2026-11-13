function fr = mc_sweep_population_fractions(D)
% MC_SWEEP_POPULATION_FRACTIONS The E and I fractions of the network a gain-by-bias result was run on.
%
%   fr = MC_SWEEP_POPULATION_FRACTIONS(D)   % D from read_mc_gain_sweep; fr = [f_E, f_I]
%
% Read, in order of precedence, from
%   1. D.settings.model_args, the overrides the stage applied to the preset
%      ('f' as a 1 x 2 row, or the alias 'f_E'), as run_mc_gain_sweep does;
%   2. the preset the run names (D.preset_name): its f and cell_type_names.
% The network must have two cell types named E and I (in either order), the
% ones rate_E and rate_I describe. With neither source the call errors,
% because a guessed fraction would silently weight the mean rate wrongly.
%
% See also: mc_sweep_quantity, read_mc_gain_sweep, srnn_param_preset

names = {'E', 'I'};
f = [];
if ~isempty(D.preset_name)
    d = srnn_param_preset(D.preset_name);
    if ~(isfield(d, 'f') && isfield(d, 'cell_type_names'))
        error('mc_sweep_population_fractions:NoFractions', ...
            'Preset %s states no f and cell_type_names.', D.preset_name);
    end
    f = d.f;
    names = d.cell_type_names;
end
args = {};
if isfield(D.settings, 'model_args'); args = D.settings.model_args; end
for k = 1:2:numel(args) - 1
    switch args{k}
        case 'f',   f = args{k + 1};
        case 'f_E', f = [args{k + 1}, 1 - args{k + 1}];
    end
end
if isempty(f)
    error('mc_sweep_population_fractions:NoFractions', ...
        ['The result names no preset and no f override, so the E and I fractions ' ...
         'are unknown (source: %s).'], D.source);
end
[is_ei, where] = ismember({'E', 'I'}, names);
if numel(f) ~= 2 || numel(names) ~= 2 || ~all(is_ei)
    error('mc_sweep_population_fractions:NotEI', ...
        'Expected two cell types named E and I; the network has %s.', strjoin(names, ', '));
end
fr = f(where) / sum(f);
end
