function names = bundle_run_conditions(cfg)
% BUNDLE_RUN_CONDITIONS The conditions a bundle simulates, by name: cfg.run_conditions, or {} for all (FD-030).
%
%   names = BUNDLE_RUN_CONDITIONS(cfg)
%
%   cfg     a bundle config (paper_config shaped).
%   names   1 x N cellstr, cfg.run_conditions as given; {} when the field is
%           absent or empty, which every stage reads as "every condition its
%           preset states", what every bundle ran before the field existed.
%
% cfg.run_conditions selects from the preset's conditions which ones the
% analyses SIMULATE; cfg.conditions selects which ones the figures DRAW, in
% their left-to-right order. The two are separate because a run may hold
% more conditions than one figure set draws (mu7chain at 'essential' runs
% four and could draw three). run_all_paper_analyses passes the list to
% every stage as its 'conditions' argument and checks it against every
% preset a stage runs on before any compute (preset_with_conditions); the
% opt-in stages read it through this function as their default.
%
% Only the type is checked here: a cellstr, a char (one name) or empty.
% Whether each name is a condition of a preset is preset_with_conditions'
% check, which knows the preset.
%
% See also: preset_with_conditions, run_all_paper_analyses, mu7chain_config

names = {};
if ~isstruct(cfg) || ~isfield(cfg, 'run_conditions') || isempty(cfg.run_conditions)
    return
end
names = cfg.run_conditions;
if ischar(names) || (isstring(names) && isscalar(names))
    names = {char(names)};
end
if ~iscellstr(names)
    error('bundle_run_conditions:BadConfig', ...
        'cfg.run_conditions must be a cell of condition names (or empty for all).');
end
names = reshape(names, 1, []);
end
