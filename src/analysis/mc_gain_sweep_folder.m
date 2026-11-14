function folder = mc_gain_sweep_folder(grid, variant)
% MC_GAIN_SWEEP_FOLDER The run-directory subfolder of one run of a bundle's mc_gain_sweep stage.
%
%   folder = MC_GAIN_SWEEP_FOLDER('bias')             % 'mc_gain_sweep_bias'
%   folder = MC_GAIN_SWEEP_FOLDER('bias', 'EI3to2')   % 'mc_gain_sweep_bias_EI3to2'
%   folder = MC_GAIN_SWEEP_FOLDER('', 'EI3to2')       % 'mc_gain_sweep_EI3to2'
%   folder = MC_GAIN_SWEEP_FOLDER('')                 % ''
%
% A bundle's mc_gain_sweep stage makes several runs (mc_gain_sweep_stage_plan):
% each GRID it lists, named by the grid's tag, and for a grid each VARIANT
% network swept along gain only, named by the variant's tag. Every run
% writes to <run_dir>/<folder>/mc_gain_sweep (run_mc_gain_sweep appends the
% last part). The stage writes there and the figures read there
% (resolve_mc_gain_sweep), so the name is formed in this one place.
%
% THE NAME is 'mc_gain_sweep' followed by the tags that are not empty, joined
% by underscores, so a run directory lists a grid and its variants together.
% A bundle with one unnamed grid (the form before a stage could hold several)
% has the empty grid tag: its main run writes to the run directory itself
% ('' here, so <run_dir>/mc_gain_sweep) and its variants to
% mc_gain_sweep_<variant>, where results made before grids had tags are.
%
% A tag must be usable as part of a folder name: letters, digits and
% underscores. Two runs of one stage must not share a folder, which the
% stage's plan checks (a grid 'a' with variant 'b_c' and a grid 'a_b' with
% variant 'c' would).
%
% See also: mc_gain_sweep_stage_plan, resolve_mc_gain_sweep

arguments
    grid    (1,:) char
    variant (1,:) char = ''
end
tags = {grid, variant};
tags = tags(~cellfun(@isempty, tags));
for k = 1:numel(tags)
    if isempty(regexp(tags{k}, '^\w+$', 'once'))
        error('mc_gain_sweep_folder:BadTag', ...
            'A grid or variant tag is letters, digits and underscores; got ''%s''.', tags{k});
    end
end
folder = '';
if ~isempty(tags)
    folder = strjoin([{'mc_gain_sweep'}, tags], '_');
end
end
