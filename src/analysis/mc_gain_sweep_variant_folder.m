function folder = mc_gain_sweep_variant_folder(tag)
% MC_GAIN_SWEEP_VARIANT_FOLDER The run-directory subfolder of a variant network's gain sweep.
%
%   folder = MC_GAIN_SWEEP_VARIANT_FOLDER(tag)     % 'EI3to2' -> 'mc_gain_sweep_EI3to2'
%
% A bundle's mc_gain_sweep stage runs the main network into
% <run_dir>/mc_gain_sweep and each variant network into
% <run_dir>/<folder>/mc_gain_sweep. The stage writes there and the figures
% read there, so the name is formed in this one place. tag is the short
% name the bundle gives the variant; it must be usable as part of a folder
% name (letters, digits and underscores).
%
% See also: mc_gain_sweep_stage_plan, resolve_mc_gain_sweep

arguments
    tag (1,:) char
end
if isempty(regexp(tag, '^\w+$', 'once'))
    error('mc_gain_sweep_variant_folder:BadTag', ...
        'A variant tag is letters, digits and underscores; got ''%s''.', tag);
end
folder = ['mc_gain_sweep_' tag];
end
