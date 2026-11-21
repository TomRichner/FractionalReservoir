function model = ospec_build_model(B)
%OSPEC_BUILD_MODEL The built model of one run of the outlier_spectra stage, from the build record a run saves (no integration).
%
%   model = ospec_build_model(B)
%
% B is the struct ospec_run writes as R.build:
%   preset_name      the bundle's noisy SRNNCellTypePairs preset
%   condition        one of the preset's conditions
%   net_args         name-value overrides that make the network (R.net_args;
%                    from ospec_network: the block means, the indegree)
%   size_args        name-value overrides of the size ({} at the preset's
%                    size; n, indegree and the weight unit's reference for
%                    a smaller network)
%   added_structure  [] or the low-rank addition (OspecNet)
%   rng_seeds        1 x 2
%   T_range, fs      the run's time grid
%   lya_args         name-value Lyapunov settings
%   ode_solver       optional: the integrator (the bundle's cfg.ode_solver,
%                    through the stage's plan). Absent or empty: the rule
%                    below. A file saved before the field existed has none
%                    and rebuilds as it was built.
%
% build_from_preset's rule (preset, then condition, then the caller's
% pairs; sra1 for a noisy preset, rk4 otherwise), on OspecNet so that
% structure can be added to the drawn W and the run can be sampled.
%
% THE GLOBAL STREAM IS SET TO THE TWISTER FIRST (rng(0, 'twister')). The
% class seeds whatever generator is current, and a pool worker's is
% threefry where the client's is the twister, so without this a seed would
% name one network in the client and another on a worker
% (run_lyapunov_spectrum does the same).
%
% See also: ospec_run, ospec_rebuild_model, OspecNet, build_from_preset

arguments
    B (1,1) struct
end

[p, model_class, conditions] = srnn_param_preset(B.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('ospec_build_model:Class', 'Preset %s builds a %s, not an SRNNCellTypePairs.', B.preset_name, model_class);
end
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
ci = find(strcmp(names, B.condition), 1);
if isempty(ci)
    error('ospec_build_model:NoSuchCondition', 'Preset ''%s'' has no condition ''%s''. Available: %s.', ...
        B.preset_name, B.condition, strjoin(names, ', '));
end
cond = rmfield(conditions{ci}, 'name');
solver = 'rk4';
if isfield(p, 'sigma_u_noise') && any(p.sigma_u_noise(:) > 0); solver = 'sra1'; end
if isfield(B, 'ode_solver') && ~isempty(B.ode_solver); solver = B.ode_solver; end

rng(0, 'twister');
args = [struct2namevalue(p), struct2namevalue(cond), {'ode_solver', solver}, ...
    B.net_args, B.size_args, ...
    {'rng_seeds', B.rng_seeds, 'T_range', B.T_range, 'fs', B.fs, ...
    'store_full_state', false, 'store_decimated_state', false, 'verbose', 'near-none'}, ...
    B.lya_args, {'added_structure', B.added_structure}];
model = OspecNet(args{:});
model.build();
end
