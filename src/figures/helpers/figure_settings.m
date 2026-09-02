function s = figure_settings(model)
% FIGURE_SETTINGS The parameters worth recording in a figure's README.
%
%   s = FIGURE_SETTINGS(model)   % model is a built SRNNModel2/SRNNCellTypePairs
%
% Reads the values OFF THE BUILT OBJECT rather than off the preset struct, so
% what the README quotes is what the simulation used -- including anything the
% class filled in itself (auto-filled tau_a, resolved F, the realized spectral
% radius) and anything the caller overrode after the preset was applied.
%
% This is the same principle write_run_parameters_md applies to a sweep
% directory, and the same mistake it exists to prevent: the memory-capacity
% script used to hardcode three settings into its saved struct and two of them
% contradicted the run.
%
% Only properties that EXIST on the given class are read, so one helper serves
% both model classes without a per-class branch.

arguments
    model
end

s = struct();
want = {'n', 'indegree', 'f', 'level_of_chaos', 'tau_d', ...
        'activation', 'S_a', 'S_c', 'c', 'c_E', ...
        'sigma_u_noise', 'ode_solver', 'fs', 'T_range', ...
        'mu_tilde_relative', 'sigma_tilde_relative', ...
        'mu_E_tilde_relative', 'mu_I_tilde_relative', ...
        'F_tracks_network', 'cell_type_names', 'n_a', 'n_a_E', 'n_b_E', ...
        'rng_seeds', 'x0_std'};

for k = 1:numel(want)
    name = want{k};
    if isprop(model, name)
        try
            v = model.(name);
        catch
            continue    % a Dependent that errors in this configuration
        end
        if ~isempty(v) || islogical(v)
            s.(name) = v;
        end
    end
end

% Derived quantities a reader would otherwise have to recompute.
if isprop(model, 'R')
    try; s.R_theoretical = model.R; catch; end %#ok<NOSEMI>
end
if isprop(model, 'N_sys_eqs')
    try; s.N_sys_eqs = model.N_sys_eqs; catch; end %#ok<NOSEMI>
end
if isprop(model, 'W') && ~isempty(model.W)
    s.spectral_radius = max(abs(eig(full(model.W))));
end
end
