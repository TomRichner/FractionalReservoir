function M = s2b_build(cname, seeds, varargin)
%S2B_BUILD The reduced network of sub-experiment C of run_numerics_verification
% (qr_trial), run once with the full state kept.
%   M = s2b_build(cname, seeds, 'n', 30, 'T_range', [-8 4], 'fs', 400)
% Same construction as the stage's build_probe + qr_trial: rng(0,'twister'),
% preset < condition < overrides, sigma_u_noise 0, ode45 (RelTol = AbsTol =
% 1e-9, the class default), indegree max(2, round(0.2 n)), F_tracks_network.
% Returns t, S (nt x N), params (with the stimulus interpolant), rhs, layout.
o = struct('n', 30, 'T_range', [-8 4], 'fs', 400, ...
    'preset', 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain');
for k = 1:2:numel(varargin); o.(varargin{k}) = varargin{k + 1}; end
rng(0, 'twister');
[preset, ~, conditions] = srnn_param_preset(o.preset);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond  = rmfield(conditions{strcmp(names, cname)}, 'name');
args  = [struct2namevalue(preset), struct2namevalue(cond), ...
    {'rng_seeds', seeds, 'verbose', 'near-none', 'n', o.n, 'indegree', max(2, round(0.2 * o.n)), ...
    'F_tracks_network', true, 'sigma_u_noise', 0, 'ode_solver', 'ode45', 'fs', o.fs, ...
    'T_range', o.T_range, 'lya_method', 'none', 'store_full_state', true}];
model = SRNNNumericsProbe(args{:});
model.build();
model.run();
params = model.get_params();
params.u_interpolant = model.u_interpolant;
M = struct('cname', cname, 'seeds', seeds, 'n', o.n, 'fs', o.fs, 'T_range', o.T_range, ...
    't', model.t_out(:), 'S', model.S_out, 'params', params, 'N', model.N_sys_eqs, ...
    'layout', params.state_layout, 'ode_opts', model.ode_opts);
M.rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, params);
end
