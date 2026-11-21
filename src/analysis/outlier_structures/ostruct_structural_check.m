function c = ostruct_structural_check(plan, key, seed)
%OSTRUCT_STRUCTURAL_CHECK Does one weight draw of a network of the outlier_structures stage have its intended eigenvalue where the name says?
%
%   c = ostruct_structural_check(plan, key, seed)
%
%   plan   outlier_structures_stage_plan(cfg)
%   key    a network's key (plan.networks(k).key)
%   seed   the draw's rng_seeds (1 x 2)
%   c      one row of the stage's check table (a struct, fields below)
%
% Builds the network's model as ospec_run does (ospec_build_model, the
% same name-value pairs and size rule, the plan's first condition: the
% weights and any added structure are drawn from rng_seeds(1) whatever the
% condition), integrates NOTHING, and finds the intended eigenvalue lam0 of
% the structural W exactly as the runs will (ospec_weight_stats,
% ospec_mode_reference). Then, with R the structural bulk radius of the
% drawn matrix (ws.R_struct) and expect the network's expected eigenvalue
% in units of R (outlier_structures_stage_plan, "expect"):
%
%   pass = |lam0 / R - expect| <= plan.check_tol |expect|
%          and Re(lam0) / R >= plan.check_min_re
%
% A network with expect NaN (no outlier intended) passes with reason
% 'no outlier intended'. Costs one build and one dense eig with both
% eigenvector sets of an n x n matrix (about a second at 500 neurons).
%
% Fields of c
%   key, network (name), seed1, seed2, target, group
%   lam0_re_R, lam0_im_R   the intended structural eigenvalue over R (Im >= 0)
%   expect_re, expect_im   where it should be (NaN: nowhere in particular)
%   rel_err      |lam0/R - expect| / |expect|
%   R_struct     the structural bulk radius of the draw (units of W)
%   R_theory     the class's theoretical radius (model.R)
%   is_rightmost the intended eigenvalue is the rightmost of W
%   coh, pr, cond   as ospec_mode_reference
%   pass         logical
%   reason       why it failed ('' when it passed)
%
% See also: run_outlier_structures_stage, outlier_structures_stage_plan,
%           ospec_mode_reference, ospec_build_model

arguments
    plan (1,1) struct
    key  (1,:) char
    seed (1,2) double
end
nw = plan.networks(strcmp({plan.networks.key}, key));
if numel(nw) ~= 1
    error('ostruct_structural_check:Key', 'The plan has no network with key ''%s''.', key);
end
p = srnn_param_preset(plan.preset_name);
net = ospec_network(nw.name, p.mu_tilde_relative);
size_args = {};
if plan.n_override > 0 && plan.n_override ~= p.n
    % ospec_run's size rule (small_network_model): alpha 0.2 and the
    % smaller network's own weight unit.
    k_in = max(2, round(0.2 * plan.n_override));
    size_args = {'n', plan.n_override, 'indegree', k_in, 'F_ref_n', plan.n_override, 'F_ref_indegree', k_in};
end
B = struct('preset_name', plan.preset_name, 'condition', plan.conditions{1}, ...
    'net_args', {net.net_args}, 'size_args', {size_args}, 'added_structure', net.added_structure, ...
    'rng_seeds', seed, 'T_range', [0 1], 'fs', plan.fs, 'lya_args', {{'lya_method', 'none'}});
model = ospec_build_model(B);
params = model.get_params();
ws = ospec_weight_stats(params);
target_value = [];
if ~isempty(model.added); target_value = model.added.target; end
ref = ospec_mode_reference(ws, net.target, target_value);
lam_R = ref.lam0 / ws.R_struct;

% Where the name says the eigenvalue is
if ~isempty(nw.expect)
    expect = nw.expect;
elseif strcmp(net.target, 'lowrank')
    expect = complex(net.added_structure.a, net.added_structure.w);
elseif strcmp(net.target, 'blockmean')
    lo = model.lambda_O(:);
    [~, q] = max(real(lo));
    expect = complex(real(lo(q)), abs(imag(lo(q)))) / model.R;
else
    expect = NaN;                             % the plan refuses this; kept for a direct call
end

if isnan(expect)
    rel_err = NaN; pass = true; reason = 'no outlier intended';
else
    rel_err = abs(lam_R - expect) / abs(expect);
    reasons = {};
    if ~(rel_err <= plan.check_tol)
        reasons{end + 1} = sprintf('the intended eigenvalue is at (%.3g%+.3gi) R, %.0f%% from the expected (%.3g%+.3gi) R (tolerance %.0f%%)', ...
            real(lam_R), imag(lam_R), 100 * rel_err, real(expect), imag(expect), 100 * plan.check_tol);
    end
    if ~(real(lam_R) >= plan.check_min_re)
        reasons{end + 1} = sprintf('its real part %.3g R is below %.3g R', real(lam_R), plan.check_min_re);
    end
    pass = isempty(reasons);
    reason = strjoin(reasons, '; ');
end

expect_im = imag(expect);
if isnan(expect); expect_im = NaN; end
c = struct('key', key, 'network', nw.name, 'seed1', seed(1), 'seed2', seed(2), 'target', net.target, ...
    'group', nw.group, 'lam0_re_R', real(lam_R), 'lam0_im_R', imag(lam_R), ...
    'expect_re', real(expect), 'expect_im', expect_im, 'rel_err', rel_err, ...
    'R_struct', ws.R_struct, 'R_theory', model.R, 'is_rightmost', ref.is_rightmost, ...
    'coh', ref.coh, 'pr', ref.pr, 'cond', ref.cond, 'pass', pass, 'reason', reason);
end
