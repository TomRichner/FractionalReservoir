function dS = bm_frozen_rhs(t, S, params, free_idx)
%BM_FROZEN_RHS The model's right-hand side with every state but free_idx held where it is.
%
%   dS = bm_frozen_rhs(t, S, params, free_idx)
%
% Calls SRNNCellTypePairs.dynamics_fast and keeps only the rows free_idx of
% its result; every other row is zero. With free_idx = the dendritic rows
% (params.state_layout.x) the dendritic state evolves exactly as in the
% model, with the adaptation, depression and facilitation variables held at
% the values they had when the integration started: the frozen-adaptation
% arm. A zero derivative row is held exactly by the fixed-step integrators
% (rk4, sra1, ...). Not model code: the model class has no clamped state;
% this wraps its right-hand side from outside.
%
% Ported unchanged from block B pilot 2 (p2_frozen_rhs, 2026-11-12), which
% checked that the arm restarts on the intact trajectory to round-off and
% that the held variables do not drift.
%
% See also: bm_run_network
full = SRNNCellTypePairs.dynamics_fast(t, S, params);
dS = zeros(size(full));
dS(free_idx) = full(free_idx);
end
