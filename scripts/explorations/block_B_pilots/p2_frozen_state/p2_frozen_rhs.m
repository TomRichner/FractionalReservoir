function dS = p2_frozen_rhs(t, S, params, free_idx)
% P2_FROZEN_RHS The model's right-hand side with every state but free_idx held.
%
%   dS = p2_frozen_rhs(t, S, params, free_idx)
%
% Calls SRNNCellTypePairs.dynamics_fast and keeps only the rows free_idx of
% its result; every other row is zero. With free_idx = params.state_layout.x
% the dendritic state x evolves exactly as in the model (the frozen a, b, g
% still enter the rate and the synaptic output through S), while SFA, every
% route's depression pools (the reserve of a chain route among them) and
% facilitation stay where they are. Nothing here knows which route model a
% route uses, so it freezes any of them.
%
% A zero derivative row is held EXACTLY by the fixed-step integrators (rk4,
% euler, heun, sra1): each forms y + h * (weighted sum of zeros), and noise
% is added to the x rows only.
%
% Block B pilot 2 (2026-11-12). Not model code: the model class never keeps
% a clamped state variable, and this does not add one. It is a wrapper for
% a rate measurement from outside the class.
full = SRNNCellTypePairs.dynamics_fast(t, S, params);
dS = zeros(size(full));
dS(free_idx) = full(free_idx);
end
