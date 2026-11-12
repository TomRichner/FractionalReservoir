function [u_ex, t_ex] = p1b_dc_step(params, T, fs, rng_seed, input_config)
% P1B_DC_STEP A DC step of one amplitude to a chosen set of neurons, from t_on to T.
%
%   [u_ex, t_ex] = P1B_DC_STEP(params, T, fs, rng_seed, input_config)
%
% The input_config.generator of pilot 1b's noise-free cases. Same signature
% and the same time grid as dc_staircase_stimulus (t_ex = 0:1/fs:T), so a
% step to every neuron with no noise is the staircase [0 amp] held for t_on
% seconds, sample for sample. The model pads t < 0 (the settling) with zero
% input.
%
% input_config fields
%   amp        DC amplitude, input units (the same units as u)
%   t_on       s; input is amp for t >= t_on, 0 before
%   target     indices of the stimulated neurons (1..n, E first); empty or
%              absent = every neuron
%   intrinsic_drive  required by the model class, unused here
%
% rng_seed is unused: the step has no random part.
%
% See also: dc_staircase_stimulus, p1b_stim_bursting

dt   = 1 / fs;
t_ex = (0:dt:T)';
nt   = numel(t_ex);

target = 1:params.n;
if isfield(input_config, 'target') && ~isempty(input_config.target)
    target = input_config.target;
end
on = (t_ex >= input_config.t_on)';

u_ex = zeros(params.n, nt);
u_ex(target, :) = repmat(input_config.amp * on, numel(target), 1);
end
