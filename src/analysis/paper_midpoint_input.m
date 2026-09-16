function [u,t]=paper_midpoint_input(params,T,fs,~,config)
% PAPER_MIDPOINT_INPUT Uniform deterministic step to every neuron.
%
% config.step_time is the absolute onset; config.step_off (optional, default
% Inf) the absolute offset, so a step that occupies the MIDDLE THIRD of the
% display window is step_time = 10, step_off = 20 for a [0 30] s window (TR,
% 2026-09-15: the mu7revised bundle's step stayed on to the end).
t=(0:1/fs:T)';
u=zeros(params.n,numel(t));
step_off=Inf; if isfield(config,'step_off') && ~isempty(config.step_off), step_off=config.step_off; end
u(:,t>=config.step_time & t<step_off)=config.amplitude;
end
