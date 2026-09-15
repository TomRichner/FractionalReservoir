function [u,t]=paper_midpoint_input(params,T,fs,~,config)
% PAPER_MIDPOINT_INPUT Uniform deterministic step; config stores absolute onset.
t=(0:1/fs:T)';
u=zeros(params.n,numel(t));
u(:,t>=config.step_time)=config.amplitude;
end
