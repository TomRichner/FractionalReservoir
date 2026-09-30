function [W,groups,info] = ied_connectivity(model,cfg)
% Fixed anatomical labels for modular runs; diagnostic labels in random runs.
groups=zeros(model.n,1);
for q=1:model.n_cellTypes
    ix=model.type_indices{q};
    groups(ix)=min(cfg.groups,ceil((1:numel(ix))'*cfg.groups/numel(ix)));
end
if strcmp(cfg.topology,'random')
    W=model.W;
else
    previous=rng; guard=onCleanup(@()rng(previous));
    rng(cfg.seed,'twister');
    within=groups==groups';
    probability=cfg.p_between*ones(model.n);
    probability(within)=cfg.p_within;
    means=cfg.wI*ones(model.n);
    means(:,model.type_indices{1})=cfg.wE;
    weights=means+cfg.w_sd*randn(model.n);
    weights(~within)=cfg.between_scale*weights(~within);
    mask=rand(model.n)<probability;
    mask(1:model.n+1:end)=false;
    W=sparse(weights.*mask);
end
within=groups==groups';
info=struct('topology',cfg.topology,'realized_indegree',nnz(W)/model.n, ...
    'within_edges',nnz(W.*within),'between_edges',nnz(W.*~within), ...
    'spectral_abscissa',max(real(eig(full(W)))), ...
    'spectral_radius',max(abs(eig(full(W)))));
end
