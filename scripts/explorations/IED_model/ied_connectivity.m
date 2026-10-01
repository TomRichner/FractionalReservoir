function [W,groups,info] = ied_connectivity(model,cfg)
% Fixed anatomical labels for modular runs; diagnostic labels in random runs.
previous=rng; guard=onCleanup(@()rng(previous));
groups=zeros(model.n,1);
if strcmp(cfg.topology,'embedded')
    rng(cfg.seed+9173,'twister'); groups(:)=1;
    for q=1:model.n_cellTypes
        ix=model.type_indices{q}; ix=ix(randperm(numel(ix)));
        count=round(cfg.focus_n*cfg.f(q));
        if q==model.n_cellTypes, count=cfg.focus_n-round(cfg.focus_n*cfg.f(1)); end
        assert((cfg.groups-1)*count<=numel(ix));
        for j=1:cfg.groups-1
            groups(ix((j-1)*count+(1:count)))=j+1;
        end
    end
else
    for q=1:model.n_cellTypes
        ix=model.type_indices{q};
        groups(ix)=min(cfg.groups,ceil((1:numel(ix))'*cfg.groups/numel(ix)));
    end
end
if strcmp(cfg.topology,'random')
    W=model.W;
else
    rng(cfg.seed,'twister');
    within=groups==groups';
    if strcmp(cfg.topology,'embedded')
        probability=cfg.p_background*ones(model.n);
        focus_pairs=groups>1 & groups'>1;
        probability(focus_pairs & ~within)=cfg.p_focus_bridge;
        probability(focus_pairs & within)=cfg.p_within;
    else
        probability=cfg.p_between*ones(model.n);
        probability(within)=cfg.p_within;
    end
    means=cfg.wI*ones(model.n);
    means(:,model.type_indices{1})=cfg.wE;
    weight_sd=cfg.w_sd*ones(model.n);
    if strcmp(cfg.topology,'embedded') && isfield(cfg,'background_E')
        core=focus_pairs & within;
        background_means=cfg.background_I*ones(model.n);
        background_means(:,model.type_indices{1})=cfg.background_E;
        means(~core)=background_means(~core);
        weight_sd(~core)=cfg.background_sd;
    end
    weights=means+weight_sd.*randn(model.n);
    if ~strcmp(cfg.topology,'embedded')
        weights(~within)=cfg.between_scale*weights(~within);
    end
    mask=rand(model.n)<probability;
    mask(1:model.n+1:end)=false;
    W=sparse(weights.*mask);
end
within=groups==groups';
info=struct('topology',cfg.topology,'realized_indegree',nnz(W)/model.n, ...
    'within_edges',nnz(W.*within),'between_edges',nnz(W.*~within), ...
    'spectral_abscissa',max(real(eig(full(W)))), ...
    'spectral_radius',max(abs(eig(full(W)))));
info.group_sizes=accumarray(groups,1)';
end
