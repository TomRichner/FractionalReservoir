function out_dir = run_paper_illustrations(cfg,run_dir)
% RUN_PAPER_ILLUSTRATIONS Compute and save the revised explanatory traces.
% Called only by the analysis master when cfg.illustrations is true.
% Network: paired reference conditions; uniform positive step at the midpoint
% of the displayed [0,30] s interval. Starts at -15 s to settle before display.
% Single neuron: within-bundle n=1/W=0 preset, no adaptation / 1TS SFA / 1TS STD.
% Does not estimate Lyapunov exponents: main 3 supplies those separately.
arguments
    cfg struct
    run_dir char
end
setup_paths();
out_dir=run_dir;
variants={'representative_dynamics','single_neuron'};
presets={cfg.preset_name,cfg.single_neuron_preset};
for j=1:2
    settings=struct('preset_name',presets{j},'display_window',cfg.illustration_display_window, ...
        'step_time',mean(cfg.illustration_display_window),'step_amplitude',cfg.illustration_step_amp, ...
        'T_range',[-15 cfg.illustration_display_window(2)],'fs',400,'seeds',[1 2], ...
        'protocol','Uniform positive input step to all neurons; fixed displayed neuron indices; no LLE computation.');
    [~,~,conds]=srnn_param_preset(presets{j}); titles=srnn_condition_titles();
    results=struct([]);
    input=struct('intrinsic_drive',0,'step_time',settings.step_time, ...
        'amplitude',settings.step_amplitude,'generator',@paper_midpoint_input);
    for c=1:numel(conds)
        model=build_from_preset(presets{j},conds{c}.name,'T_range',settings.T_range, ...
            'input_config',input,'fs',settings.fs,'rng_seeds',settings.seeds, ...
            'lya_method','none','plot_deci',2,'verbose',cfg.verbose);
        model.run(); pd=model.plot_data; pp=model.get_params();
        assert(SRNNCellTypePairs.routes_identical(pp),'Illustration requires identical outgoing routes.');
        r=struct('name',conds{c}.name,'title',titles(conds{c}.name),'t',pd.t, ...
            'u',[],'x',[],'r',[],'syn',[],'sfa',[],'std',[],'selected',{{}},'n',model.n);
        for q=1:model.n_cellTypes
            name=model.cell_type_names{q}; ix=1:min(4,size(pd.x.(name),1)); r.selected{q}=ix;
            r.u=[r.u; pd.u.(name)(ix,:)];
            r.x=[r.x; pd.x.(name)(ix,:)];
            r.r=[r.r; pd.r.(name)(ix,:)];
            post=fieldnames(pd.synaptic_output.(name));
            r.syn=[r.syn; pd.synaptic_output.(name).(post{1})(ix,:)];
            A=pd.a.(name); B=pd.b.(name).(post{1});
            sa=zeros(numel(ix),numel(pd.t)); sb=ones(numel(ix),numel(pd.t));
            if ~isempty(A), sa=pp.c_eff(q)*reshape(sum(A(ix,:,:),2),numel(ix),[]); end
            if ~isempty(B), sb=reshape(prod(B(ix,:,:),2),numel(ix),[]); end
            r.sfa=[r.sfa;sa]; r.std=[r.std;sb];
        end
        results(c)=r;
    end
    folder=fullfile(run_dir,variants{j}); if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,[variants{j} '_data.mat']),'results','settings','-v7.3');
    vprintf(cfg.verbose,'minimal','  saved %s illustrations\n',variants{j});
end
end

