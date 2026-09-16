function out_dir = run_paper_illustrations(cfg,run_dir)
% RUN_PAPER_ILLUSTRATIONS Compute and save the revised explanatory traces.
% Called only by the analysis master when cfg.illustrations is true.
% Network: paired reference conditions; uniform positive step at the midpoint
% of the displayed [0,30] s interval. Starts at -15 s to settle before display.
% Single neuron: within-bundle n=1/W=0 preset, no adaptation / 1TS SFA / 1TS STD.
% Network example also saves leading local/finite top-K rates. Accumulation
% begins at t=0 with a separate negative-time alignment interval. Single-neuron
% mechanism columns still skip Lyapunov estimation.
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
        'protocol','Uniform positive input step to all neurons; fixed displayed neuron indices.');
    lya_args={'lya_method','none'};
    settings.lya_method='none';
    if j==1
        settings.lya_method='topk'; settings.lya_K=15; settings.lya_dt=.05;
        settings.lya_T_interval=[cfg.illustration_lya_start_s settings.T_range(2)];
        settings.lya_warmup=cfg.illustration_lya_warmup_s;
        assert(settings.lya_T_interval(1)-settings.lya_warmup>=settings.T_range(1));
        lya_args={'lya_method','topk','lya_K',settings.lya_K,'lya_K_auto',false, ...
            'lya_dt',settings.lya_dt,'lya_T_interval',settings.lya_T_interval, ...
            'lya_warmup',settings.lya_warmup};
    end
    [~,~,conds]=srnn_param_preset(presets{j}); titles=srnn_condition_titles();
    results=[];
    input=struct('intrinsic_drive',0,'step_time',settings.step_time, ...
        'amplitude',settings.step_amplitude,'generator',@paper_midpoint_input);
    for c=1:numel(conds)
        model=build_from_preset(presets{j},conds{c}.name,'T_range',settings.T_range, ...
            'input_config',input,'fs',settings.fs,'rng_seeds',settings.seeds, ...
            'plot_deci',2,'verbose',cfg.verbose,lya_args{:});
        model.run(); pd=model.plot_data; pp=model.get_params();
        assert(SRNNCellTypePairs.routes_identical(pp),'Illustration requires identical outgoing routes.');
        r=struct('name',conds{c}.name,'title',titles(conds{c}.name),'t',pd.t, ...
            'u',[],'x',[],'r',[],'syn',[],'sfa',[],'std',[],'selected',{{}},'n',model.n, ...
            'cell_type_names',{model.cell_type_names},'t_lya',[],'local_lambda',[], ...
            'finite_lambda',[],'LLE',NaN);
        if j==1
            L=model.lya_results;
            % Label estimates at segment END, not before their first data exist.
            r.t_lya=L.t_lya(:)+round(settings.lya_dt*settings.fs)/settings.fs;
            r.local_lambda=L.local_LE_spectrum_t(:,1);
            r.finite_lambda=L.finite_LE_spectrum_t(:,1); r.LLE=L.LLE;
            valid=isfinite(r.finite_lambda);
            assert(any(valid) && min(r.t_lya(valid))>=settings.lya_T_interval(1)+settings.lya_dt-1e-8, ...
                'run_paper_illustrations:EarlyFiniteEstimate','Finite estimates must follow their first completed segment.');
        end
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
        if c==1, results=r; else, results(c)=r; end   % struct([]) has no fields and cannot be indexed into
    end
    folder=fullfile(run_dir,variants{j}); if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,[variants{j} '_data.mat']),'results','settings','-v7.3');
    vprintf(cfg.verbose,'minimal','  saved %s illustrations\n',variants{j});
end
end

