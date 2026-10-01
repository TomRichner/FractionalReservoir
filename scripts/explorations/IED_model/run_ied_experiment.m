function result = run_ied_experiment(selection,run_tag)
% One exploratory network only; call again only after reviewing this run.
% A config struct and a fresh tag also support replay of the frozen base.
if nargin<2, run_tag='20260930'; end
assert(ischar(run_tag) && ~isempty(regexp(run_tag,'^[a-zA-Z0-9_-]+$','once')));
setup_paths();
if isstruct(selection), cfg=selection; else, cfg=ied_run_config(selection); end
id=cfg.id; root=fileparts(which('setup_paths'));
data_dir=fullfile(root,'data','IED_model',run_tag,sprintf('run%02d',id));
fig_dir=fullfile(root,'figs','IED_model',run_tag,sprintf('run%02d',id));
assert(~isfile(fullfile(data_dir,'run.mat')),'Completed run exists; preserve it.');
if ~isfolder(data_dir), mkdir(data_dir); end
if ~isfolder(fig_dir), mkdir(fig_dir); end
base='celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium';
[args,~,conds]=srnn_param_preset(base);
cond=conds{find(cellfun(@(s)strcmp(s.name,'sfa1_std1'),conds),1)};
args.tau_a=cond.tau_a; args.synapse_config=cond.synapse_config;
args.n=cfg.n; args.f=cfg.f; args.indegree=cfg.indegree;
args.mu_tilde_relative=cfg.mu; args.sigma_tilde_relative=cfg.sigma;
args.level_of_chaos=cfg.gain; args.mu_S_c=cfg.Sc; args.sigma_S_c=cfg.Sc_sd;
args.tau_a={cfg.tau_a(1),cfg.tau_a(2)}; args.tau_a_spread=cfg.tau_spread; args.c=cfg.c;
for pre={'E','I'}
    for post={'E','I'}
        if strcmp(pre{1},'I') && ~cfg.std_on_I
            args.synapse_config.(pre{1}).(post{1})=struct();
        else
            args.synapse_config.(pre{1}).(post{1}).std= ...
                struct('tau_rec',cfg.std_rec,'tau_rel',cfg.std_rel);
        end
    end
end
args.sigma_u_noise=cfg.noise; args.noise_seed=cfg.noise_seed;
args.rng_seeds=[cfg.seed cfg.seed+1]; args.ode_solver='sra1'; args.fs=cfg.fs;
args.T_range=[-cfg.warmup cfg.duration]; args.T_plot=[0 cfg.duration]; args.plot_deci=2;
args.lya_method='topk'; args.lya_K=cfg.K; args.lya_K_auto=false;
args.lya_dt=.05; args.lya_warmup=cfg.warmup; args.lya_T_interval=[0 cfg.duration];
args.filter_local_lya=false; args.store_full_state=cfg.K>=25; args.verbose='minimal';
args.input_config.intrinsic_drive=0; args.input_config.no_stim_pattern=true(1,3);
nv=struct2namevalue(args);
if strcmp(cfg.topology,'random'), m=SRNNCellTypePairs(nv{:});
else, m=IEDExplorationNetwork(nv{:}); end
started=tic; m.build(); [W,groups,connectivity]=ied_connectivity(m,cfg);
if ~strcmp(cfg.topology,'random'), m.replace_connectivity(W); end
if strcmp(cfg.topology,'embedded')
    prior_rng=rng; rng(cfg.seed+7001,'twister');
    setpoints=m.S_c_vec; focus=groups>1;
    setpoints(focus)=cfg.focus_sc+cfg.focus_sd*randn(nnz(focus),1);
    rng(prior_rng); m.replace_setpoints(setpoints);
end
assert(all(m.u_ex==0,'all'),'External input must be zero.');
m.run();
p=m.plot_data; x=[p.x.E;p.x.I]; r=[p.r.E;p.r.I];
sfa=zeros(1,numel(p.t)); resources=[];
for q=1:2
    name=m.cell_type_names{q}; a=p.a.(name);
    sfa=sfa+cfg.c(q)*reshape(sum(a,[1 2]),1,[])/cfg.n;
    for post={'E','I'}
        b=p.b.(name).(post{1});
        if ~isempty(b), resources=[resources;reshape(prod(b,2),size(b,1),[])]; end %#ok<AGROW>
    end
end
params=m.get_params();
tau_stats=zeros(2,2);
for q=1:2, tau_stats(q,:)=[mean(params.tau_a_matrix{q}(:)) std(params.tau_a_matrix{q}(:))]; end
data=struct('t',p.t(:)','x',single(x),'r',single(r),'groups',groups, ...
    'E',m.type_indices{1},'I',m.type_indices{2},'sfa',reshape(sfa,1,[]), ...
    'resources',mean(resources,1),'lya',m.lya_results,'W',W, ...
    'tau_stats',tau_stats,'Sc_vec',m.S_c_vec);
if cfg.K>=25
    data.state_t=m.t_out; data.full_state=m.S_out; data.params=m.cached_params;
end
det=ied_detect_events(data.t,data.x,data.r,data.groups,data.lya);
files=ied_plot_run(data,det,cfg,fig_dir);
native=m.plot(); exportgraphics(native,fullfile(fig_dir,'native_model.png'),'Resolution',120); close(native);
seconds=toc(started); save(fullfile(data_dir,'run.mat'),'data','det','cfg','connectivity','seconds','-v7.3');
result=struct('cfg',cfg,'metrics',det.metrics,'connectivity',connectivity, ...
    'tau_stats',tau_stats,'seconds',seconds,'files',{files},'run_tag',run_tag);
ied_record_run(result,root);
vprintf('minimal','minimal','Run %02d: events %d (%.1f/min), width %.3f s, recruit %.3f, global %.3f, r %.3f, corr %.3f, lambda %.3f; %.1f s\n', ...
    id,det.metrics.event_count,det.metrics.events_per_min,det.metrics.median_width, ...
    det.metrics.median_recruitment,det.metrics.global_fraction,det.metrics.mean_rate, ...
    det.metrics.mean_corr,det.metrics.lambda1,seconds);
end

function ied_record_run(result,root)
file=fullfile(root,'docs','IED_model','progress_2026_09_30.md');
cfg=result.cfg; metrics=result.metrics;
fid=fopen(file,'a'); guard=onCleanup(@()fclose(fid));
fprintf(fid,'\n## Run %02d\n\n%s\n\nRun tag: `%s`.\n\n',cfg.id,cfg.rationale,result.run_tag);
fprintf(fid,'Exact selected configuration:\n\n```json\n%s\n```\n\n',jsonencode(cfg,'PrettyPrint',true));
fprintf(fid,'| Metric | Value |\n|---|---:|\n');
fields=fieldnames(metrics);
for k=1:numel(fields), fprintf(fid,'| %s | %.6g |\n',fields{k},metrics.(fields{k})); end
fprintf(fid,'\nRealized indegree %.3f; within/between edges %d/%d; structural abscissa %.4f. ', ...
    result.connectivity.realized_indegree,result.connectivity.within_edges,result.connectivity.between_edges,result.connectivity.spectral_abscissa);
fprintf(fid,'Realized E tau mean/SD %.4f/%.4f s; I %.4f/%.4f s. Runtime %.1f s.\n\n',result.tau_stats',result.seconds);
fprintf(fid,'![Overview](../../figs/IED_model/%s/run%02d/overview.png)\n\n',result.run_tag,cfg.id);
fprintf(fid,'![Event detail](../../figs/IED_model/%s/run%02d/event_zoom.png)\n\n',result.run_tag,cfg.id);
fprintf(fid,'Native plot: [model.plot()](../../figs/IED_model/%s/run%02d/native_model.png). ',result.run_tag,cfg.id);
fprintf(fid,'Saved trajectory/configuration: `data/IED_model/%s/run%02d/run.mat`.\n\n',result.run_tag,cfg.id);
fprintf(fid,'**Interpretation after visual review:** Pending.\n');
end
