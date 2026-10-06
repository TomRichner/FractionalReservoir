function [m,groups,connectivity] = ied_build_network(cfg,input_config,varargin)
% IED_BUILD_NETWORK Build (do not run) the IED network of a frozen config.
% Same construction as run_ied_experiment, which stays untouched as the record
% of the pilot; this version accepts an external-input configuration and
% name-value overrides (T_range, lya_method, ...) so the network can be
% stimulated. Default input is zero, as in the pilot.
if nargin<2 || isempty(input_config)
    input_config=struct('intrinsic_drive',0,'no_stim_pattern',true(1,3));
end
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
args.filter_local_lya=false; args.store_full_state=false; args.verbose='minimal';
fields=fieldnames(input_config);
args.input_config=struct();
for k=1:numel(fields), args.input_config.(fields{k})=input_config.(fields{k}); end
for k=1:2:numel(varargin), args.(varargin{k})=varargin{k+1}; end
nv=struct2namevalue(args);
if strcmp(cfg.topology,'random'), m=SRNNCellTypePairs(nv{:});
else, m=IEDExplorationNetwork(nv{:}); end
m.build(); [W,groups,connectivity]=ied_connectivity(m,cfg);
if ~strcmp(cfg.topology,'random'), m.replace_connectivity(W); end
if strcmp(cfg.topology,'embedded')
    prior_rng=rng; rng(cfg.seed+7001,'twister');
    setpoints=m.S_c_vec; focus=groups>1;
    setpoints(focus)=cfg.focus_sc+cfg.focus_sd*randn(nnz(focus),1);
    rng(prior_rng); m.replace_setpoints(setpoints);
end
end
