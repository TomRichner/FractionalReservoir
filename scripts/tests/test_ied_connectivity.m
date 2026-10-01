% Matrix replacement must reach the cached dynamics without changing states.
[args,~,conds]=srnn_param_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium');
args.n=20; args.indegree=8; args.tau_a=conds{2}.tau_a;
args.synapse_config=conds{2}.synapse_config; args.ode_solver='sra1';
args.lya_method='none'; args.T_range=[0 .1]; args.fs=400;
nv=struct2namevalue(args); probe=IEDExplorationNetwork(nv{:}); probe.build();
state=probe.S0; old=probe.W; W=old*.7; probe.replace_connectivity(W);
assert(isequal(probe.W,W) && isequal(probe.cached_params.W,W));
assert(isequal(probe.S0,state) && all(probe.u_ex==0,'all'));
assert(probe.sigma_u_noise==.025 && all(probe.n_a==1));
probe.replace_setpoints(.7*ones(probe.n,1));
assert(isequal(probe.cached_params.S_c_vec,.7*ones(probe.n,1)));
assert(all(probe.cached_params.activation_function(zeros(probe.n,1))==0));
assert(isequal(probe.W,W) && isequal(probe.S0,state));
disp('IED connectivity replacement: cached matrix, initial states, noise, and 1TS checks passed.');
