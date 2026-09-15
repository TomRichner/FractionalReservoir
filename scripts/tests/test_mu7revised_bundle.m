% Focused non-simulation checks: bundle independence and exact mechanisms.
cfg_test=mu7revised_config();
assert(strcmp(cfg_test.run_mode,'fast'));
assert(strcmp(cfg_test.run_dir,'data/mu7revised_fast'));
assert(strcmp(mu7revised_config('medium').run_dir,'data/mu7revised_medium'));
assert(cfg_test.transient_gain_duration_s/2+5 < cfg_test.transient_gain_duration_s-cfg_test.transient_gain_horizon_s);
[p_test,cls_test,c_test]=srnn_param_preset(cfg_test.preset_name);
[p_old,~,c_old]=srnn_param_preset(sfaEI_mu7_grouped_figures_config().preset_name);
assert(isequaln(p_test,p_old) && isequaln(c_test,c_old));
[p_mc,~,c_mc]=srnn_param_preset(cfg_test.mc_preset);
p_ref=p_test; p_ref.sigma_u_noise=0;
assert(isequaln(p_mc,p_ref) && isequaln(c_mc,c_test));
network_test=build_from_preset(cfg_test.preset_name,'sfa3_std2','lya_method','none', ...
    'T_range',[-15 30],'fs',40,'input_config', ...
    struct('intrinsic_drive',0,'generator',@paper_midpoint_input,'step_time',15,'amplitude',.5));
assert(network_test.n==500 && network_test.n_cellTypes==2);
assert(SRNNCellTypePairs.routes_identical(network_test.get_params()));
assert(all(network_test.u_ex(:,network_test.t_ex<15)==0,'all'));
assert(all(network_test.u_ex(:,network_test.t_ex>=15)==.5,'all'));
[p_one,~,c_one]=srnn_param_preset(cfg_test.single_neuron_preset);
assert(p_one.n==1 && p_one.n_cellTypes==1 && p_one.sigma_u_noise==0);
assert(isempty(c_one{1}.tau_a{1}) && isempty(fieldnames(c_one{1}.synapse_config)));
assert(isscalar(c_one{2}.tau_a{1}) && isempty(fieldnames(c_one{2}.synapse_config)));
assert(isempty(c_one{3}.tau_a{1}) && isscalar(c_one{3}.synapse_config.E.E.std.tau_rec));
for q_test=1:3
    m_test=build_from_preset(cfg_test.single_neuron_preset,c_one{q_test}.name,'lya_method','none', ...
        'T_range',[-15 30],'fs',40,'input_config', ...
        struct('intrinsic_drive',0,'generator',@paper_midpoint_input,'step_time',15,'amplitude',.5));
    assert(all(m_test.W(:)==0));
    assert(all(m_test.u_ex(:,m_test.t_ex<15)==0,'all'));
    assert(all(m_test.u_ex(:,m_test.t_ex>=15)==.5,'all'));
    assert(m_test.t_ex(1)==-15);
end
[u_test,t_test]=paper_midpoint_input(struct('n',2),30,400,1,struct('step_time',15,'amplitude',.5));
assert(all(u_test(:,t_test<15)==0,'all') && all(u_test(:,t_test>=15)==.5,'all'));
bad_test=false;
try
    run_transient_gain('run_mode','fast','duration_s',20,'horizon_s',5);
catch err_test
    bad_test=strcmp(err_test.identifier,'run_transient_gain:ShortTrajectory');
end
assert(bad_test);
disp('PASS: mu7revised presets, mode paths, horizon preflight, 1TS mechanisms and midpoint input; no simulations run.');
