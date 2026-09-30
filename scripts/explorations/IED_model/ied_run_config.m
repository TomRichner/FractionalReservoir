function cfg = ied_run_config(id)
% IED_RUN_CONFIG Frozen single-run choices for the sequential pilot.
cfg = struct('id',id,'n',350,'f',[.5 .5],'indegree',100, ...
    'mu',[7 -7;7 -7],'sigma',1.5*ones(2),'gain',1, ...
    'Sc',[.2 .2],'Sc_sd',[.1 .1],'tau_a',[.25 .25], ...
    'tau_spread',[.25 .25],'c',[.5 .5], ...
    'std_rec',2,'std_rel',.25,'std_on_I',true, ...
    'noise',.025,'fs',400,'duration',50,'warmup',10,'K',8, ...
    'seed',42,'noise_seed',224779,'topology','random','groups',5, ...
    'p_within',.30,'p_between',.003,'wE',.35,'wI',-.50, ...
    'w_sd',.08,'between_scale',.5,'rationale','', ...
    'interpretation','Pending visual review.');
switch id
    case 1
        cfg.rationale='350-neuron mu7-derived 1TS reference; frozen original weight normalization, no input, additive Wiener noise. Establish whether any local events are already present.';
    case 2
        cfg.tau_a=[.75 .75]; cfg.tau_spread=[.5 .5];
        cfg.rationale='Relative to run01, slow the 1TS SFA from 0.25 to 0.75 s and broaden log-normal spread 0.25 to 0.5. Test whether heterogeneous recovery weakens common excursions.';
    case 3
        cfg=ied_run_config(2); cfg.id=3;
        cfg.Sc=[.45 .45]; cfg.Sc_sd=[.15 .15];
        cfg.rationale='Relative to run02, raise setpoint means 0.20 to 0.45 and SD 0.10 to 0.15 to reduce tonic firing and introduce a more excitable-but-quieter heterogeneous operating point.';
    case 4
        cfg=ied_run_config(3); cfg.id=4; cfg.indegree=30;
        cfg.rationale='Relative to run03, reduce expected indegree 100 to 30 while keeping per-edge weights and original frozen normalization unchanged. Test whether sparse recurrence breaks shared pulses; expect weaker feedback.';
    case 5
        cfg=ied_run_config(4); cfg.id=5; cfg.gain=sqrt(100/30);
        cfg.rationale='Relative to run04, restore approximate random-weight variance with gain sqrt(100/30)=1.826 while retaining indegree30. This also increases mean edge weights; it is not an exact matched-feedback control.';
    case 6
        cfg=ied_run_config(5); cfg.id=6; cfg.f=[.7 .3]; cfg.K=25;
        cfg.rationale='Relative to run05, increase excitatory fraction 0.5 to 0.7 to test whether sparse excitation recruits local discharges. Expand diagnostics to top25 per user request; sums across K8 and K25 are not directly comparable.';
    case 7
        cfg=ied_run_config(6); cfg.id=7; cfg.topology='modular'; cfg.gain=1;
        cfg.rationale='Replace random recurrence with five 70-neuron E/I modules. Absolute edge means E=0.35/I=-0.50, SD0.08, within p0.30, between p0.003 and bridge weights scaled0.5. Retain run06 cellular parameters; old mu/sigma/gain do not set the replaced matrix. Aim for independent local recruitment instead of tonic global excitation.';
    case 8
        cfg=ied_run_config(7); cfg.id=8; cfg.wE=.45;
        cfg.rationale='Relative to run07, strengthen only within-module excitatory edge mean 0.35 to 0.45 (bridges use the same presynaptic mean before 0.5 scaling). Seek local amplification and sharper events while retaining weak between-module coupling.';
    otherwise
        error('ied_run_config:UnknownRun','Run %d has not yet been selected.',id);
end
end
