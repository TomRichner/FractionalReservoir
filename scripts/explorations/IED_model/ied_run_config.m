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
cfg.focus_n=35; cfg.focus_sc=.70; cfg.focus_sd=.03;
cfg.p_background=.01; cfg.p_focus_bridge=.001;
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
    case 9
        cfg=ied_run_config(8); cfg.id=9; cfg.wE=.60;
        cfg.rationale='Relative to run08, increase E mean 0.45 to 0.60 to move modules closer to discharge-generating feedback. Test whether weak broad humps become sharper cluster pulses or excessive tonic activity.';
    case 10
        cfg=ied_run_config(9); cfg.id=10; cfg.std_on_I=false;
        cfg.rationale='Relative to run09, retain depression on E outputs but remove it on I outputs. Preserve inhibitory transmission during an excitatory excursion so feedback can terminate it sharply. All other cellular and connectivity parameters stay fixed.';
    case 11
        cfg=ied_run_config(10); cfg.id=11; cfg.wI=-.40;
        cfg.rationale='Relative to run10, weaken inhibitory edge mean -0.50 to -0.40 (E stays0.60). Direct E:I weight-balance test: allow larger local recruitment while retaining non-depressing inhibitory termination.';
    case 12
        cfg=ied_run_config(10); cfg.id=12; cfg.wI=-.60;
        cfg.rationale='Relative to run10, strengthen I mean -0.50 to -0.60, bracketing the -0.40 run11. Direct opposite E:I structural-balance perturbation; E stays0.60 and inhibitory outputs remain non-depressing.';
    case 13
        cfg=ied_run_config(11); cfg.id=13;
        cfg.tau_a=[.40 .75]; cfg.tau_spread=[.70 .50];
        cfg.rationale='Relative to leading candidate run11, speed E SFA median0.75 to0.40 s and broaden its log-SD0.5 to0.7; I remains0.75/0.5. Test whether faster heterogeneous E feedback narrows local pulses without suppressing them.';
    case 14
        cfg=ied_run_config(13); cfg.id=14; cfg.noise=.05;
        cfg.rationale='Relative to run13, double input-referred Wiener amplitude0.025 to0.05 on the same noise seed. Test whether stronger independent stochastic drive gives more irregular or stronger cluster events without global recruitment.';
    case 15
        cfg=ied_run_config(13); cfg.id=15; cfg.seed=73; cfg.noise_seed=224810;
        cfg.rationale='Freeze run13 physical parameters, then change network/heterogeneity seed42 to73 and Wiener seed224779 to224810. One new-realization check before the 15-run checkpoint; no parameter retuning to this seed.';
    case 16
        cfg=ied_run_config(13); cfg.id=16;
        cfg.groups=10; cfg.Sc=[.50 .50]; cfg.p_between=.001;
        cfg.rationale='Extra targeted trial after Tom asked for sparser, more clustered and less rhythmic events: halve module size70 to35 neurons (10 modules, same within p0.30/per-edge weights), raise setpoint means0.45 to0.50, and reduce between-edge p0.003 to0.001. This intentionally weakens feedback toward a quieter subthreshold background; preserve 1TS/no-input/Wiener settings.';
    case 17
        cfg=ied_run_config(16); cfg.id=17; cfg.wE=.80;
        cfg.rationale='Final targeted adjustment to run16: raise E edge mean0.60 to0.80 while retaining 35-neuron modules, higher setpoints and p_between0.001. Seek visibly sharper local discharges on the quieter scaffold without returning to the 70-neuron rhythmic regime.';
    case 18
        cfg=ied_run_config(17); cfg.id=18;
        cfg.topology='embedded'; cfg.groups=4; cfg.p_within=.60;
        cfg.rationale='New architecture requested by Tom: three scattered dense 35-neuron E/I foci embedded in a sparse random background. Background p0.01 and Sc mean0.50/SD0.15; within-focus p0.60 and Sc mean0.70/SD0.03. Direct focus-focus p0.001. E/I means0.80/-0.40, E-only STD, heterogeneous 1TS SFA and Wiener noise remain. High focus setpoints aim for quiet intervals followed by recurrent regenerative bursts.';
    case 19
        cfg=ied_run_config(18); cfg.id=19;
        cfg.background_E=.15; cfg.background_I=-.20; cfg.background_sd=.04;
        cfg.rationale='Relative to run18, retain dense core weights0.80/-0.40 (SD0.08), but weaken every non-core route to E0.15/I-0.20 (SD0.04). This includes background, core-background and direct core-core bridges. Test whether dense high-Sc foci remain burst-capable without global cascades.';
    otherwise
        error('ied_run_config:UnknownRun','Run %d has not yet been selected.',id);
end
end
