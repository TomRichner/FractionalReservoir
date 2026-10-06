function out = ied_stim_experiment(opts)
% IED_STIM_EXPERIMENT Paired no-stimulation / stimulation protocol on the
% frozen run13 IED-like network (ied_base_config).
%
% One simulation per (network seed, amplitude): uniform input is 0 for the
% first hold s and then steps to the amplitude for another hold s. The first
% settle s of each half are discarded and the two remaining windows are
% analyzed identically: candidate events (detector v1, applied per window
% and also with the no-stimulation excursion criterion held fixed), subnet-mean
% Welch PSDs and band powers, pairwise correlation, mean rate, SFA feedback
% and synaptic resources. Results are saved per run under
% data/IED_model/<tag>/ and never overwritten unless asked. With opts.lya the
% leading Lyapunov exponent is tracked too (slower); the trajectory itself is
% unchanged, so a run can be redone with 'overwrite' to add it.
arguments
    opts.amplitudes (1,:) double = [0.05 0.10 0.20 0.30]
    opts.seeds (1,:) double = 42
    opts.hold (1,1) double = 60
    opts.settle (1,1) double = 10
    opts.tag (1,:) char = 'stim_20261005'
    opts.store_full (1,1) logical = false
    opts.overwrite (1,1) logical = false
    opts.verbose (1,1) logical = true
    opts.lya (1,1) logical = false      % track the leading exponent (top-K QR, K = 8)
end
setup_paths(); root=fileparts(which('setup_paths'));
data_dir=fullfile(root,'data','IED_model',opts.tag);
if ~isfolder(data_dir), mkdir(data_dir); end
base=ied_base_config(); rows=table();
for s=opts.seeds
    for A=opts.amplitudes
        file=fullfile(data_dir,sprintf('seed%03d_amp%s.mat',s,amp_tag(A)));
        if isfile(file) && ~opts.overwrite
            R=load(file,'summary'); rows=[rows;R.summary]; %#ok<AGROW>
            continue
        end
        cfg=base; cfg.seed=s; cfg.noise_seed=base.noise_seed+(s-base.seed);
        cfg.duration=2*opts.hold;
        input_config=struct('dc_levels',[0 A],'hold_dur',opts.hold,'ramp_dur',0, ...
            'noise_intensity',0,'intrinsic_drive',[],'generator',@dc_staircase_stimulus);
        started=tic;
        if opts.lya
            lya_args={'lya_method','topk','lya_K',8,'lya_K_auto',false,'lya_dt',.05, ...
                'lya_warmup',cfg.warmup,'lya_T_interval',[0 cfg.duration]};
        else
            lya_args={'lya_method','none'};
        end
        [m,groups,connectivity]=ied_build_network(cfg,input_config,lya_args{:}, ...
            'T_range',[-cfg.warmup cfg.duration],'T_plot',[0 cfg.duration]);
        m.run();
        run=analyze_run(m,groups,cfg,A,opts);
        run.connectivity=connectivity; run.seconds=toc(started);
        run.summary.seconds(:)=run.seconds;
        summary=run.summary; signals=run.signals; detections=run.detections; %#ok<NASGU>
        if opts.store_full || s==opts.seeds(1)
            full=run.full; %#ok<NASGU>
            save(file,'summary','signals','detections','full','cfg','connectivity','-v7.3');
        else
            save(file,'summary','signals','detections','cfg','connectivity','-v7.3');
        end
        rows=[rows;summary]; %#ok<AGROW>
        if opts.verbose
            fprintf('seed %d amp %.2f: events/min %.1f -> %.1f (fixed thr %.1f), rate %.3f -> %.3f, sfa %.3f -> %.3f, res %.3f -> %.3f; %.0f s\n', ...
                s,A,summary.events_per_min(1),summary.events_per_min(2),summary.events_per_min_fixed(2), ...
                summary.mean_rate(1),summary.mean_rate(2),summary.sfa_feedback(1),summary.sfa_feedback(2), ...
                summary.resources(1),summary.resources(2),run.seconds);
        end
    end
end
writetable(removevars(rows,{'psd','psd_f'}),fullfile(data_dir,'summary.csv'));
out=struct('summary',rows,'data_dir',data_dir,'opts',opts);
end

function tag=amp_tag(A)
tag=strrep(sprintf('%.3f',A),'.','p');
end

function run=analyze_run(m,groups,cfg,A,opts)
p=m.plot_data; t=p.t(:)'; fs=1/median(diff(t));
x=[p.x.E;p.x.I]; r=[p.r.E;p.r.I];
sfa=zeros(1,numel(t)); resources=[];
for q=1:2
    name=m.cell_type_names{q}; a=p.a.(name);
    sfa=sfa+cfg.c(q)*reshape(sum(a,[1 2]),1,[])/cfg.n;
    for post={'E','I'}
        b=p.b.(name).(post{1});
        if ~isempty(b), resources=[resources;reshape(prod(b,2),size(b,1),[])]; end %#ok<AGROW>
    end
end
resources=mean(resources,1);
G=max(groups); module_x=zeros(G,numel(t));
for g=1:G, module_x(g,:)=mean(x(groups==g,:),1); end
u=m.u_ex(1,:); tu=m.t_ex(:)';
windows={[opts.settle opts.hold],[opts.hold+opts.settle 2*opts.hold]};
labels={'no_stim','stim'}; base_excess=[];
summary=table(); detections=struct('label',{},'det',{},'det_fixed',{},'events',{},'window',{});
for w=1:2
    sel=t>=windows{w}(1) & t<windows{w}(2);
    det=ied_detect_events(t(sel),x(:,sel),r(:,sel),groups,struct());
    if w==1
        base_excess=det.threshold-det.baseline; det_fixed=det;
    else
        % Fixed-criterion control: the no-stimulation excursion criterion
        % (threshold above each neuron's median) re-centred on the
        % stimulated segment's own median, so a uniform DC offset does not
        % by itself count as an event.
        det_fixed=ied_detect_events(t(sel),x(:,sel),r(:,sel),groups,struct(),det.baseline+base_excess);
    end
    [events,spatial]=ied_event_details(struct('t',t(sel),'x',x(:,sel),'groups',groups),det);
    [psd_f,psd,band]=module_psd(module_x(:,sel),fs);
    row=table(string(labels{w}),cfg.seed,A,w==2,det.metrics.event_count,det.metrics.events_per_min, ...
        det_fixed.metrics.event_count,det_fixed.metrics.events_per_min, ...
        det.metrics.median_recruitment,det.metrics.median_width, ...
        median(events.WaveformWidth,'omitnan'),median(events.WaveformAmplitude,'omitnan'), ...
        det.metrics.mean_rate,det.metrics.saturation_fraction,mean(sfa(sel)),mean(resources(sel)), ...
        spatial.within_corr,spatial.between_corr,det.metrics.mean_corr, ...
        band(1),band(2),band(3),band(4),NaN, ...
        'VariableNames',{'segment','seed','amplitude','stim','event_count','events_per_min', ...
        'event_count_fixed','events_per_min_fixed','median_recruitment','median_width', ...
        'waveform_fwhm','waveform_amplitude','mean_rate','saturation','sfa_feedback','resources', ...
        'within_corr','between_corr','mean_corr','delta','theta','alpha','beta','seconds'});
    row.psd={psd}; row.psd_f={psd_f};
    summary=[summary;row]; %#ok<AGROW>
    detections(w)=struct('label',labels{w},'det',rmfield(det,'active'),'det_fixed',rmfield(det_fixed,'active'), ...
        'events',events,'window',windows{w});
end
signals=struct('t',t,'module_x',single(module_x),'mean_x_E',single(mean(x(m.type_indices{1},:),1)), ...
    'mean_x_I',single(mean(x(m.type_indices{2},:),1)),'mean_r',single(mean(r,1)),'sfa',single(sfa), ...
    'resources',single(resources),'u',single(interp1(tu,u,t)),'groups',groups, ...
    'E',m.type_indices{1},'I',m.type_indices{2});
lya=m.lya_results;
if ~isempty(lya) && isfield(lya,'local_LE_spectrum_t')
    % Local leading rate sampled every lya_dt, shifted by half a step as in
    % ied_detect_events; the figure accumulates it within each window.
    signals.lya_t=lya.t_lya(:)'+.025; signals.lya_local=lya.local_LE_spectrum_t(:,1)';
    signals.lya_accumulated=lya.finite_LE_spectrum_t(:,1)'; signals.lambda1=lya.LE_spectrum(1);
end
run=struct('summary',summary,'signals',signals,'detections',detections, ...
    'full',struct('x',single(x),'r',single(r)));
end

function [f,psd,band]=module_psd(channels,fs)
% Welch PSD per module-mean channel (10-s Hamming windows, 50% overlap),
% averaged arithmetically across channels, as the clinical SOZ contacts were.
win=round(10*fs); noverlap=round(win/2); nfft=max(win,2^nextpow2(win));
psd=[];
for k=1:size(channels,1)
    seg=double(channels(k,:))-mean(double(channels(k,:)));
    [pxx,f]=pwelch(seg,hamming(win),noverlap,nfft,fs);
    if isempty(psd), psd=zeros(size(pxx)); end
    psd=psd+pxx/size(channels,1);
end
f=f(:)'; psd=psd(:)';
bands=[1 4;4 8;8 13;13 30]; band=zeros(1,4);
for b=1:4
    sel=f>=bands(b,1) & f<bands(b,2);
    band(b)=trapz(f(sel),psd(sel));
end
end
