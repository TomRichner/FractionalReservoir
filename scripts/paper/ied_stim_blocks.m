function out = ied_stim_blocks(opts)
% IED_STIM_BLOCKS Block-design stimulation experiment on the frozen run13
% IED-like network (ied_base_config), analogous to the clinical recordings'
% separate no-stimulation and stimulation blocks. The simulation step of the
% ied_stim bundle (ied_stim_config, ied_stim_run; FD-031 item 4), moved from
% scripts/explorations/IED_model on 2026-11-17.
%
%   out = ied_stim_blocks('data_dir', d)                  % the defaults below
%   out = ied_stim_blocks('data_dir', d, 'seeds', [42 73], 'n_no', 1, ...)
%   (ied_stim_run passes every field of ied_stim_config(mode).sim)
%
% For each network seed, n_no no-stimulation blocks (uniform input 0) and
% n_stim stimulation blocks (uniform input = amplitude, constant through the
% block) are simulated as separate runs: same network (connectivity and
% set-points from the network seed), different noise seed per block. Each run
% discards its warmup (cfg.warmup) and its first settle s; the remaining
% block_len s are analyzed: candidate events (detector v1, ied_detect_events,
% per block), subnet-mean Welch PSD, mean rate, SFA feedback and depression
% 1 - b_1 (b_1 the first depression variable, available resources, averaged
% over all E-presynaptic synapses; I synapses carry no depression here).
% fs replaces ied_base_config's integration rate: its 400 Hz is the value
% the model was frozen with, not a property of the model; runs from
% 2026-11-17 on use 500 Hz (Tom, FD-031 item 14).
%
% Noise seed of block k (k = 1..n_no no stim, n_no+1..n_no+n_stim stim):
%   base.noise_seed + (seed - base.seed) + 1000*k
% so no block shares noise with stim_20261005 (k = 0 there).
%
% Jobs (seed x block) run in a parfor on the session's pool (pool module).
% One file per job under data_dir, skipped if present. The settings are
% recorded in data_dir/run_settings.mat, and a call with other settings on
% the same folder is refused, so no run folder mixes two runs' blocks.
arguments
    opts.data_dir (1,:) char {mustBeNonempty}
    opts.seeds (1,:) double = [42 73 101 137 211 307 409 503 601 709 811 907]
    opts.amplitude (1,1) double = 0.20
    opts.n_no (1,1) double = 4
    opts.n_stim (1,1) double = 4
    opts.block_len (1,1) double = 300
    opts.settle (1,1) double = 10
    opts.fs (1,1) double = 500
    opts.overwrite (1,1) logical = false
end
data_dir=opts.data_dir;
if ~isfolder(data_dir), mkdir(data_dir); end
check_settings(data_dir,opts);
base=ied_base_config(); base.fs=opts.fs;
nB=opts.n_no+opts.n_stim;
[S,B]=ndgrid(opts.seeds,1:nB); S=S(:); B=B(:);
files=arrayfun(@(s,b)fullfile(data_dir,sprintf('seed%03d_block%d.mat',s,b)),S,B,'UniformOutput',false);
todo=find(~cellfun(@isfile,files) | opts.overwrite);
fprintf('[ied_stim_blocks] %d of %d jobs to run into %s (fs %g Hz)\n',numel(todo),numel(files),data_dir,opts.fs);
if ~isempty(todo), pool_ensure(); end
o=opts;
parfor j=1:numel(todo)
    k=todo(j); run_block(S(k),B(k),files{k},base,o);
end
rows=table();
for k=1:numel(files)
    R=load(files{k},'summary'); rows=[rows;R.summary]; %#ok<AGROW>
end
writetable(removevars(rows,{'psd','psd_f'}),fullfile(data_dir,'summary.csv'));
save(fullfile(data_dir,'summary.mat'),'rows','opts');
out=struct('summary',rows,'data_dir',data_dir,'opts',opts);
end

function check_settings(data_dir,opts)
% Record the run's settings; refuse a folder written with other ones. A
% folder from before run_settings.mat (stim_blocks_20261117) is checked
% against the cfg saved in its first block file (fs and block duration).
f=fullfile(data_dir,'run_settings.mat');
keys={'seeds','amplitude','n_no','n_stim','block_len','settle','fs'};
new=rmfield(opts,setdiff(fieldnames(opts),keys));
if isfile(f)
    old=load(f,'settings'); old=old.settings;
    for k=1:numel(keys)
        assert(isequal(old.(keys{k}),new.(keys{k})),'ied_stim_blocks:SettingsDiffer', ...
            'Run folder %s was written with %s = %s, not %s. Delete it or use another run_dir.', ...
            data_dir,keys{k},mat2str(old.(keys{k})),mat2str(new.(keys{k})));
    end
    return
end
first=dir(fullfile(data_dir,'seed*_block*.mat'));
if ~isempty(first)
    c=load(fullfile(data_dir,first(1).name),'cfg'); c=c.cfg;
    assert(c.fs==opts.fs && abs(c.duration-(opts.settle+opts.block_len))<1e-9, ...
        'ied_stim_blocks:SettingsDiffer', ...
        'Run folder %s holds blocks at fs %g Hz, %g s; this run asks %g Hz, %g s.', ...
        data_dir,c.fs,c.duration,opts.fs,opts.settle+opts.block_len);
end
settings=new; %#ok<NASGU>
save(f,'settings');
end

function run_block(s,b,file,base,o)
started=tic;
stim=b>o.n_no; A=o.amplitude*stim;
cfg=base; cfg.seed=s; cfg.noise_seed=base.noise_seed+(s-base.seed)+1000*b;
cfg.duration=o.settle+o.block_len;
input_config=struct('dc_levels',A,'hold_dur',cfg.duration,'ramp_dur',0, ...
    'noise_intensity',0,'intrinsic_drive',[],'generator',@dc_staircase_stimulus);
[m,groups]=ied_build_network(cfg,input_config,'lya_method','none', ...
    'T_range',[-cfg.warmup cfg.duration],'T_plot',[0 cfg.duration]);
m.run();
p=m.plot_data; t=p.t(:)'; fs=1/median(diff(t));
x=[p.x.E;p.x.I]; r=[p.r.E;p.r.I];
sfa=zeros(1,numel(t)); b1=[];
for q=1:2
    name=m.cell_type_names{q}; a=p.a.(name);
    sfa=sfa+cfg.c(q)*reshape(sum(a,[1 2]),1,[])/cfg.n;
    for post={'E','I'}
        bb=p.b.(name).(post{1});
        if ~isempty(bb), b1=[b1;reshape(bb(:,1,:),size(bb,1),[])]; end %#ok<AGROW>
    end
end
depression=1-mean(b1,1);
G=max(groups); module_x=zeros(G,numel(t));
for g=1:G, module_x(g,:)=mean(x(groups==g,:),1); end
sel=t>=o.settle & t<cfg.duration;
det=ied_detect_events(t(sel),x(:,sel),r(:,sel),groups,struct());
[psd_f,psd]=module_psd(module_x(:,sel),fs);
summary=table(cfg.seed,b,stim,A,cfg.noise_seed,o.block_len,det.metrics.event_count, ...
    det.metrics.event_count/(o.block_len/60),det.metrics.mean_rate,mean(sfa(sel)), ...
    mean(depression(sel)),toc(started), ...
    'VariableNames',{'seed','block','stim','amplitude','noise_seed','block_len','event_count', ...
    'events_per_min','mean_rate','sfa_feedback','depression','seconds'});
summary.psd={psd}; summary.psd_f={psd_f};
ev=[det.events.time];
signals=struct('t',single(t(sel)),'module_x',single(module_x(:,sel)), ...
    'mean_x_E',single(mean(x(m.type_indices{1},sel),1)), ...
    'mean_x_I',single(mean(x(m.type_indices{2},sel),1)),'events',ev,'groups',groups);
save_block(file,summary,signals,cfg);
end

function save_block(file,summary,signals,cfg) %#ok<INUSD>
save(file,'summary','signals','cfg','-v7.3');
end

function [f,psd]=module_psd(channels,fs)
% as ied_stim_experiment: 10-s Hamming windows, 50% overlap, arithmetic
% mean over the subnet-mean channels
win=round(10*fs); noverlap=round(win/2); nfft=max(win,2^nextpow2(win));
psd=[];
for k=1:size(channels,1)
    seg=double(channels(k,:))-mean(double(channels(k,:)));
    [pxx,f]=pwelch(seg,hamming(win),noverlap,nfft,fs);
    if isempty(psd), psd=zeros(size(pxx)); end
    psd=psd+pxx/size(channels,1);
end
f=f(:)'; psd=psd(:)';
end
