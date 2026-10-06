function out = fig_ied_stim(opts)
% FIG_IED_STIM Model counterpart of the clinical stimulation figure: effect of
% sustained uniform input on candidate IED-like events, module-mean PSD, band
% power, pairwise correlation and the adaptation state of the frozen run13
% IED-like network. Reads the per-run files written by ied_stim_experiment;
% nothing is re-simulated.
%
% Panels mirror Fig_Tica_Main / Fig_TicaFullband_PSD_Supp:
%   (A) candidate events per minute without and with stimulation at the
%       selected amplitude, one line per network realization (the no-stimulation
%       value is the unstimulated half of the same run); the full amplitude
%       sweep and the fixed-criterion control are in the stats file;
%   (B) subnet-mean PSD, geometric mean across realizations with 95% CI,
%       no stimulation versus the selected amplitude;
%   (C) band-power change in dB, one marker per realization, red median;
%   (D) change in mean pairwise correlation within and between subnets;
%   (E) change in mean rate, SFA feedback and synaptic resources (percent);
%   (F) example traces, 30 s without and with stimulation, from the
%       realization and amplitude named in opts: mean E and mean I dendritic
%       state, the five subnet means, and beneath them the leading Lyapunov
%       exponent accumulated from the start of the window (needs a run made
%       with ied_stim_experiment(...,'lya',true)). Windows are chosen by a
%       fixed rule (event count closest to the half's median).
% The five 70-neuron modules of ied_connectivity are called subnets here.
arguments
    opts.tag (1,:) char = 'stim_20261005'
    opts.amplitude (1,1) double = 0.20     % the amplitude for panels B-F
    opts.example_seed (1,1) double = 42
    opts.out_dir (1,:) char = ''
    opts.visible (1,1) logical = false
end
setup_paths(); root=fileparts(which('setup_paths'));
data_dir=fullfile(root,'data','IED_model',opts.tag);
if isempty(opts.out_dir), opts.out_dir=fullfile(root,'figs','IED_model',opts.tag); end
if ~isfolder(opts.out_dir), mkdir(opts.out_dir); end
files=dir(fullfile(data_dir,'seed*_amp*.mat')); assert(~isempty(files),'No runs in %s',data_dir);
S=table();
for k=1:numel(files)
    R=load(fullfile(files(k).folder,files(k).name),'summary'); S=[S;R.summary]; %#ok<AGROW>
end
seeds=unique(S.seed)'; amps=unique(S.amplitude)'; nS=numel(seeds);
purple=[0.494 0.247 0.596]; green=[0 0.62 0.45]; red=[0.85 0.2 0.1]; gray=[0.6 0.6 0.6];
fs_tick=14;

fig=figure('Color','w','Visible',onoff(opts.visible),'Units','inches','Position',[0 0 20 12]);
tl=tiledlayout(fig,2,10,'TileSpacing','loose','Padding','compact');

%% (A) paired candidate-event rate, no stimulation versus the selected amplitude
% Rates at every tested level are still gathered for the stats file; the
% panel shows only the two conditions, like the human panel A.
levels=[0 amps]; rate=NaN(nS,numel(levels)); rate_fixed=rate;
for i=1:nS
    base=S(S.seed==seeds(i) & ~S.stim,:); rate(i,1)=mean(base.events_per_min); rate_fixed(i,1)=rate(i,1);
    for j=1:numel(amps)
        row=S(S.seed==seeds(i) & S.stim & S.amplitude==amps(j),:);
        rate(i,j+1)=row.events_per_min; rate_fixed(i,j+1)=row.events_per_min_fixed;
    end
end
ja=find(amps==opts.amplitude)+1;
ax=nexttile(tl,[1 2]); hold(ax,'on');
for i=1:nS, plot(ax,[1 2],[rate(i,1) rate(i,ja)],'-','Color',gray,'LineWidth',1.2); end
scatter(ax,ones(nS,1),rate(:,1),44,purple,'filled','MarkerFaceAlpha',.85);
scatter(ax,2*ones(nS,1),rate(:,ja),44,green,'filled','MarkerFaceAlpha',.85);
xlim(ax,[0.6 2.4]); xticks(ax,[1 2]); xticklabels(ax,{sprintf('No\\newlineStimulation'),'Stimulation'});
ylim(ax,[0 max(rate(:,[1 ja]),[],'all')*1.1]);
ylabel(ax,'Candidate Events per Minute');
style(ax,fs_tick); ax.XAxis.TickLabelInterpreter='tex'; ax.XTickLabelRotation=0; letter(ax,'(A)');

%% paired tables at the selected amplitude
A0=table(); A1=table();
for i=1:nS
    A0=[A0;S(S.seed==seeds(i) & ~S.stim & S.amplitude==opts.amplitude,:)]; %#ok<AGROW>
    A1=[A1;S(S.seed==seeds(i) &  S.stim & S.amplitude==opts.amplitude,:)]; %#ok<AGROW>
end

%% (B) PSD, geometric mean with 95% CI across realizations
ax=nexttile(tl,[1 2]); hold(ax,'on');
f=A0.psd_f{1}; edges=logspace(log10(0.3),log10(100),41); fc=sqrt(edges(1:end-1).*edges(2:end));
P0=binned_log_psd(A0.psd,f,edges); P1=binned_log_psd(A1.psd,f,edges);
[m0,lo0,hi0]=geo_ci(P0); [m1,lo1,hi1]=geo_ci(P1);
fill(ax,[fc fliplr(fc)],[lo0 fliplr(hi0)],purple,'FaceAlpha',.18,'EdgeColor','none');
fill(ax,[fc fliplr(fc)],[lo1 fliplr(hi1)],green,'FaceAlpha',.18,'EdgeColor','none');
h0=plot(ax,fc,m0,'-','Color',purple,'LineWidth',2.2);
h1=plot(ax,fc,m1,'--','Color',green,'LineWidth',2.2);
set(ax,'XScale','log','YScale','log'); xlim(ax,[0.3 100]); xticks(ax,[0.3 1 3 10 30 100]);
xticklabels(ax,{'0.3','1','3','10','30','100'}); ax.XTickLabelRotation=0;
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Subnet-Mean PSD (a.u.^2/Hz)');
legend(ax,[h0 h1],{'No Stimulation','Stimulation'},'Location','southwest','Box','off','FontSize',11);
style(ax,fs_tick); letter(ax,'(B)');

%% (C) band power change (dB)
ax=nexttile(tl,[1 2]); hold(ax,'on');
bands={'delta','theta','alpha','beta'}; names={'Delta','Theta','Alpha','Beta'};
dB=zeros(nS,4);
for b=1:4, dB(:,b)=10*log10(A1.(bands{b})./A0.(bands{b})); end
dot_median(ax,dB,names,red,fs_tick);
ylabel(ax,'Subnet-Mean Power Change (dB)'); letter(ax,'(C)');

%% (D) correlation change
ax=nexttile(tl,[1 2]); hold(ax,'on');
dC=[A1.within_corr-A0.within_corr, A1.between_corr-A0.between_corr];
dot_median(ax,dC,{sprintf('Within\\newlineSubnet'),sprintf('Between\\newlineSubnet')},red,fs_tick);
ylabel(ax,'Mean Correlation Change'); letter(ax,'(D)');

%% (E) adaptation state change (%)
ax=nexttile(tl,[1 2]); hold(ax,'on');
pct=100*[A1.mean_rate./A0.mean_rate-1, A1.sfa_feedback./A0.sfa_feedback-1, A1.resources./A0.resources-1];
dot_median(ax,pct,{'Rate','SFA','Resources'},red,fs_tick);
ylabel(ax,'Change with Stimulation (%)'); letter(ax,'(E)');

%% (F) example traces with the cumulative leading exponent
ex=load(fullfile(data_dir,sprintf('seed%03d_amp%s.mat',opts.example_seed,amp_tag(opts.amplitude))), ...
    'signals','detections','cfg');
sig=ex.signals; t=double(sig.t); G=size(sig.module_x,1);
has_lya=isfield(sig,'lya_t');
if ~has_lya
    warning('fig_ied_stim:NoLyapunov', ...
        'Example run has no Lyapunov tracking; rerun ied_stim_experiment with ''lya'',true,''overwrite'',true.');
end
spacing=0.35; win_len=30;
titles={'No Stimulation',sprintf('Stimulation (Uniform Input %.2g)',opts.amplitude)}; cols={purple,green};
rows_x=[double(sig.mean_x_E);double(sig.mean_x_I);double(sig.module_x)]; nR=size(rows_x,1);
row_names=[{'Mean E','Mean I'} compose('Subnet %d',1:G)];
yl=[-0.25 (nR-1)*spacing+0.45];
t0=zeros(1,2); cum=cell(1,2); cum_t=cell(1,2);
for w=1:2
    ev=[ex.detections(w).det.events.time];
    t0(w)=pick_window(ex.detections(w).window,ev,t,double(sig.module_x),win_len);
    if has_lya
        lt=double(sig.lya_t); ll=double(sig.lya_local); ks=lt>=t0(w) & lt<=t0(w)+win_len;
        cum_t{w}=lt(ks)-t0(w); cum{w}=cumsum(ll(ks))./(1:nnz(ks));
    end
end
if has_lya
    settled=[cum{1}(cum_t{1}>=2) cum{2}(cum_t{2}>=2)];
    lle_ylim=[min([settled 0])-0.1 max([settled 0])+0.1];
else
    lle_ylim=[-1 0.2];
end
for w=1:2
    inner=tiledlayout(tl,5,1,'TileSpacing','none','Padding','tight');
    inner.Layout.Tile=11+5*(w-1); inner.Layout.TileSpan=[1 5];
    ax=nexttile(inner,[4 1]); hold(ax,'on');
    ev=[ex.detections(w).det.events.time]; t1=t0(w)+win_len;
    sel=t>=t0(w) & t<=t1; tt=t(sel)-t0(w);
    for k=1:nR
        trace=rows_x(k,sel); trace=trace-median(trace);
        plot(ax,tt,trace+(nR-k)*spacing,'-','Color',cols{w},'LineWidth',0.9);
    end
    in=ev(ev>=t0(w) & ev<=t1)-t0(w);
    if ~isempty(in), plot(ax,in,((nR-1)*spacing+0.33)*ones(size(in)),'v','Color',[.15 .15 .15],'MarkerFaceColor',[.15 .15 .15],'MarkerSize',6); end
    yticks(ax,(0:nR-1)*spacing); yticklabels(ax,row_names(end:-1:1));
    xlim(ax,[0 win_len]); ylim(ax,yl);
    title(ax,titles{w},'FontWeight','normal','Color',cols{w},'FontSize',fs_tick);
    style(ax,fs_tick); ax.XAxis.Visible='off';
    if w==1, text(ax,-0.07,1.04,'(F)','Units','normalized','FontSize',16,'Clipping','off','Tag','letter'); end
    if w==2
        x0=win_len-0.6; y0=yl(1)+0.02;
        plot(ax,[x0 x0],[y0 y0+0.2],'k-','LineWidth',2);
        text(ax,x0+0.15,y0+0.1,'0.2 a.u.','FontSize',12);
    end
    ax2=nexttile(inner,[1 1]); hold(ax2,'on');
    if has_lya, plot(ax2,cum_t{w},cum{w},'-','Color',cols{w},'LineWidth',1.8); end
    yline(ax2,0,'--','Color',[.5 .5 .5],'LineWidth',1);
    xlim(ax2,[0 win_len]); ylim(ax2,lle_ylim);
    xlabel(ax2,'Time (s)'); ylabel(ax2,sprintf('Cumulative\\newline\\lambda_1 (s^{-1})'));
    style(ax2,fs_tick);
end
set(findall(fig,'-property','FontSize'),'FontSize',fs_tick);
set(findall(fig,'Type','text','Tag','letter'),'FontSize',16);

%% save and record
tag='fig_ied_stim';
exportgraphics(fig,fullfile(opts.out_dir,[tag '.png']),'Resolution',200);
savefig(fig,fullfile(opts.out_dir,[tag '.fig']));
print(fig,fullfile(opts.out_dir,[tag '.svg']),'-dsvg','-vector');
stats=struct('seeds',seeds,'amplitudes',amps,'selected_amplitude',opts.amplitude, ...
    'rate_by_level',rate,'rate_fixed_by_level',rate_fixed,'levels',levels, ...
    'band_dB',dB,'band_names',{names},'corr_change',dC,'state_pct',pct, ...
    'rate_ratio',A1.events_per_min./A0.events_per_min);
save(fullfile(opts.out_dir,[tag '_stats.mat']),'stats');
write_stats_md(fullfile(opts.out_dir,[tag '_stats.md']),stats,A0,A1);
out=struct('fig',fig,'stats',stats,'out_dir',opts.out_dir);
if ~opts.visible, close(fig); out.fig=[]; end
end

%% ---------------------------------------------------------------------------
function s=onoff(b), if b, s='on'; else, s='off'; end, end
function tag=amp_tag(A), tag=strrep(sprintf('%.3f',A),'.','p'); end
function style(ax,fs)
box(ax,'off'); set(ax,'FontSize',fs,'LineWidth',1,'TickDir','out','TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
end
function letter(ax,s)
text(ax,-0.22,1.06,s,'Units','normalized','FontSize',16,'FontWeight','normal','Clipping','off','Tag','letter');
end
function dot_median(ax,M,names,red,fs)
n=size(M,2); rng(7,'twister');
for c=1:n
    jit=0.18*(rand(size(M,1),1)-0.5);
    scatter(ax,c+jit,M(:,c),36,[.25 .25 .25],'filled','MarkerFaceAlpha',.75);
    plot(ax,[c-0.3 c+0.3],median(M(:,c))*[1 1],'-','Color',red,'LineWidth',3);
end
yline(ax,0,'-','Color',[.5 .5 .5],'LineWidth',1);
xlim(ax,[0.4 n+0.6]); xticks(ax,1:n); xticklabels(ax,names); ax.XTickLabelRotation=0;
style(ax,fs); ax.XAxis.TickLabelInterpreter='tex';
end
function P=binned_log_psd(psds,f,edges)
P=NaN(numel(psds),numel(edges)-1);
for i=1:numel(psds)
    p=psds{i};
    for b=1:numel(edges)-1
        sel=f>=edges(b) & f<edges(b+1);
        if any(sel), P(i,b)=mean(p(sel)); end
    end
end
end
function [m,lo,hi]=geo_ci(P)
L=log10(P); n=sum(~isnan(L),1); mu=mean(L,1,'omitnan'); sd=std(L,0,1,'omitnan');
tcrit=tinv(0.975,max(n-1,1)); half=tcrit.*sd./sqrt(n);
m=10.^mu; lo=10.^(mu-half); hi=10.^(mu+half);
end
function [t0,t1]=pick_window(window,ev,t,module_x,win_len)
% Fixed rule: candidate win_len-s windows start every 5 s inside the half
% (dropping its first and last 5 s); take the one whose event count is
% closest to the median count, breaking ties by the RMS closest to the
% median RMS.
starts=window(1)+5:5:window(2)-5-win_len;
counts=arrayfun(@(a)nnz(ev>=a & ev<a+win_len),starts);
rms=arrayfun(@(a)sqrt(mean(var(module_x(:,t>=a & t<a+win_len),0,2))),starts);
score=abs(counts-median(counts))+1e-3*abs(rms-median(rms));
[~,k]=min(score); t0=starts(k); t1=t0+win_len;
end
function write_stats_md(file,stats,A0,A1)
fid=fopen(file,'w'); guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# fig_ied_stim statistics\n\nSelected amplitude %.3g; %d realizations (seeds %s).\n\n', ...
    stats.selected_amplitude,numel(stats.seeds),mat2str(stats.seeds));
fprintf(fid,'## Candidate events per minute by input level (median [min, max])\n\n| Level | Per-window detector | Fixed criterion |\n|---|---|---|\n');
for j=1:numel(stats.levels)
    r=stats.rate_by_level(:,j); rf=stats.rate_fixed_by_level(:,j);
    fprintf(fid,'| %.3g | %.1f [%.1f, %.1f] | %.1f [%.1f, %.1f] |\n',stats.levels(j),median(r),min(r),max(r),median(rf),min(rf),max(rf));
end
fprintf(fid,'\n## At the selected amplitude (median [min, max] across realizations; event-shape rows omit realizations with no events)\n\n| Quantity | No stimulation | Stimulation | Change |\n|---|---|---|---|\n');
q={'events_per_min','Events/min';'mean_rate','Mean rate';'sfa_feedback','SFA feedback';'resources','Synaptic resources'; ...
   'within_corr','Within-subnet correlation';'between_corr','Between-subnet correlation';'median_recruitment','Median recruitment';'waveform_fwhm','Waveform FWHM (s)'};
for k=1:size(q,1)
    a=A0.(q{k,1}); b=A1.(q{k,1});
    fprintf(fid,'| %s | %.3g [%.3g, %.3g] | %.3g [%.3g, %.3g] | %+.3g |\n',q{k,2},median(a,'omitnan'),min(a),max(a),median(b,'omitnan'),min(b),max(b),median(b-a,'omitnan'));
end
fprintf(fid,'\nEvent-rate ratio stimulation / no stimulation: median %.2f, range [%.2f, %.2f]; below 1 in %d of %d.\n\n', ...
    median(stats.rate_ratio),min(stats.rate_ratio),max(stats.rate_ratio),nnz(stats.rate_ratio<1),numel(stats.rate_ratio));
fprintf(fid,'## Band power change (dB), median [min, max]\n\n');
for b=1:4
    fprintf(fid,'- %s: %+.2f [%+.2f, %+.2f]; below 0 in %d of %d\n',stats.band_names{b},median(stats.band_dB(:,b)),min(stats.band_dB(:,b)),max(stats.band_dB(:,b)),nnz(stats.band_dB(:,b)<0),size(stats.band_dB,1));
end
fprintf(fid,'\n## State change (%%), median [min, max]\n\n');
nm={'Mean rate','SFA feedback','Synaptic resources'};
for k=1:3
    fprintf(fid,'- %s: %+.1f [%+.1f, %+.1f]\n',nm{k},median(stats.state_pct(:,k)),min(stats.state_pct(:,k)),max(stats.state_pct(:,k)));
end
end
