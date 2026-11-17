function out = fig_ied_stim_blocks(opts)
% FIG_IED_STIM_BLOCKS Model stimulation figure from the block run written by
% ied_stim_blocks, drawn to mirror the clinical panels (Tom, 2026-11-17).
%
% Main model figure (fig_ied_stim_blocks.{png,fig,svg}), the source of panels
% D-F of fig_stim_combined's mirror layout:
%   traces  30 s of one network's no-stimulation block and stimulation block
%           (opts.example_seed, blocks opts.example_blocks): mean E, mean I and
%           the five subnet means, candidate events as triangles. Window by the
%           fixed rule of fig_ied_stim (event count closest to the block's
%           median over 30-s windows).
%   rate    candidate-event rate normalised as the clinical panel: each block's
%           rate divided by that network's pooled no-stimulation rate (all
%           no-stimulation events / all no-stimulation time), one marker per
%           block, one line per network joining its pooled no-stimulation
%           (100 %) and pooled stimulation rates; log axis with percent ticks,
%           zero-event blocks drawn on a "0" floor below the lowest tick.
%   PSD     subnet-mean PSD, per network the mean over its blocks, then median
%           and IQR across networks (as the clinical SOZ PSD).
% Also written: fig_ied_stim_blocks_rate (the rate panel alone) and
% fig_ied_stim_blocks_supp_state (supplement: % change under stimulation of
% mean rate, SFA feedback and depression 1 - b_1, one point per network from
% its block means, median line), plus _stats.{md,mat}.
arguments
    opts.data_dir (1,:) char = 'C:/Users/m218089/Desktop/github_repos/FractionalReservoir/data/IED_model/stim_blocks_20261117'
    opts.out_dir (1,:) char = 'C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/IED_model/stim_blocks_20261117'
    opts.example_seed (1,1) double = 42
    opts.example_blocks (1,2) double = [1 5]
    opts.floor_pct (1,1) double = 0.5     % where a zero-event block is drawn
    opts.visible (1,1) logical = false
end
setup_paths();
if ~isfolder(opts.out_dir), mkdir(opts.out_dir); end
L=load(fullfile(opts.data_dir,'summary.mat')); S=L.rows; run_opts=L.opts;
seeds=unique(S.seed,'stable')'; nS=numel(seeds);
purple=[0.494 0.247 0.596]; green=[0 0.62 0.45]; red=[0.85 0.2 0.1];
cols=net_colors(nS); fs=14;

% ---- per-network pooled rates and normalised block rates -----------------
pooled=zeros(nS,2); norm_blk=cell(nS,2); state=zeros(nS,3,2); psd=cell(nS,2);
for i=1:nS
    for c=1:2
        T=S(S.seed==seeds(i) & S.stim==(c==2),:);
        pooled(i,c)=sum(T.event_count)/(sum(T.block_len)/60);
        state(i,:,c)=[mean(T.mean_rate) mean(T.sfa_feedback) mean(T.depression)];
        psd{i,c}=mean(cat(1,T.psd{:}),1);
    end
    for c=1:2
        T=S(S.seed==seeds(i) & S.stim==(c==2),:);
        norm_blk{i,c}=100*T.events_per_min/pooled(i,1);
    end
end
ratio=pooled(:,2)./pooled(:,1);
f=S.psd_f{1};

fig=figure('Color','w','Visible',onoff(opts.visible),'Units','inches','Position',[0 0 16 10]);
tl=tiledlayout(fig,2,2,'TileSpacing','loose','Padding','compact');

% ---- traces --------------------------------------------------------------
titles={'No Stimulation',sprintf('Stimulation (Uniform Input %.2g)',run_opts.amplitude)}; tc={purple,green};
spacing=0.35; win_len=30; ex_info=cell(1,2);
for w=1:2
    E=load(fullfile(opts.data_dir,sprintf('seed%03d_block%d.mat',opts.example_seed,opts.example_blocks(w))),'signals');
    sig=E.signals; t=double(sig.t); G=size(sig.module_x,1);
    rows_x=[double(sig.mean_x_E);double(sig.mean_x_I);double(sig.module_x)]; nR=size(rows_x,1);
    names=[{'Mean E','Mean I'} compose('Subnet %d',1:G)];
    ev=sig.events; t0=pick_window([t(1) t(end)],ev,t,double(sig.module_x),win_len);
    ex_info{w}=sprintf('seed %d block %d, %.0f-%.0f s of the analysed block',opts.example_seed,opts.example_blocks(w),t0-t(1),t0-t(1)+win_len);
    ax=nexttile(tl,w); hold(ax,'on');
    sel=t>=t0 & t<=t0+win_len; tt=t(sel)-t0;
    for k=1:nR
        tr=rows_x(k,sel); plot(ax,tt,tr-median(tr)+(nR-k)*spacing,'-','Color',tc{w},'LineWidth',0.9);
    end
    in=ev(ev>=t0 & ev<=t0+win_len)-t0;
    if ~isempty(in), plot(ax,in,((nR-1)*spacing+0.33)*ones(size(in)),'v','Color',[.15 .15 .15],'MarkerFaceColor',[.15 .15 .15],'MarkerSize',6); end
    yticks(ax,(0:nR-1)*spacing); yticklabels(ax,names(end:-1:1));
    xlim(ax,[0 win_len]); ylim(ax,[-0.25 (nR-1)*spacing+0.45]);
    title(ax,titles{w},'FontWeight','normal','Color',tc{w});
    style(ax,fs); ax.XAxis.Visible='off';
    if w==2
        yl=ax.YLim; x0=win_len-0.6; y0=yl(1)+0.02;
        plot(ax,[x0 x0],[y0 y0+0.2],'k-','LineWidth',2);
        text(ax,x0+0.15,y0+0.1,'0.2 a.u.','FontSize',12);
    end
end

% ---- normalised rate ------------------------------------------------------
ax=nexttile(tl,3); rate_axes(ax,norm_blk,pooled,cols,seeds,opts.floor_pct,fs);
% ---- PSD -------------------------------------------------------------------
ax=nexttile(tl,4); hold(ax,'on');
keep=f>=0.3 & f<=100; fk=f(keep); hs=gobjects(1,2); cc={purple,green};
for c=1:2
    P=cat(1,psd{:,c}); P=P(:,keep); q=quantile(P,[0.25 0.5 0.75],1);
    fill(ax,[fk fliplr(fk)],[q(1,:) fliplr(q(3,:))],cc{c},'FaceAlpha',.18,'EdgeColor','none');
    hs(c)=plot(ax,fk,q(2,:),'-','Color',cc{c},'LineWidth',2);
end
set(ax,'XScale','log','YScale','log'); xlim(ax,[0.3 100]); xticks(ax,[0.3 1 3 10 30 100]);
xticklabels(ax,{'0.3','1','3','10','30','100'}); xlabel(ax,'Frequency (Hz)');
ylabel(ax,'Subnet-Mean PSD (a.u.^2/Hz)');
title(ax,sprintf('Median, IQR Over Networks (n = %d)',nS),'FontWeight','normal');
legend(ax,hs,{'No Stim.','Stim.'},'Location','southwest','Box','off');
style(ax,fs);
save_fig(fig,opts.out_dir,'fig_ied_stim_blocks');

% ---- standalone rate figure --------------------------------------------------
f2=figure('Color','w','Visible','off','Units','inches','Position',[0 0 7 4.5]);
rate_axes(axes(f2),norm_blk,pooled,cols,seeds,opts.floor_pct,fs);
save_fig(f2,opts.out_dir,'fig_ied_stim_blocks_rate'); close(f2);

% ---- supplement: state change -------------------------------------------------
pct=100*(state(:,:,2)./state(:,:,1)-1);
f3=figure('Color','w','Visible','off','Units','inches','Position',[0 0 4.5 4]);
ax=axes(f3); hold(ax,'on'); rng(7,'twister');
for c=1:3
    jit=0.18*(rand(nS,1)-0.5);
    scatter(ax,c+jit,pct(:,c),36,[.25 .25 .25],'filled','MarkerFaceAlpha',.75);
    plot(ax,[c-0.3 c+0.3],median(pct(:,c))*[1 1],'-','Color',red,'LineWidth',3);
end
yline(ax,0,'-','Color',[.5 .5 .5],'LineWidth',1);
xlim(ax,[0.4 3.6]); xticks(ax,1:3); xticklabels(ax,{'Rate','SFA','Depression'});
ylabel(ax,'Change with Stimulation (%)'); style(ax,fs);
save_fig(f3,opts.out_dir,'fig_ied_stim_blocks_supp_state'); close(f3);

% ---- stats -------------------------------------------------------------------
stats=struct('seeds',seeds,'pooled_epm',pooled,'ratio',ratio,'norm_blocks',{norm_blk}, ...
    'state_means',state,'state_pct',pct,'example',{ex_info},'run_opts',run_opts);
save(fullfile(opts.out_dir,'fig_ied_stim_blocks_stats.mat'),'stats');
write_stats(fullfile(opts.out_dir,'fig_ied_stim_blocks_stats.md'),stats,S);
out=struct('fig',fig,'stats',stats,'out_dir',opts.out_dir);
if ~opts.visible, close(fig); out.fig=[]; end
end

% ------------------------------------------------------------------------------
function rate_axes(ax,norm_blk,pooled,cols,seeds,floor_pct,fs)
hold(ax,'on'); nS=numel(seeds); h=gobjects(1,nS);
fl=@(v) max(v,floor_pct);
offs=linspace(-0.22,0.22,nS);
yline(ax,100,'-','Color',[.6 .6 .6],'LineWidth',1);
for i=1:nS
    x=[1 2]+offs(i); ps=100*pooled(i,2)/pooled(i,1);
    plot(ax,x,fl([100 ps]),'-','Color',cols(i,:),'LineWidth',1.5);
    for c=1:2
        v=norm_blk{i,c};
        hc=scatter(ax,x(c)*ones(size(v)),fl(v),40,cols(i,:),'filled'); if c==1, h(i)=hc; end
    end
end
set(ax,'YScale','log'); ylim(ax,[floor_pct*0.8 400]);
yt=[floor_pct 1 2 5 10 20 50 100 200]; yticks(ax,yt);
yticklabels(ax,[{'0'} compose('%g%%',yt(2:end))]);
xlim(ax,[0.5 2.5]); xticks(ax,[1 2]); xticklabels(ax,{'No Stimulation','Stimulation'});
ylabel(ax,'Event Rate (% of No Stim.)');
lg=legend(ax,h,compose('N%d',1:nS),'Location','eastoutside','Box','off','NumColumns',2); lg.ItemTokenSize=[8 8];
style(ax,fs);
end
function c=net_colors(n)
c=[lines(7); 0.6 0.6 0.6; 0.74 0.74 0.13; 0.09 0.75 0.81; 0 0 0; 0.55 0.34 0.29];
c=c(1:n,:);
end
function save_fig(fig,d,tag)
exportgraphics(fig,fullfile(d,[tag '.png']),'Resolution',200);
savefig(fig,fullfile(d,[tag '.fig']));
print(fig,fullfile(d,[tag '.svg']),'-dsvg','-vector');
end
function s=onoff(b), if b, s='on'; else, s='off'; end, end
function style(ax,fs)
box(ax,'off'); set(ax,'FontSize',fs,'LineWidth',1,'TickDir','out','TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
end
function t0=pick_window(window,ev,t,module_x,win_len)
starts=window(1)+5:5:window(2)-5-win_len;
counts=arrayfun(@(a)nnz(ev>=a & ev<a+win_len),starts);
rms=arrayfun(@(a)sqrt(mean(var(module_x(:,t>=a & t<a+win_len),0,2))),starts);
[~,k]=min(abs(counts-median(counts))+1e-3*abs(rms-median(rms))); t0=starts(k);
end
function write_stats(file,st,S)
fid=fopen(file,'w'); g=onCleanup(@()fclose(fid));
o=st.run_opts; nS=numel(st.seeds);
fprintf(fid,'# fig_ied_stim_blocks statistics\n\nRun: ied_stim_blocks, amplitude %.2g, %d networks (seeds %s), %d no-stimulation + %d stimulation blocks per network, %g s analysed per block after %g s settle (plus the model warmup); noise seed per block base.noise_seed + (seed - 42) + 1000*block.\n\n', ...
    o.amplitude,nS,mat2str(st.seeds),o.n_no,o.n_stim,o.block_len,o.settle);
fprintf(fid,'Example traces: no stimulation %s; stimulation %s.\n\n',st.example{1},st.example{2});
fprintf(fid,'## Candidate events per minute, pooled per condition, and normalised block rates\n\n| Network | Seed | No stim (ev/min) | Stim (ev/min) | Stim / no stim | No-stim blocks (%% of pooled) | Stim blocks (%% of pooled no stim) |\n|---|---|---|---|---|---|---|\n');
for i=1:nS
    fprintf(fid,'| N%d | %d | %.1f | %.1f | %.2f | %s | %s |\n',i,st.seeds(i),st.pooled_epm(i,1),st.pooled_epm(i,2),st.ratio(i), ...
        strjoin(compose('%.0f',st.norm_blocks{i,1}),', '),strjoin(compose('%.0f',st.norm_blocks{i,2}),', '));
end
fprintf(fid,'\nRatio stimulation / no stimulation: median %.2f, range [%.2f, %.2f]; lower under stimulation in %d of %d networks.\n', ...
    median(st.ratio),min(st.ratio),max(st.ratio),nnz(st.ratio<1),nS);
nb=S(~S.stim,:); sb=S(S.stim,:);
fprintf(fid,'Events per block: no stimulation median %d [%d, %d]; stimulation median %d [%d, %d]; zero-event blocks: %d no stim, %d stim.\n\n', ...
    median(nb.event_count),min(nb.event_count),max(nb.event_count),median(sb.event_count),min(sb.event_count),max(sb.event_count), ...
    nnz(nb.event_count==0),nnz(sb.event_count==0));
fprintf(fid,'## State change with stimulation (%%; per network from block means), median [min, max]\n\n');
nm={'Mean rate','SFA feedback','Depression 1 - b_1'};
for k=1:3
    fprintf(fid,'- %s: %+.1f [%+.1f, %+.1f]\n',nm{k},median(st.state_pct(:,k)),min(st.state_pct(:,k)),max(st.state_pct(:,k)));
end
end
