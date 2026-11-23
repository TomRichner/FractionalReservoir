function out = fig_ied_stim_blocks(opts)
% FIG_IED_STIM_BLOCKS Model stimulation figure from the block run written by
% ied_stim_blocks, drawn to mirror the clinical panels (Tom, 2026-11-17).
%
% Main model figure (fig_ied_stim_blocks.{png,fig,svg}), the source of panels
% D-F of fig_stim_combined's mirror layout:
%   traces  5 s (win_len; 30 s until 2026-11-22) of one network's no-stimulation and stimulation block
%           (opts.example_seed, blocks opts.example_blocks): mean E, mean I and
%           the five subnet means, candidate events as triangles along the
%           top and open circles on the subnets that took part (event_rows),
%           a triangle only where there is a circle. Window by rule
%           (choose_example; Tom, 2026-11-23, R11): the no-stimulation
%           window with 3 to 6 marked detections (nearest 4.5, then the
%           most typical rms),
%           the stimulation window at the same offset in the matching
%           block. Until then (window_rule 'median_count') the rule of
%           fig_ied_stim: event count closest to the block's median.
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
%
%   out = fig_ied_stim_blocks('data_dir', run_dir, 'out_dir', d)          % all three
%   out = fig_ied_stim_blocks(..., 'which', 'main' | 'rate' | 'supp_state')
%
% A figure of the ied_stim bundle (ied_stim_config registers one entry per
% 'which', each into its own folder; FD-031 item 4); moved from
% scripts/explorations/IED_model on 2026-11-17. The stats files are written
% with 'main' (and 'all'). example_seed defaults to the run's first network,
% example_blocks to its first no-stimulation and first stimulation block.
arguments
    opts.data_dir (1,:) char {mustBeNonempty}
    opts.out_dir (1,:) char {mustBeNonempty}
    opts.which (1,:) char {mustBeMember(opts.which,{'all','main','rate','supp_state'})} = 'all'
    opts.example_seed (1,1) double = NaN
    opts.example_blocks (1,2) double = [NaN NaN]
    opts.floor_pct (1,1) double = 0.5     % where a zero-event block is drawn
    opts.visible (1,1) logical = false
    % Figure 8 D (Tom, 2026-11-22; Plan 3 phase 0, question 6): 5 s of
    % signal, open circles on the subnet traces at the detections, a mark
    % along the top only where a drawn trace has a circle (as the patient
    % panel, top_marks_drawn_only), a 0.5 a.u. and 1 s scale bar.
    opts.win_len (1,1) double = 5         % s of signal drawn (30 until 2026-11-22)
    opts.amp_bar (1,1) double = 0.5       % a.u., the amplitude scale bar
    opts.time_bar (1,1) double = 1        % s, the time scale bar
    % The example window (Tom, 2026-11-23, R11: "find another example that
    % has a few more events marked for the no-stimulation case"). 'marked':
    % the rule of choose_example below; 'median_count': the rule until then
    % (pick_window in each example block, event count closest to its median).
    opts.window_rule (1,:) char {mustBeMember(opts.window_rule,{'marked','median_count'})} = 'marked'
    opts.marked_range (1,2) double = [3 6]   % marked detections wanted in the no-stimulation window
end
if ~isfolder(opts.out_dir), mkdir(opts.out_dir); end
L=load(fullfile(opts.data_dir,'summary.mat')); S=L.rows; run_opts=L.opts;
if isnan(opts.example_seed), opts.example_seed=S.seed(1); end
blocks_given=~any(isnan(opts.example_blocks));
if ~blocks_given, opts.example_blocks=[1 run_opts.n_no+1]; end
want=@(w) any(strcmp(opts.which,{'all',w}));
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
pct=100*(state(:,:,2)./state(:,:,1)-1);
fig=[]; ex_info={'',''}; rule_note='';

if want('main')
fig=figure('Color','w','Visible',onoff(opts.visible),'Units','inches','Position',[0 0 16 10]);
tl=tiledlayout(fig,2,2,'TileSpacing','loose','Padding','compact');

% ---- traces --------------------------------------------------------------
titles={'No Stimulation',sprintf('Stimulation (Uniform Input %.2g)',run_opts.amplitude)}; tc={purple,green};
spacing=0.35; win_len=opts.win_len; ex_info=cell(1,2);
load_block=@(b) getfield(load(fullfile(opts.data_dir,sprintf('seed%03d_block%d.mat',opts.example_seed,b)),'signals'),'signals');
if strcmp(opts.window_rule,'marked')
    no_blocks=1:run_opts.n_no; if blocks_given, no_blocks=opts.example_blocks(1); end
    [ex_blocks,ex_off,rule_note]=choose_example(load_block,no_blocks,opts.example_blocks(2),blocks_given, ...
        run_opts.n_no,run_opts.n_stim,win_len,opts.marked_range);
else
    ex_blocks=opts.example_blocks; ex_off=[NaN NaN];
    rule_note='each block''s window with the event count closest to its median over windows (pick_window)';
end
for w=1:2
    sig=load_block(ex_blocks(w)); t=double(sig.t); G=size(sig.module_x,1);
    rows_x=[double(sig.mean_x_E);double(sig.mean_x_I);double(sig.module_x)]; nR=size(rows_x,1);
    names=[{'Mean E','Mean I'} compose('Subnet %d',1:G)];
    ev=sig.events;
    if isnan(ex_off(w)), t0=pick_window([t(1) t(end)],ev,t,double(sig.module_x),win_len); else, t0=t(1)+ex_off(w); end
    ex_info{w}=sprintf('seed %d block %d, %.0f-%.0f s of the analysed block',opts.example_seed,ex_blocks(w),t0-t(1),t0-t(1)+win_len);
    ax=nexttile(tl,w); hold(ax,'on');
    sel=t>=t0 & t<=t0+win_len; tt=t(sel)-t0;
    off=zeros(1,nR);
    for k=1:nR
        tr=rows_x(k,sel); off(k)=median(tr);
        plot(ax,tt,tr-off(k)+(nR-k)*spacing,'-','Color',tc{w},'LineWidth',0.9);
    end
    % Circles: each subnet that took part in a detection, at its own peak
    % (event_rows); a mark along the top only for a detection with a circle.
    in=ev(ev>=t0 & ev<=t0+win_len);
    [has,pk_t,pk_x]=event_rows(t,double(sig.module_x),in);
    for j=1:numel(in)
        for g=find(has(:,j))'
            k=2+g;   % the subnet's row, under Mean E and Mean I
            plot(ax,pk_t(g,j)-t0,pk_x(g,j)-off(k)+(nR-k)*spacing,'o','Color','k', ...
                'MarkerFaceColor','none','MarkerSize',7,'LineWidth',1,'Tag','event_circle');
        end
    end
    top=in(any(has,1))-t0;
    if ~isempty(top), plot(ax,top,((nR-1)*spacing+0.33)*ones(size(top)),'v','Color',[.15 .15 .15],'MarkerFaceColor',[.15 .15 .15],'MarkerSize',6,'Tag','event_mark'); end
    ex_info{w}=sprintf('%s; %d detections in the window, %d with a circle, %d circles',ex_info{w}, ...
        numel(in),numel(top),nnz(has));
    yticks(ax,(0:nR-1)*spacing); yticklabels(ax,names(end:-1:1));
    xlim(ax,[0 win_len]); ylim(ax,[-0.25 (nR-1)*spacing+0.45]);
    title(ax,titles{w},'FontWeight','normal','Color',tc{w});
    style(ax,fs); ax.XAxis.Visible='off';
    if w==2
        % As the patient panel: an L of a time and an amplitude bar right of
        % the strip and under it, the amplitude label along its bar.
        yl=ax.YLim; x1=win_len+0.15; y0=yl(1)-0.05;
        plot(ax,[x1-opts.time_bar x1],[y0 y0],'k-','LineWidth',1.5,'Clipping','off','Tag','scalebar');
        plot(ax,[x1 x1],[y0 y0+opts.amp_bar],'k-','LineWidth',1.5,'Clipping','off','Tag','scalebar');
        text(ax,x1-opts.time_bar/2,y0-0.02,sprintf('%g s',opts.time_bar),'FontSize',12,'HorizontalAlignment','center', ...
            'VerticalAlignment','top','Clipping','off','Tag','scalebar_label');
        text(ax,x1+0.06,y0+opts.amp_bar/2,sprintf('%g a.u.',opts.amp_bar),'FontSize',12,'Rotation',90, ...
            'HorizontalAlignment','center','VerticalAlignment','top','Clipping','off','Tag','scalebar_label');
    end
end

% ---- normalised rate ------------------------------------------------------
ax=nexttile(tl,3); rate_axes(ax,norm_blk,pooled,cols,seeds,opts.floor_pct,fs,false);   % no key (Figure 8 E)
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
end

% ---- standalone rate figure --------------------------------------------------
if want('rate')
f2=figure('Color','w','Visible','off','Units','inches','Position',[0 0 7 4.5]);
rate_axes(axes(f2),norm_blk,pooled,cols,seeds,opts.floor_pct,fs,true);
save_fig(f2,opts.out_dir,'fig_ied_stim_blocks_rate'); close(f2);
end

% ---- supplement: state change -------------------------------------------------
if want('supp_state')
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
end

% ---- stats -------------------------------------------------------------------
stats=struct('seeds',seeds,'pooled_epm',pooled,'ratio',ratio,'norm_blocks',{norm_blk}, ...
    'state_means',state,'state_pct',pct,'example',{ex_info},'example_rule',rule_note,'run_opts',run_opts);
if want('main')
    save(fullfile(opts.out_dir,'fig_ied_stim_blocks_stats.mat'),'stats');
    write_stats(fullfile(opts.out_dir,'fig_ied_stim_blocks_stats.md'),stats,S);
end
out=struct('fig',fig,'stats',stats,'out_dir',opts.out_dir);
if ~opts.visible && ~isempty(fig), close(fig); out.fig=[]; end
end

% ------------------------------------------------------------------------------
function rate_axes(ax,norm_blk,pooled,cols,seeds,floor_pct,fs,with_key)
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
if with_key
    lg=legend(ax,h,compose('N%d',1:nS),'Location','eastoutside','Box','off','NumColumns',2); lg.ItemTokenSize=[8 8];
end
style(ax,fs);
end
function c=net_colors(n)
% The colors of the patient rate panel's subjects S1 to S12, in order (its
% lines, read from cscs_stim_stability's Fig_Awake_Main_PanelA_Normalized.fig
% on 2026-11-22), so that network k of Figure 8 E has subject k's color in B
% (Tom, 2026-11-22). Until then: lines(7), then five more.
c=[0.12 0.47 0.71; 1 0.5 0.05; 0.84 0.15 0.16; 0.55 0.34 0.29; 0.89 0.47 0.76; 0.5 0.5 0.5; ...
   0.09 0.75 0.81; 0.74 0.74 0.13; 0 0 0; 0.68 0.78 0.91; 1 0.73 0.47; 0.6 0.2 0.5];
c=c(1+mod(0:n-1,size(c,1)),:);
end
function [has,pk_t,pk_x]=event_rows(t,module_x,ev)
% Which subnets took part in each detection, and where to circle them.
% The block files keep the detection times only, not the groups of each
% detection (ied_detect_events works on every neuron, ied_stim_blocks
% saves the subnet means), so the rule is applied to the subnet means: a
% subnet takes part when its mean dendritic potential, smoothed over 25 ms
% as the detector smooths, rises within 0.12 s of the detection (the
% detector's merge window) above its block median plus 3 robust SD (the
% detector's per-neuron rule). The circle is at that subnet's maximum in
% the window, on the drawn (unsmoothed) trace.
fsamp=1/median(diff(t));
xs=movmean(module_x,max(1,round(.025*fsamp)),2);
base=median(xs,2); rsd=1.4826*median(abs(xs-base),2);
G=size(module_x,1); n=numel(ev);
has=false(G,n); pk_t=nan(G,n); pk_x=nan(G,n);
for j=1:n
    k=find(abs(t-ev(j))<=.12);
    for g=1:G
        if max(xs(g,k))>base(g)+3*rsd(g)
            [~,i]=max(module_x(g,k)); has(g,j)=true; pk_t(g,j)=t(k(i)); pk_x(g,j)=module_x(g,k(i));
        end
    end
end
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
function [blocks,off,note]=choose_example(load_block,no_blocks,stim_block,blocks_given,n_no,n_stim,win_len,range)
% THE EXAMPLE WINDOW of Figure 8 D, by rule (Tom, 2026-11-23, R11). Over the
% example network's no-stimulation blocks (no_blocks), every window of
% win_len s on pick_window's grid (starts 5 s apart, from 5 s into the
% block to win_len + 5 s before its end) is scored by its MARKED detections:
% the detections that get a top mark, i.e. at least one circled subnet
% (event_rows, the rule of the drawing; the detector is not touched). The
% windows with range(1) to range(2) marked detections are the candidates;
% of them the one whose marked count is closest to the middle of the range,
% then whose signal rms is closest to the median rms over all windows
% (pick_window's tie-break: a typical stretch), then the earliest. With no
% candidate, the window with the most marked detections. The stimulation window is the MATCHING one: the same offset
% in the stimulation block of the same rank (no-stimulation block k goes
% with block n_no + k; the given stimulation block when blocks were given).
cand=zeros(0,5);    % block, start offset, marked, rms, rank
for b=no_blocks
    sig=load_block(b); t=double(sig.t); mx=double(sig.module_x); ev=sig.events;
    [has]=event_rows(t,mx,ev); marked=ev(any(has,1));
    starts=t(1)+5:5:t(end)-5-win_len;
    for a=starts
        cand(end+1,:)=[b, a-t(1), nnz(marked>=a & marked<a+win_len), ...
            sqrt(mean(var(mx(:,t>=a & t<a+win_len),0,2))), size(cand,1)]; %#ok<AGROW>
    end
end
in=cand(:,3)>=range(1) & cand(:,3)<=range(2);
if any(in)
    c=cand(in,:); [~,k]=min(1e3*abs(c(:,3)-mean(range))+abs(c(:,4)-median(cand(:,4)))+1e-9*c(:,5));
    how=sprintf('of the %d windows with %d to %d marked detections, the one with its count closest to %g, then rms closest to the median over all %d windows', ...
        nnz(in),range(1),range(2),mean(range),size(cand,1));
else
    c=cand; [~,k]=max(c(:,3)-1e-9*c(:,5));
    how=sprintf('no window had %d to %d marked detections; the one with the most (%d)',range(1),range(2),c(k,3));
end
b=c(k,1); off=[c(k,2) c(k,2)];
if blocks_given, sb=stim_block; else, sb=n_no+min(find(no_blocks==b,1),n_stim); end
blocks=[b sb];
note=sprintf(['%s s windows on a 5 s grid over no-stimulation blocks %s; %s; the stimulation window at the same offset ' ...
    'in block %d. A marked detection is one with at least one circled subnet (event_rows)'], ...
    num2str(win_len),mat2str(no_blocks),how,sb);
end
function write_stats(file,st,S)
fid=fopen(file,'w'); g=onCleanup(@()fclose(fid));
o=st.run_opts; nS=numel(st.seeds);
fprintf(fid,'# fig_ied_stim_blocks statistics\n\nRun: ied_stim_blocks, amplitude %.2g, %d networks (seeds %s), %d no-stimulation + %d stimulation blocks per network, %g s analysed per block after %g s settle (plus the model warmup); noise seed per block base.noise_seed + (seed - 42) + 1000*block.\n\n', ...
    o.amplitude,nS,mat2str(st.seeds),o.n_no,o.n_stim,o.block_len,o.settle);
if isfield(o,'fs'), fprintf(fid,'Integration rate fs = %g Hz.\n\n',o.fs);
else, fprintf(fid,'Integration rate fs = 400 Hz (ied_base_config; run made before fs was a run setting).\n\n'); end
fprintf(fid,'Example traces: no stimulation %s; stimulation %s. Window rule: %s.\n\n',st.example{1},st.example{2},st.example_rule);
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
