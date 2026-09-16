function out=fig_main3_stability(run_dir)
% Saved eigen-density and switching examples. Never compute new trajectories.
f=fullfile(run_dir,'local_lyapunov','local_lyapunov_data.mat'); D=load(f);
a=fig_eig_heatmap('run_dir',run_dir,'save',false,'visible',false,'density_scale','loglog');
fig=figure('Visible','off','Color','w','Position',[40 40 1450 1150]);
old=findall(a.figs,'Type','axes');
% Tile order is the condition order; findall returns the reverse creation order.
[~,order]=sort(arrayfun(@(x)x.Layout.Tile,old)); old=old(order);
st=manuscript_style(); titles={'No adaptation','Single-timescale adaptation','Multiple-timescale adaptation'};
selected=cellfun(@(v)v(round(linspace(1,numel(v),min(12,numel(v))))), ...
    D.results(1).type_indices,'UniformOutput',false);
for c=1:3
    r=D.results(c); tk=r.topk; col=st.condition_color(r.name);
    ax=copyobj(old(c),fig); ax.Units='normalized'; ax.Position=[.105+(c-1)*.30 .745 .235 .18];
    ax.Tag=sprintf('stability_r1_c%d',c);
    texts=findall(ax,'Type','text');
    for k=1:numel(texts)
        if startsWith(string(texts(k).String),'\lambda_1 ='), delete(texts(k)); end
    end
    title(ax,titles{c},'FontSize',16,'FontWeight','normal','Color',col);
    ax.Title.Units='normalized'; ax.Title.Position=[.5 1.14 0];
    if c==1, letter(ax,'A'); end
    if c==3
        cb=colorbar(ax); cb.Units='normalized'; cb.Position=[.955 .755 .012 .16];
        cb.Label.String='log log density'; cb.Ticks=[0 .25 .5]; cb.Box='off';
        ax.Position=[.705 .745 .235 .18];
    end
    ax=newaxis(2,c,.535,.125); offset=0;
    assert(isequal(r.type_indices,D.results(1).type_indices),'Neuron indexing differs between conditions.');
    for q=1:numel(selected)
        ids=selected{q};
        if q==1, cmap=excitatory_colormap(numel(ids)); else, cmap=inhibitory_colormap(numel(ids)); end
        for k=1:numel(ids)
            plot(ax,r.t_ex,r.u_ex(ids(k),:),'Color',cmap(k,:),'LineWidth',.65);
        end
        offset=offset+numel(ids);
    end
    ax.XAxis.Visible='off'; if c==1, ylabel(ax,sprintf('Input (%d neurons)',offset)); letter(ax,'B'); end
    ax=newaxis(3,c,.315,.14); t=tk.t_lya(:); finite=tk.finite_LE_spectrum_t(:,1);
    % Saved t_lya denotes segment starts. Preserve the archived time convention.
    plot(ax,[0 60],[0 0],':','Color',[.65 .65 .65]);
    plot(ax,t,finite,'Color',col,'LineWidth',1.7,'Tag','accumulated_lambda');
    if c==1, ylim(ax,[-.5 4]); else, ylim(ax,[-.5 .5]); end
    ax.XAxis.Visible='off';
    last=find(isfinite(finite),1,'last'); value=finite(last); yl=ylim(ax);
    if c<=2, delta=-.20; else, delta=.16; end
    ypos=min(yl(2)-.1*diff(yl),max(yl(1)+.12*diff(yl),value+delta*diff(yl)));
    text(ax,55,ypos,sprintf('\\lambda_1 = %+.3f',value),'HorizontalAlignment','right','FontSize',14,'Color',col);
    if c==1, ylabel(ax,'Accumulated \lambda_1 (s^{-1})'); letter(ax,'C'); end
    ax=newaxis(4,c,.105,.135); h=sum(max(tk.local_LE_spectrum_t,0),2)/log(2);
    plot(ax,t,h,'Color',col,'LineWidth',1.1); ax.TickDir='out'; xlabel(ax,'Time (s)');
    if c==1, ylabel(ax,'Local positive sum (bit/s)'); letter(ax,'D'); end
end
close(a.figs);
sep=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'Visible','off','Tag','column_dividers'); hold(sep,'on');
for xpos=[.3775 .6775], plot(sep,[xpos xpos],[.08 .95],'Color',[.88 .88 .88],'LineWidth',2.5); end
text(sep,.51,.018,sprintf('Saved switching example: K = %d; accumulation [%g, %g] s', ...
    D.settings.K,D.settings.lya_T_interval),'HorizontalAlignment','center','FontSize',14);
uistack(sep,'bottom'); set(findall(fig,'-property','FontSize'),'FontSize',14);
for c=1:3, ax=findobj(fig,'Tag',sprintf('stability_r1_c%d',c)); ax.Title.FontSize=16; end
notes={sprintf('Current switching data: K=%d; accumulation [%g,%g] s; saved segment-start timestamps retained. No early estimates fabricated. Requested K=100 and accumulation [1,60] s are prepared in mu7revised but not run.',D.settings.K,D.settings.lya_T_interval), ...
    ['Input shows 24 fixed actual neurons (12 evenly spaced E and 12 I), identical indices across conditions. ' ...
    'Local positive sum = sum(max(saved local QR rates,0),2)/log(2). This is a provisional local h_KS proxy, not established invariant KS entropy; broader review deferred. Eigenvalue and switching stages are separate examples.']};
out=struct('figs',fig,'files',{{}},'source',{{a.source,f}},'notes',{notes});
    function ax=newaxis(row,col,y,height)
        ax=axes(fig,'Position',[.105+(col-1)*.30 y .235 height], ...
            'Tag',sprintf('stability_r%d_c%d',row,col),'FontSize',14,'LineWidth',1,'Box','off');
        hold(ax,'on'); xlim(ax,[0 60]);
    end
    function letter(ax,ch)
        text(ax,-.34,1.05,['(' ch ')'],'Units','normalized','FontSize',14, ...
            'VerticalAlignment','top','Clipping','off');
    end
end
