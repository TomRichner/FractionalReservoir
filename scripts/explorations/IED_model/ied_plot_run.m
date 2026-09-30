function files = ied_plot_run(data,det,cfg,outdir)
% Diagnostic views deliberately show individual neurons and recruitment.
t=data.t; x=double(data.x); r=double(data.r); g=data.groups;
palette=lines(cfg.groups); fig=figure('Visible','off','Color','w','Position',[40 40 1500 1100]);
tl=tiledlayout(fig,6,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); hold(ax,'on');
plot(ax,t,mean(x(data.E,:),1),'Color',[.85 .25 .1]);
plot(ax,t,mean(x(data.I,:),1),'Color',[.1 .4 .7]);
ylabel(ax,'Mean x (a.u.)'); legend(ax,{'E','I'},'Location','eastoutside');
title(ax,sprintf('Run %02d | %s | %.1f Events/min | Recruitment %.2f | lambda1 %.3f', ...
    cfg.id,cfg.topology,det.metrics.events_per_min,det.metrics.median_recruitment,det.metrics.lambda1));
ax=nexttile(tl); hold(ax,'on');
for j=1:cfg.groups, plot(ax,t,mean(x(g==j,:),1),'Color',palette(j,:)); end
ylabel(ax,'Group Mean x');
legend(ax,compose('Group %d',1:cfg.groups),'Location','eastoutside');
ax=nexttile(tl); [~,order]=sortrows([g (1:numel(g))'],[1 2]);
imagesc(ax,t,1:numel(g),x(order,:)); axis(ax,'xy');
clim(ax,[-.3 .8]); colorbar(ax); ylabel(ax,'Neuron by Group');
ax=nexttile(tl); hold(ax,'on');
for j=1:cfg.groups, plot(ax,t,det.group_recruitment(j,:),'Color',palette(j,:)); end
plot(ax,t,mean(det.active,1),'k','LineWidth',1.4);
ylabel(ax,'Recruitment'); ylim(ax,[0 1]);
ax=nexttile(tl); hold(ax,'on');
plot(ax,t,mean(r,1),'k'); plot(ax,t,data.sfa,'Color',[.8 .3 .1]);
plot(ax,t,data.resources,'Color',[0 .55 .4]);
ylabel(ax,'Rate / Feedback'); legend(ax,{'Mean Rate','SFA Feedback','Resources'},'Location','eastoutside');
ax=nexttile(tl); hold(ax,'on');
if ~isempty(det.growth.t)
    local=plot(ax,det.growth.t,det.growth.local,'Color',[.65 .65 .65]);
    accumulated=plot(ax,det.growth.t,det.growth.accumulated,'k','LineWidth',1.5);
    yline(ax,0,'--','Color',[0 .6 0],'HandleVisibility','off');
    yyaxis(ax,'right'); expansion=plot(ax,det.growth.t,det.growth.expansion,'Color',[.5 .1 .6]);
    ylabel(ax,sprintf('Positive Top-%d Growth (bits/s)',cfg.K));
    ax.YColor=[.5 .1 .6];
    yyaxis(ax,'left');
    ax.YColor=[.15 .15 .15];
    legend(ax,[local accumulated expansion],{'Local \lambda_1','Accumulated \lambda_1','Positive Growth Sum'},'Location','eastoutside');
end
ylabel(ax,'Growth (s^{-1})'); xlabel(ax,'Time (s)');
axs=findall(fig,'Type','axes'); linkaxes(axs,'x'); xlim(axs,[t(1) t(end)]);
for k=1:numel(axs)
    for e=1:numel(det.events)
        xline(axs(k),det.events(e).time,':','Color',[.75 .3 .75],'HandleVisibility','off');
    end
end
set(axs,'FontSize',11);
files={fullfile(outdir,'overview.png'),fullfile(outdir,'event_zoom.png')};
exportgraphics(fig,files{1},'Resolution',140); close(fig);
if isempty(det.events)
    center=t(round(numel(t)/2));
else
    [~,e]=max([det.events.recruitment]); center=det.events(e).time;
end
fig=figure('Visible','off','Color','w','Position',[40 40 1500 850]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
window=t>=max(t(1),center-2) & t<=min(t(end),center+2);
ax=nexttile(tl); imagesc(ax,t(window),1:numel(g),x(order,window));
axis(ax,'xy'); clim(ax,[-.3 .8]); colorbar(ax); ylabel(ax,'Neuron by Group');
title(ax,sprintf('Run %02d | Event Detail around %.3f s',cfg.id,center));
ax=nexttile(tl); hold(ax,'on');
for j=1:cfg.groups, plot(ax,t(window),mean(x(g==j,window),1),'Color',palette(j,:)); end
ylabel(ax,'Group Mean x (a.u.)');
ax=nexttile(tl); hold(ax,'on');
strength=max(x(:,abs(t-center)<.15)-det.baseline,[],2);
[~,ix]=sort(strength,'descend'); pick=ix(1:min(8,numel(ix)));
for j=1:numel(pick), plot(ax,t(window),x(pick(j),window)+.4*(j-1)); end
ylabel(ax,'Selected x, Offset 0.4'); xlabel(ax,'Time (s)');
exportgraphics(fig,files{2},'Resolution',140); close(fig);
end
