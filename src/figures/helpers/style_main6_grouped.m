function style_main6_grouped(fig)
% STYLE_MAIN6_GROUPED Style five native plots as four row-major panel groups.
% C contains the complementary cumulative-capacity and reconstruction axes;
% D is memory horizon. Preserve all plotted data, statistics and axis limits.
axs=findall(fig,'Type','axes');
assert(numel(axs)==5,'style_main6_grouped:Axes','Expected five scientific axes.');
is_tau=arrayfun(@(ax) contains(string(ax.XLabel.String),'slowest'),axs);
assert(nnz(is_tau)==2,'style_main6_grouped:Timescales','Expected two timescale axes.');
top=axs(is_tau); bottom=axs(~is_tau);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),top)); top=top(order);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),bottom)); bottom=bottom(order);
fig.Units='pixels'; fig.Position=[40 40 1500 1050];
set(axs,'Units','normalized','PositionConstraint','innerposition', ...
    'FontSize',14,'LineWidth',1,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
top(1).Position=[.085 .57 .36 .34];
top(2).Position=[.60 .57 .36 .34];
for k=1:2
    top(k).XLabel.String=strrep(top(k).XLabel.String,' (E and I)','');
end
for k=1:3
    bottom(k).Position=[.085+(k-1)*.315 .16 .235 .27];
end
% Sparse, interpretable ticks within the original limits.
top(2).YTick=[0 .5 1];
bottom(1).YTick=[0 8 16];
bottom(2).YTick=[0 1];
bottom(3).YTick=[0 6 12];
delete(findall(fig,'Tag','main6_panel_letter'));
anchors=[top(:);bottom([1 3])];
for k=1:4
    text(anchors(k),-.14,1.07,sprintf('(%c)','A'+k-1), ...
        'Units','normalized','Clipping','off','FontSize',14, ...
        'FontWeight','bold','Interpreter','none','Tag','main6_panel_letter');
end
lg=findall(fig,'Type','legend');
for k=1:numel(lg)
    lg(k).Units='normalized'; lg(k).FontSize=14;
    if any(contains(string(lg(k).String),'Timescale'))
        lg(k).NumColumns=3;
        lg(k).Box='off';
        drawnow;
        p=lg(k).Position; p(1)=.5-p(3)/2; p(2)=.035; lg(k).Position=p;
    else
        lg(k).Location='east';
    end
end
set(findall(fig,'-property','FontSize'),'FontSize',14);
end
