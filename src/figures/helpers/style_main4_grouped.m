function style_main4_grouped(fig)
% STYLE_MAIN4_GROUPED Presentation-only styling of the native grouped figure.
% Also accepts an existing saved FIG, without recomputing data or statistics.
st = manuscript_style();
% Fixed70% width/85% height of the1300-by-950 creation baseline.
fig.Units='pixels'; fig.Position=[40 40 910 807.5];
set(findall(fig,'-property','FontSize'),'FontSize',14);
axs = findall(fig,'Type','axes');
set(axs,'LineWidth',1,'XTickLabelRotation',0);
assert(numel(axs)==9,'style_main4_grouped:Axes','Expected six sensitivity and three rate panels.');
for k=1:numel(axs)
    ax=axs(k);
    if ax.Position(2)>.3
        ax.YTick=[-1 0 1];
        if ax.Position(2)>.65
            p=ax.Position; p(2)=.69; ax.Position=p;
        end
    else
        ax.XTick=[0 .5 1];
        ax.Title.String=char(strrep(strjoin(string(ax.Title.String),' '),newline,' '));
        % The plain line is the saved rate-bin median; ConstantLine objects
        % are the separate zero references and retain their original styling.
        set(findall(ax,'Type','line'),'Color',[.5 .5 .5],'LineWidth',2.5);
        names={'no_adaptation','sfa1_std1','sfa3_std2'};
        for j=1:numel(names)
            if strcmp(strrep(ax.Title.String,newline,' '),st.condition_title(names{j}))
                ax.Title.Color=st.condition_color(names{j});
            end
        end
    end
end
lg=findall(fig,'Type','legend');
names={'no_adaptation','sfa1_std1','sfa3_std2'};
colored=cell(1,3);
for j=1:3
    colored{j}=sprintf('\\color[rgb]{%.6f,%.6f,%.6f}%s',st.condition_color(names{j}),st.condition_title(names{j}));
end
for j=1:numel(lg)
    % Include the lambda_1 = 0 boundary in the shared legend. Reuse a native
    % reference line so its green dashed key matches the sensitivity panels.
    leg_ax=lg(j).Axes;
    refs=findall(leg_ax,'Type','constantline');
    zero_ref=refs(arrayfun(@(h)h.Value==0 && strcmp(h.LineStyle,'--'),refs));
    assert(~isempty(zero_ref),'style_main4_grouped:ZeroReference','Expected a zero-growth reference line.');
    zero_ref=zero_ref(1);
    handles=lg(j).PlotChildren;
    handles=handles(~arrayfun(@(h)isa(h,'matlab.graphics.chart.decoration.ConstantLine'),handles));
    zero_ref.HandleVisibility='on';
    zero_ref.DisplayName='Zero Growth';
    lg(j)=legend(leg_ax,[handles(:);zero_ref],[colored {'Zero Growth'}],'Interpreter','tex');
end
set(lg,'Orientation','vertical','NumColumns',1,'Box','off','Units','normalized');
drawnow;
% A shared vertical legend above the upper-right panel, clear of the curves.
for k=1:numel(lg)
    p=lg(k).Position; p(1)=.97-p(3); p(2)=.98-p(4);
    lg(k).Position=p;
    % Reserve space for all four keys while keeping upper-row x labels in place.
    upper=axs(arrayfun(@(ax)ax.Position(2)>.65,axs));
    for j=1:numel(upper)
        q=upper(j).Position;
        q(4)=min(q(4),p(2)-.015-q(2));
        upper(j).Position=q;
    end
end
cb=findall(fig,'Type','colorbar');
set(cb,'Box','off','FontSize',14);
for k=1:numel(cb)
    cb(k).Label.String='Weights E:I'; cb(k).Label.FontSize=14;
    % The native colorbar leaves room for ticks but its label falls outside
    % the grouped canvas. Reserve a right margin for the full label.
    p=cb(k).Position; p(1)=.90; cb(k).Position=p;
end
rate_axes=axs(arrayfun(@(ax)ax.Position(2)<.3,axs));
[~,last]=max(arrayfun(@(ax)ax.Position(1),rate_axes));
ax=rate_axes(last); p=ax.Position; ax.Title.Units='normalized';
title_x=ax.Title.Position(1)*p(3);
p(3)=min(p(3),.88-p(1)); ax.Position=p;
% Keep the long condition title centered at its original location.
ax.Title.Position(1)=title_x/p(3);
% Block A spans both sensitivity rows; block B is the rate row.
delete(findall(fig,'Tag','main4_block_A')); delete(findall(fig,'Tag','main4_block_B'));
upper=axs(arrayfun(@(ax)ax.Position(2)>.65,axs));
[~,first]=min(arrayfun(@(ax)ax.Position(1),upper)); p=upper(first).Position;
annotation(fig,'textbox',[p(1)-.015 p(2)+p(4)+.002 .045 .035],'String','(A)', ...
    'LineStyle','none','FontSize',14,'Margin',0,'VerticalAlignment','bottom','Tag','main4_block_A');
annotation(fig,'textbox',[.005 .31 .04 .035],'String','(B)', ...
    'LineStyle','none','FontSize',14,'Margin',0,'Tag','main4_block_B');

end
