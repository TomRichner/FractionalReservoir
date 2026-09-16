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
for j=1:numel(lg), lg(j).String=colored; lg(j).Interpreter='tex'; end
set(lg,'Orientation','vertical','NumColumns',1,'Box','off','Units','normalized');
drawnow;
% A shared vertical legend above the upper-right panel, clear of the curves.
for k=1:numel(lg)
    p=lg(k).Position; p(1)=.97-p(3); p(2)=.98-p(4);
    lg(k).Position=p;
end
cb=findall(fig,'Type','colorbar');
set(cb,'Box','off','FontSize',14);
for k=1:numel(cb), cb(k).Label.FontSize=14; end
% Block A spans both sensitivity rows; block B is the rate row.
delete(findall(fig,'Tag','main4_block_A')); delete(findall(fig,'Tag','main4_block_B'));
upper=axs(arrayfun(@(ax)ax.Position(2)>.65,axs));
[~,first]=min(arrayfun(@(ax)ax.Position(1),upper)); p=upper(first).Position;
annotation(fig,'textbox',[p(1)-.015 p(2)+p(4)+.002 .045 .035],'String','(A)', ...
    'LineStyle','none','FontSize',14,'Margin',0,'VerticalAlignment','bottom','Tag','main4_block_A');
annotation(fig,'textbox',[.005 .31 .04 .035],'String','(B)', ...
    'LineStyle','none','FontSize',14,'Margin',0,'Tag','main4_block_B');

end
