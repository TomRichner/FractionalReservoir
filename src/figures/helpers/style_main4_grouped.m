function style_main4_grouped(fig)
% STYLE_MAIN4_GROUPED Presentation-only styling of the native grouped figure.
% Also accepts an existing saved FIG, without recomputing data or statistics.
st = manuscript_style();
set(findall(fig,'-property','FontSize'),'FontSize',14);
axs = findall(fig,'Type','axes');
assert(numel(axs)==9,'style_main4_grouped:Axes','Expected six sensitivity and three rate panels.');
for k=1:numel(axs)
    ax=axs(k);
    if ax.Position(2)>.3
        ax.YTick=[-1 0 1];
        if ax.Position(2)>.65
            p=ax.Position; p(2)=.69; ax.Position=p;
        end
    else
        % The plain line is the saved rate-bin median; ConstantLine objects
        % are the separate zero references and retain their original styling.
        set(findall(ax,'Type','line'),'Color',[.5 .5 .5],'LineWidth',2.5);
        names={'no_adaptation','sfa1_std1','sfa3_std2'};
        for j=1:numel(names)
            if strcmp(ax.Title.String,st.condition_title(names{j}))
                ax.Title.Color=st.condition_color(names{j});
            end
        end
    end
end
lg=findall(fig,'Type','legend');
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
annotation(fig,'textbox',[.005 .955 .04 .035],'String','(A)', ...
    'LineStyle','none','FontSize',14,'Margin',0,'Tag','main4_block_A');
annotation(fig,'textbox',[.005 .31 .04 .035],'String','(B)', ...
    'LineStyle','none','FontSize',14,'Margin',0,'Tag','main4_block_B');

end
