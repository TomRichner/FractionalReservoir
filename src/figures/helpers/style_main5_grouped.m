function style_main5_grouped(fig)
% STYLE_MAIN5_GROUPED Arrange native saved margin and gain graphics in one row.
% Presentation only: preserve all curve data, styles, horizons and gain limits.
st=manuscript_style();
axs=findall(fig,'Type','axes');
assert(numel(axs)==4,'style_main5_grouped:Axes','Expected one margin and three gain axes.');
is_margin=arrayfun(@(ax) any(contains(string(ax.Title.String),'non-normal')),axs);
assert(nnz(is_margin)==1,'style_main5_grouped:Margin','Expected one margin panel.');
margin=axs(is_margin); gains=axs(~is_margin);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),gains)); gains=gains(order);
fig.Units='pixels'; fig.Position=[40 60 1295 406];
set(findall(fig,'-property','FontSize'),'FontSize',14);
set(axs,'LineWidth',1,'Units','normalized','PositionConstraint','innerposition', ...
    'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
ordered=[margin;gains(:)];
for k=1:4
    ordered(k).Position=[.055+(k-1)*.242 .35 .19 .48];
    ordered(k).Tag=sprintf('main5_panel%d',k);
end
margin.YLim=[0 100]; margin.YTick=[0 25 50 75];
margin.Title.String={'non-normal margin per state','(5-95%, IQR, median)'};
names={'no_adaptation','sfa1_std1','sfa3_std2'};
labels={'No adaptation','Single-timescale','Multiple-timescale'};
colored=cell(1,3);
for k=1:3
    c=st.condition_color(names{k});
    colored{k}=sprintf('\\color[rgb]{%.6f,%.6f,%.6f}%s',c,labels{k});
    gains(k).Title.Color=c;
end
margin.TickLabelInterpreter='tex'; margin.XTick=1:3;
margin.XTickLabel=colored; margin.XTickLabelRotation=30;
gains(1).Title.String='No Adaptation';
gains(2).Title.String={'Single-Timescale','Adaptation'};
gains(3).Title.String={'Multiple-Timescale','Adaptation'};
lg=findall(fig,'Type','legend');
set(lg,'Units','normalized','NumColumns',3,'Box','off','FontSize',14);
drawnow;
for k=1:numel(lg)
    center=ordered(3).Position(1)+ordered(3).Position(3)/2;
    ordered(3).XLabel.Units='normalized'; label_extent=ordered(3).XLabel.Extent;
    label_bottom=ordered(3).Position(2)+label_extent(2)*ordered(3).Position(4);
    p=lg(k).Position; p(1)=center-p(3)/2; p(2)=label_bottom-.025-p(4); lg(k).Position=p;
end
delete(findall(fig,'Tag','main5_block_A')); delete(findall(fig,'Tag','main5_block_B'));
for k=1:2
    annotation(fig,'textbox',[ordered(k).Position(1)-.03 .88 .025 .06], ...
        'String',sprintf('(%c)',char('A'+k-1)),'LineStyle','none','FontSize',14, ...
        'Margin',0,'Tag',sprintf('main5_block_%c',char('A'+k-1)));
end

end
