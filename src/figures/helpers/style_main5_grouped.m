function style_main5_grouped(fig,margin_names,gain_names)
% STYLE_MAIN5_GROUPED Arrange native saved margin and gain graphics in one row.
% Presentation only: preserve all gain curve data, styles, horizons and limits.
%
%   STYLE_MAIN5_GROUPED(fig, margin_names, gain_names)
%   STYLE_MAIN5_GROUPED(fig)          % restyle a saved FIG: names read from it
%
% CONDITIONS (FD-006). One margin panel (one box per condition, left to
% right in margin_names' order) and one gain panel per condition (gain_names,
% left to right). Tick labels, tick colours and gain-panel titles are looked
% up from those NAMES (grouped_condition_label: 'tick' for the margin ticks,
% 'two_line' for the gain titles; colours from manuscript_style). Without
% names the gain names are read back from the gain panels' titles and the
% margin is assumed to hold the same conditions. Layout: 1 + N panel slots;
% with N = 3 the 1295 x 406 canvas and the slot geometry of before, with
% more the canvas widens by a slot per condition and the slots keep their
% pixel size.
st=manuscript_style();
axs=findall(fig,'Type','axes');
is_margin=arrayfun(@(ax) any(contains(string(ax.Title.String),'non-normal')),axs);
assert(nnz(is_margin)==1,'style_main5_grouped:Margin','Expected one margin panel.');
margin=axs(is_margin); gains=axs(~is_margin);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),gains)); gains=gains(order);
n=numel(gains);
if nargin<3 || isempty(gain_names)
    gain_names=grouped_condition_label(arrayfun(@(ax) {ax.Title.String},gains),'name_of_title');
end
if nargin<2 || isempty(margin_names), margin_names=gain_names; end
assert(numel(gain_names)==n,'style_main5_grouped:Axes', ...
    'Expected one gain panel per condition (%d names, %d gain panels).',numel(gain_names),n);
assert(numel(margin.XTick)==numel(margin_names) || isempty(margin.XTick), ...
    'style_main5_grouped:MarginTicks','The margin panel has %d ticks for %d conditions.', ...
    numel(margin.XTick),numel(margin_names));
% Slot geometry of the four-panel row, in units of its 1295-px width; a
% wider row keeps each slot's pixel size (u scales normalized x).
u=4/(n+1);
fig.Units='pixels'; fig.Position=[40 60 1295/u 406];
set(findall(fig,'-property','FontSize'),'FontSize',14);
set(axs,'LineWidth',1,'Units','normalized','PositionConstraint','innerposition', ...
    'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
ordered=[margin;gains(:)];
for k=1:n+1
    ordered(k).Position=[(.055+(k-1)*.242)*u .35 .19*u .48];
    ordered(k).Tag=sprintf('main5_panel%d',k);
end
% Reserve room for the descriptive two-line margin label at export.
margin.Position([1 3])=[.08 .165]*u;
margin.YLim=[0 100]; margin.YTick=[0 25 50 75];
margin.YLabel.String={'Non-Normal Amplification','Margin (s^{-1})'};
margin.Title.String={'non-normal margin per state','(5-95%, IQR, median)'};
labels=grouped_condition_label(margin_names,'tick');
colored=cell(1,numel(margin_names));
for k=1:numel(margin_names)
    colored{k}=sprintf('\\color[rgb]{%.6f,%.6f,%.6f}%s',st.condition_color(margin_names{k}),labels{k});
end
margin.TickLabelInterpreter='tex'; margin.XTick=1:numel(margin_names);
margin.XTickLabel=colored; margin.XTickLabelRotation=30;
titles=grouped_condition_label(gain_names,'two_line');
for k=1:n
    gains(k).Title.String=titles{k};
    gains(k).Title.Color=st.condition_color(gain_names{k});
end
lg=findall(fig,'Type','legend');
set(lg,'Units','normalized','NumColumns',3,'Box','off','FontSize',14);
drawnow;
for k=1:numel(lg)
    % Centred under the gain row; the label height is the first gain panel's.
    center=(gains(1).Position(1)+gains(end).Position(1)+gains(end).Position(3))/2;
    gains(1).XLabel.Units='normalized'; label_extent=gains(1).XLabel.Extent;
    label_bottom=gains(1).Position(2)+label_extent(2)*gains(1).Position(4);
    p=lg(k).Position; p(1)=center-p(3)/2; p(2)=label_bottom-.025-p(4); lg(k).Position=p;
end
delete(findall(fig,'Tag','main5_block_A')); delete(findall(fig,'Tag','main5_block_B'));
for k=1:2
    annotation(fig,'textbox',[ordered(k).Position(1)-.03*u .88 .025*u .06], ...
        'String',sprintf('(%c)',char('A'+k-1)),'LineStyle','none','FontSize',14, ...
        'Margin',0,'Tag',sprintf('main5_block_%c',char('A'+k-1)));
end

end
