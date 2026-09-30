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
fig.Units='pixels'; fig.Position=[40 40 1100 780];
set(axs,'Units','normalized','PositionConstraint','innerposition', ...
    'FontSize',14,'LineWidth',1,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
top(1).Position=[.085 .57 .36 .34];
top(2).Position=[.60 .57 .36 .34];
for k=1:2
    top(k).XLabel.String='Longest SFA Time Constant (s)';
    top(k).XLabel.Interpreter='none';
    title(top(k),'Multiple-Timescale Adaptation','FontWeight','normal', ...
        'Interpreter','none');
end
top(1).YLabel.String='Lyapunov Exp., \lambda_1 (s^{-1})';
top(1).YLabel.Interpreter='tex';
% Match the zero-growth boundary used in the other stability figures.
st=manuscript_style();
zero_ref=findall(top(1),'Type','constantline');
zero_ref=zero_ref(arrayfun(@(h) h.Value==0,zero_ref));
assert(isscalar(zero_ref),'style_main6_grouped:ZeroReference', ...
    'Expected one zero-growth reference line.');
set(zero_ref,'Color',st.zeroline_color,'LineStyle','--', ...
    'LineWidth',st.zeroline_lw);
top(2).YLabel.String={'Leading-Vector','Squared-Norm Fraction'};
top(2).YLabel.Interpreter='none';
% State coordinates have their own palette; condition colors encode a
% different quantity in A and C-D. Markers provide a second state cue.
state_names={'SFA','STD','x'};
state_colors=[213 94 0; 0 158 115; 142 68 173]/255;
state_markers={'o','s','^'};
for k=1:numel(state_names)
    h=findall(top(2),'Type','line','DisplayName',state_names{k});
    assert(numel(h)==1,'style_main6_grouped:States','Expected one line per state block.');
    set(h,'Color',state_colors(k,:),'MarkerFaceColor',state_colors(k,:), ...
        'Marker',state_markers{k});
    if k==3, h.DisplayName='Dendritic State'; end
end
for k=1:3
    bottom(k).Position=[.085+(k-1)*.315 .16 .235 .27];
end
bottom(2).YLabel.String={'Input Reconstruction','Score (R^2)'};
bottom(2).YLabel.Interpreter='tex';
bottom(3).XTickLabel={'No Adaptation','1TS','MTS'};
bottom(3).XLabel.String='Adaptation Condition';
bottom(3).XLabel.Interpreter='none';
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
        'FontWeight','normal','Interpreter','none','Tag','main6_panel_letter');
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
        lg(k).String={'SFA','STD','Dendritic State'};
        lg(k).Location='east';
    end
end
set(findall(fig,'-property','FontSize'),'FontSize',14);
end
