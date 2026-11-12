function style_main6_grouped(fig,tau_names)
% STYLE_MAIN6_GROUPED Style the native plots as four row-major panel groups.
% (A) timescale LLE, (B) leading-vector fractions, (C) the complementary
% cumulative-capacity and reconstruction axes, (D) memory horizon. Preserve
% all plotted data, statistics and axis limits.
%
%   STYLE_MAIN6_GROUPED(fig, tau_names)
%   STYLE_MAIN6_GROUPED(fig)        % restyle a saved FIG: names read from it
%
% CONDITIONS (FD-006). The tau sweep holds the fully-adapted conditions
% only, so the top row has one LLE panel and one leading-vector panel per
% condition the sweep drew (tau_names, left to right): the chain alone in
% the paper's three-condition set, chain and product in the review set.
% (A) is the LLE panels side by side in the left half, (B) the
% leading-vector panels in the right half; with one condition each is one
% panel, at the geometry of before. Every top panel is titled by the NAME of
% its condition (srnn_condition_titles), so the product MTS is never called
% "Multiple-Timescale Adaptation". The memory-horizon ticks (D) are the
% conditions fig_memory_capacity drew, read from its own tick labels
% (condition titles) and relabelled by name (grouped_condition_label
% 'abbrev': No Adaptation / 1TS / MTS / prodMTS); a literal three-label list
% here cycled silently with four (inventory R5). Without tau_names the top
% panels' names are read from their titles, which a single-panel sweep
% figure does not carry: then the names must be passed.
axs=findall(fig,'Type','axes');
is_tau=arrayfun(@(ax) contains(strjoin(cellstr(string(ax.XLabel.String)),' '),{'slowest','Longest SFA'}),axs);
top=axs(is_tau); bottom=axs(~is_tau);
assert(numel(bottom)==3,'style_main6_grouped:Axes','Expected three memory-capacity axes, found %d.',numel(bottom));
% Leading-vector panels carry one line per state block, DisplayName 'SFA'.
is_block=arrayfun(@(ax) ~isempty(findall(ax,'Type','line','DisplayName','SFA')),top);
lle=top(~is_block); blk=top(is_block);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),lle)); lle=lle(order);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),blk)); blk=blk(order);
[~,order]=sort(arrayfun(@(ax) ax.Position(1),bottom)); bottom=bottom(order);
n=numel(lle);
assert(n>=1 && numel(blk)==n,'style_main6_grouped:Timescales', ...
    'Expected one LLE and one leading-vector panel per tau condition; found %d and %d.',n,numel(blk));
if nargin<2 || isempty(tau_names)
    tau_names=grouped_condition_label(arrayfun(@(ax) {ax.Title.String},lle),'name_of_title');
end
tau_names=cellstr(tau_names);
assert(numel(tau_names)==n,'style_main6_grouped:TauNames', ...
    '%d tau condition names for %d LLE panels.',numel(tau_names),n);
fig.Units='pixels'; fig.Position=[40 40 1100 780];
set(axs,'Units','normalized','PositionConstraint','innerposition', ...
    'FontSize',14,'LineWidth',1,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
% Each half holds n panels; with n = 1 these are the two panels of before.
gap=.06; w=(.36-(n-1)*gap)/n;
for k=1:n
    lle(k).Position=[.085+(k-1)*(w+gap) .57 w .34];
    blk(k).Position=[.60+(k-1)*(w+gap) .57 w .34];
end
titles=grouped_condition_label(tau_names,'title');
st=manuscript_style();
for k=1:n
    for ax=[lle(k) blk(k)]
        ax.XLabel.String='Longest SFA Time Constant (s)';
        ax.XLabel.Interpreter='none';
        title(ax,titles{k},'FontWeight','normal','Interpreter','tex');
    end
    if k>1
        lle(k).YLabel.String=''; blk(k).YLabel.String='';
    end
    % Match the zero-growth boundary used in the other stability figures.
    zero_ref=findall(lle(k),'Type','constantline');
    zero_ref=zero_ref(arrayfun(@(h) h.Value==0,zero_ref));
    assert(isscalar(zero_ref),'style_main6_grouped:ZeroReference', ...
        'Expected one zero-growth reference line.');
    set(zero_ref,'Color',st.zeroline_color,'LineStyle','--', ...
        'LineWidth',st.zeroline_lw);
end
lle(1).YLabel.String='Lyapunov Exp., \lambda_1 (s^{-1})';
lle(1).YLabel.Interpreter='tex';
blk(1).YLabel.String={'Leading-Vector','Squared-Norm Fraction'};
blk(1).YLabel.Interpreter='none';
% State coordinates have their own palette; condition colors encode a
% different quantity in A and C-D. Markers provide a second state cue.
state_names={'SFA','STD','x'};
state_colors=[213 94 0; 0 158 115; 142 68 173]/255;
state_markers={'o','s','^'};
for b=1:n
    for k=1:numel(state_names)
        h=findall(blk(b),'Type','line','DisplayName',state_names{k});
        if k==3 && isempty(h), h=findall(blk(b),'Type','line','DisplayName','Dendritic State'); end
        assert(isscalar(h),'style_main6_grouped:States','Expected one line per state block.');
        set(h,'Color',state_colors(k,:),'MarkerFaceColor',state_colors(k,:), ...
            'Marker',state_markers{k});
        if k==3, h.DisplayName='Dendritic State'; end
    end
end
for k=1:3
    bottom(k).Position=[.085+(k-1)*.315 .16 .235 .27];
end
bottom(2).YLabel.String={'Input Reconstruction','Score (R^2)'};
bottom(2).YLabel.Interpreter='tex';
% The horizon panel's conditions, from its own tick labels (titles, by name).
% A FIG restyled a second time already carries the abbreviations; it keeps them.
try
    mc_names=grouped_condition_label(cellstr(bottom(3).XTickLabel),'name_of_title');
    bottom(3).XTickLabel=grouped_condition_label(mc_names,'abbrev');
catch ME
    if ~strcmp(ME.identifier,'grouped_condition_label:UnknownTitle'), rethrow(ME); end
end
bottom(3).XTickLabelRotation=0;
% Three ticks fit on one line (as before); with more, 'No Adaptation' collides
% with its neighbour in the narrow panel, so multi-word labels break.
if numel(bottom(3).XTick)>3
    bottom(3).TickLabelInterpreter='tex';
    bottom(3).XTickLabel=strrep(cellstr(bottom(3).XTickLabel),' ','\newline');
end
bottom(3).XLabel.String='Adaptation Condition';
bottom(3).XLabel.Interpreter='none';
% Sparse, interpretable ticks within the original limits.
for b=1:n, blk(b).YTick=[0 .5 1]; end
bottom(1).YTick=[0 8 16];
bottom(2).YTick=[0 1];
bottom(3).YTick=[0 6 12];
delete(findall(fig,'Tag','main6_panel_letter'));
anchors=[lle(1);blk(1);bottom([1 3])];
% The letter sits .14 of a full-width top panel to the left of its panel.
x_letter=[-.14*.36/w -.14*.36/w -.14 -.14];
for k=1:4
    text(anchors(k),x_letter(k),1.07,sprintf('(%c)','A'+k-1), ...
        'Units','normalized','Clipping','off','FontSize',14, ...
        'FontWeight','normal','Interpreter','none','Tag','main6_panel_letter');
end
lg=findall(fig,'Type','legend');
for k=1:numel(lg)
    lg(k).Units='normalized'; lg(k).FontSize=14;
    if ~all(ismember(lg(k).String,{'SFA','STD','x','Dendritic State'}))
        % The memory-capacity condition legend: one row of up to three keys.
        lg(k).NumColumns=min(3,numel(lg(k).String));
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
