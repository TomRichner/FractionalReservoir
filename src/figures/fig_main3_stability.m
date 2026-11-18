function out=fig_main3_stability(run_dir,conditions,opts)
% Saved eigen-density and switching examples. Never compute new trajectories.
%
% CONDITIONS (FD-006): one column per condition, in the order of
% `conditions` (a cell of names; empty = every condition the local_lyapunov
% stage saved, in its order), resolved by select_conditions against that
% stage's results. The SAME names, in the same order, are passed to
% fig_eig_heatmap as its selection, so row (A) and rows (B)-(D) of a column
% are always one condition even though they come from two stages; each
% heat-map tile is then taken by condition (its tile index in the selected
% order, checked against its title), never by position in a 1 x 3 grid
% (fig_eig_heatmap is 2 x 2 at four conditions). Titles and colours from the
% condition NAME. N columns keep the three-column pitch and the canvas
% widens with them (out.width_scale, 1 for three columns, whose geometry is
% unchanged). Row (C)'s y range is [-.5 4] for no_adaptation (by name) and
% [-.5 .5] otherwise, and its lambda_1 label sits above the accumulated
% curve when that curve ends within a tenth of the axis range of zero, below
% it otherwise (by value; for the manuscript's three columns that is the
% placement they had).
%
% SHOW_LAMBDA_ROW (default false since 2026-11-17, FD-031 item 7: the paper's
% layout, approved by Tom). false removes row (C) of the layout above, the
% lambda_1 traces, and reorders the remaining rows top to bottom as (A)
% switching input, (B) expansion, (C) eigenvalue density heat map, on a
% shorter canvas of the same width. Condition titles then sit over the input
% row, and each column's accumulated lambda_1 (the value the removed row
% printed) is kept as a small label in the top right corner of its expansion
% panel. true draws the old four-row layout described above, kept reachable.
%
% K, the number of QR directions the local_lyapunov stage tracked (top-K), is
% read from the saved settings and written into the notes (provenance): the
% expansion row sums the positive local rates of those K directions only, so
% when more than K directions expand it is a lower bound.
arguments
    run_dir char
    conditions = {}
    opts.show_lambda_row (1,1) logical = false
end
show_lambda=opts.show_lambda_row;
% Row geometry, [y height] in normalized figure units, and the canvas height.
if show_lambda
    fig_h=1150; g.heat=[.745 .18]; g.input=[.535 .125]; g.lambda=[.315 .14]; g.exp=[.105 .135];
    g.cbar=[.755 .16]; g.div=[.08 .95];
else
    fig_h=900; g.input=[.755 .155]; g.exp=[.45 .2]; g.heat=[.065 .25];
    g.cbar=[.08 .22]; g.div=[.04 .965];
end
f=fullfile(run_dir,'local_lyapunov','local_lyapunov_data.mat'); D=load(f);
[names,cidx]=select_conditions({D.results.name},conditions);
D.results=D.results(cidx);
n_col=numel(names); width_scale=1+(n_col-3)*.30; s=1/width_scale;
T=D.settings.T;
a=fig_eig_heatmap('run_dir',run_dir,'save',false,'visible',false,'density_scale','loglog', ...
    'conditions',names);
fig=figure('Visible','off','Color','w','Position',[40 40 1450*width_scale fig_h]);
old=findall(a.figs,'Type','axes');
% Tile k is the k-th selected condition (fig_eig_heatmap draws them in the
% selection's order, row-major); findall returns the reverse creation order.
[~,order]=sort(arrayfun(@(x)x.Layout.Tile,old)); old=old(order);
assert(numel(old)==n_col,'fig_main3_stability:Tiles', ...
    'fig_eig_heatmap drew %d tiles for %d conditions.',numel(old),n_col);
st=manuscript_style(); titles=grouped_condition_label(names,'sentence');
heat_titles=titles_by_name(names);
selected=cellfun(@(v)v(round(linspace(1,numel(v),min(12,numel(v))))), ...
    D.results(1).type_indices,'UniformOutput',false);
for c=1:n_col
    r=D.results(c); tk=r.topk; col=st.condition_color(r.name);
    src=old(c);
    assert(strcmp(char(strjoin(cellstr(src.Title.String),' ')),heat_titles{c}), ...
        'fig_main3_stability:TileMismatch','Heat-map tile %d is titled ''%s'', expected ''%s''.', ...
        c,char(strjoin(cellstr(src.Title.String),' ')),heat_titles{c});
    heat_pos=[(.105+(c-1)*.30)*s g.heat(1) .235*s g.heat(2)];
    ax=copyobj(src,fig); ax.Units='normalized'; ax.Position=heat_pos;
    if show_lambda, ax.Tag=sprintf('stability_r1_c%d',c); else, ax.Tag=sprintf('stability_r3_c%d',c); end
    % Recolor only this figure's copied Re(lambda)=0 reference.
    refs=findall(ax,'Type','line');
    for k=1:numel(refs)
        if numel(refs(k).XData)==2 && all(refs(k).XData==0)
            set(refs(k),'Color',st.zeroline_color,'LineStyle','--', ...
                'LineWidth',st.zeroline_lw,'Tag','zero_eigenvalue_reference');
        end
    end
    texts=findall(ax,'Type','text');
    for k=1:numel(texts)
        if startsWith(string(texts(k).String),'\lambda_1 ='), delete(texts(k)); end
    end
    if show_lambda
        title(ax,titles{c},'FontSize',16,'FontWeight','normal','Color',col);
        ax.Title.Units='normalized'; ax.Title.Position=[.5 1.14 0];
        if c==1, letter(ax,'A'); end
    else
        title(ax,'');
    end
    if c==n_col
        cb=colorbar(ax); cb.Units='normalized'; cb.Position=[(.955+(n_col-3)*.30)*s g.cbar(1) .012*s g.cbar(2)];
        cb.Label.String='log log density'; cb.Ticks=[0 .25 .5]; cb.Box='off';
        ax.Position=heat_pos;
    end
    if show_lambda, ax=newaxis(2,c,g.input(1),g.input(2)); else, ax=newaxis(1,c,g.input(1),g.input(2)); end
    offset=0;
    assert(isequal(r.type_indices,D.results(1).type_indices),'Neuron indexing differs between conditions.');
    for q=1:numel(selected)
        ids=selected{q};
        if q==1, cmap=excitatory_colormap(numel(ids)); else, cmap=inhibitory_colormap(numel(ids)); end
        for k=1:numel(ids)
            plot(ax,r.t_ex,r.u_ex(ids(k),:),'Color',cmap(k,:),'LineWidth',.65);
        end
        offset=offset+numel(ids);
    end
    ax.XAxis.Visible='off';
    if show_lambda
        if c==1, ylabel(ax,'Input'); letter(ax,'B'); end
    else
        title(ax,titles{c},'FontSize',16,'FontWeight','normal','Color',col);
        ax.Title.Units='normalized'; ax.Title.Position=[.5 1.22 0];
        if c==1, ylabel(ax,'Input'); end
    end
    t=tk.t_lya(:); finite=tk.finite_LE_spectrum_t(:,1);
    last=find(isfinite(finite),1,'last'); value=finite(last);
    if show_lambda
    ax=newaxis(3,c,g.lambda(1),g.lambda(2));
    % Saved t_lya denotes segment starts. Preserve the archived time convention.
    plot(ax,[0 T],[0 0],'--','Color',st.zeroline_color, ...
        'LineWidth',st.zeroline_lw,'Tag','zero_growth_reference');
    local_line=plot(ax,t,tk.local_LE_spectrum_t(:,1),'Color',.75+.25*col,'LineWidth',.5,'Tag','local_lambda');
    finite_line=plot(ax,t,finite,'Color',col,'LineWidth',1.7,'Tag','accumulated_lambda');
    if strcmp(r.name,'no_adaptation'), ylim(ax,[-.5 4]); else, ylim(ax,[-.5 .5]); end
    ax.XAxis.Visible='off';
    yl=ylim(ax);
    if abs(value)>=.1*diff(yl), delta=-.20; else, delta=.16; end
    ypos=min(yl(2)-.1*diff(yl),max(yl(1)+.12*diff(yl),value+delta*diff(yl)));
    text(ax,T-.08*T,ypos,sprintf('\\lambda_1 = %+.3f',value),'HorizontalAlignment','right','FontSize',14,'Color',col,'BackgroundColor','w','Margin',1);
    if c==1
        ylabel(ax,'$\lambda_1\;(s^{-1})$','Interpreter','latex'); letter(ax,'C');
        legend(ax,[local_line finite_line],{'Local','Accumulated'},'Location','northoutside', ...
            'Orientation','horizontal','Box','off','FontSize',14);
    end
    end
    if show_lambda, ax=newaxis(4,c,g.exp(1),g.exp(2)); else, ax=newaxis(2,c,g.exp(1),g.exp(2)); end
    h=sum(max(tk.local_LE_spectrum_t,0),2)/log(2);
    plot(ax,t,h,'Color',col,'LineWidth',1.1); ax.TickDir='out'; xlabel(ax,'Time (s)');
    if c==1, ylabel(ax,'Expansion (bits/s)'); if show_lambda, letter(ax,'D'); end, end
    if ~show_lambda
        % The accumulated lambda_1 the removed row printed, kept unobtrusively.
        % Just above the panel's top right corner, clear of the trace.
        text(ax,1,1.02,sprintf('\\lambda_1 = %+.3f s^{-1}',value),'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12,'Clipping','off', ...
            'Color',col,'Tag','accumulated_lambda_label');
    end
end
close(a.figs);
sep=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'Visible','off','Tag','column_dividers'); hold(sep,'on');
for xpos=(.3775+(0:n_col-2)*.30)*s, plot(sep,[xpos xpos],g.div,'Color',[.88 .88 .88],'LineWidth',2.5); end
if ~show_lambda
    % Row letters in figure coordinates, one column, where letter() puts them
    % on a full-width axis (the heat map's plot box is narrower than its axis).
    rows={'A',g.input;'B',g.exp;'C',g.heat};
    for k=1:size(rows,1)
        y=rows{k,2}(1)+1.05*rows{k,2}(2);
        text(sep,(.105-.42*.235)*s,y,['(' rows{k,1} ')'],'FontSize',14,'VerticalAlignment','top', ...
            'Clipping','off','Tag','row_letter');
    end
end
uistack(sep,'bottom'); set(findall(fig,'-property','FontSize'),'FontSize',14);
for c=1:n_col, ax=findobj(fig,'Tag',sprintf('stability_r1_c%d',c)); ax.Title.FontSize=16; end
set(findall(fig,'Tag','accumulated_lambda_label'),'FontSize',12);
if show_lambda
notes={sprintf('Current switching data: K=%d; accumulation [%g,%g] s; saved segment-start timestamps retained. Settings are read from the completed saved analysis.',D.settings.K,D.settings.lya_T_interval), ...
    ['Input shows 24 fixed actual neurons (12 evenly spaced E and 12 I), identical indices across conditions. ' ...
    'Row C displays saved local leading rates (thin/light) and accumulated rates (prominent); fixed requested y limits clip local excursions. ' ...
    'Row D shows expansive transients: sum(max(saved local QR rates,0),2)/log(2), the positive local growth summed over the tracked QR directions. Division by log(2) expresses this growth in bits/s; it is not a measurement of KS entropy or information transmission. Eigenvalue and switching stages are separate examples.'], ...
    sprintf('Columns, left to right: %s (titled and coloured by condition name; heat-map tiles taken by condition).',strjoin(names,', '))};
else
notes={sprintf('Paper layout, without the lambda_1 row (show_lambda_row false, the default since FD-031). Switching data: top-K QR, K=%d tracked directions; accumulation [%g,%g] s; saved segment-start timestamps retained. Settings are read from the completed saved analysis.',D.settings.K,D.settings.lya_T_interval), ...
    ['Rows, top to bottom: (A) input, 24 fixed actual neurons (12 evenly spaced E and 12 I), identical indices across conditions; ' ...
    sprintf('(B) expansive transients: sum(max(saved local QR rates,0),2)/log(2), the positive local growth summed over the K = %d tracked QR directions, in bits/s; a lower bound when more than K directions expand; it is not a measurement of KS entropy or information transmission; ',D.settings.K) ...
    'the small label in each (B) panel is the accumulated leading Lyapunov exponent lambda_1 at the end of the accumulation window (the value the removed row printed); ' ...
    '(C) eigenvalue density heat map. Eigenvalue and switching stages are separate examples.'], ...
    sprintf('Columns, left to right: %s (titled and coloured by condition name; heat-map tiles taken by condition).',strjoin(names,', '))};
end
out=struct('figs',fig,'files',{{}},'source',{{a.source,f}},'notes',{notes},'width_scale',width_scale);
    function ax=newaxis(row,col,y,height)
        ax=axes(fig,'Position',[(.105+(col-1)*.30)*s y .235*s height], ...
            'Tag',sprintf('stability_r%d_c%d',row,col),'FontSize',14,'LineWidth',1,'Box','off');
        hold(ax,'on'); xlim(ax,[0 T]);
    end
    function letter(ax,ch)
        text(ax,-.42,1.05,['(' ch ')'],'Units','normalized','FontSize',14, ...
            'VerticalAlignment','top','Clipping','off');
    end
end
