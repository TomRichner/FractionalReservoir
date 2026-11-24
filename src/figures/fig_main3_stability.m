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
% widens with them (out.width_scale, 1 for three columns). Row (C)'s y range is [-.5 4] for no_adaptation (by name) and
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
%
% RITZ_OVERLAY (FD-033, default false: the figure as it was). true passes
% fig_eig_heatmap its 'ritz_overlay' with ritz_k, ritz_scaling, ritz_style
% and ritz_inset, so each heat-map panel also carries the Ritz values of the
% leading Lyapunov subspace (and its inset near zero when asked); the
% points and ticks are children of the copied axes and come with them, and
% each inset is put back in the same place within its panel's plot box. The
% eig_heatmap stage must have been run with the Ritz pass.
%
% SPECTRA_SOURCE (FD-034, default 'eig_heatmap': the figure as it was).
% 'local_lyapunov' takes the heat map, and its Ritz overlay, from the states
% run_local_lyapunov sampled along the switching run itself
% ('spectra_n_samples'; fig_eig_heatmap's option of the same name), so the
% expansion row and the heat-map row are one trajectory and one top-K pass:
% the Ritz values are the instantaneous rates, with rotation, of the
% directions whose positive local rates the expansion row sums.
%
% REAL_HISTOGRAM_ROW (Tom, 2026-11-20; default false: the figure as it was;
% needs show_lambda_row false). true adds a fourth row, (D), under the heat
% map: fig_eig_real_histogram for the same conditions and spectra_source
% (the eigenvalues of row (C) projected on the real axis, log counts per
% snapshot), each panel as wide as its column. The numbers and a caption
% sentence are in out.table_text (fig_grouped_main saves it as
% <tag>_table.md). As Tom asked in his review of 2026-11-22, the row has no
% dotted curve, no key, no labels and no leader lines, runs from -11 to +3
% s^-1 and is named "Eigenvalues per snapshot". What he took out stays
% reachable:
%   HISTOGRAM_XLIM (default [-11 3], s^-1)   the row's x limits; [] puts it
%       on the heat map's x limits (the bins are always those of the heat
%       map's window).
%   HISTOGRAM_REFERENCE (default false)   the dotted projection of a uniform
%       disk (fig_eig_real_histogram's 'reference'). It needs the frozen
%       gains in the stage's file (disk_by_cond); without them it is the
%       unit-gain disk.
%   HISTOGRAM_LABELS (default false)   fig_eig_real_histogram's 'labels'
%       ("neurons at floor or ceiling" at the spike at -1/tau_d, "adaptive
%       states" at the pile near zero, "network modes" at the plateau, each
%       with its leader; with the reference, its key). They are placed on
%       the heat map's real axis, so they need histogram_xlim [].
%
% THE CANVAS (2026-11-22) is manuscript_style's canvas_w (9 in) for three
% columns, with the type at manuscript_style's sizes by role (canvas_type)
% and axis lines at axis_lw; the height is what the rows need (8.6 in with
% row (D), 6.6 in without, 11 in with the lambda_1 row), all less than the
% workstation's screen: MATLAB shrinks a window larger than the screen, at
% creation or at export and not every time (seen on the earlier 1152 px
% canvas), so the size is checked after the figure is made. On 9 in the
% column headings are on two lines, the y axes are named in the first
% column only, and the row letters sit above each row's top left corner.
%
% ROW (A) draws the input to 12 E and 12 I neurons (rows of the saved u_ex,
% evenly spaced within each type's indices), each trace in the color of the
% type of the neuron it drives: the E and I line color maps of
% fig_main2_dynamics (excitatory_colormap, inhibitory_colormap). ROW (B) has
% one y axis in every column, 0 to 75 bits/s; larger peaks are clipped.
%
% THE PAPER'S FIGURE 3 (figure plan of 2026-11-21) is spectra_source
% 'local_lyapunov', real_histogram_row true, ritz_overlay false: the heat
% map of the switching run without the Ritz points, their exponent ticks or
% the inset. The overlay stays an option (ritz_overlay, ritz_inset) for the
% appendix version. The lambda_1 label above row (B) is at the figure's
% text size; its tag is 'accumulated_lambda_value'. Rows (C) and (D) are
% labeled as rates, in s^-1.
arguments
    run_dir char
    conditions = {}
    opts.show_lambda_row (1,1) logical = false
    opts.ritz_overlay (1,1) logical = false
    opts.ritz_k (1,1) double = 0
    opts.ritz_scaling (1,:) char = 'raw'
    opts.ritz_inset (1,1) logical = false
    opts.ritz_style (1,:) char = 'density'
    opts.spectra_source (1,:) char {mustBeMember(opts.spectra_source, {'eig_heatmap', 'local_lyapunov'})} = 'eig_heatmap'
    opts.real_histogram_row (1,1) logical = false
    opts.histogram_labels (1,1) logical = false
    opts.histogram_reference (1,1) logical = false
    opts.histogram_xlim double = [-11 3]
end
show_lambda=opts.show_lambda_row;
EXPANSION_YLIM=[0 75];   % bits/s, row (B) of every column
assert(~(opts.real_histogram_row && show_lambda),'fig_main3_stability:HistogramRow', ...
    'real_histogram_row is a row under the heat map of the layout without the lambda_1 row: set show_lambda_row false.');
assert(~opts.histogram_labels || isempty(opts.histogram_xlim),'fig_main3_stability:HistogramLabels', ...
    'The labels of row (D) are placed on the heat map''s real axis: pass ''histogram_xlim'', [] with them.');
assert(isempty(opts.histogram_xlim) || (numel(opts.histogram_xlim)==2 && diff(opts.histogram_xlim)>0), ...
    'fig_main3_stability:HistogramXLim','histogram_xlim is [lo hi] in s^-1, or [] for the heat map''s limits.');
st=manuscript_style();
% THE CANVAS (2026-11-22): manuscript_style's canvas_w (9 in) for three
% columns; the height is what the rows need, in inches below. Row geometry,
% [y height], and the columns (left edge of the first, width, gap) in
% normalized figure units. The columns leave room on the left for a
% two-line y label beside four-character tick labels, between them for the
% tick labels of the next column and the divider, and on the right for the
% color bar with its label.
col_x0=.10; col_w=.205; col_gap=.089; col_pitch=col_w+col_gap;
if show_lambda
    fig_h_in=11; g.heat=[.745 .18]; g.input=[.535 .125]; g.lambda=[.315 .14]; g.exp=[.105 .135];
    g.cbar=[.755 .16]; g.div=[.08 .93];
elseif opts.real_histogram_row
    fig_h_in=8.6; g.input=[.855 .085]; g.exp=[.665 .135]; g.heat=[.335 .225]; g.hist=[.06 .17];
    g.cbar=[.345 .205]; g.div=[.03 .945];
else
    fig_h_in=6.6; g.input=[.80 .12]; g.exp=[.52 .19]; g.heat=[.06 .31];
    g.cbar=[.075 .28]; g.div=[.03 .945];
end
f=fullfile(run_dir,'local_lyapunov','local_lyapunov_data.mat'); D=load(f);
[names,cidx]=select_conditions({D.results.name},conditions);
D.results=D.results(cidx);
n_col=numel(names); width_scale=1+(n_col-3)*col_pitch; s=1/width_scale;
T=D.settings.T;
ritz_args={};
if opts.ritz_overlay
    ritz_args={'ritz_overlay',true,'ritz_k',opts.ritz_k,'ritz_scaling',opts.ritz_scaling,'ritz_inset',opts.ritz_inset,'ritz_style',opts.ritz_style};
end
a=fig_eig_heatmap('run_dir',run_dir,'save',false,'visible',false,'density_scale','loglog', ...
    'conditions',names,ritz_args{:},'spectra_source',opts.spectra_source);
ppi=get(groot,'ScreenPixelsPerInch');
fig=figure('Visible','off','Color','w','Position',[40 40 st.canvas_w*ppi*width_scale fig_h_in*ppi]);
drawnow;
% MATLAB shrinks a window larger than the screen; the layout below is made for this size.
assert(max(abs(fig.Position(3:4)-[st.canvas_w*width_scale fig_h_in]*ppi))<1,'fig_main3_stability:CanvasShrunk', ...
    'The figure is %g x %g px, not the %g x %g px asked for (larger than the screen?).', ...
    fig.Position(3:4),[st.canvas_w*width_scale fig_h_in]*ppi);
old=findall(a.figs,'Type','axes');
% The Ritz insets are axes of their own, not tiles; they follow their panel below.
is_inset=strcmp(get(old,'Tag'),'ritz_inset'); insets=old(is_inset); old=old(~is_inset);
% Tile k is the k-th selected condition (fig_eig_heatmap draws them in the
% selection's order, row-major); findall returns the reverse creation order.
[~,order]=sort(arrayfun(@(x)x.Layout.Tile,old)); old=old(order);
assert(numel(old)==n_col,'fig_main3_stability:Tiles', ...
    'fig_eig_heatmap drew %d tiles for %d conditions.',numel(old),n_col);
titles=grouped_condition_label(names,'sentence');   % one line: what the inner figures title their tiles
col_titles=grouped_condition_label(names,'two_line');   % the column headings: two lines on the 9 in canvas
heat_titles=titles_by_name(names);
if opts.real_histogram_row
    % Row (D): fig_eig_real_histogram for the same selection and source.
    hh=fig_eig_real_histogram('run_dir',run_dir,'save',false,'visible',false,'verbose','near-none', ...
        'conditions',names,'spectra_source',opts.spectra_source,'labels',opts.histogram_labels, ...
        'reference',opts.histogram_reference);
    hist_old=findall(hh.figs,'Type','axes');
    [~,order]=sort(arrayfun(@(x)x.Layout.Tile,hist_old)); hist_old=hist_old(order);
    assert(numel(hist_old)==n_col,'fig_main3_stability:Tiles', ...
        'fig_eig_real_histogram drew %d tiles for %d conditions.',numel(hist_old),n_col);
end
selected=cellfun(@(v)v(round(linspace(1,numel(v),min(12,numel(v))))), ...
    D.results(1).type_indices,'UniformOutput',false);
for c=1:n_col
    r=D.results(c); tk=r.topk; col=st.condition_color(r.name);
    src=old(c);
    assert(strcmp(char(strjoin(cellstr(src.Title.String),' ')),heat_titles{c}), ...
        'fig_main3_stability:TileMismatch','Heat-map tile %d is titled ''%s'', expected ''%s''.', ...
        c,char(strjoin(cellstr(src.Title.String),' ')),heat_titles{c});
    heat_pos=[(col_x0+(c-1)*col_pitch)*s g.heat(1) col_w*s g.heat(2)];
    ax=copyobj(src,fig); ax.Units='normalized'; ax.Position=heat_pos; ax.LineWidth=st.axis_lw;
    if show_lambda, ax.Tag=sprintf('stability_r1_c%d',c); else, ax.Tag=sprintf('stability_r3_c%d',c); end
    % Recolor only this figure's copied Re(lambda)=0 reference.
    refs=findall(ax,'Type','line');
    for k=1:numel(refs)
        if numel(refs(k).XData)==2 && all(refs(k).XData==0)
            set(refs(k),'Color',st.zeroline_color,'LineStyle','--', ...
                'LineWidth',st.zeroline_lw,'Tag','zero_eigenvalue_reference');
        end
    end
    % Rates, in s^-1; ticks every 10 on the real axis, fixed so that a later
    % re-layout (fig_grouped_main) cannot change them or rotate their labels.
    xt=10*(ceil(ax.XLim(1)/10):floor(ax.XLim(2)/10));
    if numel(xt)>=2, ax.XTick=xt; end
    ax.XTickLabelRotation=0; xlabel(ax,'Re(\lambda) (s^{-1})');
    % Named in the first column only, as the other rows are: the gap between
    % the columns holds the tick labels and the divider.
    if c==1, ylabel(ax,'Im(\lambda) (s^{-1})'); else, ylabel(ax,''); end
    texts=findall(ax,'Type','text');
    for k=1:numel(texts)
        if startsWith(string(texts(k).String),'\lambda_1 ='), delete(texts(k)); end
    end
    if show_lambda
        title(ax,col_titles{c},'FontSize',st.title_fs,'FontWeight','normal','Color',col);
        ax.Title.Units='normalized'; ax.Title.Position=[.5 1.14 0];
        if c==1, letter(ax,'A'); end
    else
        title(ax,'');
    end
    if c==n_col
        cb=colorbar(ax); cb.Units='normalized'; cb.Position=[heat_pos(1)+(col_w+.014)*s g.cbar(1) .012*s g.cbar(2)];
        cb.Label.String='Log log density'; cb.Ticks=[0 .25 .5]; cb.Box='off'; cb.LineWidth=st.axis_lw;   % label in sentence case (Tom, 2026-11-23, G1)
        ax.Position=heat_pos;
    end
    for k=1:numel(insets)
        if insets(k).UserData.panel~=c, continue; end
        box=plot_box_position(ax); rel=insets(k).UserData.rel;
        ins=copyobj(insets(k),fig); ins.Units='normalized';
        ins.Position=[box(1)+rel(1)*box(3) box(2)+rel(2)*box(4) rel(3)*box(3) rel(4)*box(4)];
    end
    if opts.real_histogram_row
        % The eigenvalues of the heat map above, projected on the real axis:
        % the column's width, on histogram_xlim (Tom, 2026-11-22: -11 to +3
        % s^-1, the part between the leak rate and zero; the bins are those
        % of the heat map's window), or on the heat map's x limits when empty.
        assert(strcmp(char(strjoin(cellstr(hist_old(c).Title.String),' ')),titles{c}), ...
            'fig_main3_stability:TileMismatch','Histogram tile %d is not ''%s''.',c,titles{c});
        assert(max(abs(hh.edges([1 end])-ax.XLim))<=1e-9*diff(ax.XLim),'fig_main3_stability:HistogramWindow', ...
            'The histogram''s bins do not span the heat map''s x limits.');
        hx=copyobj(hist_old(c),fig); hx.Units='normalized'; hx.Position=[heat_pos(1) g.hist(1) heat_pos(3) g.hist(2)];
        hx.Tag=sprintf('stability_r4_c%d',c); hx.LineWidth=st.axis_lw; hx.XTickLabelRotation=0;
        if isempty(opts.histogram_xlim)
            hx.XLim=ax.XLim; hx.XTick=ax.XTick;
        else
            hx.XLim=opts.histogram_xlim; hx.XTick=5*(ceil(hx.XLim(1)/5):floor(hx.XLim(2)/5));
        end
        title(hx,''); xlabel(hx,'Re(\lambda) (s^{-1})');
        if strcmp(hx.YScale,'log')
            % Every second decade, labeled plainly: the default leaves one
            % label on this short axis, and every decade is too many for it.
            dec=ceil(log10(hx.YLim(1))):floor(log10(hx.YLim(2))); dec=dec(mod(dec,2)==1);
            hx.YTick=10.^dec; hx.YTickLabel=arrayfun(@(v)sprintf('%g',v),10.^dec,'UniformOutput',false);
        end
        if c==1, ylabel(hx,{'Eigenvalues','per snapshot'}); else, ylabel(hx,''); end
    end
    if show_lambda, ax=newaxis(2,c,g.input(1),g.input(2)); else, ax=newaxis(1,c,g.input(1),g.input(2)); end
    assert(isequal(r.type_indices,D.results(1).type_indices),'Neuron indexing differs between conditions.');
    assert(isequal(reshape(r.cell_type_names,1,[]),{'E','I'}),'fig_main3_stability:CellTypes', ...
        'The input row colors by E and I; this run''s cell types are %s.',strjoin(r.cell_type_names,', '));
    % Each trace is the input to ONE neuron (a row of the saved u_ex), in
    % the color of that neuron's type: the E and I line color maps and the
    % drawing order of fig_main2_dynamics (I first, E on top).
    for q=numel(selected):-1:1
        ids=selected{q};
        if q==1, cmap=excitatory_colormap(numel(ids)); else, cmap=inhibitory_colormap(numel(ids)); end
        for k=numel(ids):-1:1
            plot(ax,r.t_ex,r.u_ex(ids(k),:),'Color',cmap(k,:),'LineWidth',.5,'Tag',['input_' r.cell_type_names{q}]);
        end
    end
    ax.XAxis.Visible='off';
    if show_lambda
        if c==1, ylabel(ax,'Input'); letter(ax,'B'); end
    else
        title(ax,col_titles{c},'FontSize',st.title_fs,'FontWeight','normal','Color',col);
        ax.Title.Units='normalized'; ax.Title.Position=[.5 1.12 0];
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
    text(ax,T-.08*T,ypos,sprintf('\\lambda_1 = %+.3f',value),'HorizontalAlignment','right','FontSize',st.tick_fs,'Color',col,'BackgroundColor','w','Margin',1);
    if c==1
        ylabel(ax,'$\lambda_1\;(s^{-1})$','Interpreter','latex'); letter(ax,'C');
        legend(ax,[local_line finite_line],{'Local','Accumulated'},'Location','northoutside', ...
            'Orientation','horizontal','Box','off','FontSize',st.tick_fs);
    end
    end
    if show_lambda, ax=newaxis(4,c,g.exp(1),g.exp(2)); else, ax=newaxis(2,c,g.exp(1),g.exp(2)); end
    h=sum(max(tk.local_LE_spectrum_t,0),2)/log(2);
    plot(ax,t,h,'Color',col,'LineWidth',1.1); ax.TickDir='out'; xlabel(ax,'Time (s)');
    % One y axis in every column (Tom, 2026-11-22); peaks above it are clipped.
    ylim(ax,EXPANSION_YLIM); ax.YTick=EXPANSION_YLIM(1):25:EXPANSION_YLIM(2);
    if c==1, ylabel(ax,{'Expansion','(bits/s)'}); if show_lambda, letter(ax,'D'); end, end
    if ~show_lambda
        % The accumulated lambda_1 the removed row printed, kept unobtrusively.
        % Just above the panel's top right corner, clear of the trace.
        text(ax,1,1.02,sprintf('\\lambda_1 = %+.3f s^{-1}',value),'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',st.tick_fs,'Clipping','off', ...
            'Color',col,'Tag','accumulated_lambda_value');
    end
end
close(a.figs);
if opts.real_histogram_row, close(hh.figs); end
sep=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'Visible','off','Tag','column_dividers'); hold(sep,'on');
% Between the columns, clear of the tick labels of every row: the place is
% measured on the drawn axes (column_dividers), not typed in; again below,
% once the type has its final size, and on the composite's canvas.
column_dividers(sep,'stability',n_col,g.div,'Color',[.88 .88 .88],'LineWidth',2.5);
if ~show_lambda
    % Row letters in figure coordinates, one column at the canvas's left
    % edge, each just above its row's top: beside the row the margin holds
    % the y label, which is as long as the row is tall.
    rows={'A',g.input;'B',g.exp;'C',g.heat};
    if opts.real_histogram_row, rows(end+1,:)={'D',g.hist}; end
    for k=1:size(rows,1)
        y=rows{k,2}(1)+rows{k,2}(2)+.006;
        text(sep,.004*s,y,['(' rows{k,1} ')'],'FontSize',st.letter_fs,'VerticalAlignment','bottom', ...
            'Clipping','off','Tag','row_letter');
    end
end
uistack(sep,'bottom');
canvas_type(fig);   % every piece of type at its role's size (manuscript_style)
set(findall(fig,'Tag','ritz_inset'),'FontSize',9); set(findall(fig,'-regexp','Tag','^ritz_(annotation|tick_label)$'),'FontSize',9);
column_dividers(fig);   % the type is final now: place the dividers on it
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
if opts.ritz_overlay
    notes{end+1}=sprintf(['Ritz overlay on the eigenvalue density (ritz_k %d, 0 = all the run tracked; %s norm; inset %d): points are the ' ...
        'eigenvalues of the Jacobian projected on the leading Lyapunov directions at the same sampled states as the density, in the ' ...
        'condition''s color; ticks on the bottom edge are that pass''s Lyapunov exponents. In a stable adapting network these directions ' ...
        'are the band of slow single-neuron adaptation modes, not the network''s collective modes. Numbers and how to read them: the ' ...
        'table fig_eig_heatmap writes beside its own figure (fig_eig_heatmap_ritz*_table.md).'], ...
        opts.ritz_k,opts.ritz_scaling,opts.ritz_inset);
end
if strcmp(opts.spectra_source,'local_lyapunov')
    % Said last, because it replaces the "separate examples" statement above.
    notes{end+1}=sprintf(['spectra_source local_lyapunov (FD-034): the heat map and its Ritz values are NOT a separate example. ' ...
        'They are taken at states of the switching run of the expansion row, and the Ritz values come from the top-%d pass ' ...
        'whose positive local rates that row sums. That pass: alignment over [%g, %g] s, accumulation over [%g, %g] s.'], ...
        D.settings.K,D.settings.lya_T_interval(1)-D.settings.lya_warmup,D.settings.lya_T_interval(1),D.settings.lya_T_interval);
end
out=struct('figs',fig,'files',{{}},'source',{{a.source,f}},'notes',{notes},'width_scale',width_scale);
if opts.real_histogram_row
    if isempty(opts.histogram_xlim), hist_span='(C)''s x limits';
    else, hist_span=sprintf('[%g, %g] s^-1',opts.histogram_xlim); end
    if opts.histogram_reference, dotted='Dotted: the projection of that uniform disk on the real axis (a semicircle). ';
    else, dotted=''; end
    row_d=['Row (D): for each condition, the real parts of the eigenvalues of the Jacobian at the sampled states of row (C), ' ...
        'as a histogram (bars: eigenvalues per snapshot in each bin of (C)''s real axis, logarithmic count axis), drawn over ' ...
        hist_span '. ' dotted ...
        'Dashed: zero. Eigenvalues right of zero per snapshot against the expectation of a uniform disk (the network''s n ' ...
        'eigenvalues, centered at -1/tau_d, its radius the random-matrix radius of the connectivity times the ' ...
        'root-mean-square frozen gain at those snapshots): ' ...
        strjoin(arrayfun(@(k)sprintf('%s %.2f against %.2f',names{k},hh.reference(k).pos_per_state,hh.reference(k).pos_expected),1:n_col, ...
        'UniformOutput',false),'; ') ' (the table''s last column; its other columns give the disk and the rest).'];
    out.notes{end+1}=row_d;
    if opts.histogram_labels
        out.notes{end+1}='Labels on row (D), each once: "neurons at floor or ceiling" left of the spike at -1/tau_d in the first column; "network modes" above the plateau, left of the spike, in the second; "adaptive states" beside the pile near zero in the last column that has adaptive variables; with the dotted curve, its key right of the last panel, outside the axes ("uniform disk, same radius").';
    end
    out.table_text=sprintf('%s\n\n## Row (D) of the composite\n\n%s\n',hh.table_text,row_d);
end
    function ax=newaxis(row,col,y,height)
        ax=axes(fig,'Position',[(col_x0+(col-1)*col_pitch)*s y col_w*s height], ...
            'Tag',sprintf('stability_r%d_c%d',row,col),'FontSize',st.tick_fs,'LineWidth',st.axis_lw,'Box','off');
        hold(ax,'on'); xlim(ax,[0 T]);
    end
    function letter(ax,ch)
        text(ax,-.42,1.05,['(' ch ')'],'Units','normalized','FontSize',st.letter_fs, ...
            'VerticalAlignment','top','Clipping','off');
    end
end
