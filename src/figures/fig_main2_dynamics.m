function out=fig_main2_dynamics(run_dir,source_fig_root,conditions)
% FIG_MAIN2_DYNAMICS Seven-row dynamics sheet (input step first); no simulation in the plotter.
% New runs use saved numeric representative_dynamics. Current mu7 archive has
% only raster artwork: retain its trace pixels under calibrated native axes.
% Never substitute a different network's native example or invent early LLEs.
%
% CONDITIONS (FD-006): one column per condition, in the order of
% `conditions` (a cell of names; empty = every condition the stage saved, in
% its order), resolved by select_conditions against the saved results' names,
% so a name the run lacks is an error. Column titles and colours come from
% the condition NAME (grouped_condition_label 'sentence', manuscript_style),
% not from the title the stage saved. N columns keep the three-column pitch
% (one column and one gap of the three-column sheet) and the canvas widens
% with them; out.width_scale is the canvas width relative to the common
% canvas (manuscript_style canvas_w; 1 for three columns). The sheet itself
% is that canvas (see THE LAYOUT below). The archived raster
% (no numeric stage) holds exactly the three columns of the mu7 run it was
% calibrated on, no_adaptation / sfa1_std1 / sfa3_std2, and is selected from
% by the same rule.
arguments
    run_dir char
    source_fig_root char
    conditions = {}
end
root=fileparts(which('setup_paths'));
if ~startsWith(source_fig_root,'/') && isempty(regexp(source_fig_root,'^[A-Za-z]:','once'))
    source_fig_root=fullfile(root,source_fig_root);
end
f=fullfile(run_dir,'representative_dynamics','representative_dynamics_data.mat');
numeric=isfile(f);
if numeric
    D=load(f); window=D.settings.display_window; source=f;
    [names,cidx]=select_conditions({D.results.name},conditions);
    D.results=D.results(cidx);
    notes={'Native saved representative trajectories; neuron identity and E/I palettes preserved across state rows.'};
else
    source=fullfile(source_fig_root,'fig_example_timeseries','fig_example_timeseries.png');
    assert(isfile(source),'fig_main2_dynamics:MissingSource','No saved representative data or selected archived artwork.');
    fid=fopen(source,'rb'); guard=onCleanup(@()fclose(fid)); bytes=fread(fid,Inf,'*uint8');
    md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(bytes,'int8'));
    digest=reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]); clear guard
    assert(strcmpi(digest,'22c1371f18ae2938a0e6d7250dc7ee98eeb379e551e3cae4cf6e57a71380e2b9'), ...
        'fig_main2_dynamics:UncalibratedRaster','Raster calibration is valid only for the explicitly verified mu7 source.');
    artwork=imread(source); window=[0 20];
    % The archived artwork's own columns, left to right (it predates the chain).
    [names,cidx]=select_conditions({'no_adaptation','sfa1_std1','sfa3_std2'},conditions);
    notes={['Current mu7 source has no native FIG or saved x/SFA/STD trajectories. ' ...
        'Data traces remain raster artwork, calibrated to the original axes; labels/axes/legend/dividers are native. ' ...
        'No data values are digitized, interpolated into new estimates, or borrowed from another stage.'], ...
        ['Original E/I neuron shades are retained. Exact individual-neuron recoloring cannot be recovered from ' ...
        'flattened antialiased overlaps. Future numeric plots use a wider lightness/hue spread within reddish E and bluish I palettes.'], ...
        ['Old first-row neutral text pixels and colored legend strokes are masked. ' ...
        'Only pixels overwritten by source legend glyphs/strokes are missing; underlying trajectories cannot be recovered there. Finite LLE still begins at 5 s; ' ...
        'the current archive has no midpoint input step. Fixed requested y limits clip outside-range source artwork.'], ...
        ['Calibration (pixels): x origins[297,2897,5496], widths[2214,2214,2214] for[0,20] s; ' ...
        'row bounds[101,898],[1105,1902],[2108,2906],[3112,3909],[4115,4912],[5119,5916]. ' ...
        'x row maps[-10,10]; rate/synaptic/SFA rows[0,1]; depression[0,1.02]. ' ...
        'LLE tick calibration: pixel5690.5=0,34.6 pixels per inverse second. Raster-coordinate precision is about one source pixel.']};
end
% THE LAYOUT (Plan 3, 2026-11-22: the common canvas, manuscript_style's
% canvas_w, with the type at its role's size in points). The sheet is drawn
% at the full canvas width, and fig_grouped_main copies it onto the whole
% composite: what is measured here in inches is then the same on the
% composite. (Until the second round it was 0.95 of the canvas, copied into
% the middle 0.95; with the letters moved right the drawing no longer
% filled the canvas, so the margin went to the columns.) Everything below is in inches and turned into
% normalized units; nothing is scaled from a wider sheet. N columns keep
% the three-column pitch and the canvas widens with them.
st=manuscript_style();
n_col=numel(names);
% Second round (Tom, 2026-11-22): the panel letters moved right and up
% (letter_dx, letter_dy, inches from the sheet's left edge and the old
% height), the E/I key moved left into row (B) of the last column, and the
% three columns widened into the room this gives (G.left and G.right were
% 1.52 and .50 in). On the left the three-line label of row (E) fills its
% row and the gaps beside it, so the letters can move right only as far as
% its left edge allows: little is gained there. The key's place gives the
% rest.
G.left=1.50;    % panel letter, y label (up to three lines) and tick labels of column 1
G.right=.16;    % what the last column's title and tick marks need beyond the axes
G.letter_dx=.06; G.letter_dy=.08;
G.gap=.56;      % between columns: the divider, and row (G)'s own tick labels
G.top=.62;      % the column titles (two lines)
G.bottom=.40;   % the scale bar's label under row (G)
G.row_h=1.02; G.row_gap=.20;
sheet_w3=st.canvas_w;                                  % three columns
G.ax_w=(sheet_w3-G.left-G.right-2*G.gap)/3;
sheet_w=G.left+G.right+n_col*G.ax_w+(n_col-1)*G.gap;
sheet_h=G.top+G.bottom+7*G.row_h+6*G.row_gap;
width_scale=sheet_w/st.canvas_w;                       % 1 for three columns
ppi=get(groot,'ScreenPixelsPerInch');
fig=figure('Visible','off','Color','w','Position',[40 40 sheet_w*ppi sheet_h*ppi]);
% SEVEN rows, (A)-(G): the input step first, the local lambda(t) smoothed
% by a 2-Hz low-pass, half of the saved neurons in the model's E/I colours
% (excitatory_colormap / inhibitory_colormap in src/plotting, the palettes SRNNModel2.plot uses; TR, 2026-09-16).
ax_all=gobjects(7,n_col);
% Row (F) is b1, the FIRST depression variable (available resources, 0 to 1)
% of every condition, as run_paper_illustrations saves it (FD-031 item 12;
% Tom, 2026-11-17: b1 is the synaptic depression). Not the route's depression
% factor `std` (prod(b) on a product route), which the stage also saves.
rows={'u','x','r','syn','sfa','b1','lambda'};
% Words in the figure's own typeface, only the symbols in italics (Tom,
% 2026-11-22, note 6: the labels need not be LaTeX apart from the symbols).
labels={'Input, {\itu_i}',{'Dendritic','potential, {\itx_i}'},'Spike rate, {\itr_i}', ...
    {'Synaptic','output, \theta_{\iti}'},{'Spike-','frequency','adaptation'}, ...
    {'Synaptic','depression'},[regexp(st.label_growth,'^(\S+) (.*)$','tokens','once') {'(s^{-1})'}]};
% Sentence case (Tom, 2026-11-23, G1). Row (G) is manuscript_style's
% label_growth, "Growth rate, \lambda_1" (G3; it was "Lyapunov Exp."), on
% the three lines the row had: 'Growth' / 'rate, \lambda_1' / its unit
% (Plan 3 fix X5), as row (E) has three lines.
% (The label stays "Synaptic depression": with ", b_1" its long second line
% runs into the labels of rows E and G. The provenance names b1.)
titles=grouped_condition_label(names,'sentence');
% Row (G), Tom 2026-11-22: no adaptation on plus or minus 5, every adapting
% condition on plus or minus 1.
lambda_lim=ones(1,n_col); lambda_lim(strcmp(names,'no_adaptation'))=5;
for c=1:n_col
    for j=1:7
        x0=G.left+(c-1)*(G.ax_w+G.gap); y0=G.bottom+(7-j)*(G.row_h+G.row_gap);
        ax=axes(fig,'Position',[x0/sheet_w y0/sheet_h G.ax_w/sheet_w G.row_h/sheet_h], ...
            'Tag',sprintf('dynamics_r%d_c%d',j,c)); ax_all(j,c)=ax; hold(ax,'on');
        if numeric
            draw_numeric(ax,D.results(c),rows{j});
        else
            draw_artwork(ax,artwork,j,cidx(c));
        end
        set(ax,'FontSize',st.tick_fs,'LineWidth',st.axis_lw,'Box','off','YDir','normal','XLim',window);
        ax.XAxis.Visible='off';
        if c==1
            ylabel(ax,labels{j},'Interpreter','tex','FontSize',st.label_fs);
            if j==7
                % Three lines: the first ("Growth") would touch the
                % letter (G) above it; the label sits .06 in lower.
                ax.YLabel.Units='normalized';
                ax.YLabel.Position(2)=.5-.06/G.row_h;
            end
            % Near the sheet's left edge, a little above the top of the row
            % (moved right and up, second round).
            text(ax,(G.letter_dx-G.left)/G.ax_w,1.05+G.letter_dy/G.row_h,sprintf('(%c)',char('A'+j-1)),'Units','normalized', ...
                'FontSize',st.letter_fs,'HorizontalAlignment','left','VerticalAlignment','top','Clipping','off', ...
                'Tag','dynamics_letter');
        end
        % Few ticks (Tom, 2026-11-22): an axis from 0 to 1 has 0 and 1 only;
        % the input row runs 0 to .35 with 0 and .25.
        switch j
            case 1, ylim(ax,[0 .35]); yticks(ax,[0 .25]);
            case 2, ylim(ax,[-6 6]); yticks(ax,[-5 0 5]);
            case {3,4}, ylim(ax,[0 1]); yticks(ax,[0 1]);
            case 5, ylim(ax,[0 .6]); yticks(ax,[0 .5]);
            case 6
                % b1 lies in [0,1]; the archived raster keeps its calibrated [0,1.02].
                if numeric, ylim(ax,[0 1]); else, ylim(ax,[0 1.02]); end
                yticks(ax,[0 1]);
            case 7, ylim(ax,lambda_lim(c)*[-1 1]); yticks(ax,lambda_lim(c)*[-1 0 1]);
        end
        % Tick labels on the left column only; in row (G) also on a column
        % whose axis is not the left column's (the tick marks stay).
        if c>1 && ~(j==7 && lambda_lim(c)~=lambda_lim(1)), ax.YTickLabel={}; end
        if j==7
            yline(ax,0,'--','Color',[0 .55 0],'LineWidth',.5,'Tag','lambda_zero');
        end
        if j==1
            % "Single-timescale adaptation" is wider than a column of the
            % 9 in canvas: the last word goes on a second line.
            t=titles{c}; k=find(t==' ',1,'last');
            if ~isempty(k) && numel(t)>16, t={t(1:k-1),t(k+1:end)}; end
            title(ax,t,'FontSize',st.title_fs,'FontWeight','normal','Color',st.condition_color(names{c}));
            ax.Title.Units='normalized'; ax.Title.Position=[.5 1+.10/G.row_h 0];
        end
    end
end
% One E/I key. Current raster shades remain original; numeric data use the
% wider within-type color spread (no model simulation is involved).
if numeric
    ce=excitatory_colormap(1); ci=inhibitory_colormap(1);
else
    tc=SRNNCellTypePairs.type_colors(2); ce=tc(1,:); ci=tc(2,:);
end
le=plot(ax_all(2,end),NaN,NaN,'Color',ce,'LineWidth',.5);
li=plot(ax_all(2,end),NaN,NaN,'Color',ci,'LineWidth',.5);
% In row (B) of the last column, one line ("E  I") along the top right of
% the axes (second round, Tom 2026-11-22: the key further left, so that the
% columns can take the room it had beside them). The adapting columns' x
% traces stay near 0 on a [-6 6] axis, so the top of the axes is empty
% there; the sidecar says where it is. Before: right of the axes, outside.
% The axes keep their place: the legend is positioned by hand.
row_b=ax_all(2,end).Position;
lg=legend(ax_all(2,end),[le li],{'E','I'},'FontSize',st.tick_fs,'Box','off','AutoUpdate','off','Tag','EI_key', ...
    'Orientation','horizontal');
lg.ItemTokenSize=[14 18]; lg.Units='normalized'; drawnow;
q=lg.Position;
lg.Position=[row_b(1)+row_b(3)-q(3), row_b(2)+row_b(4)-q(4), q(3), q(4)];
ax_all(2,end).Position=row_b;
% Lower-left lambda panel: the time bar UNDER the axes, clear of the traces
% and of the bottom tick label, and its label under it, "5 s" as the scale
% bars of Figures 1 and 8 (Plan 3 fix X4; it was drawn inside the axes, on
% the -5 level, labelled "5 seconds"). Offsets in inches, turned into data
% units of the row. A quarter of a 20 s window or a third of a 30 s one: 5 s or 10 s.
yl=ax_all(7,1).YLim; per_in=diff(yl)/G.row_h;
bar_s=10; if diff(window)<30, bar_s=5; end
bar_x=window(1)+diff(window)/6+[0 bar_s];
bar_y=yl(1)-.08*per_in;
plot(ax_all(7,1),bar_x,bar_y*[1 1],'k','LineWidth',4,'Clipping','off','Tag','time_bar');
text(ax_all(7,1),mean(bar_x),bar_y-.05*per_in,sprintf('%d s',bar_s),'FontSize',st.tick_fs,'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Clipping','off','Tag','time_bar_label');
sep=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'Visible','off','Tag','column_dividers'); hold(sep,'on');
% Between the columns, clear of the tick labels of every row: the place is
% measured on the drawn axes (column_dividers), not typed in, and measured
% again on the composite's canvas (fig_grouped_main). From under row (G) to
% the top of row (A), below the titles.
column_dividers(sep,'dynamics',n_col,[(G.bottom-.25)/sheet_h 1-(G.top-.04)/sheet_h],'Color',[.88 .88 .88],'LineWidth',2.5);
uistack(sep,'bottom');
notes{end+1}=sprintf('Time shown: %g to %g s; the scale bar under row (G) is %d s.',window,bar_s);
notes{end+1}=['Numeric rendering draws the first 30 percent of the saved neurons of each type (round(0.3*25): 25 saved -> 8 displayed; ' ...
    'Tom, 2026-11-22: about 60 percent of the 13 drawn before), ' ...
    'the same indices in every state row, in the model''s own E/I palettes (excitatory_colormap / inhibitory_colormap, ' ...
    'the colours SRNNModel2.plot uses). Row (A) is the uniform input step. Row (G): the local rate is low-pass filtered at 2 Hz ' ...
    '(2nd-order Butterworth, zero phase) for display only and drawn in a lighter tint of the column''s condition color; ' ...
    'the finite-time (accumulated) lambda_1 is unfiltered, in the condition color. The E/I key is one line at the top right of row (B) of the last column, inside the axes. ' ...
    'The green dashed 0.5-point zero reference and the y labels (words in the figure typeface, symbols in italics) apply to both sources. Y ticks: 0 and 1 on the unit axes, 0 and 0.25 on the input row (limits 0 to 0.35); tick labels on the left column only, except row (G), where no adaptation is on plus or minus 5 and the adapting conditions are on plus or minus 1 with their own tick labels (the local rate is clipped there).'];
notes{end+1}=sprintf('Columns, left to right: %s (titled and coloured by condition name).',strjoin(names,', '));
if numeric
    notes{end+1}=['Row (F), Synaptic depression: b1, the first depression variable (available resources, 0 to 1) ' ...
        'of each neuron''s first outgoing route, field `b1` of the saved results, for every condition: the only b ' ...
        'of a one-variable route, the releasable pool of a chained route, the first factor of a product route. ' ...
        'Not the route''s depression factor `std`. Blank where the condition has no depression.'];
end
out=struct('figs',fig,'files',{{}},'source',{{source}},'notes',{notes},'width_scale',width_scale);
end
function draw_numeric(ax,r,field)
if strcmp(field,'lambda')
    if ~isfield(r,'finite_lambda') || isempty(r.finite_lambda)
        text(ax,.5,.5,'finite estimate not saved','Units','normalized','HorizontalAlignment','center','FontSize',14);
        return
    end
    % The accumulated lambda_1 in the column's condition color, the local
    % rate a lighter tint of it (polish list, 2026-11-22; they were black
    % and gray, and no condition is drawn in black or gray).
    st=manuscript_style(); col=[0 0 0];
    if st.condition_color.isKey(r.name), col=st.condition_color(r.name); end
    plot(ax,r.t_lya,lowpass_2hz(r.t_lya,r.local_lambda),'Color',.6+.4*col,'LineWidth',.8,'Tag','local_lambda');
    plot(ax,r.t_lya,r.finite_lambda,'Color',col,'LineWidth',1.6,'Tag','accumulated_lambda');
    return
end
if strcmp(field,'u')
    % Every neuron receives the same step; one black trace.
    plot(ax,r.t,r.u(1,:),'k','LineWidth',1.25); return
end
if strcmp(r.name,'no_adaptation') && ismember(field,{'sfa','std','b1'})
    return
end
if strcmp(field,'b1')
    % A condition without depression (SFA only) leaves the row blank.
    if isfield(r,'std_model') && strcmp(r.std_model,'none'), return; end
    assert(isfield(r,'b1') && ~isempty(r.b1),'fig_main2_dynamics:NoB1', ...
        ['The saved representative_dynamics results have no b1 (the first depression variable). ' ...
        'They predate FD-031 item 12: rerun the paper_illustrations stage for this run.']);
end
Y=r.(field); counts=cellfun(@numel,r.selected); offsets=[0 cumsum(counts)];
for q=numel(counts):-1:1
    % EVERY saved neuron, same indices for every state row (TR, 2026-09-15:
    % half of four per type was too few; the stage now saves 25 per type).
    % Half of the saved neurons (TR, 2026-09-16), then about 60 percent of
    % those (Tom, 2026-11-22, second round): the first 30 percent, 8 of 25.
    % The model's E / I colours.
    n=max(1,round(.3*counts(q))); ids=offsets(q)+(1:n);
    if q==1, cmap=excitatory_colormap(n); else, cmap=inhibitory_colormap(n); end
    lw=.8; if n>8, lw=.5; end
    for k=n:-1:1, plot(ax,r.t,Y(ids(k),:),'Color',cmap(k,:),'LineWidth',lw); end
end
end
function y=lowpass_2hz(t,y)
% Zero-phase 2nd-order Butterworth low-pass at 2 Hz over the finite samples;
% NaNs (before the accumulation window) are left in place.
ok=isfinite(y); if nnz(ok)<12, return; end
fs=1/median(diff(t(ok))); if fs<=4, return; end
[b,a]=butter(2,2/(fs/2));
y(ok)=filtfilt(b,a,y(ok));
end
function draw_artwork(ax,I,row,col)
if row==1
    text(ax,.5,.5,'input not archived','Units','normalized','HorizontalAlignment','center','FontSize',14); return
end
row=row-1;   % the archived raster has six rows, x .. lambda
% Calibrated native-axis placement of unmodified trace pixels. Original axes,
% titles, legends and labels are outside the retained artwork except legend
% glyphs/strokes, which are already occluded in the only surviving raster.
xl=[297 2897 5496]; xr=xl+2214;
yt=[101 1105 2108 3112 4115 5119]; yb=[898 1902 2906 3909 4912 5916];
if col==1 && ismember(row,[4 5])
    return
end
cols=(xl(col)+4):(xr(col)-4); rr=(yt(row)+4):(yb(row)-4);
C=I(rr,cols,:);
if row<=5
    % These rows contain only colored scientific traces: neutral pixels are
    % embedded axis/tick/legend text artwork. Do not erase a large legend box,
    % since the original legend background is transparent.
    neutral=max(C,[],3)-min(C,[],3)<3 & mean(double(C),3)<245;
    for ch=1:3, channel=C(:,:,ch); channel(neutral)=255; C(:,:,ch)=channel; end
end
if row==1
    % Only the two opaque horizontal legend samples, not a rectangular field
    % around the legend: retain every still-visible underlying trace pixel.
    maskrows=(rr>=210 & rr<=230) | (rr>=285 & rr<=307);
    maskcols=cols>=xl(col)+1458 & cols<=xl(col)+1714;
    C(maskrows,maskcols,:)=255;
end
if row==6
    ydata=-(double(rr([1 end]))-5690.5)/34.6;
else
    hi=[10 1 1 1 1.02]; lo=[-10 0 0 0 0];
    ydata=hi(row)-(double(rr([1 end]))-yt(row))/(yb(row)-yt(row))*(hi(row)-lo(row));
end
xdata=(double(cols([1 end]))-xl(col))/2214*20;
image(ax,'XData',xdata,'YData',ydata,'CData',C,'Tag','archived_trace_artwork');
end

function cmap=dynamics_palette(type,n)
% Greater hue/lightness separation while keeping E reddish and I bluish.
if type==1
    base=[.70 .02 .06;1 .32 .18;.38 .015 .12;.85 .18 .43];
else
    base=[.015 .16 .48;0 .65 .85;.12 .35 .85;0 .42 .48];
end
if n<=size(base,1), cmap=base(1:n,:); return; end
cmap=interp1(linspace(0,1,size(base,1)),base,linspace(0,1,n),'linear');
end
