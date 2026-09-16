function out=fig_main2_dynamics(run_dir,source_fig_root)
% FIG_MAIN2_DYNAMICS Seven-row dynamics sheet (input step first); no simulation in the plotter.
% New runs use saved numeric representative_dynamics. Current mu7 archive has
% only raster artwork: retain its trace pixels under calibrated native axes.
% Never substitute a different network's native example or invent early LLEs.
arguments
    run_dir char
    source_fig_root char
end
root=fileparts(which('setup_paths'));
if ~startsWith(source_fig_root,'/') && isempty(regexp(source_fig_root,'^[A-Za-z]:','once'))
    source_fig_root=fullfile(root,source_fig_root);
end
f=fullfile(run_dir,'representative_dynamics','representative_dynamics_data.mat');
numeric=isfile(f);
if numeric
    D=load(f); window=D.settings.display_window; source=f;
    assert(numel(D.results)==3,'fig_main2_dynamics:Conditions','Expected three saved conditions.');
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
fig=figure('Visible','off','Color','w','Position',[40 40 1450 1200]);
% SEVEN rows, (A)-(G): the input step first, the local lambda(t) smoothed
% by a 2-Hz low-pass, half of the saved neurons in the model's E/I colours
% (excitatory_colormap / inhibitory_colormap in src/plotting, the palettes SRNNModel2.plot uses; TR, 2026-09-16).
ax_all=gobjects(7,3);
rows={'u','x','r','syn','sfa','std','lambda'};
labels={'Input, $u_i$',{'Dendritic','potential, $x_i$'},'Spike rate, $r_i$', ...
    {'Synaptic','output, $\theta_i$'},{'SFA','$\frac{c_i}{K}\sum_k a_{ik}$'}, ...
    {'STD','$\prod_m b_{im}$'},'$\lambda_1$'};
titles={'No adaptation','Single-timescale adaptation','Multiple-timescale adaptation'};
st=manuscript_style(); names={'no_adaptation','sfa1_std1','sfa3_std2'};
for c=1:3
    for j=1:7
        ax=axes(fig,'Position',[.105+(c-1)*.30 .085+(7-j)*.1225 .245 .092], ...
            'Tag',sprintf('dynamics_r%d_c%d',j,c)); ax_all(j,c)=ax; hold(ax,'on');
        if numeric
            draw_numeric(ax,D.results(c),rows{j});
        else
            draw_artwork(ax,artwork,j,c);
        end
        set(ax,'FontSize',14,'LineWidth',1,'Box','off','YDir','normal','XLim',window);
        ax.XAxis.Visible='off';
        if c==1
            ylabel(ax,labels{j},'Interpreter','latex','FontSize',14);
            text(ax,-.32,1.05,sprintf('(%c)',char('A'+j-1)),'Units','normalized', ...
                'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top','Clipping','off');
        end
        switch j
            case 1, ylim(ax,[-.05 .55]); yticks(ax,[0 .25 .5]);
            case 2, ylim(ax,[-6 6]); yticks(ax,[-5 0 5]);
            case {3,4}, ylim(ax,[0 1]); yticks(ax,[0 .5 1]);
            case 5, ylim(ax,[0 .6]); yticks(ax,[0 .5]);
            case 6, ylim(ax,[0 1.02]); yticks(ax,[0 .5 1]);
            case 7, ylim(ax,[-5 5]); yticks(ax,[-5 0 5]);
        end
        if j==7
            yline(ax,0,'--','Color',[0 .55 0],'LineWidth',.5,'Tag','lambda_zero');
        end
        if j==1
            title(ax,titles{c},'FontSize',16,'FontWeight','normal','Color',st.condition_color(names{c}));
            ax.Title.Units='normalized'; ax.Title.Position=[.5 1.20 0];
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
le=plot(ax_all(2,3),NaN,NaN,'Color',ce,'LineWidth',.5);
li=plot(ax_all(2,3),NaN,NaN,'Color',ci,'LineWidth',.5);
legend(ax_all(2,3),[le li],{'E','I'},'FontSize',14,'Box','off','Location','northeast');
% Requested lower-left lambda panel, with explicit data coordinates.
plot(ax_all(7,1),[5 15],[-4.8 -4.8],'k','LineWidth',4,'Tag','ten_second_bar');
text(ax_all(7,1),10,-5.6,'10 seconds','FontSize',14,'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Clipping','off','Tag','ten_second_label');
sep=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'Visible','off','Tag','column_dividers'); hold(sep,'on');
for xpos=[.3775 .6775]
    plot(sep,[xpos xpos],[.07 .945],'Color',[.88 .88 .88],'LineWidth',2.5);
end
uistack(sep,'bottom');
notes{end+1}=['Numeric rendering draws the first half of the saved neurons of each type (25 saved -> 12 displayed), ' ...
    'the same indices in every state row, in the model''s own E/I palettes (excitatory_colormap / inhibitory_colormap, ' ...
    'the colours SRNNModel2.plot uses). Row (A) is the uniform input step. Row (G): the local rate is low-pass filtered at 2 Hz ' ...
    '(2nd-order Butterworth, zero phase) for display only; the finite-time lambda_1 is unfiltered. ' ...
    'The green dashed 0.5-point zero reference and descriptive LaTeX labels apply to both sources.'];
out=struct('figs',fig,'files',{{}},'source',{{source}},'notes',{notes});
end
function draw_numeric(ax,r,field)
if strcmp(field,'lambda')
    if ~isfield(r,'finite_lambda') || isempty(r.finite_lambda)
        text(ax,.5,.5,'finite estimate not saved','Units','normalized','HorizontalAlignment','center','FontSize',14);
        return
    end
    plot(ax,r.t_lya,lowpass_2hz(r.t_lya,r.local_lambda),'Color',[.6 .6 .6],'LineWidth',.8);
    plot(ax,r.t_lya,r.finite_lambda,'k','LineWidth',1.25);
    return
end
if strcmp(field,'u')
    % Every neuron receives the same step; one black trace.
    plot(ax,r.t,r.u(1,:),'k','LineWidth',1.25); return
end
if strcmp(r.name,'no_adaptation') && ismember(field,{'sfa','std'})
    msg='no SFA'; if strcmp(field,'std'), msg='no STD'; end
    text(ax,.5,.5,msg,'Units','normalized','HorizontalAlignment','center','FontSize',14); return
end
Y=r.(field); counts=cellfun(@numel,r.selected); offsets=[0 cumsum(counts)];
for q=numel(counts):-1:1
    % EVERY saved neuron, same indices for every state row (TR, 2026-09-15:
    % half of four per type was too few; the stage now saves 25 per type).
    % Half of the saved neurons (TR, 2026-09-16), the model's E / I colours.
    n=max(1,round(counts(q)/2)); ids=offsets(q)+(1:n);
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
    msg='no SFA'; if row==5, msg='no STD'; end
    text(ax,.5,.5,msg,'Units','normalized','HorizontalAlignment','center','FontSize',14); return
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
