function out = fig_stim_combined(clinical_dir, model_fig, opts)
% FIG_STIM_COMBINED The paper's stimulation figure: clinical stimulation on
% top, the IED-model stimulation experiment below, one continuous lettering.
%
%   out = FIG_STIM_COMBINED(clinical_dir, model_fig, 'out_dir', d)
%
% clinical_dir  folder holding the clinical panels saved one per file, as
%               cscs_stim_stability writes them under
%               results/manuscript_stim/tica_awake/panels/ (read only):
%                 Fig_Awake_Main_PanelD.fig             example traces
%                 Fig_Awake_Main_PanelA_Normalized.fig  IED rate, % of no stim
%                 Fig_Awake_Main_PanelB_PSD.fig         SOZ PSD, median + IQR
%               (file names overridable through opts.clinical_files).
%               'clinical_psd_fig', f takes the PSD panel from the .fig at the
%               full path f instead (one axes, drawn like PanelB_PSD), e.g.
%               the no-harmonic-removal version under psd_2hz_check/.
% model_fig     the model figure written by fig_ied_stim (its .fig): panels
%               A event rate, B PSD, C band power, D correlation, E state
%               change, F example traces with cumulative lambda_1. A newer
%               model version is swapped in by passing its .fig; the axes are
%               found by ROLE (y label, title, line colour), not by order, so
%               any figure fig_ied_stim draws will do.
%
% Nothing is recomputed: every panel is the axes of its source figure copied
% into a new canvas, re-positioned, and restyled to the grouped main figures'
% look (Helvetica, one font size, "(A)" letters top left, no box, ticks out).
% Figure-level header text of the sources (the clinical "S9, SOZ Contacts"
% line, the sources' own panel letters) is dropped; scale bars, event
% triangles, stimulation ticks and legends are kept.
%
% Layout (letters run left to right, top to bottom):
%   A  clinical example traces, no stimulation | stimulation
%   B  clinical IED rate per participant      C  clinical SOZ PSD
%   D  model event rate  E  model band power  F  model PSD   (D under B, F under C)
%   G  model example traces, no stimulation over stimulation  | H correlation change
%                                                              | I state change
%
% Layout 'mirror' (default since 2026-11-17, Tom): the model half mirrors the
% clinical half panel for panel, and only the matching panels are kept:
%   A  clinical traces, no stimulation | stimulation
%   B  clinical IED rate                C  clinical SOZ PSD
%   D  model traces, no stimulation | stimulation (no lambda_1 strip; time
%      axis replaced by a scale bar, as in A)
%   E  model event rate                 F  model PSD
% Its text prints at the grouped main figures' size (9 pt on a 7.0 x 7.65 in
% canvas; fig_grouped_main's 14 pt on an 11-in PNG is 8.9 pt at 7 in), with
% the panel letters larger (12 pt, 'letter_size'; Tom, 2026-11-17), and
% the rows are packed so the gap between the halves matches those within.
% The model rate axes is found by a y label containing 'Events per Minute' or
% 'Event Rate'; a legend attached to it in the source (one entry per network)
% is carried over like the clinical one. Layout 'v1' is the nine-panel A-I
% figure of 2026-11-17 (figs/stim_combined_20261117).
%
% Output: <out_dir>/<tag>.{png,fig,svg} and <tag>_provenance.md.
%
% See also: fig_ied_stim, fig_grouped_main, manuscript_style
arguments
    clinical_dir (1,:) char
    model_fig (1,:) char
    opts.out_dir (1,:) char = ''
    opts.tag (1,:) char = 'Fig_Stim_Combined'
    opts.clinical_files struct = struct('traces','Fig_Awake_Main_PanelD.fig', ...
        'rate','Fig_Awake_Main_PanelA_Normalized.fig','psd','Fig_Awake_Main_PanelB_PSD.fig')
    opts.clinical_psd_fig (1,:) char = ''   % full path of the clinical PSD panel; '' = clinical_dir/clinical_files.psd
    opts.clinical_note (1,:) char = ''   % provenance line, e.g. the source repo commit
    opts.model_note (1,:) char = ''
    opts.drop_text (1,:) cell = {'SOZ Contacts$'}   % regexps of header texts to delete
    % NaN = the layout's default. mirror: 9 / 8 pt on 7.0 x 7.65 in, which
    % prints the letters, labels and ticks at the size of the grouped main
    % figures (14 pt on their 11.04-in-wide PNG = 8.9 pt at 7 in; 2026-11-17).
    % v1: 7.5 / 6.5 pt on 7.0 x 9.4 in, as it was drawn.
    opts.font_size (1,1) double = NaN   % points, at the printed size
    opts.annot_size (1,1) double = NaN  % channel names, scale bars, legends
    opts.letter_size (1,1) double = NaN % panel letters; NaN = 12 pt (mirror), font_size (v1)
    opts.size_in (1,2) double = [NaN NaN]   % PNAS full page width is 7 in
    opts.resolution (1,1) double = 450
    opts.save (1,1) logical = true
    opts.visible (1,1) logical = false
    opts.layout (1,:) char {mustBeMember(opts.layout,{'mirror','v1'})} = 'mirror'
    opts.time_bar (1,1) double = 5   % s, model traces scale bar (mirror layout)
end
assert(isfolder(clinical_dir),'fig_stim_combined:NoClinical','No folder %s',clinical_dir);
assert(isfile(model_fig),'fig_stim_combined:NoModel','No model figure %s',model_fig);
cf = opts.clinical_files;
src = struct('traces',fullfile(clinical_dir,cf.traces),'rate',fullfile(clinical_dir,cf.rate), ...
    'psd',fullfile(clinical_dir,cf.psd),'model',model_fig);
if ~isempty(opts.clinical_psd_fig), src.psd = opts.clinical_psd_fig; end
for fn = fieldnames(src)'
    assert(isfile(src.(fn{1})),'fig_stim_combined:MissingSource','Missing %s',src.(fn{1}));
end
purple = [0.494 0.247 0.596]; green = [0 0.62 0.45];
is_mirror = strcmp(opts.layout,'mirror');
if isnan(opts.font_size),  opts.font_size  = pick_default(is_mirror,9,7.5); end
if isnan(opts.annot_size), opts.annot_size = pick_default(is_mirror,8,6.5); end
if isnan(opts.letter_size), opts.letter_size = pick_default(is_mirror,12,opts.font_size); end
if any(isnan(opts.size_in))
    if is_mirror, opts.size_in = [7.0 7.65]; else, opts.size_in = [7.0 9.4]; end
end

% The canvas is the printed size (PNAS full page is 7 in wide), so font sizes
% are true point sizes. It must also fit the screen: exportgraphics renders
% a figure taller than the monitor clipped to the monitor height (seen
% 2026-11-17 on a 1073-px-high display), whatever fig.Position says.
fig = figure('Visible',onoff(opts.visible),'Color','w','Units','inches', ...
    'Position',[0.2 0.2 opts.size_in]);
drawnow;
mon = get(0,'MonitorPositions'); px = getpixelposition(fig);
if px(4) > max(mon(:,4)) - 40
    warning('fig_stim_combined:TooTall', ...
        'Canvas %d px high on a %d-px screen; the PNG will be clipped.',round(px(4)),max(mon(:,4)));
end

% ---- target rectangles (normalized inner positions) ----------------------
L = 0.085;                         % left edge of left-column axes
R = 0.975;                         % right edge
yA = 0.815; hA = 0.160;            % clinical traces
yB = 0.630; hB = 0.130;            % clinical rate / PSD
yD = 0.430; hD = 0.130;            % model rate / band power / PSD
xC = 0.660;                        % PSD column (C and F)
pos.cl_trace_no  = [L+0.020 yA 0.405 hA];
pos.cl_trace_st  = [0.555 yA 0.405 hA];
pos.cl_rate      = [L yB 0.415 hB];
pos.cl_psd       = [xC yB R-xC hB];
pos.m_rate       = [L yD 0.140 hD];
pos.m_band       = [0.330 yD 0.230 hD];
pos.m_psd        = [xC yD R-xC hD];
% model traces: no stimulation above stimulation, each with its lambda_1 strip
xT = 0.115; wT = 0.505; hTr = 0.100; hLy = 0.035; gap = 0.045; y0 = 0.040;
pos.m_ly_st  = [xT y0 wT hLy];
pos.m_tr_st  = [xT y0+hLy wT hTr];
pos.m_ly_no  = [xT y0+hLy+hTr+gap wT hLy];
pos.m_tr_no  = [xT y0+2*hLy+hTr+gap wT hTr];
yTop = y0+2*hLy+2*hTr+gap;
xH = 0.770; wH = R-xH; hH = (yTop-y0-0.095)/2;
pos.m_state = [xH y0+0.025 wH hH];
pos.m_corr  = [xH y0+0.025+hH+0.070 wH hH];
if strcmp(opts.layout,'mirror')
    % Rows in inches from the bottom of the canvas, so each gap holds the
    % text it has to (tick labels, x label, title, letter) and no more: the
    % gap between the halves (C's x label, D's title) is drawn as tight as
    % the gaps within them.
    H = opts.size_in(2);
    hR = 1.25/H;                       % rate and PSD rows
    yE = 0.45/H;                       % model rate / PSD
    yD = 2.20/H; hD = 1.30/H;          % model traces (time bar hangs below)
    yB = 4.17/H;                       % clinical rate / PSD
    yA = 5.90/H; hA = 1.52/H;          % clinical traces
    L = 0.115; xC = 0.665; wR = 0.325; wT = 0.385;
    pos.cl_trace_no = [0.125 yA wT hA];
    pos.cl_trace_st = [0.965-wT yA wT hA];
    pos.cl_rate = [L yB wR hR];
    pos.cl_psd  = [xC yB R-xC hR];
    pos.m_tr_no = [0.135 yD wT-0.010 hD];
    pos.m_tr_st = [0.925-(wT-0.010) yD wT-0.010 hD];   % room for the amplitude bar's label
    pos.m_rate  = [L yE wR hR];
    pos.m_psd   = [xC yE R-xC hR];
end
% letter x per panel: left column at the canvas edge, others left of the y label
xl = struct('A',0.005,'B',0.005,'C',xC-0.085,'D',0.005,'E',0.330-0.085,'F',xC-0.085, ...
    'G',0.005,'H',xH-0.10,'I',xH-0.10);

% ---- clinical ---------------------------------------------------------------
sf = openfig(src.traces,'invisible'); guard1 = onCleanup(@()close_if(sf));
drop_texts(sf,opts.drop_text);
ax = findall(sf,'Type','axes');
a_no = place(pick(ax,@(a)strcmp(strtrim(char(a.Title.String)),'No Stimulation')),pos.cl_trace_no);
a_st = place(pick(ax,@(a)strcmp(strtrim(char(a.Title.String)),'Stimulation')),pos.cl_trace_st);
close_if(sf);

sf = openfig(src.rate,'invisible');
[a_rate,lg_rate] = place_with_legend(sf,pos.cl_rate);
close_if(sf);

sf = openfig(src.psd,'invisible');
[a_psd,lg_psd] = place_with_legend(sf,pos.cl_psd);
close_if(sf);

% ---- model ------------------------------------------------------------------
sf = openfig(src.model,'invisible'); guard2 = onCleanup(@()close_if(sf));
delete(findall(sf,'Type','text','Tag','letter'));
ax = findall(sf,'Type','axes');
ylab = @(a) strjoin(cellstr(a.YLabel.String),' ');
title_of = @(a) strtrim(char(a.Title.String));
lg_mrate = [];
if strcmp(opts.layout,'mirror')
    [m_rate,lg_mrate] = place(pick(ax,@(a)contains(ylab(a),'Events per Minute') || contains(ylab(a),'Event Rate')),pos.m_rate,sf);
    [m_psd,lg_mpsd] = place(pick(ax,@(a)contains(ylab(a),'PSD')),pos.m_psd,sf);
    m_tr_no = place(pick(ax,@(a)strcmp(title_of(a),'No Stimulation')),pos.m_tr_no);
    m_tr_st = place(pick(ax,@(a)startsWith(title_of(a),'Stimulation')),pos.m_tr_st);
    close_if(sf);
    mirror_finish();
else
m_rate  = place(pick(ax,@(a)contains(ylab(a),'Events per Minute')),pos.m_rate);
[m_psd,lg_mpsd] = place(pick(ax,@(a)contains(ylab(a),'PSD')),pos.m_psd,sf);
m_band  = place(pick(ax,@(a)contains(ylab(a),'Power Change')),pos.m_band);
m_corr  = place(pick(ax,@(a)contains(ylab(a),'Correlation Change')),pos.m_corr);
m_state = place(pick(ax,@(a)contains(ylab(a),'Change with Stimulation')),pos.m_state);
title_of = @(a) strtrim(char(a.Title.String));
m_tr_no = place(pick(ax,@(a)strcmp(title_of(a),'No Stimulation')),pos.m_tr_no);
m_tr_st = place(pick(ax,@(a)startsWith(title_of(a),'Stimulation')),pos.m_tr_st);
is_ly = @(a) contains(ylab(a),'lambda_1');
m_ly_no = place(pick(ax,@(a)is_ly(a) && line_near(a,purple)),pos.m_ly_no);
m_ly_st = place(pick(ax,@(a)is_ly(a) && line_near(a,green)),pos.m_ly_st);
close_if(sf);

% ---- shorter y labels: the combined canvas gives each panel less height
relabel = {a_rate,'IED Rate (% of No Stim.)'; m_rate,'Events per Minute'; ...
    m_band,'Power Change (dB)'; m_psd,'PSD (a.u.^2/Hz)'; ...
    m_corr,'Correlation Change'; m_state,'Change (%)'};
for k = 1:size(relabel,1)
    h = relabel{k,1};   % through a variable: c{k}.Prop = v would overwrite the cell
    h.YLabel.String = relabel{k,2};
end

% ---- house style ------------------------------------------------------------
all_ax = findall(fig,'Type','axes');
set(all_ax,'FontName','Helvetica','FontSize',opts.font_size,'TickDir','out', ...
    'LineWidth',1,'Box','off','TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
set(findall(fig,'Type','text'),'FontName','Helvetica');
set(all_ax,'LineWidth',0.6);
% the sources were drawn on canvases two to three times larger: thin the
% marks in proportion so they keep their weight relative to the text
for h = findall(fig,'Type','line')'
    h.LineWidth = max(0.5,0.6*h.LineWidth); h.MarkerSize = max(2.5,0.6*h.MarkerSize);
end
for h = findall(fig,'Type','scatter')', h.SizeData = 0.35*h.SizeData; end
for h = findall(fig,'Type','constantline')', h.LineWidth = max(0.5,0.6*h.LineWidth); end
for a = all_ax'
    a.Title.FontWeight = 'normal';
end
% annotation texts inside axes (channel names, scale bars) a size smaller
for t = findall(fig,'Type','text')'
    if ~isempty(t.String) && ~is_label_of(t)
        t.FontSize = opts.annot_size;
    end
end
% the clinical condition labels sit in the traces axes as text; keep them title size
for t = findall([a_no a_st],'Type','text')'
    if any(strcmp(t.String,{'No Stimulation','Stimulation'})), t.FontSize = opts.font_size; end
end
set(findall(fig,'Type','legend'),'FontSize',opts.annot_size,'FontName','Helvetica','Box','off');
% the clinical PSD axes title is a statistic note, not a header: keep it small
a_psd.Title.FontSize = opts.annot_size;
% its saved limits cut the top of the IQR band at 1-2 Hz: let it fit the data
a_psd.YLimMode = 'auto';
% model trace row names and lambda_1 strips
m_tr_no.Title.String = 'Model, No Stimulation'; m_tr_st.Title.String = strrep(title_of(m_tr_st),'Stimulation','Model, Stimulation');
m_ly_no.XLabel.String = ''; m_ly_no.XTickLabel = {};
% one-line lambda_1 label: the two-line source label runs into the row names
set([m_ly_no.YLabel m_ly_st.YLabel],'String','\lambda_1 (s^{-1})');
% the narrow model rate panel: short condition names on one line
m_rate.XTickLabel = {'No Stim.','Stim.'};
% keep the model PSD legend inside, bottom left
if ~isempty(lg_mpsd), lg_mpsd.Location = 'southwest'; end
if ~isempty(lg_psd),  lg_psd.Location  = 'southwest'; end
% clinical rate legend: participants in two columns just right of the axes
if ~isempty(lg_rate)
    lg_rate.Units = 'normalized'; lg_rate.Location = 'none';
    lg_rate.NumColumns = 2;
    p = a_rate.Position; w = lg_rate.Position(3); h = lg_rate.Position(4);
    lg_rate.Position = [p(1)+p(3)+0.008 p(2)+(p(4)-h)/2 w h];
end

% ---- letters ------------------------------------------------------------------
order = {a_no,'A';a_rate,'B';a_psd,'C';m_rate,'D';m_band,'E';m_psd,'F';m_tr_no,'G';m_corr,'H';m_state,'I'};
for k = 1:size(order,1)
    letter(fig,order{k,1},order{k,2},xl.(order{k,2}),opts.letter_size);
end
end   % layout v1
drawnow;

out = struct('fig',fig,'files',{{}},'sources',src,'letters',{panel_map(opts.layout,src.psd)});
if opts.save
    outdir = opts.out_dir;
    assert(~isempty(outdir),'fig_stim_combined:NoOutDir','Give out_dir.');
    if ~isfolder(outdir), mkdir(outdir); end
    for ext = {'.png','.fig','.svg'}
        f = fullfile(outdir,[opts.tag ext{1}]); if isfile(f), delete(f); end
    end
    keep = fig.Position;
    for attempt = 1:3   % see fig_grouped_main: exportgraphics can resize the canvas
        exportgraphics(fig,fullfile(outdir,[opts.tag '.png']),'Resolution',opts.resolution,'Padding',10);
        if isequal(fig.Position(3:4),keep(3:4)), break; end
        fig.Position = keep; drawnow;
    end
    savefig(fig,fullfile(outdir,[opts.tag '.fig']));
    print(fig,fullfile(outdir,[opts.tag '.svg']),'-dsvg','-vector');
    out.files = existing_outputs(outdir,opts.tag);
    write_provenance(fullfile(outdir,[opts.tag '_provenance.md']),src,opts);
end
if ~opts.visible && ~opts.save, return; end
if ~opts.visible, close(fig); out.fig = []; end

% ============================================================================
    function mirror_finish()
        % house style and letters for the mirror layout
        a_rate.YLabel.String = {'IED Rate','(% of No Stim.)'};   % two lines: one runs past the panel
        m_psd.YLabel.String = 'PSD (a.u.^2/Hz)';
        if contains(ylab(m_rate),'Events per Minute'), m_rate.YLabel.String = 'Events per Minute'; end
        axs_all = findall(fig,'Type','axes');
        set(axs_all,'FontName','Helvetica','FontSize',opts.font_size,'TickDir','out', ...
            'Box','off','TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1,'LineWidth',0.6);
        set(findall(fig,'Type','text'),'FontName','Helvetica');
        for hh = findall(fig,'Type','line')'
            hh.LineWidth = max(0.5,0.6*hh.LineWidth); hh.MarkerSize = max(2.5,0.6*hh.MarkerSize);
        end
        for hh = findall(fig,'Type','scatter')', hh.SizeData = 0.35*hh.SizeData; end
        for hh = findall(fig,'Type','constantline')', hh.LineWidth = max(0.5,0.6*hh.LineWidth); end
        for aa = axs_all', aa.Title.FontWeight = 'normal'; end
        for tt = findall(fig,'Type','text')'
            if ~isempty(tt.String) && ~is_label_of(tt), tt.FontSize = opts.annot_size; end
        end
        for tt = findall([a_no a_st],'Type','text')'
            if any(strcmp(tt.String,{'No Stimulation','Stimulation'})), tt.FontSize = opts.font_size; end
        end
        set(findall(fig,'Type','legend'),'FontSize',opts.annot_size,'FontName','Helvetica','Box','off');
        a_psd.Title.FontSize = opts.annot_size; a_psd.YLimMode = 'auto';
        m_psd.Title.FontSize = opts.annot_size;
        s_rate = ylab(m_rate);   % the normalised rate label, split like the clinical one
        if contains(s_rate,' (% of'), m_rate.YLabel.String = strsplit(strrep(s_rate,' (% of','|(% of'),'|'); end
        % model traces drawn like the clinical ones: no axis lines, row names
        % only on the left, time given by a scale bar on the stimulation side
        m_tr_no.Title.String = 'Model, No Stimulation';
        m_tr_st.Title.String = strrep(title_of(m_tr_st),'Stimulation','Model, Stimulation');
        for aa = [m_tr_no m_tr_st]
            aa.XAxis.Visible = 'off'; aa.YAxis.Color = 'none';
            aa.YAxis.TickLabelColor = [0 0 0]; aa.TickLength = [0 0];
        end
        m_tr_st.YTickLabel = {};
        xlm = m_tr_st.XLim; ylm = m_tr_st.YLim; x1 = xlm(2); x0 = x1-opts.time_bar;
        yb = ylm(1)-0.04*diff(ylm);
        line(m_tr_st,[x0 x1],[yb yb],'Color','k','LineWidth',1,'Clipping','off');
        text(m_tr_st,(x0+x1)/2,yb,sprintf('%g s',opts.time_bar),'HorizontalAlignment','center', ...
            'VerticalAlignment','top','FontSize',opts.annot_size,'FontName','Helvetica','Clipping','off');
        if numel(m_rate.XTickLabel)==2 && contains(ylab(m_rate),'Events per Minute')
            m_rate.XTickLabel = {'No Stimulation','Stimulation'};
        end
        if ~isempty(lg_mpsd), lg_mpsd.Location = 'southwest'; end
        if ~isempty(lg_psd),  lg_psd.Location  = 'southwest'; end
        pairs = {a_rate,lg_rate; m_rate,lg_mrate};
        for kk = 1:2
            aa = pairs{kk,1}; lg = pairs{kk,2};
            if isempty(lg), continue; end
            lg.Units = 'normalized'; lg.Location = 'none'; lg.NumColumns = 2;
            lg.ItemTokenSize = [8 8];   % one marker per entry: a long token has no room
            pp = aa.Position; w = lg.Position(3); hh = lg.Position(4);
            lg.Position = [pp(1)+pp(3)+0.008 pp(2)+(pp(4)-hh)/2 w hh];
        end
        ord = {a_no,'A',0.005; a_rate,'B',0.005; a_psd,'C',xC-0.095; ...
            m_tr_no,'D',0.005; m_rate,'E',0.005; m_psd,'F',xC-0.095};
        for kk = 1:size(ord,1), letter(fig,ord{kk,1},ord{kk,2},ord{kk,3},opts.letter_size); end
    end
    function [axn,lgn] = place(axs,p,srcfig)
        % copy one source axes (with its legend, if any) into fig at p;
        % freeze its limits first so the copy keeps the source's ranges
        axs.XLim = axs.XLim; axs.YLim = axs.YLim;
        lg0 = [];
        if nargin > 2
            lg0 = findall(srcfig,'Type','legend');
            lg0 = lg0(arrayfun(@(l) isequal(l.Axes,axs),lg0));
        end
        if isempty(lg0)
            axn = copyobj(axs,fig); lgn = [];
        else
            hc = copyobj([lg0(1) axs],fig); lgn = hc(1); axn = hc(2);
        end
        axn.Units = 'normalized'; axn.PositionConstraint = 'innerposition'; axn.Position = p;
    end
    function [axn,lgn] = place_with_legend(srcfig,p)
        axs = findall(srcfig,'Type','axes');
        assert(isscalar(axs),'fig_stim_combined:Axes','Expected one axes in %s',srcfig.FileName);
        [axn,lgn] = place(axs,p,srcfig);
    end
end

% ---------------------------------------------------------------------------
function a = pick(ax,test)
ok = false(numel(ax),1);
for k = 1:numel(ax), ok(k) = test(ax(k)); end
assert(nnz(ok)==1,'fig_stim_combined:Pick','Expected exactly one matching axes, found %d.',nnz(ok));
a = ax(ok);
end
function tf = line_near(a,col)
L = findobj(a,'Type','line'); tf = false;
for k = 1:numel(L)
    if norm(L(k).Color-col) < 0.05, tf = true; return; end
end
end
function drop_texts(f,patterns)
t = findall(f,'Type','text');
for k = 1:numel(t)
    s = t(k).String; if iscell(s), s = strjoin(s,' '); end
    if any(cellfun(@(p) ~isempty(regexp(char(s),p,'once')),patterns)), delete(t(k)); end
end
end
function tf = is_label_of(t)
a = t.Parent; tf = false;
if ~isa(a,'matlab.graphics.axis.Axes'), return; end
tf = any(t == [a.XLabel a.YLabel a.Title a.Subtitle]);
end
function letter(fig,ax,ch,x,fs)
% "(A)" at the top left of the panel, just above the axes, as the grouped figures do
ip = ax.Position; y = ip(2)+ip(4)+0.006;
w = 0.04*max(1,fs/7.5);   % wide enough that a larger letter does not wrap
annotation(fig,'textbox',[x y w 0.02],'String',['(' ch ')'],'FontSize',fs, ...
    'FontName','Helvetica','LineStyle','none','VerticalAlignment','bottom', ...
    'HorizontalAlignment','left','Margin',0,'Tag','letter');
end
function v = pick_default(tf,a,b)
if tf, v = a; else, v = b; end
end
function close_if(f)
if ~isempty(f) && isvalid(f), close(f); end
end
function s = onoff(b)
if b, s = 'on'; else, s = 'off'; end
end
function m = panel_map(layout,psd_file)
psd_src = 'panels/Fig_Awake_Main_PanelB_PSD';
if nargin > 1
    [d,n] = fileparts(psd_file); [~,d] = fileparts(d); psd_src = [d '/' n];
end
if nargin > 0 && strcmp(layout,'mirror')
    m = {'A','clinical example traces (S9, SOZ), no stimulation and stimulation','clinical composite panel A (panels/Fig_Awake_Main_PanelD)'; ...
         'B','clinical IED rate per participant, % of no stimulation','clinical composite panel B (panels/Fig_Awake_Main_PanelA_Normalized)'; ...
         'C','clinical SOZ PSD, median and IQR, n = 12',['clinical composite panel C (' psd_src ')']; ...
         'D','model example traces, no stimulation and stimulation','model figure, traces axes (titles No Stimulation / Stimulation)'; ...
         'E','model event rate','model figure, rate axes'; ...
         'F','model subnet-mean PSD','model figure, PSD axes'};
    return
end
m = {'A','clinical example traces (S9, SOZ), no stimulation and stimulation','clinical composite panel A (panels/Fig_Awake_Main_PanelD)'; ...
     'B','clinical IED rate per participant, % of no stimulation','clinical composite panel B (panels/Fig_Awake_Main_PanelA_Normalized)'; ...
     'C','clinical SOZ PSD, median and IQR, n = 12',['clinical composite panel C (' psd_src ')']; ...
     'D','model candidate events per minute','fig_ied_stim panel A'; ...
     'E','model band-power change (dB)','fig_ied_stim panel C'; ...
     'F','model subnet-mean PSD','fig_ied_stim panel B'; ...
     'G','model example traces with cumulative lambda_1','fig_ied_stim panel F'; ...
     'H','model correlation change','fig_ied_stim panel D'; ...
     'I','model rate, SFA and resources change (%)','fig_ied_stim panel E'};
end
function write_provenance(file,src,opts)
fid = fopen(file,'w'); guard = onCleanup(@()fclose(fid));
fprintf(fid,'# %s provenance\n\nWritten by `src/figures/fig_stim_combined.m` on %s. Nothing recomputed: each panel is the axes of its source `.fig`, copied, re-positioned and restyled.\n\n', ...
    opts.tag,char(datetime('now','Format','yyyy-MM-dd HH:mm')));
fprintf(fid,'Contains clinical recordings: do not add to git (`git add -f`) or publish outside the manuscript.\n\n');
fprintf(fid,'## Sources\n\n');
fprintf(fid,'- Clinical traces: `%s`\n- Clinical rate: `%s`\n- Clinical PSD: `%s`\n- Model: `%s`\n',src.traces,src.rate,src.psd,src.model);
if ~isempty(opts.clinical_note), fprintf(fid,'- Clinical note: %s\n',opts.clinical_note); end
if ~isempty(opts.model_note),    fprintf(fid,'- Model note: %s\n',opts.model_note); end
fprintf(fid,'\n## Panel letters\n\n| New | Shows | Source panel |\n|---|---|---|\n');
m = panel_map(opts.layout,src.psd);
for k = 1:size(m,1), fprintf(fid,'| (%s) | %s | %s |\n',m{k,:}); end
fprintf(fid,['\n## Canvas\n\nLayout `%s`, %.2f x %.2f in (the printed size), exported at %d dpi. ' ...
    'Panel letters %g pt; axis labels, tick labels and trace titles %g pt; channel names, scale bars, ' ...
    'legends and the PSD notes %g pt.\n'],opts.layout,opts.size_in,opts.resolution,opts.letter_size,opts.font_size,opts.annot_size);
fprintf(fid,'\nDropped: the clinical header text "S9, SOZ Contacts" (the caption carries it) and the sources'' own panel letters.\n');
end
