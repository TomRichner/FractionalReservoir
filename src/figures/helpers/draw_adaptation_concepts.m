function A = draw_adaptation_concepts(fig, preset_name, pos, opts)
%DRAW_ADAPTATION_CONCEPTS The four concept diagrams of adaptation (Figure 1): SFA-shifted sigmoids and discs, STD-scaled sigmoids and discs.
%
%   A = draw_adaptation_concepts(fig, preset_name)
%   A = draw_adaptation_concepts(fig, preset_name, pos)
%   A = draw_adaptation_concepts(fig, preset_name, pos, 'arrows', true, 'font_size', 14)
%
%   fig          the figure to draw into
%   preset_name  the preset whose fully adapted condition sets the
%                depression endpoint (full_adaptation_condition)
%   pos          4 x 4, normalized [x y w h] of the SFA sigmoid, the SFA
%                discs, the STD sigmoid and the STD discs; default: the
%                bottom row of fig_main1_concepts
%   arrows       (false) 'Increasing SFA' over the SFA curves with an arrow
%                pointing right, 'Increasing STD' over the STD curves with
%                an arrow pointing down, as in the first concept figure
%                (Fig_FI_curve, commit 1d3e9b4); Figure 1 asks for them
%   font_size    (14) of the arrows' words
%
%   A.axes       1 x 4 axes handles, tags sfa_sigmoid, sfa_discs,
%                std_sigmoid, std_discs
%   A.sfa_colors, A.std_colors, A.b_levels, A.is_chain
%   A.std_note, A.frozen_note   the provenance sentences for the caller's notes
%   A.arrows     the two arrows (figure annotations) and their words (text
%                in the curves' axes); empty without 'arrows'
%
% Moved out of fig_main1_concepts (Plan 2, package B) so that the new
% Figure 1 (fig_main1_concept_networks) draws the same panel without
% copying it; the drawing is unchanged. The depression endpoint at the
% raw reference rate r_ref = 0.25 is the ROUTE MODEL's (steady_state and
% depression handles through the class's route compiler): a product route
% gives prod_m 1/(1 + r_ref tau_rec_m / tau_rel_m), a chain route b1 of
% (C4), the releasable pool. The intermediate levels are illustrative
% frozen factors, not simulated steady states. Font sizes are the
% caller's to set.
%
% See also: fig_main1_concepts, fig_main1_concept_networks, full_adaptation_condition

arguments
    fig
    preset_name (1,:) char
    pos double = [.065 .10 .18 .235; .295 .10 .19 .235; .555 .10 .18 .235; .77 .10 .18 .235]
    opts.arrows (1,1) logical = false     % 'Increasing SFA' and 'Increasing STD' with their arrows
    opts.font_size (1,1) double = 14      % of the arrows' text
end

[~,~,conditions]=srnn_param_preset(preset_name);
full=full_adaptation_condition(conditions);
cond=conditions{find(cellfun(@(c)strcmp(c.name,full),conditions),1)};
std=cond.synapse_config.E.E.std;
r_ref=.25;
stf=[]; if isfield(cond.synapse_config.E.E,'stf'), stf=cond.synapse_config.E.E.stf; end
route=SRNNCellTypePairs.compile_route(struct('std',std,'stf',stf), ...
    struct('npre',1,'label','E->E','err_id_prefix','fig_main1_concepts'));
ss=route.model.steady_state(r_ref,route.P);
b_each=ss.s(:,1:route.P.n_b);   % per-variable steady state: pools (product) or [b1 b2] (chain)
b_min=route.model.depression(b_each,route.P);
is_chain=strcmp(route.P.model,'pool');
assert(b_min>0 && b_min<1);
a_levels=linspace(0,1,5);
b_levels=linspace(1,b_min,5); % intermediate frozen total factors, not trajectories
sfa_colors=a_levels(:)*[.95 .40 .04];
% Paired black -> dark blue -> teal, with the tallest/largest level black.
std_colors=[0 0 0;.04 .22 .65;.08 .40 .85;0 .55 .70;0 .65 .55];
x=linspace(-.6,1.8,450);
phi=@(z)SRNNCellTypePairs.logisticSigmoid(z,.4);
ax1=axes(fig,'Position',pos(1,:),'Tag','sfa_sigmoid'); hold(ax1,'on');
ax2=axes(fig,'Position',pos(2,:),'Tag','sfa_discs'); hold(ax2,'on');
ax3=axes(fig,'Position',pos(3,:),'Tag','std_sigmoid'); hold(ax3,'on');
ax4=axes(fig,'Position',pos(4,:),'Tag','std_discs'); hold(ax4,'on');
theta=linspace(0,2*pi,400);
for k=1:5
    plot(ax1,x,phi(x-.6*a_levels(k)),'Color',sfa_colors(k,:),'LineWidth',2,'Tag',sprintf('level_%d',k));
    center=-.25-1.2*a_levels(k);
    plot(ax2,center+.8*cos(theta),.8*sin(theta),'Color',sfa_colors(k,:),'LineWidth',2,'Tag',sprintf('level_%d',k));
    plot(ax3,x,b_levels(k)*phi(x),'Color',std_colors(k,:),'LineWidth',2,'Tag',sprintf('level_%d',k));
    radius=1.6*b_levels(k);
    plot(ax4,-1+radius*cos(theta),radius*sin(theta),'Color',std_colors(k,:),'LineWidth',2,'Tag',sprintf('level_%d',k));
end
for ax=[ax1 ax3]
    set(ax,'XLim',[x(1) x(end)],'YLim',[0 1.02],'XTick',[0 1],'YTick',[0 1], ...
        'LineWidth',1.0,'FontSize',14,'Box','off');
    ax.XAxis.LineWidth=1.0; ax.YAxis.LineWidth=1.0;
    % Sentence case (polish list, 2026-11-22; they were in title case).
    xlabel(ax,'Dendritic potential','FontSize',14); ylabel(ax,'Synaptic output','FontSize',14);
end
concept_axes(ax2,[-2.5 1],[-1.25 1.25]);
concept_axes(ax4,[-2.9 1],[-1.9 1.9]);
if is_chain
    std_note=sprintf(['STD endpoint (%s, chained vesicle pools): tau_rec=%s s (releasable-pool refill, ' ...
        'reserve refill), tau_rel=%s s, n2=%g, raw reference rate r_ref=%.2f. ' ...
        'Steady state (C4) of route_pool: releasable pool b1=%.9g, reserve b2=%.9g; the depression ' ...
        'factor is b1 alone, B_min=%.9g. ' ...
        'Five displayed frozen factors=%s. Intermediate factors are illustrative, not simulated steady-state points.'], ...
        full,mat2str(std.tau_rec),mat2str(std.tau_rel),std.n2,r_ref,b_each(1),b_each(2),b_min,mat2str(b_levels));
    frozen_note=['not fits to the reference piecewise activation. STD uses a frozen factor B (the releasable ' ...
        'pool b1) multiplying the entire curve, not the self-consistent rate-dependent steady-state ' ...
        'input-output relation. The reserve b2 acts on the output only through b1.'];
else
    std_note=sprintf(['STD endpoint: tau_rec=%s s, tau_rel=%s s, raw reference rate r_ref=%.2f. ' ...
        'Each b_m,ss=1/(1+r_ref*tau_rec,m/tau_rel,m)=%s; total product B_min=%.9g. ' ...
        'Five displayed frozen total factors=%s. Intermediate factors are illustrative, not simulated steady-state points.'], ...
        mat2str(std.tau_rec),mat2str(std.tau_rel),r_ref,mat2str(b_each),b_min,mat2str(b_levels));
    frozen_note=['not fits to the reference piecewise activation. STD uses a frozen product B multiplying the entire curve, ' ...
        'not the self-consistent rate-dependent steady-state input-output relation. ' ...
        'The two-timescale product must not be identified with a single depression variable.'];
end
arrows=gobjects(0);
if opts.arrows
    arrows=draw_arrows(fig,ax1,ax3,opts.font_size);
end
A=struct('axes',[ax1 ax2 ax3 ax4],'sfa_colors',sfa_colors,'std_colors',std_colors,'b_levels',b_levels, ...
    'is_chain',is_chain,'std_note',std_note,'frozen_note',frozen_note,'arrows',arrows);
end

function h=draw_arrows(fig,ax_sfa,ax_std,fs)
% The arrows of the first concept figure (Fig_FI_curve, 2026-07-28,
% commit 1d3e9b4), brought back as they were (Tom, 2026-11-22): over the
% SFA curves 'Increasing SFA' with a brown arrow under it pointing right
% (the curve shifts right), over the STD curves 'Increasing STD' with a
% dark-blue arrow pointing down (the curve is scaled down); black text,
% 'vback2' heads. Placed in data units of each axes. Annotations live in
% the figure: the caller sets the axes' final positions first.
% On the 9 in canvas the open corner left of the curves is narrower than
% the words (they crossed the y axis there), so the words sit on one line
% just above each axes (the caller leaves about 0.4 in for them): the SFA
% arrow under its words, over the tops of the curves; the STD arrow down
% from its words, across the plateaus it scales.
h=gobjects(1,4);
p=d2f(ax_sfa,[-.15 .85],[1.08 1.08]);
h(1)=annotation(fig,'arrow',p.x,p.y,'Color',[.35 .2 0],'HeadStyle','vback2','LineWidth',1,'Tag','arrow_sfa');
h(2)=text(ax_sfa,.35,1.13,'Increasing SFA','FontSize',fs,'Color',[0 0 0],'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','Clipping','off','Tag','text_sfa');
p=d2f(ax_std,[1.25 1.25],[1.08 .27]);
h(3)=annotation(fig,'arrow',p.x,p.y,'Color',[0 .25 .5],'HeadStyle','vback2','LineWidth',1,'Tag','arrow_std');
h(4)=text(ax_std,1.25,1.13,'Increasing STD','FontSize',fs,'Color',[0 0 0],'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','Clipping','off','Tag','text_std');
end

function p=d2f(ax,x,y)
% Data coordinates of ax to normalized figure coordinates (linear axes,
% limits set).
u=ax.Units; ax.Units='normalized'; pos=ax.Position; ax.Units=u;
xl=ax.XLim; yl=ax.YLim;
p.x=pos(1)+(x-xl(1))/diff(xl)*pos(3);
p.y=pos(2)+(y-yl(1))/diff(yl)*pos(4);
end

function concept_axes(ax,xr,yr)
axis(ax,'off'); axis(ax,'equal');
hx=plot(ax,[xr(1) xr(2)-.15],[0 0],'k','LineWidth',1);
hy=plot(ax,[0 0],[yr(1)*.85 yr(2)*.85],'k','LineWidth',1);
uistack([hx hy],'bottom');
text(ax,xr(2)-.10,0,'Re','FontSize',14,'VerticalAlignment','middle');
text(ax,0,yr(2)*.85,'Im','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','bottom');
xlim(ax,xr); ylim(ax,yr);
end
