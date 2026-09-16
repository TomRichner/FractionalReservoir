function out = fig_main1_concepts(preset_name, intro_root)
% FIG_MAIN1_CONCEPTS Native archived intro plus four analytic mechanism panels.
% Reads saved MATLAB graphics only: never integrates or estimates an exponent.
% intro_root names the source explicitly. Older split .figs and newer combined
% .figs hold the same Sompolinsky illustration, independently of the E/I run.
arguments
    preset_name char
    intro_root char
end
root=fileparts(which('setup_paths'));
if ~startsWith(intro_root,'/') && isempty(regexp(intro_root,'^[A-Za-z]:','once'))
    intro_root=fullfile(root,intro_root);
end
fig=figure('Visible','off','Color','w','Position',[40 40 910 665]);
combined=fullfile(intro_root,'Fig_Intro_Concepts.fig');
if isfile(combined)
    old=openfig(combined,'invisible'); guard=onCleanup(@()close(old));
    aa=findall(old,'Type','axes'); assert(numel(aa)==6);
    pos=vertcat(aa.Position); [~,ix]=sort(pos(:,2),'descend');
    eig_ax=left_to_right(aa(ix(1:3))); trace_ax=left_to_right(aa(ix(4:6)));
    sources={combined};
else
    fe=fullfile(intro_root,'eigenspectra','panelA_eigenspectrum.fig');
    ft=fullfile(intro_root,'statetraces','panelA_bottom_traces.fig');
    assert(isfile(fe)&&isfile(ft),'fig_main1_concepts:MissingNativeIntro', ...
        'The explicitly selected native intro source is missing: %s',intro_root);
    olde=openfig(fe,'invisible'); oldt=openfig(ft,'invisible');
    guard=onCleanup(@()close([olde oldt]));
    eig_ax=left_to_right(findall(olde,'Type','axes'));
    trace_ax=left_to_right(findall(oldt,'Type','axes'));
    assert(numel(eig_ax)==3 && numel(trace_ax)==3);
    sources={fe,ft};
end
for k=1:3
    a=copyobj(eig_ax(k),fig); a.Units='normalized';
    a.Position=[.055+(k-1)*.315 .705 .265 .255];
    a.PositionConstraint='innerposition'; a.Tag=sprintf('intro_eigen_%d',k);
    a=copyobj(trace_ax(k),fig); a.Units='normalized';
    a.Position=[.065+(k-1)*.315 .465 .265 .225];
    a.PositionConstraint='innerposition'; a.Tag=sprintf('intro_trace_%d',k);
    a.LineWidth=1.0; a.YAxis.LineWidth=1.0;
    % Remove only the archived 10-unit scale annotation, never a trajectory.
    tt=findall(a,'Type','text');
    for j=1:numel(tt)
        if strcmp(strtrim(string(tt(j).String)),'10 s'), delete(tt(j)); end
    end
    ll=findall(a,'Type','line');
    for j=1:numel(ll)
        if numel(ll(j).XData)==2 && ll(j).LineWidth==4 && ...
                abs(diff(ll(j).XData)-10)<1e-10 && diff(ll(j).YData)==0
            delete(ll(j));
        end
    end
end
clear guard
% Source illustration values are untouched: only copied graphics are styled.
% STD endpoint from the actual reference MTS route at raw rate r_ref=0.25.
[~,~,conditions]=srnn_param_preset(preset_name);
full=full_adaptation_condition(conditions);
cond=conditions{find(cellfun(@(c)strcmp(c.name,full),conditions),1)};
std=cond.synapse_config.E.E.std;
r_ref=.25;
b_each=1./(1+r_ref.*std.tau_rec./std.tau_rel);
b_min=prod(b_each);
assert(b_min>0 && b_min<1);
a_levels=linspace(0,1,5);
b_levels=linspace(1,b_min,5); % intermediate frozen total factors, not trajectories
sfa_colors=a_levels(:)*[.95 .40 .04];
% Paired black -> dark blue -> teal, with the tallest/largest level black.
std_colors=[0 0 0;.04 .22 .65;.08 .40 .85;0 .55 .70;0 .65 .55];
x=linspace(-.6,1.8,450);
phi=@(z)SRNNCellTypePairs.logisticSigmoid(z,.4);
ax1=axes(fig,'Position',[.065 .10 .18 .235],'Tag','sfa_sigmoid'); hold(ax1,'on');
ax2=axes(fig,'Position',[.295 .10 .19 .235],'Tag','sfa_discs'); hold(ax2,'on');
ax3=axes(fig,'Position',[.555 .10 .18 .235],'Tag','std_sigmoid'); hold(ax3,'on');
ax4=axes(fig,'Position',[.77 .10 .18 .235],'Tag','std_discs'); hold(ax4,'on');
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
    xlabel(ax,'Dendritic potential','FontSize',14); ylabel(ax,'synaptic output','FontSize',14);
end
concept_axes(ax2,[-2.5 1],[-1.25 1.25]);
concept_axes(ax4,[-2.9 1],[-1.9 1.9]);
row_labels=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1], ...
    'Visible','off','Tag','intro_row_labels');
for j=1:3
    row_y=[.955 .685 .39];
    text(row_labels,.015,row_y(j),sprintf('(%c)',char('A'+j-1)), ...
        'FontSize',14,'VerticalAlignment','top','Tag',sprintf('row_label_%d',j));
end
set(findall(fig,'-property','FontSize'),'FontSize',14);
notes={['For the current mu7 figure-only configuration, the native introductory source is explicitly ' ...
    'figs/sfaEI_fast/fig_introductory_concepts (clean source commit4406b56). It uses the same ' ...
    'sompolinsky_pairs preset, gammas[0.9,1.6,2.5], seed0, 15 traces and [0,60]s as the mu7 intro. ' ...
    'The intervening fig_introductory_concepts code change only combined the two rows; simulation code is unchanged.'], ...
    sprintf(['STD endpoint: tau_rec=%s s, tau_rel=%s s, raw reference rate r_ref=%.2f. ' ...
    'Each b_m,ss=1/(1+r_ref*tau_rec,m/tau_rel,m)=%s; total product B_min=%.9g. ' ...
    'Five displayed frozen total factors=%s. Intermediate factors are illustrative, not simulated steady-state points.'], ...
    mat2str(std.tau_rec),mat2str(std.tau_rel),r_ref,mat2str(b_each),b_min,mat2str(b_levels)), ...
    ['The sigmoid panels are conceptual logistic curves (center0.4; enlarged total SFA shift0.6), ' ...
    'not fits to the reference piecewise activation. STD uses a frozen product B multiplying the entire curve, ' ...
    'not the self-consistent rate-dependent steady-state input-output relation. ' ...
    'The two-timescale product must not be identified with a single depression variable.'], ...
    ['Each SFA curve and shifted disk shares the exact black-to-orange color; each STD curve and shrinking disk ' ...
    'shares the exact black-through-dark-blue-to-teal color. Disk shifts/radii are effective-connectivity intuition, not exact ' ...
    'transformations of the full active Jacobian.'], ...
    ['Top panels are copied native saved graphics; all scientific XData/YData, neuron selections and gains are unchanged. The10-unit scale bar and its text are omitted because the example time is arbitrary up to rescaling. ' ...
    'Only position, row labels, 14-point fonts and trace y-axis width1.0 are changed. No simulation or Lyapunov calculation is run.']};
out=struct('figs',fig,'files',{{}},'source',{sources},'notes',{notes}, ...
    'sfa_colors',sfa_colors,'std_colors',std_colors,'b_levels',b_levels);
end
function aa=left_to_right(aa)
pos=vertcat(aa.Position); [~,ix]=sort(pos(:,1)); aa=aa(ix);
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
