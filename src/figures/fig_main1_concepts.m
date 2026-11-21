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
% The bottom row (SFA and STD sigmoids and discs, the STD endpoint from the
% route model of the fully adapted condition at raw rate r_ref = 0.25) is
% drawn by draw_adaptation_concepts, shared with fig_main1_concept_networks.
A=draw_adaptation_concepts(fig,preset_name);
sfa_colors=A.sfa_colors; std_colors=A.std_colors; b_levels=A.b_levels;
std_note=A.std_note; frozen_note=A.frozen_note;
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
    std_note, ...
    ['The sigmoid panels are conceptual logistic curves (center0.4; enlarged total SFA shift0.6), ' ...
    frozen_note], ...
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
