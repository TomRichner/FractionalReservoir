function [h_digraph, dgA] = plot_network_graph_widthRange_color_R(A,max_weight,EI_vec)
%PLOT_NETWORK_GRAPH Summary of this function goes here
%   Detailed explanation goes here

    dgA = digraph(A');
%     LWidths = 5*abs(dgA.Edges.Weight)/max(abs(dgA.Edges.Weight))+1;
    LWidths = 2.5*abs(dgA.Edges.Weight)/max_weight+1;
%     ASize = 15*0.1/max_weight;
    ASize=12;

    BR = lines(2); % get blue and red 
    B = BR(1,:);
    R = BR(2,:); 
    n_E = sum(EI_vec==1);
    n_I = sum(EI_vec==-1);
    % NodeColor_mat = [repmat(B,n_E,1); repmat(R,n_I,1)];

    inhibitory_node_color = [0.7 0 0]; % dark red
    NodeColor_mat = [lines(n_E); repmat(inhibitory_node_color,n_I,1)];

% %     h_digraph = plot(dgA,'EdgeLabel',round(dgA.Edges.Weight*100)/100,'LineWidth',LWidths,'ArrowSize',15,'ArrowPosition',0.92,'NodeColor', lines(length(A)),'MarkerSize',18,'NodeLabel',{});
%     h_digraph = plot(dgA,'Layout','circle','EdgeLabel',{},'LineWidth',LWidths,'ArrowSize',ASize,'ArrowPosition',0.94,'NodeColor', lines(length(A)),'MarkerSize',18,'NodeLabel',{},'LineStyle','-');
    h_digraph = plot(dgA,'Layout','circle','EdgeLabel',{},'LineWidth',LWidths,'ArrowSize',ASize,'ArrowPosition',0.94,'NodeColor', NodeColor_mat,'MarkerSize',15,'NodeLabel',{},'LineStyle','-');

%     h_digraph = plot(dgA,'Layout','circle','EdgeLabel',{},'LineWidth',LWidths,'ArrowSize',15,'ArrowPosition',0.94,'NodeColor', lines(length(A)),'MarkerSize',18,'NodeLabel',{});


    h_digraph.NodeFontSize = 11;
    h_digraph.EdgeFontSize = 9;
    
    inhibitory_graph = rmedge(dgA,find(dgA.Edges.Weight>0));
    excitatory_graph = rmedge(dgA,find(dgA.Edges.Weight<0));
    highlight(h_digraph,inhibitory_graph,'EdgeColor','r','LineStyle','-')
%     highlight(h_digraph,excitatory_graph,'EdgeColor',[0 0.65 0.35])
    
%     title('Network directed graph')

    set (gca, 'xdir', 'reverse' )

end

