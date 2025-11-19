function h_plot = plot_J_eff_graph(J, max_weight, color_limits)
%PLOT_J_EFF_GRAPH Plots a directed graph of J_eff with weighted edges and color-coded connections.
%   J: The Jacobian matrix (or J_eff matrix).
%   max_weight: The maximum absolute weight for scaling edge widths.
%   color_limits: [min_val, max_val] for the colormap.

    % Create directed graph
    % Transpose J because digraph(A) creates edges from row i to col j if A(i,j) is non-zero,
    % but typically connectivity matrices are defined such that W(i,j) is connection from j to i.
    % However, let's check standard MATLAB digraph behavior vs the project convention.
    % In create_W_matrix (implied), usually W(i,j) is j -> i. 
    % digraph(A) takes A as adjacency matrix where A(i,j) is edge i -> j.
    % So if J(i,j) represents influence of j on i, we want an edge from j to i.
    % Thus, if we use digraph(J'), then (i,j) entry of J' is J(j,i), which is influence of i on j.
    % Wait.
    % If J(i,j) is partial derivative of f_i wrt x_j (influence of j on i).
    % We want edge j -> i.
    % digraph(A) draws edge i -> j with weight A(i,j).
    % So we want the matrix passed to digraph to have entry (j,i) = J(i,j).
    % That is exactly the transpose of J. 
    % So digraph(J') is correct assuming J(i,j) is influence of j on i.
    
    dgJ = digraph(J');

    % Remove self-loops (diagonal elements)
    dgJ = rmedge(dgJ, 1:numnodes(dgJ), 1:numnodes(dgJ));
    
    % Calculate line widths
    % Scale widths based on max_weight. 
    % Ensure max_weight is not zero to avoid division by zero.
    if max_weight == 0
        max_weight = 1;
    end
    
    % Base width + scaled width
    % Similar to: LWidths = 2.5*abs(dgA.Edges.Weight)/max_weight+1;
    weights = dgJ.Edges.Weight;
    abs_weights = abs(weights);
    LWidths = 1 * abs_weights / max_weight; % Slightly thinner base than before might be good for dense graphs
    
    % Arrow size
    ASize = 10; 

    % Plot graph
    % Nodes are grey. Edges colored by weight.
    h_plot = plot(dgJ, 'Layout', 'circle', ...
        'EdgeLabel', {}, ...
        'LineWidth', LWidths, ...
        'ArrowSize', ASize, ...
        'ArrowPosition', 0.94, ...
        'NodeColor', [0.6 0.6 0.6], ... % Grey nodes
        'MarkerSize', 4, ...
        'NodeLabel', {}, ...
        'LineStyle', '-');
    
    % Color edges by weight
    h_plot.EdgeCData = weights;
    
    % Apply colormap and limits
    colormap(gca, bluewhitered_colormap(256));
    clim(color_limits);
    
    % Clean up appearance
    set(gca, 'Visible', 'off'); % Hide axes
    set(gca, 'XTick', [], 'YTick', []); % Remove ticks
    axis equal;

end

