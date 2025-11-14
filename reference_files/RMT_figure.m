close all
clear all
clc

rng(3) % reproducibility


set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 20);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
make_imagesc_nans_transparent = false; % if false, NaNs are treated as zeros
n = 250;
scale_n = 1; % scale factor for larger matrices
assert(mod(scale_n * n, 1) == 0, 'scale_n * n must be an integer');
b = 1;
mu = 0;
f = 0.66667;
mu_E = 1; % relative to b
mean_indegree = round(0.5*n); % desired mean number of connections in and out if sparse


G = struct();
%% 1 Base Ginibre ensemble
g = 1;
G(g).rmt = RMT(n, b, mu);
G(g).rmt.description = ['G(0,1)'];

%% Rajan
g = 2;
G(g).rmt = G(1).rmt.copy();
G(g).rmt.set_rajan_means(f, mu_E);
G(g).rmt.description = {['Rajan, f = 2/3, mu_E = 1, mu_I = -2']};

%% make network sparse
g = 3;
G(g).rmt = G(2).rmt.copy();
G(g).rmt.apply_sparsity(mean_indegree);
G(g).rmt.description = ['Sparse connections'];

activation_one = rand(n,1)<0.5;
activation_two = rand(n,1)<0.5;

%% apply sparse activation one
g = 4;
G(g).rmt = G(3).rmt.copy();
G(g).rmt.zero_columns(activation_one);
G(g).rmt.description = ['Sparse Activation One'];

%% homeostasis row zero
g = 5;
G(g).rmt = G(4).rmt.copy();
G(g).rmt.row_sum_to_zero();
row_sum_Z_diff = G(g).rmt.row_sum_Z_diff; % these were the modified weights 
abscissa5 = max(real(eig(G(g).rmt.A))); % abscissa of this matrix
G(g).rmt.description = ['Rows summed to zero'];

%% shifted to EOS, This will be subplot(1,2,1)
g = 6;
G(g).rmt = G(5).rmt.copy();
G(g).rmt.shift_diagonal(-abscissa5)
G(g).rmt.compute_eigenvalues();
G(g).rmt.set_plot_circle(abscissa5, -abscissa5);
G(g).rmt.description = ['Activation One, row-zeroed, and shifted to EOS'];

%% apply sparse activation two, This will be subplot(1,2,2)
g = 7;
G(g).rmt = G(3).rmt.copy(); % reset to 3, before activation one was applied
G(g).rmt.A = G(g).rmt.A + row_sum_Z_diff; % apply the previous homeostatic row-zero matrix
G(g).rmt.zero_columns(activation_two); % zero according to the new activation
G(g).rmt.shift_diagonal(-abscissa5) % shift accordin previous homeostatic shift
G(g).rmt.compute_eigenvalues();
G(g).rmt.set_plot_circle(abscissa5, -abscissa5);
G(g).rmt.description = ['Activation Two'];

%% figure 1, eigenvalue distribution
f1 = figure(1);
% set(f1, 'Position', [100 100 710 1060], 'Color', 'white')
set(f1,'Color', 'white')

% tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
tiledlayout(1, 2, 'TileSpacing', 'compact');

g_to_plot = [6, 7]
ax = gobjects(2, 1);
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g)
    ax(i_g) = nexttile;
    G(g).rmt.plot_circle(ax(i_g));
    hold on
    G(g).rmt.plot_eigenvalue_distribution(ax(i_g));
    hold off
end

linkaxes(ax,'xy')

% Customize axes after linking
for i_g = 1:length(ax)
    axes(ax(i_g));
    x_lim = xlim;
    y_lim = ylim;
    axis off;
    
    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    text(1.35*x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    % text(x_lim(1), y_lim(1), G(i_g).rmt.description, ...
    %     'HorizontalAlignment', 'left', ...
    %     'VerticalAlignment', 'top', ...
    %     'FontWeight', 'normal', ...
    %     'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off
end
axis equal
set(f1, 'Position', [200   200   713   434])

%% figure 2, modified imagesc plot of A
% Collect all A matrix values to determine clims
all_A_values = [];
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    all_A_values = [all_A_values; G(g).rmt.A(:)];
end

% Calculate symmetric clims based on 95% percentile
max_abs_vals = abs(all_A_values);
clim_val = prctile(max_abs_vals, 95);
clims = [-clim_val, clim_val];

% Create custom colormap: blue (negative) -> black (zero) -> red (positive)
n_colors = 256;
half = n_colors / 2;
% Blue to black for negative values
blues = [linspace(0, 0, half)', linspace(0, 0, half)', linspace(1, 0, half)'];
% Black to red for positive values
reds = [linspace(0, 1, half)', linspace(0, 0, half)', linspace(0, 0, half)'];
custom_cmap = [blues; reds];

%% Create concatenated matrix with constant pixel size (1x2 layout)
% Parameters for layout
n_cols = 2;
n_rows = 1;
col_spacing = 50;  % NaN pixels between columns

% Find maximum matrix size
max_size = 0;
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    max_size = max(max_size, size(G(g).rmt.A, 1));
end

% Create cell array to hold padded matrices and track padding amounts
padded_matrices = cell(length(g_to_plot), 1);
pad_before_horiz_array = zeros(length(g_to_plot), 1);
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    A_plot = G(g).rmt.A;
    if make_imagesc_nans_transparent
        A_plot(~G(g).rmt.dense_mask) = NaN;
    else
        A_plot(~G(g).rmt.dense_mask) = 0;
    end
    
    % Pad matrix to max_size
    % Horizontal: centered (left-right)
    % Vertical: aligned to top (no padding above, all padding below)
    current_size = size(A_plot, 1);
    pad_amount = max_size - current_size;
    
    % Horizontal padding (centered)
    pad_before_horiz = floor(pad_amount / 2);
    pad_after_horiz = pad_amount - pad_before_horiz;
    
    % Vertical padding (aligned to top, all padding below)
    pad_before_vert = 0;
    pad_after_vert = pad_amount;
    
    % Store padding amounts for label positioning
    pad_before_horiz_array(i_g) = pad_before_horiz;
    
    padded_matrices{i_g} = NaN(max_size, max_size);
    padded_matrices{i_g}(pad_before_vert+1:pad_before_vert+current_size, ...
                         pad_before_horiz+1:pad_before_horiz+current_size) = A_plot;
end

% Build the concatenated matrix (1 row with 2 columns)
concat_matrix = [];
for i_col = 1:n_cols
    concat_matrix = [concat_matrix, padded_matrices{i_col}];
    
    % Add column spacing (except after last column)
    if i_col < n_cols
        concat_matrix = [concat_matrix, NaN(max_size, col_spacing)];
    end
end

% Plot the concatenated matrix
f2 = figure(2);
set(f2, 'Position', [100 100 2400 1000], 'Color', 'white')
ax2 = axes('Parent', f2);

h = imagesc(ax2, concat_matrix);
set(h, 'AlphaData', ~isnan(concat_matrix));
set(ax2, 'Color', [1 1 1]);

colormap(ax2, custom_cmap);
caxis(ax2, clims);
axis(ax2, 'equal', 'tight');
box off

% % Add ylabel-style labels at appropriate positions
% ylabel_offset = -15;  % Position labels this many pixels from the left edge of actual data
% for i_g = 1:length(g_to_plot)
%     g = g_to_plot(i_g);
%     i_col = i_g;
% 
%     % Calculate position
%     % x position: at the left edge of the non-NaN submatrix
%     x_pos = (i_col - 1) * (max_size + col_spacing) + pad_before_horiz_array(i_g) + ylabel_offset;
%     % y position: at the top of the matrix
%     y_pos = 1;
% 
%     % Add text annotation as ylabel
%     desc = G(g).rmt.description;
%     text(ax2, x_pos, y_pos, desc, ...
%         'HorizontalAlignment', 'right', ...
%         'VerticalAlignment', 'bottom', ...
%         'FontWeight', 'normal', ...
%         'Rotation', 90);
% end

set(ax2, 'Visible', 'off');

%% figure 3 colorbar
% Create separate colorbar figure
% Create modified colormap with thin white line at center for colorbar only
custom_cmap_white = custom_cmap;
% % Replace center 2-3 pixels with white
if make_imagesc_nans_transparent % otherwise don't add the white line
    center_idx = round(size(custom_cmap, 1) / 2);
    white_line_width = 2; % number of pixels for white line
    white_indices = center_idx + (-floor(white_line_width/2):floor(white_line_width/2));
    custom_cmap_white(white_indices, :) = repmat([1 1 1], length(white_indices), 1);
end

show_colorbar_ticks = false; % set to false to hide all ticks

f3 = figure(3);
set(f3, 'Position', [100 100 291 213], 'Color', 'white')
ax3 = axes('Parent', f3);
colormap(ax3, custom_cmap_white);
cb = colorbar(ax3, 'Location', 'west');
caxis(ax3, clims);
set(ax3, 'Visible', 'off');
if show_colorbar_ticks
    set(cb, 'Box', 'off', 'TickLength', 0, 'Ticks', [-1, 0, 1]);
else
    set(cb, 'Box', 'off', 'TickLength', 0, 'Ticks', []);
end
ylabel(cb, 'Connection Strength');

%% Figure 4: Histograms of A matrices (ignoring NaNs)
f4 = figure(4);
set(f4, 'Position', [100 100 800 400], 'Color', 'white')
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax4 = gobjects(length(g_to_plot), 1);
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    ax4(i_g) = nexttile;
    A_plot = G(g).rmt.A;
    A_plot(~G(g).rmt.dense_mask) = NaN;  % Set non-dense elements to NaN
    A_plot(logical(eye(size(A_plot)))) = NaN;  % Exclude diagonal entries
    
    % Check if this RMT has E/I structure defined
    if ~isempty(G(g).rmt.E) && ~isempty(G(g).rmt.I)
        % Separate E and I column values
        A_E_cols = A_plot(:, G(g).rmt.E);  % Excitatory columns
        A_I_cols = A_plot(:, G(g).rmt.I);  % Inhibitory columns
        
        A_E_vals = A_E_cols(:);
        A_E_vals = A_E_vals(~isnan(A_E_vals));  % Remove NaNs
        
        A_I_vals = A_I_cols(:);
        A_I_vals = A_I_vals(~isnan(A_I_vals));  % Remove NaNs
        
        % Plot overlapping histograms with transparency
        hold on;
        histogram(A_E_vals, 50, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % Red for E
        histogram(A_I_vals, 50, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % Blue for I
        hold off;
    else
        % Original gray histogram when no E/I structure
        A_vals = A_plot(:);
        A_vals = A_vals(~isnan(A_vals));  % Remove NaNs
        hold on
        histogram(A_vals, 50, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % overlapping identical histograms to match the colors above
        histogram(A_vals, 50, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        hold off
    end
    
    desc = G(g).rmt.description;
    % title(desc);
    grid off;
    % axis off;
end

% Link x and y axes across all subplots
linkaxes(ax4, 'xy');

%% Figure 5: Row sums as vertical plots
f5 = figure(5);
set(f5, 'Renderer', 'painters');
set(f5, 'Position', [588         197        1200        1040], 'Color', 'white')
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax5 = gobjects(length(g_to_plot), 1);
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    ax5(i_g) = nexttile;
    A_curr = G(g).rmt.A;
    row_sums = sum(A_curr, 2, 'omitnan');  % Sum across columns, ignoring NaNs
    
    % Plot faint red-grey line at zero
    if g == 6
        % For g=6, plot red line first so row_sum appears on top
        % plot([0 0], [1 n], 'Color', [1 0 0], 'LineWidth', 3);
        % hold on
        % plot(row_sums, 1:n, 'Color',[0 0 0]', 'LineWidth', 1.5);  
        stairs([row_sums; row_sums(end)], 1:n+1, 'Color',[0 0 0]', 'LineWidth', 1.5);
        % hold off;
    else
        % plot(row_sums, 1:n, 'Color',[0 0 0], 'LineWidth', 1.5);
        stairs([row_sums; row_sums(end)], 1:n+1, 'Color',[0 0 0], 'LineWidth', 1.5);
        % hold on
        % plot([0 0], [1 n], 'Color', [1 0 0], 'LineWidth', 3);
        % hold off;
    end
    
    xlabel('Row Sum');
    ylabel('Row Index');
    desc = G(g).rmt.description;
    % title(desc);
    grid off;
    set(gca, 'YDir', 'reverse');  % Make row 1 at top
    axis off;
end

% Link x and y axes across all subplots
linkaxes(ax5, 'xy');
% Set the xlims of ax5(1) to 5x its current xlim
xlim(ax5(1), 5*xlim(ax5(1)));

%% Figure 6: horizontal concatenated imagesc plot of [double(activation_one), nan_pad, double(activation_two)] similar to figure 2
% Create concatenated vector with constant spacing (1x2 layout)
n_cols_6 = 2;
col_spacing_6 = 250;  % NaN pixels between columns

% Reshape activation vectors to column vectors and convert to double
act_one_vec = double(not(activation_one(:)))';
act_two_vec = double(not(activation_two(:)))';

% Build the concatenated matrix (n rows, 2 columns with spacing)
concat_activation = [act_one_vec, ones(1, col_spacing_6), act_two_vec];

% Create black/white colormap
bw_cmap = [0 0 0; 1 1 1];  % [black; white]

% Plot the concatenated activation vectors
f6 = figure(6);
set(f6, 'Renderer', 'painters');
set(f6, 'Position', [100 100 1200 65], 'Color', 'white')
ax6 = axes('Parent', f6);

h6 = imagesc(ax6, concat_activation);
% set(h6, 'AlphaData', ~isnan(concat_activation));
set(ax6, 'Color', [1 1 1]);

colormap(ax6, bw_cmap);
caxis(ax6, [0 1]);
axis(ax6, 'equal', 'tight');
box off
set(ax6, 'Visible', 'off');

% save_some_figs_to_folder_2('RMT_figs', 'RMT_example', [], [])