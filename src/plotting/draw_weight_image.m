function h = draw_weight_image(ax, W, type, clim_val, opts)
%DRAW_WEIGHT_IMAGE A small weight matrix as an image: excitatory (positive) entries red, inhibitory (negative) blue, absent white.
%
%   h = draw_weight_image(ax, W, type, clim_val)
%   h = draw_weight_image(ax, W, type, clim_val, 'divider', true, 'line_width', 1)
%
%   W         m x m, (post, pre): columns are the presynaptic neurons, so
%             under Dale's law a column is all red (E) or all blue (I)
%   type      m x 1, 1 E and 2 I (weight_image_sample's S.type), used for
%             the optional divider between the E and the I columns
%   clim_val  the color range [-clim_val, clim_val]; pass the same value to
%             every panel that is compared. A larger |w| saturates.
%
% The colormap is redwhiteblue_colormap (as SRNNCellTypePairs.plot_W and
% the first paper's Figure 1, ConnectivityAdaptation's
% Fig_1_RMT_examples), so zero is exactly white. Square pixels, no ticks,
% a thin frame. Options: divider (false) draws a thin gray line between
% the last E column and the first I column; line_width (1) of the frame.
%
% See also: weight_image_sample, redwhiteblue_colormap, SRNNCellTypePairs.plot_W

arguments
    ax
    W (:,:) double
    type (:,1) double
    clim_val (1,1) double {mustBePositive}
    opts.divider (1,1) logical = false
    opts.line_width (1,1) double = 1
end
m = size(W, 1);
h = imagesc(ax, W);
colormap(ax, redwhiteblue_colormap(256));
clim(ax, [-clim_val, clim_val]);
axis(ax, 'image');
set(ax, 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', opts.line_width, 'XColor', 'k', 'YColor', 'k', ...
    'TickLength', [0 0]);
if opts.divider
    nE = nnz(type == 1);
    hold(ax, 'on');
    plot(ax, [nE nE] + 0.5, [0.5 m + 0.5], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.75);
    hold(ax, 'off');
end
end
