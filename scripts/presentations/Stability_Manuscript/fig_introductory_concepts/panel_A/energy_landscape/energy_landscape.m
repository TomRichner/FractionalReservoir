%% energy_landscape.m
% Figure 1, Panel A (top): effective potential ("energy landscape") for the
% 1-D reduction of the network model, at the same three gains used in Panel A
% bottom (traces) and the eigenspectra.
%
% For a single self-coupled unit with a tanh nonlinearity the dynamics are a
% gradient flow:
%
%     dx/dt = -x + gamma*tanh(x) = -dU/dx,
%     U(x)  = x^2/2 - gamma*log(cosh(x)).
%
% The curvature at the origin is U''(0) = 1 - gamma, so the fixed point x = 0
% is a MINIMUM (stable) for gamma < 1 and a MAXIMUM (unstable) for gamma > 1.
% This is the same gamma = 1 transition shown by the eigenspectra (the disk
% crossing Re = 0) and by the time series (relaxation vs chaos).
%
% Arrows show the flow direction dx/dt = -dU/dx (always downhill):
%   gamma < 1 : arrows converge on x = 0   -> stable, ball settles
%   gamma > 1 : arrows diverge from x = 0  -> unstable, ball rolls away
%
% The ball is drawn at x = 0 in every panel -- the same fixed point that the
% eigenspectrum linearizes about -- so the story is told by the change in the
% shape of the landscape underneath it, not by moving the ball.

%% Setup
close all; clear; clc;
setup_paths();

this_dir = fileparts(mfilename('fullpath'));

%% ---- Editable parameters -------------------------------------------------
gammas      = [0.9, 1.6, 2.5];   % must match panelA_bottom_traces.m
x_lim       = 2.6;               % landscape drawn over x in [-x_lim, x_lim]
                                 % 2.6 keeps all three curves in frame under
                                 % common y-limits and just clears the
                                 % gamma = 2.5 minima at x* ~ 2.46
ylim_mode   = 'common_span';     % 'common_span' -> same y-SPAN in every panel (so well
                                 %                  depths stay directly comparable and
                                 %                  gamma=1.6 stays visibly shallow), but
                                 %                  each panel gets its own vertical
                                 %                  CENTER. The ball therefore sits low in
                                 %                  the frame when trapped (gamma<1) and
                                 %                  high when perched on a hilltop
                                 %                  (gamma>1) -- left shifts down, right
                                 %                  shifts up.
                                 % 'common'      -> identical y-limits (ball at the same
                                 %                  height in all panels, since U(0)=0)
                                 % 'auto'        -> per-panel autoscale (fills each tile,
                                 %                  but exaggerates the gamma=1.6 well)
ylim_pad    = 0.12;              % fractional padding added to the common span
show_titles = false;             % gamma titles above each panel
show_ulabel = false;             % "U(x)" label on the left-most panel
show_ghosts = true;              % faint balls at the stable minima (gamma > 1)
show_arrows = true;              % flow-direction arrows
ball_size   = 15;                % marker size (diameter, points) for the ball
head_len    = 9;                 % arrowhead length (points)
ball_gap    = 4;                 % clear space between an arrow tip and a ball (points)
curve_lw    = 2.5;               % landscape line width
fig_width   = 1200;              % match panelA_bottom_traces.m for column alignment
fig_height  = 300;
% --------------------------------------------------------------------------

n_cases = numel(gammas);
U   = @(x, g) x.^2 / 2 - g * log(cosh(x));   % potential
xq  = linspace(-x_lim, x_lim, 1200);

% Per-panel extent of the landscape, and a single common span
U_lo   = arrayfun(@(g) min(U(xq, g)), gammas);
U_hi   = arrayfun(@(g) max(U(xq, g)), gammas);
span   = max(U_hi - U_lo) * (1 + 2 * ylim_pad);   % same vertical scale everywhere
y_lims = [min(U_lo) - ylim_pad * span, max(U_hi) + ylim_pad * span];  % for 'common'

%% Figure
fig = figure('Color', 'white', 'Position', [100, 300, fig_width, fig_height]);
tl  = tiledlayout(1, n_cases, 'TileSpacing', 'compact', 'Padding', 'compact'); %#ok<NASGU>
ax      = gobjects(1, n_cases);
x_stars = zeros(1, n_cases);   % nonzero stable equilibria (0 when gamma <= 1)

for k = 1:n_cases
    g = gammas(k);
    ax(k) = nexttile; hold(ax(k), 'on');

    % --- landscape ---
    plot(ax(k), xq, U(xq, g), 'k-', 'LineWidth', curve_lw);

    % --- stable minima (nonzero roots of x = gamma*tanh(x), exist for gamma>1) ---
    if g > 1
        x_star = fzero(@(x) x - g * tanh(x), [1e-3, 10 * g]);
    else
        x_star = 0;
    end

    x_stars(k) = x_star;   % arrows and balls are drawn in a second pass below

    xlim(ax(k), [-x_lim, x_lim]);
    switch lower(ylim_mode)
        case 'common_span'
            % same span, per-panel center -> ball drifts down (trapped) to up (perched)
            ctr = (U_hi(k) + U_lo(k)) / 2;
            ylim(ax(k), ctr + [-span, span] / 2);
        case 'common'
            ylim(ax(k), y_lims);
        case 'auto'
            padk = 0.18 * (U_hi(k) - U_lo(k));
            ylim(ax(k), [U_lo(k) - padk, U_hi(k) + padk]);
        otherwise
            error('Unknown ylim_mode: %s', ylim_mode);
    end
    axis(ax(k), 'off');

    if show_titles
        title(ax(k), sprintf('\\gamma = %.2f', g), 'FontWeight', 'normal', ...
            'FontSize', 13, 'Visible', 'on');
    end
end

%% Arrows and balls -- both ride one ball-radius above the surface
% Every ball sits at a critical point of U (the red ball at x = 0, the ghosts at
% the minima), so the surface is locally horizontal there and "tangent" simply
% means lifting the ball straight up by one radius. The arrows are offset along
% the surface NORMAL by the same radius, so they clear the landscape by a
% constant distance even where it is steep. MarkerSize is a diameter in POINTS,
% so all of this needs the final axes geometry -- hence this second pass.
drawnow;   % let tiledlayout settle so the axes heights are final

R_pt = ball_size / 2;                        % ball radius, in points
dU   = @(x, g) x - g * tanh(x);              % dU/dx (arrow rides its slope)

for k = 1:n_cases
    g = gammas(k);
    [sx, sy] = data_per_point(ax(k));
    R = R_pt * sy;                            % ball radius in data-y units

    % Horizontal clearance a shaft must stop short of a ball centre, so that the
    % ARROWHEAD TIP (which extends head_len beyond the shaft) still leaves
    % ball_gap of clear space. Expressed in points, so it adapts per panel --
    % this is what keeps the closely-spaced gamma = 1.6 minima from being hit.
    clear_x = (R_pt + head_len + ball_gap) * sx;
    hold(ax(k), 'on');

    % --- flow-direction arrows: dx/dt = -dU/dx (always downhill) ---
    if show_arrows
        if g > 1
            % unstable origin: flow diverges away from x = 0 toward the minima.
            % Tail clears the red ball at x = 0 (no head at that end, so no
            % head_len term); tip clears the ghost ball at x = x_star.
            x_beg = max(0.30 * x_stars(k), (R_pt + ball_gap) * sx);
            x_end = x_stars(k) - clear_x;
            draw_flow_arrow(ax(k),  x_beg,  x_end, U, dU, g, R_pt, head_len);
            draw_flow_arrow(ax(k), -x_beg, -x_end, U, dU, g, R_pt, head_len);
        else
            % stable origin: flow converges on x = 0
            x_end = max(0.28 * x_lim, clear_x);                % clear the red ball
            draw_flow_arrow(ax(k),  0.90 * x_lim,  x_end, U, dU, g, R_pt, head_len);
            draw_flow_arrow(ax(k), -0.90 * x_lim, -x_end, U, dU, g, R_pt, head_len);
        end
    end

    % --- ghost balls at the stable minima (where the ball ends up) ---
    if show_ghosts && g > 1
        xs = [-x_stars(k), x_stars(k)];
        plot(ax(k), xs, U(xs, g) + R, 'o', ...
            'MarkerSize', ball_size, 'MarkerFaceColor', [0.85 0.85 0.85], ...
            'MarkerEdgeColor', [0.45 0.45 0.45], 'LineWidth', 1.0);
    end

    % --- the ball, resting on the fixed point x = 0 ---
    plot(ax(k), 0, U(0, g) + R, 'o', 'MarkerSize', ball_size, ...
        'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);

    hold(ax(k), 'off');
end

% "U(x)" label on the left-most panel
if show_ulabel
    yl1 = ylim(ax(1));
    text(ax(1), -x_lim, yl1(2), 'U(x)', 'FontSize', 13, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

%% Save (fig, png, svg via repo helper)
save_some_figs_to_folder_2(this_dir, 'panelA_energy_landscape', fig.Number, {'fig', 'png', 'svg'});

%% Local functions
function [sx, sy] = data_per_point(ax)
% DATA_PER_POINT Data units per typographic point, separately in x and y.
% Reading the axes size in POINTS -- the same unit MarkerSize and LineWidth use
% -- converts between screen and data units with no DPI assumption, and the
% ratio survives uniform export scaling.
old_units = ax.Units;
ax.Units  = 'points';
pos       = ax.Position;
ax.Units  = old_units;
sx = diff(xlim(ax)) / pos(3);
sy = diff(ylim(ax)) / pos(4);
end

function draw_flow_arrow(ax, x_tail, x_head, U, dU, g, off_pt, head_pt)
% DRAW_FLOW_ARROW Arrow riding parallel to the landscape from x_tail to x_head,
% showing the gradient flow dx/dt = -dU/dx (always downhill).
%
% The shaft is offset from the surface along the local NORMAL by off_pt points,
% so it clears the curve by a constant distance even where the slope is steep.
% The head is a triangular patch rotated to the tangent direction at the tip,
% with the angle computed in SCREEN space so it looks correct despite the x and
% y axes having different scales.
col = [0.20 0.35 0.75];
[sx, sy] = data_per_point(ax);

xs = linspace(x_tail, x_head, 60);
ys = U(xs, g);

% Unit normal to the curve, computed in screen (point) space, then mapped back
m_screen = dU(xs, g) * sx / sy;              % slope in points-per-point
nrm      = sqrt(1 + m_screen.^2);
nx_pt    = -m_screen ./ nrm;                 % upward-ish unit normal, in points
ny_pt    =  1        ./ nrm;
xs_off   = xs + off_pt * nx_pt * sx;
ys_off   = ys + off_pt * ny_pt * sy;

% Truncate the shaft so it terminates at the MIDDLE of the arrowhead rather
% than at the tip. Running it to the tip leaves the full line width visible
% where the triangle tapers to a point, which reads as a blunt, thickened tip;
% ending it half a head back puts the join where the triangle is wide enough to
% hide it. Arc length is measured in points, so the trim is scale-correct.
seg_pt = hypot(diff(xs_off) / sx, diff(ys_off) / sy);
s_cum  = [0, cumsum(seg_pt)];
s_stop = max(s_cum(1), s_cum(end) - head_pt / 2);
keep   = s_cum <= s_stop;
x_shaft = [xs_off(keep), interp1(s_cum, xs_off, s_stop)];
y_shaft = [ys_off(keep), interp1(s_cum, ys_off, s_stop)];
plot(ax, x_shaft, y_shaft, '-', 'Color', col, 'LineWidth', 1.8);

% --- arrowhead: triangle rotated to the direction of travel ---
s      = sign(x_head - x_tail);
dir_pt = [s / sx, s * dU(x_head, g) / sy];   % tangent direction, in points
theta  = atan2(dir_pt(2), dir_pt(1));

L = head_pt;                                  % head length (points)
W = head_pt * 0.60;                           % head width  (points)
vx = [0, -L, -L];                             % triangle pointing along +x
vy = [0,  W/2, -W/2];
rx = vx * cos(theta) - vy * sin(theta);       % rotate in point space
ry = vx * sin(theta) + vy * cos(theta);

patch(ax, 'XData', xs_off(end) + rx * sx, 'YData', ys_off(end) + ry * sy, ...
    'FaceColor', col, 'EdgeColor', 'none');
end
