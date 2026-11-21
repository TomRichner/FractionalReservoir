% TEST_COMPLEX_PLANE_AXES The parts of complex_plane_axes that can be checked without looking.
%
% complex_plane_axes draws the axes of a complex-plane eigenvalue plot through
% the origin, at 1:1, with an optional break in the real axis, and returns the
% mapping the caller must plot through. What is tested here:
%
%   - the mapping from a real part to the plotted x, and its inverse
%   - that a value inside the removed stretch is an ERROR (it must not vanish,
%     and must not be drawn somewhere it is not), and is counted by cp.check
%   - that no tick falls inside the interval kept clear for the bulk, nor at
%     the origin, nor in the removed stretch
%   - that the plotted aspect is 1:1 (measured on the drawn plot box, in
%     pixels), and that two panels given the same options have the same
%     limits and the same plot box size
%
% How it LOOKS (collisions of labels, whether the break reads as a break) is
% not tested: look at figs/explorations/outliers/outliers_complex_plane_v2.png.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
% Opens one invisible figure and closes it.
%
% See also: complex_plane_axes, plane_fig_outliers

all_passed = true;
fprintf('\n=== test_complex_plane_axes ===\n');
args = {'re_lim', [-32 12], 'im_lim', [-22 22], 'far_lim', [28 42], 'gap', 3, ...
    'bulk_center', -10, 'bulk_radius', 17.5};

%% ------------------------------------------------------------------------
fprintf('\n-- the mapping and its inverse --\n');
cp = complex_plane_axes([], args{:});
re = [-32 -10 0 12 28 35 42];
all_passed = check('main window: plotted x is the real part', ...
    isequal(cp.x([-32 -10 0 12]), [-32 -10 0 12])) && all_passed;
all_passed = check('far segment: starts one gap after the main window, same scale', ...
    isequal(cp.x([28 35 42]), [15 22 29])) && all_passed;
all_passed = check('the inverse undoes the mapping', ...
    isequal(cp.x_inv(cp.x(re)), re)) && all_passed;
all_passed = check('the inverse is NaN inside the gap', ...
    all(isnan(cp.x_inv([12.5 14.9]))) && ~any(isnan(cp.x_inv([12 15])))) && all_passed;
all_passed = check('shape is kept and NaN passes through', ...
    isequaln(cp.x([1 NaN; 30 2]), [1 NaN; 17 2])) && all_passed;
all_passed = check('plotted limits are main window + gap + far segment', ...
    isequal(cp.x_lim, [-32 29]) && isequal(cp.y_lim, [-22 22])) && all_passed;
nb = complex_plane_axes([], 're_lim', [-32 12], 'im_lim', [-22 22]);
all_passed = check('no break: the mapping is the identity', ...
    isequal(nb.x([-40 0 100]), [-40 0 100]) && isequal(nb.x_lim, [-32 12])) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- a point inside the removed stretch --\n');
all_passed = check('cp.x errors with complex_plane_axes:InBreak', ...
    strcmp(error_id(@() cp.x([0 20 35])), 'complex_plane_axes:InBreak')) && all_passed;
all_passed = check('the ends of the removed stretch are not inside it', ...
    isequal(cp.x([12 28]), [12 15])) && all_passed;
all_passed = check('in_break marks it', ...
    isequal(cp.in_break([0 12 12.1 27.9 28 40]), [false false true true false false])) && all_passed;
c = cp.check([-40 0 20 21 35 50 NaN]);
all_passed = check('check counts in the break, left of the window and right of the far segment', ...
    c.n_in_break == 2 && c.n_left == 1 && c.n_right == 1) && all_passed;
all_passed = check('a far segment that starts inside the main window is refused', ...
    strcmp(error_id(@() complex_plane_axes([], 're_lim', [-32 12], 'im_lim', [-22 22], 'far_lim', [10 42])), ...
    'complex_plane_axes:BadFarLimits')) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- ticks --\n');
inside = @(t, iv) any(t > iv(1) & t < iv(2));
all_passed = check('there are ticks on both axes', ...
    ~isempty(cp.re_ticks) && ~isempty(cp.im_ticks)) && all_passed;
all_passed = check('no real tick inside the bulk interval [-27.5, 7.5]', ...
    ~inside(cp.re_ticks, [-27.5 7.5])) && all_passed;
all_passed = check('no imaginary tick inside the bulk interval [-17.5, 17.5]', ...
    ~inside(cp.im_ticks, [-17.5 17.5])) && all_passed;
all_passed = check('no tick at the origin, with or without a bulk', ...
    ~any(cp.re_ticks == 0) && ~any(cp.im_ticks == 0) && ~any(nb.re_ticks == 0) && ~any(nb.im_ticks == 0)) && all_passed;
all_passed = check('no real tick in the removed stretch, and some beyond the break', ...
    ~any(cp.in_break(cp.re_ticks)) && any(cp.re_ticks >= 28)) && all_passed;
all_passed = check('every tick is inside the limits', ...
    all(cp.re_ticks >= -32 & cp.re_ticks <= 42) && all(cp.im_ticks >= -22 & cp.im_ticks <= 22)) && all_passed;
all_passed = check('neighboring real ticks are at least tick_min_sep apart as plotted', ...
    all(diff(cp.x(cp.re_ticks)) >= 0.6 * cp.tick_step - 1e-12)) && all_passed;
sweep_ok = true;                                        % other bulks, other windows: the rule, not one case
for r = [0 3 9 17.5 25]
    for ctr = [-10, -4 + 2i, 0]
        q = complex_plane_axes([], 're_lim', [-32 12], 'im_lim', [-22 22], 'far_lim', [28 42], ...
            'bulk_center', ctr, 'bulk_radius', r);
        sweep_ok = sweep_ok && ~inside(q.re_ticks, real(ctr) + [-r r]) && ~inside(q.im_ticks, imag(ctr) + [-r r]) && ...
            ~any(q.re_ticks == 0) && ~any(q.im_ticks == 0) && ~any(q.in_break(q.re_ticks));
    end
end
all_passed = check('the same holds over a sweep of bulk centers and radii', sweep_ok) && all_passed;
all_passed = check('an explicit tick inside the bulk is an error, not a silent drop', ...
    strcmp(error_id(@() complex_plane_axes([], args{:}, 're_ticks', [-30 -10])), 'complex_plane_axes:BadTick') && ...
    strcmp(error_id(@() complex_plane_axes([], args{:}, 'im_ticks', 5)), 'complex_plane_axes:BadTick')) && all_passed;
all_passed = check('an explicit tick in the removed stretch is an error', ...
    strcmp(error_id(@() complex_plane_axes([], args{:}, 're_ticks', 20)), 'complex_plane_axes:BadTick')) && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n-- as drawn: 1:1, and the same in two panels --\n');
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 900 400]);
closer = onCleanup(@() close(fig));
ax1 = axes(fig, 'Units', 'pixels', 'Position', [40 40 380 320]);
ax2 = axes(fig, 'Units', 'pixels', 'Position', [480 40 380 320]);
d1 = complex_plane_axes(ax1, args{:});
plot(ax1, d1.x([-10 39]), [0 0], 'o');                  % data in one panel only, some of it beyond the break
d2 = complex_plane_axes(ax2, args{:}, 'draw_far', false);
drawnow;
p1 = pixel_box(ax1); p2 = pixel_box(ax2);
per_re = p1(3) / diff(ax1.XLim); per_im = p1(4) / diff(ax1.YLim);
all_passed = check(sprintf('pixels per unit of Re and of Im agree (%.3f, %.3f)', per_re, per_im), ...
    abs(per_re / per_im - 1) < 0.01) && all_passed;
all_passed = check('DataAspectRatio is [1 1 1]', ...
    isequal(ax1.DataAspectRatio, [1 1 1]) && isequal(ax2.DataAspectRatio, [1 1 1])) && all_passed;
all_passed = check('the limits are the function''s, in both panels, after plotting', ...
    isequal(ax1.XLim, d1.x_lim) && isequal(ax1.YLim, d1.y_lim) && ...
    isequal(ax2.XLim, d1.x_lim) && isequal(ax2.YLim, d1.y_lim)) && all_passed;
all_passed = check('the two plot boxes are the same size, with and without the far segment drawn', ...
    max(abs(p1(3:4) - p2(3:4))) < 0.5) && all_passed;
all_passed = check('draw_far = false leaves out the far ticks and nothing else', ...
    isequal(d2.re_ticks, d1.re_ticks(d1.re_ticks <= 12)) && isequal(d2.im_ticks, d1.im_ticks)) && all_passed;
all_passed = check('MATLAB''s own axes are off', strcmp(ax1.Visible, 'off')) && all_passed;
n_before = numel(findall(ax1, 'Tag', 'complex_plane_axes'));
complex_plane_axes(ax1, args{:});
all_passed = check('a second call replaces what the first drew', ...
    numel(findall(ax1, 'Tag', 'complex_plane_axes')) == n_before && n_before > 0) && all_passed;
clear closer

%% ------------------------------------------------------------------------
fprintf('\n');
if all_passed
    fprintf('test_complex_plane_axes: ALL TESTS PASSED\n');
else
    fprintf(2, 'test_complex_plane_axes: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function passed = check(name, condition)
if condition
    fprintf('  PASS  %s\n', name);
    passed = true;
else
    fprintf(2, '  FAIL  %s\n', name);
    passed = false;
end
end

function id = error_id(fn)
% The identifier of the error fn throws; '' when it does not throw.
try
    fn();
    id = '';
catch err
    id = err.identifier;
end
end

function p = pixel_box(ax)
% The plot box as drawn, in pixels: [x y w h].
fig = ancestor(ax, 'figure');
p = hgconvertunits(fig, tightPosition(ax), ax.Units, 'pixels', ax.Parent);
end
