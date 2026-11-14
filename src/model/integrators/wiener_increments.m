function [xi1, xi2] = wiener_increments(n_rows, n_steps, seed, schedule, t0, fs, err_id_prefix)
%WIENER_INCREMENTS The standard normals behind a run's Wiener path, with seed switches.
%
%   [xi1, xi2] = WIENER_INCREMENTS(n_rows, n_steps, seed)
%   [xi1, xi2] = WIENER_INCREMENTS(n_rows, n_steps, seed, schedule, t0, fs, err_id_prefix)
%
% xi1 and xi2 are n_rows x n_steps, the two tensors sde_fixed_step reads
% (noise.xi1, noise.xi2): column j belongs to the step that starts at
% t0 + (j - 1) / fs.
%
% WITHOUT A SCHEDULE this is the draw the model classes have always made:
% rng(seed, 'twister'), xi1 = randn(n_rows, n_steps), then xi2 the same.
%
% WITH A SCHEDULE, rows of [t_start, seed] (a model's noise_seed_schedule),
% the path is SPLICED AT THE SAME ABSOLUTE TIME INDEX. Every seed named,
% the base one and each scheduled one, defines a whole path over the run,
% drawn exactly as above. Row k of the schedule takes over at the step that
% starts at t_start(k): from that column on, both tensors hold the columns
% of seed(k)'s path, up to the next row's t_start. Nothing is re-indexed,
% so
%   - two runs with the same base seed are bit-identical before the first
%     switch whatever their schedules say afterwards (the state at
%     t_start(1) included: the step that ends there used a base column);
%   - a schedule naming the base seed reproduces the unscheduled path;
%   - two runs with the same schedule are bit-identical throughout.
% It costs one extra full draw per schedule row, which is small next to the
% integration.
%
% t_start must sit on the run's time grid, strictly inside it: after the
% first sample (the base seed owns the first step; to change the whole path
% set the base seed) and before the last (a switch there would change
% nothing). Violations raise <err_id_prefix>:InvalidNoiseSchedule.
%
% The global random stream is saved and restored, so the caller's W,
% stimulus and initial condition are untouched.
%
% See also: check_noise_seed_schedule, sde_fixed_step,
%           SRNNCellTypePairs.noise_seed_schedule

if nargin < 4, schedule = []; end
if nargin < 7 || isempty(err_id_prefix), err_id_prefix = 'wiener_increments'; end

first_col = zeros(1, 0);
if ~isempty(schedule)
    check_noise_seed_schedule(schedule, err_id_prefix);
    first_col = switch_columns(schedule(:, 1), n_steps, t0, fs, err_id_prefix);
end

stream_state = rng;
restore = onCleanup(@() rng(stream_state));

[xi1, xi2] = draw_path(n_rows, n_steps, seed);
for k = 1:numel(first_col)
    last_col = n_steps;
    if k < numel(first_col), last_col = first_col(k + 1) - 1; end
    cols = first_col(k):last_col;
    [s1, s2] = draw_path(n_rows, n_steps, schedule(k, 2));
    xi1(:, cols) = s1(:, cols);
    xi2(:, cols) = s2(:, cols);
end
end

% -------------------------------------------------------------------------
function [xi1, xi2] = draw_path(n_rows, n_steps, seed)
% One seed's whole path. xi1 first, then xi2: the order the classes have
% always drawn in, which existing results depend on.
rng(seed, 'twister');
xi1 = randn(n_rows, n_steps);
xi2 = randn(n_rows, n_steps);
end

% -------------------------------------------------------------------------
function cols = switch_columns(t_start, n_steps, t0, fs, err_id_prefix)
% The column (step index) at which each schedule row takes over.
pos  = (t_start(:)' - t0) * fs;          % steps after the first sample
cols = round(pos) + 1;
off_grid = abs(pos - round(pos)) > 1e-6;
if any(off_grid)
    error([err_id_prefix ':InvalidNoiseSchedule'], ...
        ['noise_seed_schedule: t_start = %s is not on the run''s time grid ' ...
         '(t0 = %g s, fs = %g Hz).'], mat2str(t_start(off_grid)'), t0, fs);
end
outside = cols < 2 | cols > n_steps;
if any(outside)
    error([err_id_prefix ':InvalidNoiseSchedule'], ...
        ['noise_seed_schedule: t_start = %s is not strictly inside the run ' ...
         '(%g, %g) s. The base seed owns the first step; a switch at or ' ...
         'after the last sample would change nothing.'], ...
        mat2str(t_start(outside)'), t0, t0 + n_steps / fs);
end
end
