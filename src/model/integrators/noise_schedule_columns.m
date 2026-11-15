function cols = noise_schedule_columns(t_start, n_steps, t0, fs, err_id_prefix)
%NOISE_SCHEDULE_COLUMNS The step at which each row of a noise seed schedule takes over.
%
%   cols = NOISE_SCHEDULE_COLUMNS(t_start, n_steps, t0, fs, err_id_prefix)
%
% t_start is column 1 of a noise_seed_schedule; the run has n_steps steps
% and step j starts at t0 + (j - 1) / fs. cols(k) is the step that starts
% at t_start(k): from it on, the increments come from row k's seed.
%
% One definition for every noise source (wiener_increments,
% release_increments), so that a schedule splices all of them at the same
% absolute step. It was a local function of wiener_increments until the
% release noise needed it too.
%
% t_start must sit on the run's time grid, strictly inside it: after the
% first sample (the base seed owns the first step; to change the whole
% path set the base seed) and before the last (a switch there would change
% nothing). Violations raise <err_id_prefix>:InvalidNoiseSchedule.
%
% See also: wiener_increments, release_increments, check_noise_seed_schedule

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
