function check_noise_seed_schedule(schedule, err_id_prefix)
%CHECK_NOISE_SEED_SCHEDULE Validate the shape and values of a noise seed schedule.
%
%   CHECK_NOISE_SEED_SCHEDULE(schedule, err_id_prefix)
%
% A schedule is empty (no switch) or a K x 2 real matrix, rows of
% [t_start, seed]: finite times in strictly increasing order, and seeds
% that rng accepts (integers in [0, 2^32)). Whether the times fall on a
% run's grid and inside it is checked where the grid is known, in
% wiener_increments.
%
% Called from a model's validate(), so a malformed schedule fails at
% construction or build() and not at the top of run().
%
% See also: wiener_increments, check_noise_settings

if nargin < 2 || isempty(err_id_prefix), err_id_prefix = 'check_noise_seed_schedule'; end
if isempty(schedule), return; end

id = [err_id_prefix ':InvalidNoiseSchedule'];
if ~isnumeric(schedule) || ~isreal(schedule) || ~ismatrix(schedule) || size(schedule, 2) ~= 2
    error(id, 'noise_seed_schedule must be empty or a K x 2 real matrix, rows of [t_start, seed].');
end
t_start = schedule(:, 1);
seeds   = schedule(:, 2);
if any(~isfinite(t_start)) || any(diff(t_start) <= 0)
    error(id, 'noise_seed_schedule: the switch times (column 1) must be finite and strictly increasing.');
end
if any(~isfinite(seeds)) || any(seeds ~= round(seeds)) || any(seeds < 0) || any(seeds >= 2^32)
    error(id, 'noise_seed_schedule: the seeds (column 2) must be integers in [0, 2^32).');
end
end
