function dR = release_increments(rate, sigma, fs, steps, seed, schedule, t0, n_steps, err_id_prefix)
%RELEASE_INCREMENTS The gamma increments behind a run's release noise.
%
%   dR = RELEASE_INCREMENTS(rate, sigma, fs, steps, seed)
%   dR = RELEASE_INCREMENTS(rate, sigma, fs, steps, seed, schedule, t0, n_steps, err_id_prefix)
%
% The spontaneous release drive of each presynaptic neuron over each step
% of a run (route_model, "RELEASE NOISE"). Row i is a gamma process with
% independent increments
%
%   dR_i ~ Gamma(shape = rate_i^2 dt / sigma_i^2, scale = sigma_i^2 / rate_i)
%
% over a step dt = 1/fs: mean rate_i dt, variance sigma_i^2 dt, never
% negative, independent across neurons and steps and of the model's state.
%
%   rate, sigma  n x 1 (or scalars): r_s and sigma_gamma of each neuron.
%                A row with sigma = 0 is the EXACT LEAK, dR = rate dt, and
%                nothing is drawn for it. sigma > 0 needs rate > 0.
%   fs           the run's sampling rate; a step is 1/fs.
%   steps        [first, last]: the ABSOLUTE steps wanted, counted from the
%                start of the run (step j starts at t0 + (j - 1)/fs).
%   seed         the run's noise seed (a model's noise_seed, or what it
%                derives one from). The generator adds its own offset, so
%                the release stream is not the Wiener path's.
%   schedule, t0, n_steps, err_id_prefix
%                a noise_seed_schedule, rows of [t_start, seed], with the
%                run's first sample time and its number of steps: the
%                increments of every step that starts at or after
%                t_start(k) come from seed(k)'s stream, spliced at the
%                same absolute step as wiener_increments splices the
%                Wiener path (noise_schedule_columns).
%   dR           n x (last - first + 1). Column c belongs to step
%                first + c - 1.
%
% THE INCREMENT OF A STEP IS A FUNCTION OF (seed, absolute step) AND
% NOTHING ELSE. Step j of a seed is drawn from substream j of that seed's
% 'mrg32k3a' stream: one uniform per neuron, turned into the gamma variate
% by the inverse of its distribution function. So
%   - any range of steps can be asked for, in any order and any number of
%     pieces, and gives the numbers the whole run would give: a run
%     integrated in segments consumes exactly what a single pass consumes;
%   - two runs with the same base seed are bit-identical before a
%     schedule's first switch, and a schedule naming the base seed
%     reproduces the unscheduled increments;
%   - the uniform of (neuron, step) does not depend on rate or sigma, and
%     the inverse is increasing in it: runs that differ only in sigma get
%     the same events, scaled. A sweep over sigma is not also a sweep over
%     noise realizations, and sigma -> 0 goes to the leak continuously.
% Substreams, and not seed + step or seed + block: repetition seeds are
% consecutive integers, so offsets of the seed would hand one repetition
% another's numbers; the substreams of one seed are independent by the
% generator's construction. Every neuron's uniform is drawn whether or not
% its row has sigma > 0, so a neuron's increments do not depend on which
% other neurons release.
%
% SMALL SHAPES. At the paper's values the shape per step is about 0.01:
% most increments are far below the mean (many underflow to exactly 0)
% and a few are large. Nothing here or downstream divides by a draw.
%
% The global random stream is not touched.
%
% Cost: about 25 ms per 1e5 increments (the inverse), 10 microseconds per
% step (the substream), measured 2026-11-14.
%
% See also: wiener_increments, noise_schedule_columns, sde_fixed_step,
%           SRNNCellTypePairs.release_sigma, test_release_noise

if nargin < 6, schedule = []; end
if nargin < 9 || isempty(err_id_prefix), err_id_prefix = 'release_increments'; end
id = [err_id_prefix ':InvalidParams'];

rate = double(rate(:));
sigma = double(sigma(:));
n_rows = max(numel(rate), numel(sigma));
if isscalar(rate), rate = repmat(rate, n_rows, 1); end
if isscalar(sigma), sigma = repmat(sigma, n_rows, 1); end
if numel(rate) ~= numel(sigma) || any(~isfinite(rate)) || any(~isfinite(sigma)) || ...
        any(rate < 0) || any(sigma < 0)
    error(id, 'rate and sigma must be finite, nonnegative and of one length (or scalars).');
end
if any(sigma > 0 & rate == 0)
    error(id, ['A row with sigma > 0 needs rate > 0: the gamma increment has mean rate dt ' ...
        'and variance sigma^2 dt, and no such variate exists at mean zero.']);
end
if ~(isnumeric(steps) && numel(steps) == 2 && all(steps == round(steps)) && ...
        steps(1) >= 1 && steps(2) >= steps(1) - 1)
    error(id, 'steps must be [first, last], integers with 1 <= first and first - 1 <= last.');
end
cols = steps(1):steps(2);
dt = 1 / fs;

% The leak. Rows without noise keep it exactly.
dR = repmat(rate .* dt, 1, numel(cols));
noisy = find(sigma > 0);
if isempty(noisy) || isempty(cols)
    return;
end

% Which seed owns each step: the base seed, then each schedule row from
% its own first step on.
owner = repmat(double(seed), 1, numel(cols));
if ~isempty(schedule)
    check_noise_seed_schedule(schedule, err_id_prefix);
    first_col = noise_schedule_columns(schedule(:, 1), n_steps, t0, fs, err_id_prefix);
    for k = 1:numel(first_col)
        owner(cols >= first_col(k)) = schedule(k, 2);
    end
end

% One uniform per neuron and step, from substream (step) of the owner's
% stream.
U = zeros(n_rows, numel(cols));
seeds = unique(owner);
for k = 1:numel(seeds)
    stream = RandStream('mrg32k3a', 'Seed', stream_seed(seeds(k)));
    for c = find(owner == seeds(k))
        stream.Substream = cols(c);
        U(:, c) = rand(stream, n_rows, 1);
    end
end

% The gamma variate by inversion: scale * P^{-1}(shape, U), P the
% regularized lower incomplete gamma function.
shape = rate(noisy) .^ 2 .* dt ./ sigma(noisy) .^ 2;
scale = sigma(noisy) .^ 2 ./ rate(noisy);
dR(noisy, :) = scale .* gammaincinv(U(noisy, :), repmat(shape, 1, numel(cols)));
end

% -------------------------------------------------------------------------
function s = stream_seed(noise_seed)
% The release stream's seed from a noise seed: its own offset (a prime, as
% the class's other derived seeds use), wrapped into what RandStream takes.
s = mod(noise_seed + 86028121, 2 ^ 32);
end
