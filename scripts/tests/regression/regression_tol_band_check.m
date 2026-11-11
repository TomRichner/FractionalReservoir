function [ok, diffs, info] = regression_tol_band_check(tier, opts)
% REGRESSION_TOL_BAND_CHECK Layer 3 of the tolerance tier: a chaotic tier ran in full and its statistics sit in a loose band.
%
%   ok = regression_tol_band_check('rapid')        % about 20 s, no pool needed
%   ok = regression_tol_band_check('stages', 'ignore', {'settings\.n_workers$'})   % minutes, pool up
%   ok = regression_tol_band_check('full', 'tag', 'after_refactor', ...)           % about 50 min
%   [ok, diffs, info] = regression_tol_band_check(tier, 'band', 0.1)
%
% WHAT IT IS FOR (FD-020; Plan A, section 2.2, third row). After a rewrite
% of the model's arithmetic the byte-exact tiers are expected to differ:
% their networks are chaotic, so a last-bit change grows to order one and
% no tolerance on their numbers means anything. What CAN still be asked of
% them is structural: every stage ran, every file and figure was written,
% every saved array is still there with the same size and the same count
% of non-finite values. This runs one of the existing checks
% (regression_rapid_check, regression_stages_check, regression_check)
% against its usual reference and re-reads the outcome that way:
%
%   FAIL    rows of kind missing, extra, shape, file, figure, and text
%           (a changed text value is structural: a name, a solver, a label)
%   REPORT  rows of kind value. They are summarised as a band: how many
%           leaves moved at all, how many by more than BAND relative to
%           the leaf's scale (regression_compare's measure, on the digest
%           statistics), the quantiles, and the largest few by name.
%
% A value row never fails this check. Whether the numbers are RIGHT is
% the business of the first two layers (regression_tol_check), which
% compare element by element where a tolerance means something. Read the
% band as a smoke alarm: summary statistics of a chaotic run should move
% by a sampling fluctuation, not by an order of magnitude.
%
% It is run on demand by the lead, not by run_all_tests('quick'): 'stages'
% and 'full' take minutes and need a pool up (serial results differ from
% pool results at round-off, and those two references were captured on a
% pool). The n_workers leaf of three stages differs with the pool size;
% pass 'ignore', {'settings\.n_workers$'} until the references store it.
%
% Options: band (default regression_tol_settings().rtol.band), ignore
% (cell of regular expressions, added to the reference's own), keep, and
% for 'full' tag and figures, as regression_check takes them. The
% underlying check keeps its output folder whenever anything differs,
% which after a refactor is always; unless KEEP is set this removes the
% rapid and stages folders again when the structural verdict is a pass.
%
% INFO: info.check (the underlying check's info), info.n_moved,
% info.n_outside_band, info.rel (the relative size of every value row).
%
% See also: regression_tol_check, regression_tol_settings,
%           regression_compare, regression_rapid_check,
%           regression_stages_check, regression_check
arguments
    tier (1,:) char {mustBeMember(tier, {'rapid', 'stages', 'full'})}
    opts.band double {mustBeScalarOrEmpty, mustBeNonnegative} = []
    opts.ignore cell = {}
    opts.keep (1,1) logical = false
    opts.tag (1,:) char = ''
    opts.figures (1,1) logical = true
end
band = opts.band;
if isempty(band)
    settings = regression_tol_settings();
    band = settings.rtol.band;
end
project_root = fileparts(which('setup_paths'));
fixtures = fullfile(project_root, 'scripts', 'tests', 'fixtures');
ref_names = struct('rapid', 'regression_rapid_digest.mat', 'stages', 'regression_stages_digest.mat', ...
    'full', 'regression_fast_digest.mat');

% The underlying check, in EXACT mode and kept: its printout lists every
% leaf that moved, and its digest is compared again below.
switch tier
    case 'rapid'
        [~, ~, check] = regression_rapid_check('ignore', opts.ignore, 'keep', true);
        out_dir = check.out_dir;
    case 'stages'
        [~, ~, check] = regression_stages_check('ignore', opts.ignore, 'keep', true);
        out_dir = check.out_dir;
    case 'full'
        [~, ~, check] = regression_check('tag', opts.tag, 'figures', opts.figures, ...
            'ignore', opts.ignore, 'keep', true);
        out_dir = '';      % a 50-minute run is left for the caller to remove
end

R = load(fullfile(fixtures, ref_names.(tier)));
[~, exact] = regression_compare(R.digest, check.digest, 'ignore', opts.ignore, 'verbose', false);
[~, diffs] = regression_compare(R.digest, check.digest, 'rtol', band, 'ignore', opts.ignore, 'verbose', false);
live = diffs(~diffs.ignored, :);
structural = live(~strcmp(live.kind, 'value'), :);
outside = live(strcmp(live.kind, 'value'), :);
moved = exact(~exact.ignored & strcmp(exact.kind, 'value'), :);
ok = isempty(structural);

fprintf('\nregression_tol_band_check(''%s''): %d leaves in the reference\n', tier, numel(R.digest.keys));
fprintf('  structural differences (these FAIL): %d\n', height(structural));
for k = 1:min(height(structural), 40)
    fprintf('    %-7s %s\n            %s\n', structural.kind{k}, structural.key{k}, structural.detail{k});
end
if height(structural) > 40
    fprintf('    ... and %d more (see the returned table)\n', height(structural) - 40);
end
fprintf('  value leaves that moved at all: %d\n', height(moved));
rel = moved.rel(isfinite(moved.rel));
if ~isempty(rel)
    q = quantile(rel, [0.5 0.9 1]);
    fprintf('    relative to the leaf''s scale: median %.1e, 90th percentile %.1e, largest %.1e\n', q(1), q(2), q(3));
end
fprintf('  value leaves outside the band of %g (REPORTED, not failed): %d\n', band, height(outside));
for k = 1:min(height(outside), 20)
    fprintf('    %s\n            %s  (rel %.2e)\n', outside.key{k}, outside.detail{k}, outside.rel(k));
end
if height(outside) > 20
    fprintf('    ... and %d more (see the returned table)\n', height(outside) - 20);
end
if ok
    fprintf('  PASS: everything the reference holds was produced again, with the same shapes\n');
else
    fprintf('  FAIL: the run is structurally different from the reference\n');
end

info = struct('check', check, 'n_moved', height(moved), 'n_outside_band', height(outside), ...
    'rel', moved.rel, 'band', band, 'structural', structural);
if ok && ~opts.keep && ~isempty(out_dir) && isfolder(out_dir)
    rmdir(out_dir, 's');
elseif ~isempty(out_dir) && isfolder(out_dir)
    fprintf('  kept %s\n', out_dir);
end
end
