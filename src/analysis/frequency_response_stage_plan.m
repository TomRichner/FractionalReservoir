function plan = frequency_response_stage_plan(cfg)
% FREQUENCY_RESPONSE_STAGE_PLAN What a bundle's frequency_response stage computes, validated before any compute.
%
%   plan = FREQUENCY_RESPONSE_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'frequency_response'
% (paper_opt_in_stages; FD-041, Tom, 2026-11-20: "how much does it
% attenuate modes with different oscillatory frequencies? ... could we make
% something of a frequency response for this?"). The network frequency
% response of the linearization at the calibration states that the stage
% outlier_spectra saved in the same run directory: per neuron, the complex
% gain with its adaptation and depression variables eliminated; from the
% gains the bulk radius, the rightmost eigenvalue and the loop gain against
% frequency; the checks that the elimination is exact and local; and
% single-neuron responses. Config in, plan out; nothing is computed and no
% file is touched. run_frequency_response_stage carries the plan out.
%
% IT SIMULATES NOTHING. Its input is <run_dir>/outlier_spectra/ (FD-040):
% it runs on every condition and seed that stage saved for the networks
% named here, so the number of seeds is that stage's, not this one's, and
% the settings below are the same in every run mode.
%
% THE CONFIG
%   cfg.run_mode   checked against run_mode_names; sizes nothing here
%   cfg.stages     when the bundle names its stages, 'outlier_spectra' must
%                  be among them (the input would otherwise be missing only
%                  hours into the run), unless
%                  cfg.frequency_response.input_exists is true
%   cfg.frequency_response   optional struct of overrides. A field the
%                  stage does not know is refused, so a typo cannot be
%                  ignored.
%
% THE SETTINGS (defaults; the field of cfg.frequency_response that changes each)
%   networks     {'reference', 'EEx1p5'}: the paper's network and E-to-E
%                mean x 1.5, by the outlier_spectra stage's names. Every
%                one must be in that stage's index.
%   f            logspace(-3, 2, 61) Hz: the grid of the network curves
%   f_single     logspace(-3, 2, 201) Hz: the grid of the single-neuron
%                responses
%   exact_states []: the calibration states where the dense full Jacobian
%                is taken (its spectrum, exactness, locality, the radius
%                over the complex plane). Empty = the first and the middle
%                state (1 and 7 of 12).
%   shown_state  []: the state the plane figure draws and the single-neuron
%                responses are taken at. Empty = the last of exact_states.
%   n_lead       20, n_sample 200: eigenvalues of the full Jacobian tested
%                for exactness, the rightmost and a random draw
%   exact_f      [0.001 0.03 0.5 3 30] Hz: where the closed-form gains are
%                compared with the Jacobian's blocks
%   sigma, w     -12:0.1:8 1/s and 0:0.25:30 rad/s: the grid of complex z on
%                which the predicted radius is evaluated for the bulk's edge
%   plane_all_states  true: the predicted abscissa at every calibration
%                state (false: at the check states only)
%   r0           0.25: the rate of the matched operating point (the
%                phase-lead figure)
%   input_exists false: true says the outlier_spectra output is already in
%                the run directory, so the bundle need not name that stage
%
% THE PLAN
%   plan.folder        where it writes, relative to the run directory
%                      ('frequency_response')
%   plan.file          the .mat it writes there ('frequency_response_data.mat')
%   plan.input_folder  what it reads ('outlier_spectra')
%   plan.networks      the networks
%   plan.opts          the settings above, as fresp_run_one takes them
%
% See also: run_frequency_response_stage, fresp_run_one, fresp_load_states,
%           paper_opt_in_stages

arguments
    cfg (1,1) struct
end

DEFAULTS = struct( ...
    'networks', {{'reference', 'EEx1p5'}}, ...
    'f', logspace(-3, 2, 61), ...
    'f_single', logspace(-3, 2, 201), ...
    'exact_states', [], ...
    'shown_state', [], ...
    'n_lead', 20, ...
    'n_sample', 200, ...
    'exact_f', [0.001 0.03 0.5 3 30], ...
    'sigma', -12:0.1:8, ...
    'w', 0:0.25:30, ...
    'plane_all_states', true, ...
    'r0', 0.25, ...
    'input_exists', false);
INPUT_STAGE = 'outlier_spectra';

modes = run_mode_names();
if ~isfield(cfg, 'run_mode') || ~(ischar(cfg.run_mode) || isstring(cfg.run_mode)) || ...
        ~ismember(char(cfg.run_mode), modes)
    bad_config('The frequency_response stage needs cfg.run_mode, one of %s.', strjoin(modes, ', '));
end

over = struct();
if isfield(cfg, 'frequency_response') && ~isempty(cfg.frequency_response)
    over = cfg.frequency_response;
end
if ~(isstruct(over) && isscalar(over))
    bad_config('cfg.frequency_response must be a scalar struct.');
end
unknown = setdiff(fieldnames(over), fieldnames(DEFAULTS));
if ~isempty(unknown)
    bad_config('cfg.frequency_response has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(fieldnames(DEFAULTS)', ', '));
end
o = DEFAULTS;
names = fieldnames(over);
for k = 1:numel(names)
    o.(names{k}) = over.(names{k});
end

%% The settings, one by one
if ischar(o.networks) || isstring(o.networks); o.networks = cellstr(o.networks); end
if ~(iscellstr(o.networks) && ~isempty(o.networks) && numel(unique(o.networks)) == numel(o.networks))
    bad_config('cfg.frequency_response.networks must be a non-empty cell of distinct network names.');
end
o.networks = o.networks(:)';
for name = {'f', 'f_single', 'exact_f'}
    v = o.(name{1});
    if ~(isnumeric(v) && isvector(v) && all(isfinite(v)) && all(v > 0))
        bad_config('cfg.frequency_response.%s must be a vector of positive frequencies in Hz.', name{1});
    end
    o.(name{1}) = v(:)';
end
for name = {'f', 'f_single'}
    if ~(numel(o.(name{1})) >= 2 && all(diff(o.(name{1})) > 0))
        bad_config('cfg.frequency_response.%s must be ascending, with at least two points.', name{1});
    end
end
if ~(isempty(o.exact_states) || (isnumeric(o.exact_states) && isvector(o.exact_states) && ...
        all(is_whole(o.exact_states)) && all(o.exact_states >= 1)))
    bad_config('cfg.frequency_response.exact_states must be empty or a vector of state numbers (1, 2, ...).');
end
o.exact_states = o.exact_states(:)';
if ~(isempty(o.shown_state) || (isscalar(o.shown_state) && is_whole(o.shown_state) && o.shown_state >= 1))
    bad_config('cfg.frequency_response.shown_state must be empty or one state number.');
end
for name = {'n_lead', 'n_sample'}
    v = o.(name{1});
    if ~(isscalar(v) && is_whole(v) && v >= 1)
        bad_config('cfg.frequency_response.%s must be a whole number of at least 1.', name{1});
    end
end
for name = {'sigma', 'w'}
    v = o.(name{1});
    if ~(isnumeric(v) && isvector(v) && numel(v) >= 2 && all(isfinite(v)) && all(diff(v) > 0))
        bad_config('cfg.frequency_response.%s must be an ascending grid of at least two points.', name{1});
    end
    o.(name{1}) = v(:)';
end
if o.w(1) < 0
    bad_config('cfg.frequency_response.w is the imaginary part in the upper half plane: it must start at 0 or above.');
end
if ~(isscalar(o.r0) && isnumeric(o.r0) && o.r0 > 0 && o.r0 < 1)
    bad_config('cfg.frequency_response.r0 must be a rate in (0, 1).');
end
for name = {'plane_all_states', 'input_exists'}
    v = o.(name{1});
    if ~(isscalar(v) && (islogical(v) || isnumeric(v)))
        bad_config('cfg.frequency_response.%s must be true or false.', name{1});
    end
    o.(name{1}) = logical(v);
end

%% The input stage
if isfield(cfg, 'stages') && ~isempty(cfg.stages) && ~o.input_exists && ...
        ~ismember(INPUT_STAGE, cellstr(cfg.stages))
    bad_config(['The frequency_response stage reads what the stage %s saves in the same run ' ...
        'directory, and this bundle''s cfg.stages does not name it. Name it, or set ' ...
        'cfg.frequency_response.input_exists = true when its output is already there.'], INPUT_STAGE);
end

plan = struct();
plan.folder       = 'frequency_response';
plan.file         = 'frequency_response_data.mat';
plan.input_folder = INPUT_STAGE;
plan.networks     = o.networks;
plan.opts         = rmfield(o, {'networks', 'input_exists'});
end

%% ------------------------------------------------------------------------
function tf = is_whole(v)
tf = isnumeric(v) & isfinite(v) & v == round(v);
end

function bad_config(varargin)
error('frequency_response_stage_plan:BadConfig', varargin{:});
end
