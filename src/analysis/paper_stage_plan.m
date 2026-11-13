function plan = paper_stage_plan(cfg)
% PAPER_STAGE_PLAN Which stages and sweeps a bundle's config asks run_all_paper_analyses to run.
%
%   plan = PAPER_STAGE_PLAN(cfg)
%
% Config in, validated plan out; no compute and no file access, so the entry
% point calls it in its pre-flight (before the run directory exists) and a
% test calls it directly (FD-023).
%
% WHY IT EXISTS. "How big" and "which stages" are different knobs. The run
% mode (run_mode_names) sets sizes and every stage accepts every mode; this
% selects stages. A bundle that wants a quick look at a network (the
% 'essential' set of mu7chain_config) names the few stages its figures read
% rather than paying for all ten.
%
% THE CONFIG FIELDS
%   cfg.stages   cell of stage names (paper_stage_names). Absent or empty:
%                every stage, exactly as before the field existed. An
%                OPT-IN stage (paper_stage_names' third output) runs only
%                when named here; "every stage" never includes it, and a
%                plan that does not name it does not mention it (not in
%                stages_skipped, no field in runs).
%   cfg.sweeps   cell of sweep names, the analyses the 'sweeps' stage runs.
%                Absent or empty: all three. Ignored (plan.sweeps is empty)
%                when 'sweeps' itself is not selected.
%   cfg.illustrations   as before: 'paper_illustrations' exists only when
%                this is true. Naming that stage in cfg.stages without it is
%                an error rather than a stage that silently does not run.
% A char or string is accepted for a single name. The ORDER in the config is
% ignored: stages always run in the entry point's order (dc_lle after the
% stages that feed in-paper figures, and so on).
%
% THE PLAN
%   plan.stages          the stages to run, in run order
%   plan.stages_skipped  the stages this bundle has and does not run
%   plan.sweeps          the sweeps to run, in run order ({} when the
%                        'sweeps' stage is not selected)
%   plan.sweeps_skipped  the sweeps not run
%   plan.is_selection    true when anything is skipped
%   plan.runs            a struct of logicals, one field per stage name of
%                        the first list, plus a true field for each opt-in
%                        stage selected
%
% ERRORS, all listing the valid names:
%   paper_stage_plan:UnknownStage   a name in cfg.stages that is no stage
%   paper_stage_plan:UnknownSweep   a name in cfg.sweeps that is no sweep
%   paper_stage_plan:BadSelection   not a cell of names / char / string
%   paper_stage_plan:IllustrationsOff   'paper_illustrations' selected while
%                                   cfg.illustrations is absent or false
%
% See also: paper_stage_names, run_all_paper_analyses, run_all_analyses

arguments
    cfg (1,1) struct
end

[all_stages, all_sweeps, opt_in] = paper_stage_names();

has_illustrations = isfield(cfg, 'illustrations') && ~isempty(cfg.illustrations) && ...
    all(logical(cfg.illustrations));
available = all_stages;
if ~has_illustrations
    available = setdiff(all_stages, {'paper_illustrations'}, 'stable');
end

%% Stages
% An opt-in stage joins `available` only when cfg.stages names it, so a
% bundle that names none has a plan that does not mention them.
want = selection(cfg, 'stages');
if isempty(want)
    run_stage = true(size(available));
else
    unknown = setdiff(want, [all_stages, opt_in], 'stable');
    if ~isempty(unknown)
        error('paper_stage_plan:UnknownStage', ...
            'cfg.stages names %s, which %s. The stages are: %s; opt-in: %s.', ...
            quoted(unknown), is_not(unknown, 'a stage'), strjoin(all_stages, ', '), ...
            strjoin(opt_in, ', '));
    end
    available = [available, opt_in(ismember(opt_in, want))];
    if ismember('paper_illustrations', want) && ~has_illustrations
        error('paper_stage_plan:IllustrationsOff', ...
            ['cfg.stages names ''paper_illustrations'', but cfg.illustrations is absent or false, ' ...
             'so that stage would not run. Set cfg.illustrations = true or drop the name.']);
    end
    run_stage = ismember(available, want);
end
plan = struct();
plan.stages         = available(run_stage);
plan.stages_skipped = available(~run_stage);

%% Sweeps
want = selection(cfg, 'sweeps');
unknown = setdiff(want, all_sweeps, 'stable');
if ~isempty(unknown)
    error('paper_stage_plan:UnknownSweep', ...
        'cfg.sweeps names %s, which %s. The sweeps are: %s.', ...
        quoted(unknown), is_not(unknown, 'a sweep'), strjoin(all_sweeps, ', '));
end
if ~ismember('sweeps', plan.stages)
    run_sweep = false(size(all_sweeps));
elseif isempty(want)
    run_sweep = true(size(all_sweeps));
else
    run_sweep = ismember(all_sweeps, want);
end
plan.sweeps         = all_sweeps(run_sweep);
plan.sweeps_skipped = all_sweeps(~run_sweep);

plan.is_selection = ~isempty(plan.stages_skipped) || ~isempty(plan.sweeps_skipped);
plan.runs = struct();
for k = 1:numel(all_stages)
    plan.runs.(all_stages{k}) = ismember(all_stages{k}, plan.stages);
end
selected_opt_in = opt_in(ismember(opt_in, plan.stages));
for k = 1:numel(selected_opt_in)
    plan.runs.(selected_opt_in{k}) = true;
end
end

%% ------------------------------------------------------------------------
function names = selection(cfg, field)
% cfg.(field) as a row cellstr; {} when absent or empty.
names = {};
if ~isfield(cfg, field); return; end
v = cfg.(field);
if isempty(v) || ((ischar(v) || isstring(v)) && all(strlength(v) == 0)); return; end
if ~(ischar(v) || isstring(v) || iscellstr(v))
    error('paper_stage_plan:BadSelection', ...
        'cfg.%s must be a cell of names (or one name as char or string); got a %s.', field, class(v));
end
names = unique(reshape(cellstr(v), 1, []), 'stable');
end

function s = quoted(names)
s = strjoin(cellfun(@(n) ['''' n ''''], names, 'UniformOutput', false), ', ');
end

function s = is_not(names, what)
if isscalar(names); s = ['is not ' what]; else; s = ['are not ' what 's']; end
end
