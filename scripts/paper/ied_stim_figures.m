function results = ied_stim_figures(cfg, names)
% IED_STIM_FIGURES Draw the ied_stim bundle's registered figures from a complete run directory (the figure step of ied_stim_run).
%
%   results = ied_stim_figures(ied_stim_config('fast'))
%   results = ied_stim_figures(cfg, {'fig_ied_stim_blocks', 'Fig_Stim_Combined'})
%
% Each entry of cfg.figures is called as
%   fn('data_dir', cfg.run_dir, 'out_dir', <cfg.fig_root>/<name>, 'visible', cfg.visible_figures, args{:})
% in registry order. An entry that throws is recorded and the rest still
% run. names, optional: draw only these entries. Nothing is simulated, so
% this redraws a finished run (for example the 2026-11-17 data:
% ied_stim_figures(ied_stim_config('stim_blocks_20261117'))).
%
% RESULTS is a struct array: name, in_paper, ok, seconds, out_dir, message.
arguments
    cfg (1,1) struct
    names = {}
end
if ischar(names) || isstring(names), names = cellstr(names); end
assert(isfile(fullfile(cfg.run_dir, 'summary.mat')), 'ied_stim_figures:NoRun', ...
    'No summary.mat in %s: simulate first (ied_stim_run).', cfg.run_dir);
% The run directory must be the run the config describes.
L = load(fullfile(cfg.run_dir, 'summary.mat'), 'opts'); o = L.opts;
if ~isfield(o, 'fs'), o.fs = 400; end   % runs before fs was a setting used ied_base_config's 400 Hz
for key = {'seeds','amplitude','n_no','n_stim','block_len','settle','fs'}
    assert(isequal(o.(key{1}), cfg.sim.(key{1})), 'ied_stim_figures:WrongRun', ...
        '%s was run with %s = %s; the config (%s) says %s.', cfg.run_dir, key{1}, ...
        mat2str(o.(key{1})), cfg.run_mode, mat2str(cfg.sim.(key{1})));
end
F = [cfg.figures{:}];
if ~isempty(names)
    unknown = setdiff(names, {F.name});
    assert(isempty(unknown), 'ied_stim_figures:Unknown', 'Not registered: %s', strjoin(unknown, ', '));
    F = F(ismember({F.name}, names));
end
results = struct('name', {}, 'in_paper', {}, 'ok', {}, 'seconds', {}, 'out_dir', {}, 'message', {});
for k = 1:numel(F)
    e = F(k); out_dir = fullfile(cfg.fig_root, e.name); t = tic; msg = '';
    try
        e.fn('data_dir', cfg.run_dir, 'out_dir', out_dir, 'visible', cfg.visible_figures, e.args{:});
        ok = true;
    catch err
        ok = false; msg = getReport(err, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ied_stim_figures] %s FAILED: %s\n', e.name, msg);
    end
    results(end+1) = struct('name', e.name, 'in_paper', e.in_paper, 'ok', ok, ...
        'seconds', toc(t), 'out_dir', out_dir, 'message', msg); %#ok<AGROW>
    fprintf('[ied_stim_figures] %-32s %s  %.1f s\n', e.name, pass_fail(ok), results(end).seconds);
end
end

function s = pass_fail(ok)
if ok, s = 'ok'; else, s = 'FAILED'; end
end
