function [D, data_file] = fresp_fig_load(data_file, run_dir)
% FRESP_FIG_LOAD What the frequency_response stage saved, for its figures: from one run directory, or the file given.
%
%   [D, data_file] = FRESP_FIG_LOAD(data_file, run_dir)
%
% The four frequency-response figures (fig_freqresp_network, _plane,
% _phase, _single_neuron) read
% <run_dir>/frequency_response/frequency_response_data.mat, or the file
% given explicitly, through resolve_data_file: one run directory or an
% error, no fallback. D has the file's variables: settings, runs, matched
% (run_frequency_response_stage).
%
% See also: resolve_data_file, run_frequency_response_stage, fresp_fig_save

data_file = resolve_data_file(data_file, run_dir, ...
    {fullfile(run_dir, 'frequency_response')}, 'frequency_response_data.mat', ...
    'Run the stage frequency_response (run_frequency_response_stage) first');
D = load(data_file);
missing = setdiff({'settings', 'runs', 'matched'}, fieldnames(D));
if ~isempty(missing)
    error('fresp_fig_load:BadFile', '%s lacks the variable(s) %s.', data_file, strjoin(missing, ', '));
end
end
