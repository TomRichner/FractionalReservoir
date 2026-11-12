function titles = titles_by_name(names, saved)
% TITLES_BY_NAME Display titles for condition NAMES, read from srnn_condition_titles.
%
%   titles = TITLES_BY_NAME(names)
%   titles = TITLES_BY_NAME(names, saved)
%
%   names  cell of condition names (snake_case keys), in drawing order.
%   saved  (optional) cell of the titles a stage saved with its data, in the
%          same order as names; used only for a name the map does not know.
%   titles 1 x N cellstr.
%
% WHY. Seven analysis stages save a display title beside each condition at
% compute time (inventory R3: run_eig_heatmap, run_paper_illustrations,
% run_local_lyapunov, run_lyapunov_spectrum, run_numerics_verification,
% run_transient_gain, run_stim_evoked). A figure that draws the SAVED title
% keeps whatever the condition was called when the run was made: after
% sfa3_std2 was retitled "(product)" on 2026-11-12, every figure regenerated
% from an older run would still call the product "Multiple-Timescale
% Adaptation". The map is keyed by name precisely so a saved run can be
% retitled without recomputing it, so the figures title from the name.
%
% A name the map lacks (a run directory can name anything) falls back to the
% saved title when one is given, else to the name itself; never an error.
%
% See also: srnn_condition_titles, select_conditions, mc_display_names

if nargin < 2; saved = {}; end
if ischar(names); names = {names}; end
map = srnn_condition_titles();
titles = cell(1, numel(names));
for k = 1:numel(names)
    n = names{k};
    if ischar(n) && isKey(map, n)
        titles{k} = map(n);
    elseif numel(saved) >= k && ischar(saved{k}) && ~isempty(saved{k})
        titles{k} = saved{k};
    else
        titles{k} = n;
    end
end
end
