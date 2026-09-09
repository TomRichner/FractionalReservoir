function colors = mc_condition_colors(keys)
% MC_CONDITION_COLORS One colour row per condition, keyed by NAME.
%
%   colors = MC_CONDITION_COLORS(results_all.conditions)   % n_cond x 3
%
% The memory-capacity plotters used to carry a POSITIONAL palette -- four rows
% in mc_pairs_dualStd's regime order (no_adaptation, sfa3_std0, sfa0_std2,
% sfa3_std2) -- and apply them to conditions 1..n. That is only right for that
% one preset. On the paper's 3-condition presets it gave sfa1_std1 the SFA-only
% orange and sfa3_std2 the STD-only blue: plausible-looking, wrong hues.
%
% This reads manuscript_style's condition_color map, which is keyed by the
% snake_case condition name that saved runs carry, so every figure in the
% manuscript colours a regime identically. An unknown key -- a run directory can
% name anything -- falls back to a lines() row rather than erroring, the same
% guard mc_display_names applies to titles.
%
% See also: mc_display_names, manuscript_style, srnn_condition_titles

st = manuscript_style();
fallback = lines(max(numel(keys), 1));
colors = zeros(numel(keys), 3);
for k = 1:numel(keys)
    key = keys{k};
    if ischar(key) && isKey(st.condition_color, key)
        colors(k, :) = st.condition_color(key);
    else
        colors(k, :) = fallback(k, :);
    end
end
end
