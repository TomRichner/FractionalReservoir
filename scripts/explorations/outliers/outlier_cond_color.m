function col = outlier_cond_color(st, name)
%OUTLIER_COND_COLOR The color the outlier exploration figures draw a condition in.
%
%   col = outlier_cond_color(manuscript_style(), 'sfa1_std1')
%
% manuscript_style's condition colors, except that no adaptation is the
% dark red of the Figure 3 Ritz overlay (FR-004, Tom, 2026-11-20; the value
% fig_eig_heatmap's draw_ritz substitutes for a colorless condition color),
% not manuscript_style's black, so that it cannot be confused with the
% marks that are black for another reason (the full Jacobian, the
% uniform-gain diamonds, zero lines). The shared map in manuscript_style is
% deliberately not changed here: the paper's figures change together on the
% next full figures run, and this helper goes then.
%
% See also: fig_outliers, fig_outlier_structures, manuscript_style
if strcmp(name, 'no_adaptation')
    col = [0.75 0.10 0.15];
else
    col = st.condition_color(name);
end
end
