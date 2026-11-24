function [titles, ostruct_keys] = cx_gallery_names(X, stage_seed)
% CX_GALLERY_NAMES What the figures call the structures of the gallery S3, and the key of each in the outlier_structures stage.
%
%   titles = CX_GALLERY_NAMES(C.gallery)
%   [titles, ostruct_keys] = CX_GALLERY_NAMES(C.gallery, C.settings.seed)
%
%   X            the gallery of the concept_examples stage (cx_load: 1 x m
%                struct with key, name, label, seed)
%   stage_seed   the stage's weight seed (settings.seed), the seed of a
%                gallery row that names none
%
%   titles        1 x m cellstr: the tile titles of
%                 fig_supp3_connectivity_gallery ('Reference', 'E-to-E × 1.5',
%                 ...), by the row's key; a key not listed here is titled
%                 by its network's label
%   ostruct_keys  1 x m cellstr: the key the same network has in the
%                 outlier_structures stage (outlier_structures_stage_plan):
%                 its name when it is drawn on the stage's seed,
%                 '<name>_w<a>' when the row fixes the weight seed [a b]
%                 ('reference_w27'). Whether that stage ran the network is
%                 for the caller to check against its keys.
%
% One place for both, so that S3's tiles and S7 panel B (fig_supp7_grouped,
% which draws the gallery's networks under these titles) cannot drift
% apart.
%
% See also: cx_load, concept_examples_stage_plan, fig_supp3_connectivity_gallery,
%           fig_supp7_grouped, outlier_structures_stage_plan

% 'Reference' and 'E-to-E × 1.5' (Tom, 2026-11-23: "Reference and E-to-E × 1.5
% everywhere"; until then 'Balanced' and 'E-to-E up'), as Figures 1, 5, S6
% and S7 panel A name the two networks.
TITLE = struct('balanced', 'Reference', 'EE_up', ['E-to-E ' char(215) ' 1.5'], 'EE_down', 'E-to-E down', 'IE_down', 'I-to-E down', ...
    'II_up', 'I-to-I up', 'EI_loop', 'E-I loop', 'local_outlier', 'Local outlier', ...
    'lowrank_real', 'Real low-rank mode', 'lowrank_rot', 'Rotating low-rank mode');
m = numel(X);
titles = cell(1, m);
ostruct_keys = cell(1, m);
for k = 1:m
    if isfield(TITLE, X(k).key); titles{k} = TITLE.(X(k).key); else; titles{k} = X(k).label; end
    if nargin > 1
        ostruct_keys{k} = X(k).name;
        if ~isequal(double(X(k).seed(:)'), double(stage_seed(:)'))
            ostruct_keys{k} = sprintf('%s_w%d', X(k).name, X(k).seed(1));
        end
    end
end
end
