function [ranges, factors] = mu_block_from_preset(preset_defaults)
% MU_BLOCK_FROM_PRESET Sweep ranges for the four connectivity blocks.
%
%   [ranges, factors] = MU_BLOCK_FROM_PRESET(preset_defaults)
%
% Returns an N x 2 cell {block_name, [lo hi]} covering mu_EE / mu_EI / mu_IE /
% mu_II, and the multiplier pair the ranges were built from.
%
% Each block is swept from a QUARTER to TRIPLE whatever the PRESET operates at
% rather than over fixed absolute numbers -- i.e. -75% to +200% of the default.
% mu_*_relative is a multiplier of F = 1/sqrt(n*alpha*(2-alpha)), so "the default
% level" is only meaningful relative to the preset -- and mu_tilde_relative is a
% REQUIRED constructor property with no class default to fall back on, which is
% why this reads the preset rather than ParamSpaceAnalysis2.class_default.
%
% Widened from [0.5, 2.0] (-50% to +100%). Note the default does NOT sit at the
% centre of the resulting linear axis and is not meant to: 0.25x-3x is roughly
% symmetric in RATIO, so the preset sits about a fifth of the way along. The
% percent ruler apply_percent_axis draws is what makes that readable.
%
% WAS a local subfunction of run_sensitivity_analysis, returning one block at a
% time. Promoted to its own file -- returning all four, WITH the multipliers --
% when run_param_space_analysis gained the same four axes: the two analyses must
% sweep them over IDENTICAL ranges for a grid point and a sensitivity level to
% describe the same network, and two copies of the loop would be two chances to
% drift apart. This function is the single place those ranges are decided.
%
% See also: run_sensitivity_analysis, run_param_space_analysis

factors = [0.25, 3.0];
block_names = {'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'};

ranges = cell(numel(block_names), 2);
for b = 1:numel(block_names)
    base = block_base(preset_defaults, block_names{b});
    % SORT, because add_grid_parameter rejects a descending range and the
    % inhibitory blocks are negative: 0.25x and 3x of -3 is [-0.75, -9]. After
    % sorting, an inhibitory axis runs from STRONGEST to weakest inhibition,
    % i.e. the 3x end is on the left. Same multipliers either way.
    ranges{b, 1} = block_names{b};
    ranges{b, 2} = sort(factors * base);
end
end

%% ------------------------------------------------------------------------
function v = block_base(preset_defaults, block_name)
% One entry of a preset's mu_tilde_relative, named the way the scalar aliases
% are: mu_EI is "E receives from I", i.e. (post, pre) = (1, 2).
%
% Handles either shape the block may be given in. SRNNCellTypePairs.expand_block
% accepts a full C x C matrix or a 1 x C PRESYNAPTIC row broadcast down the
% columns -- so for a row, the entry depends only on the presynaptic index and
% mu_EE == mu_IE. Errors rather than guessing.
if ~isfield(preset_defaults, 'mu_tilde_relative') || isempty(preset_defaults.mu_tilde_relative)
    error('mu_block_from_preset:NoMuBlock', ...
        ['Sweeping %s relative to its default needs the preset to set ' ...
        'mu_tilde_relative; SRNNCellTypePairs has no class default for it.'], ...
        block_name);
end
idx = struct('mu_EE_relative', [1 1], 'mu_EI_relative', [1 2], ...
             'mu_IE_relative', [2 1], 'mu_II_relative', [2 2]);
if ~isfield(idx, block_name)
    error('mu_block_from_preset:UnknownMuBlock', ...
        'Unknown block ''%s''.', block_name);
end
ij = idx.(block_name);
M = preset_defaults.mu_tilde_relative;
if isvector(M)
    M = repmat(reshape(M, 1, []), numel(M), 1);   % presynaptic row -> full block
end
v = M(ij(1), ij(2));
end
