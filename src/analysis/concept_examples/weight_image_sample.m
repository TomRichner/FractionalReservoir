function S = weight_image_sample(Wf, type_of, n_show, seed)
%WEIGHT_IMAGE_SAMPLE A representative n_show x n_show piece of a weight matrix: the submatrix on n_show/2 E and n_show/2 I neurons.
%
%   S = weight_image_sample(Wf, type_of, n_show, seed)
%
%   Wf       n x n weight matrix, (post, pre), the network's own (full)
%   type_of  n x 1, 1 for E and 2 for I
%   n_show   even number of neurons shown (30 in Figure 1)
%   seed     seed of the private stream the neurons are drawn from (the
%            global stream is not touched)
%
%   S.W        n_show x n_show: Wf(idx, idx), E neurons first. Rows are the
%              postsynaptic neurons, columns the presynaptic ones, so the
%              first n_show/2 COLUMNS are excitatory (positive, red in
%              draw_weight_image) and the rest inhibitory (blue); an absent
%              connection is 0 (white).
%   S.idx      the neurons, in the network's numbering
%   S.type     n_show x 1 (1 E, 2 I)
%   S.density  realized fraction of nonzero entries, per (post, pre) block
%              (2 x 2) and overall (S.density_all)
%   S.block_mean, S.block_sd   per block, of all entries, zeros included
%
% WHY A SUBMATRIX. Its entries are the network's own, so the connection
% probability and the block statistics are the network's in expectation:
% the same E-to-E mean raised by half, the same sparsity, the same
% distribution of nonzero weights. Anything else (resampling, a smaller
% network of the ensemble) would be a different matrix to explain. What
% the piece cannot show is structure spread thinly over all n x n
% entries: an added low-rank pattern of norm comparable to the bulk adds
% entries of order R |a + iw| / n to every pair, a few percent of a sparse
% weight, and is nearly invisible at the color scale of the sparse part.
%
% See also: draw_weight_image, run_concept_examples_stage

arguments
    Wf (:,:) double
    type_of (:,1) double
    n_show (1,1) double {mustBeInteger, mustBePositive}
    seed (1,1) double = 1
end
if mod(n_show, 2) ~= 0
    error('weight_image_sample:Odd', 'n_show must be even (half E, half I).');
end
iE = find(type_of == 1); iI = find(type_of == 2);
h = n_show / 2;
if numel(iE) < h || numel(iI) < h
    error('weight_image_sample:TooFew', 'The network has %d E and %d I neurons; %d of each are asked for.', numel(iE), numel(iI), h);
end
rs = RandStream('mt19937ar', 'Seed', seed);
pE = sort(iE(rs.randperm(numel(iE), h)));
pI = sort(iI(rs.randperm(numel(iI), h)));
idx = [pE; pI];
W = Wf(idx, idx);
type = [ones(h, 1); 2 * ones(h, 1)];
dens = zeros(2); bm = zeros(2); bs = zeros(2);
for a = 1:2
    for b = 1:2
        blk = W(type == a, type == b);
        dens(a, b) = nnz(blk) / numel(blk);
        bm(a, b) = mean(blk(:));
        bs(a, b) = std(blk(:), 1);
    end
end
S = struct('W', W, 'idx', idx, 'type', type, 'density', dens, 'density_all', nnz(W) / numel(W), ...
    'block_mean', bm, 'block_sd', bs, 'seed', seed);
end
