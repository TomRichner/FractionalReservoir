function [edge, band, by_blocks] = ospec_kinetic_edge(J, ws)
%OSPEC_KINETIC_EDGE The rightmost eigenvalue of the slow block of a full Jacobian: what the slow variables would do with the recurrent weights removed.
%
%   [edge, band, by_blocks] = ospec_kinetic_edge(J, ws)
%
%   J     the sparse full Jacobian at one state (compute_Jacobian_fast)
%   ws    ospec_weight_stats(params): is_x, neuron_of_state
%   edge  largest real part of the eigenvalues of J(slow, slow), 1/s
%         (-Inf when the model has no slow variables)
%   band  1 x 3: the 5th, 50th and 95th percentile of those real parts
%         (NaN without slow variables)
%   by_blocks   true when the block structure below was used
%
% THE SLOW BLOCK IS BLOCK DIAGONAL BY NEURON. A neuron's adaptation and
% depression variables are driven by that neuron's rate alone, and its
% rate depends, among the slow variables, on its own adaptation only. So
% J(slow, slow) has no entry between the variables of two neurons, and its
% spectrum is the union of the spectra of one small matrix per neuron
% (about 7 x 7 with multiple timescales): a few ms for 500 neurons,
% against about 10 s for a dense eig of the 3500 x 3500 block. The
% structure is CHECKED on the matrix handed in; if an entry couples two
% neurons (a route model this function was not written for) the dense eig
% is taken instead and by_blocks is false.
%
% See also: ospec_full_spectrum, ospec_dense_spectrum, ospec_weight_stats

slow = find(~ws.is_x);
band = nan(1, 3);
by_blocks = true;
if isempty(slow)
    edge = -Inf;
    return
end
Jss = J(slow, slow);
owner = ws.neuron_of_state(slow);
[i, j, v] = find(Jss);
if any(owner(i) ~= owner(j))
    by_blocks = false;
    es = real(eig(full(Jss)));
else
    % Position of every slow variable within its own neuron's block.
    [owner_sorted, order] = sort(owner);
    first = [true; diff(owner_sorted) ~= 0];
    start = find(first);
    block_of = cumsum(first);                          % block number of each sorted variable
    pos_sorted = (1:numel(owner_sorted))' - start(block_of) + 1;
    pos = zeros(numel(slow), 1); blk = zeros(numel(slow), 1);
    pos(order) = pos_sorted; blk(order) = block_of;
    sizes = accumarray(block_of, 1);
    m = max(sizes);
    n_blk = numel(sizes);
    % Blocks padded to m x m with zeros: a padded row and column add an
    % eigenvalue 0, removed below by taking each block's own size.
    B = zeros(m, m, n_blk);
    B(sub2ind(size(B), pos(i), pos(j), blk(i))) = v;
    es = zeros(numel(slow), 1);
    at = 0;
    for q = 1:n_blk
        s = sizes(q);
        es(at + 1:at + s) = real(eig(B(1:s, 1:s, q)));
        at = at + s;
    end
end
es = sort(es);
edge = es(end);
band = es(max(1, round([0.05 0.5 0.95] * numel(es))))';
end
