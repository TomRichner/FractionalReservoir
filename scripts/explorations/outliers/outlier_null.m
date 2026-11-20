function c = outlier_null(S, params, ws, n_null, n_shuffle)
%OUTLIER_NULL Where the bulk ends at one state: null draws with the observed gains, and the gain controls.
%
%   c = outlier_null(S, params, ws, n_null, n_shuffle)
%
% THE NULL FOR THE BULK. The bulk is what a matrix with no mean structure
% and these gains would have. Each draw takes the realized W, removes the
% realized block means (after sparsity, so the variance the mask makes from
% the block means stays in: 95% of the bulk variance at the reference
% network), permutes the entries within each block (same entry
% distribution, same number of synapses, no relation left between a
% neuron's weights and its gain, no low-rank part), and applies the
% observed gains to the columns. The rightmost real part and the largest
% modulus of each draw's spectrum are kept; their 95th percentiles are the
% calibrated right edge and radius at this state, to set beside the
% predicted R_eff.
%
% THE GAIN CONTROLS, on the actual W (the note's controls, and the sibling
% report's shuffle):
%   shuffle   gains permuted among the neurons of each type (keeps the mean
%             and the spread of the gains by type, breaks which neuron has
%             which gain)
%   typemean  every neuron given its type's mean gain (no spread)
%   uniform   every neuron given the grand mean gain (no E-I difference:
%             the mean structure back at its structural ratio to the bulk)
%
% Fields
%   R_eff            predicted radius at this state
%   null_re, null_abs    n_null x 1: rightmost real part, largest modulus
%   q95_re, q95_abs, med_re, med_abs    their summaries
%   obs_re           rightmost real part of the actual M
%   shuf_re          n_shuffle x 1: rightmost real part under the shuffle
%   obs_pct_in_shuffle   fraction of shuffles whose rightmost is BELOW the
%                    observed one (near 0: the actual assignment of gains
%                    to neurons keeps the edge further in than chance)
%   typemean_re, uniform_re   rightmost real part under those controls
%
% See also: outlier_effective_matrix, outlier_measures

if nargin < 4 || isempty(n_null);    n_null = 40;    end
if nargin < 5 || isempty(n_shuffle); n_shuffle = 20; end

[Mf, d, ~, pred] = outlier_effective_matrix(S, params, ws);
C  = params.n_cellTypes;
n  = params.n;
ti = ws.ti;

% Gains as an n x n mask: row i of type a, column j gets d(a, j).
Dmat = zeros(n);
for a = 1:C
    Dmat(ti{a}, :) = repmat(d(a, :), numel(ti{a}), 1);
end

% Centered W, block by block.
Wc = ws.Wf;
for a = 1:C
    for b = 1:C
        Wc(ti{a}, ti{b}) = Wc(ti{a}, ti{b}) - ws.mu(a, b);
    end
end

null_re = zeros(n_null, 1); null_abs = zeros(n_null, 1);
for k = 1:n_null
    Wn = Wc;
    for a = 1:C
        for b = 1:C
            blk = Wc(ti{a}, ti{b});
            blk(:) = blk(randperm(numel(blk)));
            Wn(ti{a}, ti{b}) = blk;
        end
    end
    ev = eig(Wn .* Dmat);
    null_re(k)  = max(real(ev));
    null_abs(k) = max(abs(ev));
end

obs_re = max(real(eig(Mf)));

shuf_re = zeros(n_shuffle, 1);
for k = 1:n_shuffle
    Dsh = Dmat;
    for b = 1:C
        cols = ti{b};
        Dsh(:, cols) = Dmat(:, cols(randperm(numel(cols))));
    end
    shuf_re(k) = max(real(eig(ws.Wf .* Dsh)));
end

Dtm = Dmat;
for a = 1:C
    for b = 1:C
        Dtm(ti{a}, ti{b}) = mean(d(a, ti{b}));
    end
end
typemean_re = max(real(eig(ws.Wf .* Dtm)));
uniform_re  = max(real(eig(ws.Wf .* mean(d(:)))));

c = struct('R_eff', pred.R_eff, 'null_re', null_re, 'null_abs', null_abs, ...
    'q95_re', quantile_(null_re, 0.95), 'q95_abs', quantile_(null_abs, 0.95), ...
    'med_re', median(null_re), 'med_abs', median(null_abs), ...
    'obs_re', obs_re, 'shuf_re', shuf_re, ...
    'obs_pct_in_shuffle', mean(shuf_re < obs_re), ...
    'typemean_re', typemean_re, 'uniform_re', uniform_re);
end

function q = quantile_(x, p)
% Quantile by linear interpolation of the sorted sample, without the
% Statistics Toolbox.
x = sort(x(:));
m = numel(x);
if m == 0; q = NaN; return; end
pos = min(max(p * m + 0.5, 1), m);
lo = floor(pos); hi = ceil(pos);
q = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
