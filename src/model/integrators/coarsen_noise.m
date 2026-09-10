function [xi1_c, xi2_c] = coarsen_noise(xi1, xi2, h, m)
%COARSEN_NOISE Rebuild the SAME Brownian path on a step m times coarser.
%   [xi1_c, xi2_c] = coarsen_noise(xi1, xi2, h, m)
%
%   xi1, xi2 are the unit-variance normals a noise struct for SDE_FIXED_STEP
%   carries (one row per noise-driven state, one column per fine step of
%   length h). The result is the pair for the SAME path on step H = m*h, with
%   m a power of 2 and size(xi1, 2) a multiple of m. Together with a new
%   noise struct at fs/m (same t0, same sigma, same idx) it lets a coarse-step
%   run consume exactly the Brownian motion a fine-step run consumed, so the
%   difference between the two is pure discretisation error -- the strong
%   error a convergence study measures.
%
%   Do not DECIMATE the columns instead: a subsampled path is a different
%   path with the wrong variance.
%
%   Increments simply add: dW_c = sum of the m fine increments. The second
%   integral does not -- over two consecutive fine steps,
%     I_c = I_a + h*dW_a + I_b
%   because the second sub-interval's area is measured from a base that has
%   already risen by dW_a. Applying that pairwise, log2(m) times, is exact.
%   The pair is then converted back to unit-variance normals by inverting the
%   Kloeden-Platen identity I = (H/2)*(dW + dZ/sqrt(3)).
%
%   Moved here from a local function of scripts/tests/test_sde_integrators.m
%   (2026-09-10) so run_numerics_verification can share it.
%
%   See also SDE_FIXED_STEP.

    assert(mod(log2(m), 1) == 0, 'coarsen_noise:NotPowerOfTwo', ...
        'coarsen_noise expects m to be a power of 2, got %g.', m);
    assert(mod(size(xi1, 2), m) == 0, 'coarsen_noise:ColumnsNotDivisible', ...
        'coarsen_noise needs size(xi1, 2) = %d to be a multiple of m = %d.', ...
        size(xi1, 2), m);
    dW = sqrt(h) * xi1;
    I10 = (h / 2) * (dW + sqrt(h) * xi2 / sqrt(3));
    hc = h;
    while size(dW, 2) > size(xi1, 2) / m
        a = 1:2:size(dW, 2);
        b = 2:2:size(dW, 2);
        I10 = I10(:, a) + hc * dW(:, a) + I10(:, b);
        dW = dW(:, a) + dW(:, b);
        hc = 2 * hc;
    end
    H = hc;
    xi1_c = dW / sqrt(H);
    xi2_c = sqrt(3) * (2 * I10 / H - dW) / sqrt(H);
end
