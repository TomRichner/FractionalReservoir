function [h_local, h_finite] = plot_local_and_finite_lle(ax, lya_results, colour, opts)
% PLOT_LOCAL_AND_FINITE_LLE Leading local rate and accumulating finite-time lambda_1.
%
%   [h_local, h_finite] = PLOT_LOCAL_AND_FINITE_LLE(ax, lya_results, colour)
%   [h_local, h_finite] = PLOT_LOCAL_AND_FINITE_LLE(ax, lya_results, colour, opts)
%
% Draws, on axes ax, the two quantities the paper keeps apart: the LOCAL
% expansion rate of the leading direction (thin grey; positive stretches are
% "transient expansion") and the ACCUMULATING finite-time lambda_1 (thick, in
% colour; its end value is the exponent the paper reports, "stable" or
% "chaotic" only over the stated window). A zero line, the y label, and the
% final lambda_1 printed in the corner.
%
% Accepts either estimator's result struct:
%   top-K / QR : local_LE_spectrum_t(:, 1), finite_LE_spectrum_t(:, 1), t_lya
%                (finite is NaN before accumulation starts; those samples
%                are not drawn)
%   Benettin   : local_lya, finite_lya, t_lya. finite_lya is the class's own
%                running estimate over the accumulation window (NaN before
%                it); if a result lacks it, the cumulative mean of local_lya
%                is drawn instead and labelled as such in the header comment
%                of the calling figure.
%
% opts.ylim        [lo hi], optional; default from the data with the zero
%                  line kept in view.
% opts.local_color default [0.55 0.55 0.55]
% opts.label       logical, default true: print lambda_1 = ... in the corner.
%
% See also: fig_example_timeseries, lyapunov_topk, plot_lyapunov

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'local_color'), opts.local_color = [0.55 0.55 0.55]; end
if ~isfield(opts, 'label'),       opts.label = true; end

t = lya_results.t_lya(:);
if isfield(lya_results, 'local_LE_spectrum_t')
    local  = lya_results.local_LE_spectrum_t(:, 1);
    finite = lya_results.finite_LE_spectrum_t(:, 1);
elseif isfield(lya_results, 'local_lya')
    local = lya_results.local_lya(:);
    if isfield(lya_results, 'finite_lya')
        finite = lya_results.finite_lya(:);
    else
        finite = cumsum(local) ./ (1:numel(local))';   % cumulative mean fallback
    end
else
    error('plot_local_and_finite_lle:NoSeries', ...
        'lya_results carries neither local_LE_spectrum_t nor local_lya.');
end
n = min([numel(t), numel(local), numel(finite)]);
t = t(1:n); local = local(1:n); finite = finite(1:n);

hold(ax, 'on');
h_local  = plot(ax, t, local, '-', 'Color', opts.local_color, 'LineWidth', 0.6);
ok = isfinite(finite);
h_finite = plot(ax, t(ok), finite(ok), '-', 'Color', colour, 'LineWidth', 2.2);
yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
hold(ax, 'off');
ylabel(ax, '\lambda (s^{-1})');
if isfield(opts, 'ylim') && ~isempty(opts.ylim)
    ylim(ax, opts.ylim);
else
    lo = min([local; 0]); hi = max([local; 0]);
    pad = 0.05 * max(hi - lo, eps);
    ylim(ax, [lo - pad, hi + pad]);
end
if opts.label && isfield(lya_results, 'LLE')
    text(ax, 0.98, 0.95, sprintf('\\lambda_1 = %+.3f s^{-1}', lya_results.LLE), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Color', colour, 'FontSize', 9);
end
end
