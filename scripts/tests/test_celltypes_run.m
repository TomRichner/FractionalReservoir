% test_celltypes_run.m
% Build / run / stability / STF-sanity / plot smoke test for SRNNModelCellTypes.
% Prints CELLTYPES_RUN_PASS / CELLTYPES_RUN_FAIL.

setup_paths();
set(0, 'DefaultFigureVisible', 'off');
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

%% 1. Build + run a data-parameterized model (Benettin Lyapunov)
model = SRNNModelCellTypes('n', 40, 'n_a', 1, 'n_b', 1, 'n_u', 1, ...
    'T_range', [-2, 8], 'fs', 200, 'rng_seeds', [3 4], ...
    'level_of_chaos', 1.5, 'rescale_by_abscissa', true, ...
    'lya_method', 'benettin', 'store_full_state', true);
model.build();
model.run();

pd = model.plot_data;
rvals = pd.r(:);
rstd  = std(rvals); rmean = mean(rvals);
alive = rstd > 1e-3 && rmean > 0.01 && rmean < 0.99 && all(isfinite(rvals));
report(sprintf('run completes, rates alive (mean=%.3f std=%.3f)', rmean, rstd), alive);
ok = ok && alive;

haveLLE = isfield(model.lya_results, 'LLE') && isfinite(model.lya_results.LLE);
report(sprintf('Benettin LLE finite (%.4f)', model.lya_results.LLE), haveLLE);
ok = ok && haveLLE;

%% 2. STF/STD sanity: under activity, b depresses (<1) and p facilitates (>p0)
p = model.cached_params;
b_end = pd.b(:, :, end);       % n x K
p_end = pd.p(:, :, end);       % n x K
b_dep = mean(b_end(:)) < 1 - 1e-3;
p_fac = mean(p_end(:)) > mean(p.p0_mat(:)) + 1e-3;
report(sprintf('STD depresses (mean b_end=%.3f < 1)', mean(b_end(:))), b_dep);
report(sprintf('STF facilitates (mean p_end=%.3f > mean p0=%.3f)', mean(p_end(:)), mean(p.p0_mat(:))), p_fac);
ok = ok && b_dep && p_fac;

%% 3. Plot smoke test
try
    fh = model.plot();
    if ~isempty(fh) && all(isgraphics(fh)), close(fh); end
    report('plot() renders', true);
catch ME
    report(['plot ERROR: ' ME.message], false); ok = false;
end

%% 4. QR Lyapunov path (small n) — exercises the analytic Jacobian end-to-end
try
    m2 = SRNNModelCellTypes('n', 12, 'n_a', 1, 'n_b', 1, 'n_u', 1, ...
        'T_range', [-2, 6], 'fs', 100, 'rng_seeds', [5 6], ...
        'level_of_chaos', 1.5, 'rescale_by_abscissa', true, 'lya_method', 'qr');
    m2.build(); m2.run();
    qr_ok = isfield(m2.lya_results, 'LE_spectrum') && all(isfinite(m2.lya_results.LE_spectrum));
    report(sprintf('QR spectrum finite (LLE=%.4f)', max(m2.lya_results.LE_spectrum)), qr_ok);
    ok = ok && qr_ok;
catch ME
    report(['QR ERROR: ' ME.message], false); ok = false;
end

if ok
    disp('CELLTYPES_RUN_PASS');
else
    disp('CELLTYPES_RUN_FAIL');
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end
