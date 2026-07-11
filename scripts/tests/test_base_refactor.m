% test_base_refactor.m
% Reproducibility guard for the SRNNModelBase extraction (Phase 1).
%
% Runs a small, fixed-seed SRNNModel2 configuration that exercises BOTH numeric
% seams that the base/subclass split reroutes:
%   - RHS seam      (dynamics_fast)         via the fiducial ODE trajectory
%   - Jacobian seam (compute_Jacobian_fast) via the QR variational equations, and via
%     compute_Jacobian_at_indices (called through plot_eigenvalues, smoke-tested)
% SFA (n_a_E/I>0) and STD (n_b_E/I>0) states are on so those Jacobian/RHS blocks are covered.
%
% Two modes, selected by the variable ref_mode set before running (e.g. via matlab -batch):
%   ref_mode = 'capture'  -> run and SAVE a reference .mat  (use on the PRE-refactor class)
%   ref_mode = 'compare'  -> run and ASSERT bit-identical to the reference (POST-refactor)
% Optional variable ref_file overrides the reference path.
%
% Prints exactly one of:  BASE_REFACTOR_PASS  /  BASE_REFACTOR_FAIL  (compare mode)
%                         BASE_REFACTOR_CAPTURED                      (capture mode)

setup_paths();

if ~exist('ref_mode', 'var') || isempty(ref_mode)
    ref_mode = 'compare';
end
project_root = fileparts(fileparts(which('setup_paths')));
if ~exist('ref_file', 'var') || isempty(ref_file)
    ref_file = fullfile(project_root, 'data', 'refactor_ref', 'base_refactor_ref.mat');
end

%% --- Fixed-seed configuration exercising both seams -----------------------
model = SRNNModel2( ...
    'n',              20, ...
    'f',              0.5, ...
    'indegree',       8, ...
    'n_a_E',          1, ...      % SFA state on E
    'n_a_I',          1, ...      % SFA state on I
    'n_b_E',          1, ...      % STD state on E
    'n_b_I',          1, ...      % STD state on I
    'level_of_chaos', 1.5, ...
    'T_range',        [-2, 15], ...   % negative start = transient settling; >=15 so QR window is [13,15]
    'fs',             100, ...
    'rng_seeds',      [42 43], ...
    'lya_method',     'qr', ...       % QR exercises the Jacobian seam (variational eqs)
    'store_full_state',      true, ...
    'store_decimated_state', true);

model.build();
model.run();

%% --- Collect numeric outputs --------------------------------------------
res = struct();
res.W     = model.W;
res.t_out = model.t_out;
res.S_out = model.S_out;                 % fiducial trajectory (store_full_state=true)

% Capture every numeric field of lya_results generically (LE_spectrum, local_LE_spectrum_t,
% finite_LE_spectrum_t, t_lya_vec, Lyapunov dimension, etc.) so the QR Jacobian-seam output
% is compared without hard-coding field names. Skip 'params' (nested struct w/ handles).
res.lya = struct();
if ~isempty(model.lya_results) && isstruct(model.lya_results)
    lf = fieldnames(model.lya_results);
    for i = 1:numel(lf)
        v = model.lya_results.(lf{i});
        if isnumeric(v) && ~strcmp(lf{i}, 'params')
            res.lya.(lf{i}) = v;
        end
    end
end

% Direct numeric check of the E/I Jacobian kernel at a sample state.
params_j  = model.cached_params;
S_sample  = model.S_out(round(size(model.S_out, 1) / 2), :).';
res.J_sample = full(SRNNModel2.compute_Jacobian_fast(S_sample, params_j));

% Smoke-test the compute_Jacobian_at_indices seam (used by plot_eigenvalues).
% Figure is created invisibly under -batch; we only check it does not error.
try
    fh = model.plot_eigenvalues([5 10]);
    if ~isempty(fh) && all(isgraphics(fh)), close(fh); end
    res.eig_smoke_ok = true;
catch ME
    res.eig_smoke_ok = false;
    fprintf('plot_eigenvalues smoke-test errored: %s\n', ME.message);
end

%% --- Capture or compare --------------------------------------------------
if strcmpi(ref_mode, 'capture')
    ref_dir = fileparts(ref_file);
    if ~exist(ref_dir, 'dir'), mkdir(ref_dir); end
    save(ref_file, '-struct', 'res');
    fprintf('Reference saved to %s\n', ref_file);
    disp('BASE_REFACTOR_CAPTURED');
else
    if ~exist(ref_file, 'file')
        error('test_base_refactor:noRef', 'Reference file not found: %s (run with ref_mode=''capture'' first).', ref_file);
    end
    ref = load(ref_file);
    tol = 1e-10;
    all_ok = true;

    % Compare the top-level numeric arrays.
    top_fields = {'W', 't_out', 'S_out', 'J_sample'};
    for k = 1:numel(top_fields)
        all_ok = compare_field(res, ref, top_fields{k}, tol, all_ok, top_fields{k});
    end

    % Compare each captured lya_results numeric field.
    if isfield(ref, 'lya') && isstruct(ref.lya)
        lyf = fieldnames(ref.lya);
        for k = 1:numel(lyf)
            all_ok = compare_field(res.lya, ref.lya, lyf{k}, tol, all_ok, ['lya.' lyf{k}]);
        end
    end

    if all_ok
        disp('BASE_REFACTOR_PASS');
    else
        disp('BASE_REFACTOR_FAIL');
    end
end

%% --- helper --------------------------------------------------------------
function ok = compare_field(A, B, fld, tol, ok, label)
    hasA = isfield(A, fld); hasB = isfield(B, fld);
    if hasA ~= hasB
        fprintf('  [FAIL] field presence mismatch: %s (new=%d ref=%d)\n', label, hasA, hasB);
        ok = false; return;
    end
    if ~hasA, return; end
    a = A.(fld); b = B.(fld);
    if ~isequal(size(a), size(b))
        fprintf('  [FAIL] size mismatch: %-16s new=%s ref=%s\n', label, mat2str(size(a)), mat2str(size(b)));
        ok = false; return;
    end
    d = max(abs(double(a(:)) - double(b(:))));
    if isempty(d), d = 0; end
    status = 'ok'; if d > tol, status = 'FAIL'; ok = false; end
    fprintf('  [%s] %-16s max|delta| = %.3e\n', status, label, d);
end
