function [model, params, ws] = ospec_rebuild_model(R)
%OSPEC_REBUILD_MODEL The model of one run of the outlier_spectra stage, built again from what the run saved, without integrating.
%
%   [model, params, ws] = ospec_rebuild_model(R)      R = ospec_load_run(run_dir, row)
%
%   model    the built OspecNet (an SRNNCellTypePairs): the same weights,
%            setpoints, time constants and initial state as the run's
%   params   model.get_params(): what the class's static functions take,
%            for instance the full Jacobian at a saved calibration state,
%              J = SRNNCellTypePairs.compute_Jacobian_fast(R.cal.S(k, :)', params)
%            (the rows and columns of each kind of variable are in
%            params.state_layout: x, a{type}, b{pre, post}, g{pre, post})
%   ws       ospec_weight_stats(params): the full W, the type indices, the
%            realized block statistics; with it,
%              M = ospec_effective_matrix(R.cal.S(k, :)', params, ws)
%            is the frozen-gain effective connectivity whose eigenvalues
%            are R.cal.ev_M(:, k)
%
% It takes a fraction of a second at 500 neurons. Nothing is simulated:
% the states a consumer needs are the saved ones (R.cal.S).
%
% WHAT IS CHECKED. That the rebuilt model has the run's size (R.n, R.N). The
% test of the stage (test_outlier_spectra_stage) asserts that the spectra
% recomputed from the saved states on the rebuilt model are the saved ones.
%
% See also: ospec_load_run, ospec_build_model, ospec_effective_matrix, ospec_run

arguments
    R (1,1) struct
end
model = ospec_build_model(R.build);
if model.n ~= R.n || model.N_sys_eqs ~= R.N
    error('ospec_rebuild_model:Size', ...
        'The rebuilt model has n = %d and %d state variables; the run saved n = %d and %d. Has the preset %s changed?', ...
        model.n, model.N_sys_eqs, R.n, R.N, R.preset_name);
end
if nargout >= 2; params = model.get_params(); end
if nargout >= 3; ws = ospec_weight_stats(params); end
end
