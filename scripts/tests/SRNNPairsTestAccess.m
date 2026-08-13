classdef SRNNPairsTestAccess < SRNNCellTypePairs
    % SRNNPairsTestAccess Test-only subclass exposing SRNNCellTypePairs's
    % protected statics, the counterpart of SRNN2TestAccess.
    %
    % Only lyapunov_sample_grid is needed so far: test_lya_window checks that
    % the two duck-typed sibling classes lay out identical Lyapunov segment
    % grids, which is the property that was broken before the window fix.
    % (dynamics_fast is already public on SRNNCellTypePairs, so it needs no
    % passthrough.)
    %
    % Never constructed -- the statics are reached through the class name.
    methods (Static)
        function [idx, t_lya, acc_start] = lya_grid(t, dt, deci, tau, T, warmup)
            [idx, t_lya, acc_start] = SRNNCellTypePairs.lyapunov_sample_grid( ...
                t, dt, deci, tau, T, warmup);
        end
    end
end
