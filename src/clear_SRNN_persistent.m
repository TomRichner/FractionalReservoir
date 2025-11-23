function clear_SRNN_persistent()
% CLEAR_SRNN_PERSISTENT Clears the persistent variables in SRNN_reservoir
%
% Usage:
%   clear_SRNN_persistent()
%
% This function clears the SRNN_reservoir function from memory, which
% resets its persistent variables (u_interpolant and t_ex_last).
% This is necessary when running multiple simulations where the input u_ex
% changes but the time vector t_ex might remain the same, as SRNN_reservoir
% only checks for changes in t_ex to decide whether to rebuild the interpolant.

    clear SRNN_reservoir;
    clear SRNN_reservoir_DDE;

end

