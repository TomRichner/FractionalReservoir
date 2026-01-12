function [t_plot, S_plot, indices] = decimate_states(t_out, S_out, deci)
% DECIMATE_STATES Decimates state trajectory for plotting
%
% Inputs:
%   t_out - Time vector
%   S_out - State matrix (nt x N)
%   deci  - Decimation factor (integer)
%
% Outputs:
%   t_plot  - Decimated time vector
%   S_plot  - Decimated state matrix
%   indices - Indices used for decimation

    indices = 1:deci:length(t_out);
    t_plot = t_out(indices);
    S_plot = S_out(indices, :);
end





