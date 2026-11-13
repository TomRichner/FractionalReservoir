function [delay_s, hold_len, d_max] = mc_delay_axis(P)
% MC_DELAY_AXIS The delays a memory-capacity protocol scores, in seconds.
%
%   [delay_s, hold_len, d_max] = MC_DELAY_AXIS(P)
%
%   P         a protocol struct with mc_protocol's fields fs, d_max_sec,
%             input_type and T_hold
%   delay_s   1 x n_delay, the delay of each scored readout (s)
%   hold_len  samples per held input value (1 unless input_type is
%             'sample_hold')
%   d_max     the longest delay in samples, SRNN_ESN_reservoir's d_max
%
% Sample-and-hold input is scored in HOLD units, one delay per hold; every
% other input type at sample resolution. This is the grid
% SRNN_ESN_reservoir.mc_readout_plan builds once a reservoir exists. The
% stages need it earlier (to size their outputs and to check a result's
% delay count), and a reader of job files needs it with no reservoir at
% all, so it is written once here for run_memory_capacity,
% run_mc_gain_sweep and read_mc_gain_sweep_jobs.
%
% See also: mc_protocol, SRNN_ESN_reservoir.mc_readout_plan

arguments
    P (1,1) struct
end

d_max = round(P.d_max_sec * P.fs);
hold_len = 1;
if strcmpi(P.input_type, 'sample_hold')
    hold_len = max(1, round(P.T_hold * P.fs));
end
n_delay = max(1, floor(d_max / hold_len));
delay_s = (1:n_delay) * (hold_len / P.fs);
end
