function name = ospec_run_name(network, condition, seed)
%OSPEC_RUN_NAME The file name of one run of the outlier_spectra stage: ospec_<network>_<condition>_s<seed1>-<seed2>.mat
%
% See also: ospec_run, ospec_read_index
name = sprintf('ospec_%s_%s_s%d-%d.mat', network, condition, seed(1), seed(2));
end
