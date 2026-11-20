function name = outlier_run_name(network, condition, seed, tag)
%OUTLIER_RUN_NAME The file name of one outlier run: run_<network>_<condition>_s<seed1>-<seed2><tag>.mat
%
% See also: outlier_run, explore_outliers
if nargin < 4; tag = ''; end
name = sprintf('run_%s_%s_s%d-%d%s.mat', network, condition, seed(1), seed(2), tag);
end
