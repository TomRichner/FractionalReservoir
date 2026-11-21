function name = plane_states_name(network, condition, seed)
%PLANE_STATES_NAME The file name of one run's per-neuron operating points (plane_states), beside outlier_run_name's.
name = sprintf('states_%s_%s_s%d-%d.mat', network, condition, seed(1), seed(2));
end
