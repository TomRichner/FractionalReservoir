function s = p1a_signals(m, info, groups)
%P1A_SIGNALS Population signals of one run of the run13 network (pilot 1a).
%
%   s = p1a_signals(m, info, groups)
%
% m is a run IEDExplorationNetwork from run13_condition_model, groups and info its second and third
% outputs. Returns, on the plot grid (200 Hz, 0..duration):
%   t           1 x nt time (s)
%   x, r        n x nt dendritic state and rate, E rows first, then I
%   groups      n x 1 module label (1..5)
%   E, I        row indices of each type in x and r
%   rate_E, rate_I     1 x nt population mean rate
%   sfa_E, sfa_I       1 x nt mean SFA term (c/K) sum_k a_k of each type
%                      (pilot 0a's definition; K = 1 for 'original')
%   sfa_doc            1 x nt sum_q c_q sum a / n, the "SFA feedback" of
%                      ied_stim_experiment (equal to the n-weighted mean of
%                      sfa_E and sfa_I when K = 1; K times it otherwise).
%                      Used only for the regression check of 'original'.
%   dep                1 x nt mean depression factor of the E->E route
%                      (= E->I), from the model's std_route_factor
%   u                  1 x nt external input to neuron 1 (uniform input)
%   module_x           G x nt module-mean x

p = m.plot_data; params = m.get_params();
s.t = p.t(:)';
s.x = [p.x.E; p.x.I]; s.r = [p.r.E; p.r.I];
s.E = (1:size(p.x.E, 1))'; s.I = size(p.x.E, 1) + (1:size(p.x.I, 1))';
% Rows of x follow type_indices order: check the modules line up.
assert(isequal(m.type_indices{1}(:), s.E) && isequal(m.type_indices{2}(:), s.I), ...
    'Type indices are not E-first contiguous; the groups vector would not line up.');
s.groups = groups(:);
s.rate_E = mean(p.r.E, 1); s.rate_I = mean(p.r.I, 1);
nt = numel(s.t); sfa = zeros(2, nt); s.sfa_doc = zeros(1, nt);
for q = 1:2
    a = p.a.(m.cell_type_names{q});                 % n_q x K x nt
    K = size(a, 2);
    total = reshape(sum(a, [1 2]), 1, []);
    sfa(q, :) = info.params.c(q) / K * total / size(a, 1);
    s.sfa_doc = s.sfa_doc + info.params.c(q) * total / m.n;
end
s.sfa_E = sfa(1, :); s.sfa_I = sfa(2, :);
D = SRNNCellTypePairs.std_route_factor(p.b.E.E, params, 1, 1);   % n_E x 1 x nt
s.dep = reshape(mean(D, 1), 1, []);
DI = SRNNCellTypePairs.std_route_factor(p.b.E.I, params, 1, 2);
assert(isequal(D, DI), 'E->E and E->I depression differ; the E->E factor would not stand for both.');
s.u = interp1(m.t_ex(:)', m.u_ex(1, :), s.t);
G = max(s.groups); s.module_x = zeros(G, nt);
for g = 1:G, s.module_x(g, :) = mean(s.x(s.groups == g, :), 1); end
end
