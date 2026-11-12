function model = route_model(name)
%ROUTE_MODEL The registry of route models: the one place a model name is resolved.
%
%   model = ROUTE_MODEL('product')    the model, a struct of function handles
%   model = ROUTE_MODEL('pool')
%   names = ROUTE_MODEL()             the valid names, a cell row
%
% A ROUTE is one (presynaptic type, postsynaptic type) pair of a
% SRNNCellTypePairs network. Its synapses carry a small local dynamical
% system, one copy per presynaptic neuron, driven by that neuron's rate r,
% and what the route hands the network is an output factor F: the synaptic
% output is theta = r .* F. A ROUTE MODEL is one such system (its state,
% its equations, its output factor and their derivatives). The model class
% knows nothing about any of them: dynamics_fast, compute_Jacobian_fast,
% jacobian_times and unpack_and_compute_states call the handles below once
% per route, identically for every model (docs/Plan_A_chained_vessicles.md,
% section 4.0).
%
% THE INTERFACE. Every model file returns a struct with the field `name`
% and these handles. npre is the number of presynaptic neurons, ns the
% number of state variables per neuron; a state s is npre x ns, ONE ROW PER
% NEURON; the rate r is npre x 1.
%
%   P  = compile(std_cfg, stf_cfg, opts)
%          Validates one route's configuration (the `std` and `stf` structs
%          of synapse_config.<pre>.<post>; [] or a struct with no fields
%          means absent) and returns the route's parameters. opts may hold
%          npre (default 1), label (the route's name in error messages,
%          'E->I'), err_id_prefix (the class named in error identifiers)
%          and model options. P always has
%            model      the model's name
%            n_b, n_g   columns of s that live in the state vector's
%                       b-block and g-block, in that order (ns = n_b + n_g)
%            ns, npre
%            nominal    the configured values, as rows. It always has
%                       tau_rec, tau_rel, tau_dec, tau_fac and G ([] where
%                       the route has none): SRNNCellTypePairs copies them
%                       into the params fields code outside the class
%                       reads. A model adds its own (the product's
%                       combine; the pool's n2 and std_model)
%            signature  a struct that is isequal for two routes exactly
%                       when they carry the same dynamics
%            coupling   ns x ns logical: (m, k) is true where ds_m/dt can
%                       depend on s_k, that is, where A(:, m, k) can be
%                       nonzero. It must not understate: an entry stated
%                       false is never read. Where it is diagonal (no
%                       variable's derivative depends on another), the
%                       class applies A in jacobian_times as ns products
%                       and not ns^2.
%          and the parameters the equations read, PER NEURON (npre x k),
%          the nominal row replicated: a per-neuron spread fills in other
%          values and needs no other code path.
%   s0 = initial(P)            1 x ns, the state every neuron starts in
%   [ds, F] = rhs(s, r, P)     ds npre x ns, the state's time derivative;
%                              F npre x 1, the output factor
%   F  = output(s, P)          the output factor alone. s may be
%                              npre x ns x nt (nt time points); F is then
%                              npre x 1 x nt
%   [A, B, F, Fs] = linearize(s, r, P)
%                              A   npre x ns x ns   d(ds_m)/d(s_k) at (i, m, k)
%                              B   npre x ns        d(ds_m)/dr
%                              F   npre x 1         the output factor
%                              Fs  npre x ns        dF/d(s_k)
%   ss = steady_state(r, P)    closed form at a constant rate r (npre x 1,
%                              or a scalar): ss.s npre x ns, ss.F npre x 1
%   D  = depression(b, P)      the depression factor from the route's
%                              b-block alone, for display code that holds a
%                              b array and not the whole state (b npre x
%                              n_b, or npre x n_b x nt, n_b >= 1; D has the
%                              shape of prod(b, 2)). The product's is
%                              prod(b) after its combine and zero floor;
%                              the pool's is b1. The one caller is
%                              SRNNCellTypePairs.std_route_factor, which
%                              calls it for every model and returns []
%                              itself for a route without depression.
%
% (Plan A section 4.0 has linearize return one struct L with these four as
% fields. They are four outputs because building that struct costs about
% 1.5 microseconds, which is a third of the whole call for the product
% model, once per route in every jacobian_times. `depression` is not in
% that section's table: it came with the chain, whose factor is not
% prod(b).)
%
% F DEPENDS ON THE STATE ONLY, never on the rate: the class carries the
% rate dependence of theta = r .* F.
%
% THE MODELS, and how a route gets one (SRNNCellTypePairs.compile_route):
%   'product'  route_product: product depression times output-only
%              facilitation; also the route with no synapse dynamics. What
%              a route gets unless it asks otherwise.
%   'pool'     route_pool: vesicle pools in series. A route asks for it
%              with std.model = 'chain'.
%
% ADDING A MODEL: write route_<name>.m returning such a struct, add its
% case here and its entry to scripts/tests/test_route_models.m, and say in
% SRNNCellTypePairs.compile_route which configuration selects it. Nothing
% else in SRNNCellTypePairs changes.
%
% See also: route_product, route_pool, SRNNCellTypePairs, test_route_models

names = {'product', 'pool'};
if nargin == 0
    model = names;
    return;
end
switch char(name)
    case 'product'
        model = route_product();
    case 'pool'
        model = route_pool();
    otherwise
        error('route_model:UnknownModel', 'Unknown route model ''%s''. Valid: %s.', ...
            char(string(name)), strjoin(names, ', '));
end
end
