# SRNN Cell-Type Generalization

`SRNNCellTypes` extends the SRNN equations to an arbitrary ordered set of cell
types. Let (p(i)\in\{1,\ldots,C\}) identify the type of neuron (i). Each type
has its own SFA and STD counts, time constants, and connectivity statistics.

## System Equations

For a neuron (i) of type (q=p(i)),

$$
\frac{dx_i}{dt}=\frac{-x_i+u_i+\sum_{j=1}^{N}w_{ij}
\left(\prod_{m=1}^{M_{p(j)}}b_{jm}\right)r_j}{\tau_d},
$$

$$
r_i=\phi\left(x_i-c_q\sum_{k=1}^{K_q}a_{ik}\right),
\qquad
\frac{da_{ik}}{dt}=\frac{-a_{ik}+r_i}{\tau_{a,qk}},
$$

$$
\frac{db_{im}}{dt}=\frac{1-b_{im}}{\tau_{rec,qm}}
-\frac{b_{im}r_i}{\tau_{rel,q}}.
$$

Disabled mechanisms have (K_q=0) or (M_q=0) and contribute no state
variables. The packed state is ordered as all type-specific SFA blocks, all
type-specific STD blocks, and finally (x).

## Cell-Type Parameters

| MATLAB property | Shape | Meaning |
|---|---:|---|
| `n_cellTypes` | scalar | Number of cell types (C) |
| `cell_type_names` | 1-by-C cellstr | Unique names and ordering of types |
| `f` | 1-by-C | Population fractions, positive and summing to one |
| `mu_tilde` | 1-by-C | Mean potential weight by presynaptic type |
| `sigma_tilde` | 1-by-C | Weight standard deviation by presynaptic type |
| `n_a`, `n_b` | 1-by-C | SFA and STD timescale counts (K_q,M_q) |
| `tau_a` | 1-by-C cell | Per-type SFA time-constant vectors |
| `c` | 1-by-C | SFA coupling strength by type |
| `tau_b_rec` | 1-by-C cell | Per-type STD recovery vectors |
| `tau_b_rel` | 1-by-C | Shared-within-type STD release constants |

`RMTCellTypes` assigns contiguous neuron ranges using the largest-remainder
method. For a presynaptic neuron (j) of type (q), a connection is retained
with probability `alpha`, and its potential nonzero weight is

$$w_{ij}=\widetilde\mu_q+\widetilde\sigma_q z_{ij},\qquad z_{ij}\sim\mathcal N(0,1).$$

Weight sign is therefore controlled by `mu_tilde`; the implementation does not
hard-code excitatory or inhibitory semantics.
