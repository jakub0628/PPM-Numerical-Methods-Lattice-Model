# $n$-vector model simulation

## Description
$n$-vector model is the generalisation of the Ising lattice model describing the interaction of spins $\mathbf{s}_{i}$ represented by $n$-dimensional unit vectors placed on the sites $i$ of a $d$-dimensional lattice $\Lambda$. The system Hamiltonian is given by

```math
\mathcal{H} = - \sum_{\langle i,j \rangle} J_{ij} \mathbf{s}_{i} \cdot \mathbf{s}_{j} - \mu \sum_{i} \mathbf{h}_{i} \mathbf{s}_{i}
```

where
- $`i`$, $`j`$ are the sites of a $d$-dimensional lattice $\Lambda$, i.e. $`i, j \in \Lambda`$
- $`\mathbf{s}_{i}`$, $`\mathbf{s}_{j}`$ are the spins at their corresponding sites represented by $`n`$-dimensional unit vectors
- $`\langle i,j \rangle`$ denotes all pairs of neighboring sites
- $`J_{ij}`$ is the coupling (interaction) strength between the spins $`\mathbf{s}_{i}`$, $`\mathbf{s}_{j}`$
- $`\mathbf{h}_{i} = (h_i, 0, ...)`$ is the external magnetic field vector aligned with the first axis
- $`\mu`$ is the magnetic moment, describing the strength of interaction between the spins and the external magnetic field


If we ignore the external field (i.e. $`\forall_{i \in \Lambda} ~ h_i = 0`$) and assume the same coupling strength, equal to unity, for all neighbours (i.e. $`\forall_{\langle i, j \rangle \in \Lambda} ~ J_{ij} = J = 1`$) this expression simplifies to

```math
\mathcal{H} = - \sum_{\langle i,j \rangle}\mathbf{s}_{i} \cdot \mathbf{s}_{j}
```

Hence, the energy of each site $i$ with $2d$ neighbours ${\langle j \rangle}$ can be expressed as

```math
E_i = - \sum_{\langle j \rangle}\mathbf{s}_{i} \cdot \mathbf{s}_{j}
```

## Metropolis–Hastings algorithm

1. For a randomly site $`i`$ calculate the energy difference after the spin $`\mathbf{s}_{i}`$ is flipped ($`\mathbf{s}_{i} \to -\mathbf{s}_{i}`$):
 $`\Delta E_i = E_i(-\mathbf{s}_{i}) - E_i(\mathbf{s}_{i})`$.
2. If $`\Delta E_i < 0`$ (i.e. the flip is energetically favourable) flip the spin.
3. Otherwise, even if the flip is energetically unfavourable, it can still happen because of a random thermal fluctuation with the probability $`p_{\text{B}} = \exp(-\beta E_i).`$ Hence, if $`p_{\text{B}} > p`$, where $`p \in [0, 1]`$ is a random number chosen from an uniform distribution, also flip the spin.
4. Repeat the procedure.

## Resources
- [Wikipedia article on the $n$-vector model](https://en.wikipedia.org/wiki/N-vector_model)
- [XY (2-vector) model simulation](https://kjslag.github.io/XY/)
- Stanley, H. E. (1968). "Dependence of Critical Properties upon Dimensionality of Spins". *Phys. Rev. Lett*. **20** (12): 589–592. [doi:10.1103/PhysRevLett.20.589](https://doi.org/10.1103%2FPhysRevLett.20.589)
## Libraries
- https://github.com/genkuroki/Ising2D.jl
- https://juliapackages.com/p/isingmodels

# TODO

- [x] Establish general aims of the project 
- [x] Find out how to implement the model
- [ ] Divide the tasks
- [ ] Implement the basic version of the algorithm
- [ ] drink coffee
