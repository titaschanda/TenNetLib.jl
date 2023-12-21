# Example: DMRG

## Ground state search

A straight-forward DMRG can be performed as follows.

```
using ITensors
using TeNLib

let
    N = 32
    sites = siteinds("S=1/2",N)
    os = OpStrings()
    
    for j=1:N-1
        os += 1, "Sz" => j, "Sz" => j+1
        os += 0.5, "S+" => j, "S-" => j+1
        os += 0.5, "S-" => j, "S+" => j+1
    end
    
    H = CouplingModel(os,sites)
    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    params = DMRGParams(;nsweeps = [10, 10], maxdim = [20, 50],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 5)

    # dmrg2 for two-site DMRG
    en, psi = dmrg2(psi0, H, params)

    # or dmrg1 for single-site DMRG
    # en, psi = dmrg1(psi0, H, params)
end
```

Here, we have used `OpStrings` and `CouplingModel`. Alternatively, standard `OpSum` and `MPO` can
be used.

Instead of using such higher-level code, one can also use lower-level functions for a better
control.
```
sysenv = StateEnvs(psi0, H)
nsite = 2 # two-site update

swdata = dmrg!(sysenv, params, nsite)

# Get energy from `Sweepdata`
energy = swdata.energy[end]

# take a shallow copy of the MPS
# if the `StateEnvs` will be updated later again
psi = getpsi(sysenv)

# Alternatively, take the psi from `StateEnvs` itself.
# NOTE: This can crash the simulation, if the MPS is modified (e.g., in measurements)
# and `StateEnvs` is going to be updated later on.
# psi = sysenv.psi
```
Most often, it is better to do a single-site DMRG (without any noise) after standard two-site update for better
convergence. Such lowe-level function using `StateEnvs` is useful for that.
```
sysenv = StateEnvs(psi0, H)

params2 = DMRGParams(;nsweeps = [10, 10], maxdim = [20, 50],
                     cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                     disable_noise_after = 5)
dmrg!(sysenv, params2, 2)

params1 = DMRGParams(;nsweeps = [10], maxdim = [50], cutoff = 0.0)
dmrg!(sysenv, params1, 1)
```

Using such lower-level function, one can also restart the simulation at later times
(e.g., by saving the `StateEnvs` using `Serialization.jl`), or one can perform
[`updateH!`](@ref updateH!(sysenv::StateEnvs{ProjMPO}, H::MPO; recalcEnv::Bool = true)) to slowly
change the Hamiltonian during DMRG simulations. Most often it is better to perform a 


Global Subspace Expansion can also be used to get rid of nasty local minimas (if needed), if the `StateEnvs` is
built from a single MPO.
```
krylov_extend!(sysenv; extension_applyH_maxdim = 40)
```
**Note**: Global Subspace Expansion can result into huge MPS bond dimension. That is why
the named input parameters of [`krylov_extend!`](@ref krylov_extend!(sysenv::StateEnvs{ProjMPO}; kwargs...)) should be chosen carefully.

## Excited state DMRG

Excited state DMRG is also straightforword.

```
# Given a ground state `psi_gr`, initial MPS `psi0`,
# and a Hamiltonian `H`

# dmrg2 for two-site DMRG
en, psi = dmrg2(psi0, H, [psi_gr], params; weight = 10.0)

# or dmrg1 for single-site DMRG
# en, psi = dmrg1(psi0, H, [psi_gr], params; weight = 10.0)
```

Similarly, using `StateEnvs`:
```
sysenv_ex = StateEnvs(psi0, H, [psi_gr]; weight = 10.0)
nsite = 2 # two-site update

swdata_ex = dmrg!(sysenv_ex, params, nsite)

# Get energy from `Sweepdata`
energy1 = swdata_ex.energy[end]

# take a shallow copy of the MPS
# if the `StateEnvs` will be updated later again
psi1 = getpsi(sysenv_ex)

# Alternatively, take the psi from `StateEnvs` itself.
# NOTE: This can crash the simulation, if the MPS is modified (e.g., in measurements)
# and `StateEnvs` is going to be updated later.
# psi1 = sysenv_ex.psi
```