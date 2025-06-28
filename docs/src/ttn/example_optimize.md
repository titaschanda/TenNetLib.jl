# Example: Optimizing TTN

## Ground state search

A straight-forward TTN optimization can be performed as follows.

```
using ITensors
using ITensorMPS
using TenNetLib

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

    psi0 = default_randomTTN(sites, 12, QN("Sz", 0))
    sweeppath = default_sweeppath(psi0)

    params = OptimizeParamsTTN(; nsweeps = [10, 10], maxdim = [20, 50],
                               cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                               disable_noise_after = 5)

    en, psi = optimize(psi0, H, params, sweeppath)
end
```

!!! warning
    For TTN optimzaion, direct `MPO` input is not supported. One can, however, prepare a
    `CouplingModel` from `MPO` using
    [`CouplingModel(mpos::MPO...)`](@ref CouplingModel(mpos::MPO...)). 
    

Instead of using such higher-level code, one can also use lower-level functions for a better
control.
```
sysenv = StateEnvsTTN(psi0, H)

swdata = optimize!(sysenv, params, sweeppath)

# Get energy from `Sweepdata`
energy = swdata.energy[end]

# take a shallow copy of the TTN
# if the `StateEnvsTTN` will be updated later again
psi = getpsi(sysenv)

# Alternatively, take the psi from `StateEnvsTTN` itself.
# NOTE: This can crash the simulation, if the TTN is modified (e.g., in measurements)
# and `StateEnvsTTN` is going to be updated later on.
# psi = sysenv.psi
```

## Excited state search

Excited state search is also straightforword.

```
# Given a ground state `psi_gr`, initial TTN `psi0`,
# and a Hamiltonian `H`

en, psi = optimize(psi0, H, [psi_gr], params, sweeppath; weight = 10.0)
```

Similarly, using `StateEnvsTTN`:
```
sysenv_ex = StateEnvsTTN(psi0, H, [psi_gr]; weight = 10.0)
swdata_ex = optimize!(sysenv_ex, params, sweeppath)

# Get energy from `Sweepdata`
energy1 = swdata_ex.energy[end]

# take a shallow copy of the TTN
# if the `StateEnvsTTN` will be updated later again
psi1 = getpsi(sysenv_ex)

# Alternatively, take the psi from `StateEnvsTTN` itself.
# NOTE: This can crash the simulation, if the TTN is modified (e.g., in measurements)
# and `StateEnvsTTN` is going to be updated later.
# psi1 = sysenv_ex.psi
```