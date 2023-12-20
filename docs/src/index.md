# TeNLib

A Tensor Network Library (TeNLib) built on top of [ITensors.jl](https://github.com/ITensor/ITensors.jl) for quantum many-body problems.

| **Build Status** | **Documentation** |
|:----------------:|:-----------------:|
| [![Build Status](https://github.com/titaschanda/TeNLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/titaschanda/TeNLib.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![Build Status](https://github.com/titaschanda/TeNLib.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/titaschanda/TeNLib.jl/actions/workflows/documentation.yml?query=branch%3Amain) |

The source code for TeNLib can be found on [GitHub](https://github.com/titaschanda/TeNLib.jl)

The documentation for TeNLib can be found [**here**](https://titaschanda.github.io/TeNLib.jl/dev/).

## Overview

TeNLib features popular Tensor Network (TN) codes with multi-layered abstraction, that provides varyling level of control to the user. Currently, TeNLib contains codes for
* *(a)* Finite-size Matrix-Product States (MPS): Different varaints of DMRG and TDVP (including subspace expansion).
* *(b)* Tree Tensor Network (TTN): Variational search for the ground state and first few excited states.


## Installation

Currently, TeNLib.jl is not registered on Julia General Registry. To install the library (along with ITensors.jl), you can use the following steps:

```
$ julia

julia> ]

pkg> add ITensors

pkg> add https://github.com/titaschanda/TeNLib.jl
```

## Future functionality?

Here is a list for future additions in the decreasing order of priority. Any help / suggestions will be helpful.
* Augmented Tree Tensor Network (aTTN) for variational ground state search for 2D problems.
* Infinite DMRG (iDMRG) and/or Variational Uniform Matrix Product States (VUMPS) to tackle 1D / quasi-1D problems directly at the thermodynamic limit.
* Projected Entangled Pair States (PEPS) for 2D problems.
* Real-time evolution method using PEPS and TTN.

## Example: A simple DMRG code

The following code is for a simple DMRG run at the highest level of abstraction without any additional control.

```
using ITensors
using TeNLib

let
    N = 32
    sites = siteinds("S=1/2",N; conserve_qns = qn)
    os = OpStrings()
    
    for j=1:N-1
        os += 1, "Sz" => j,"Sz" => j+1
        os += 0.5, "S+" =>j, "S-" => j+1
        os += 0.5, "S-"=>j, "S+" => j+1
    end
    
    H = CouplingModel(os,sites)
    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    params = DMRGParams(;nsweeps = [5, 5], maxdim = [20, 50],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 3)

    # dmrg2 for two-site DMRG
    en, psi = dmrg2(psi0, H, params)
end
```
