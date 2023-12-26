# TeNLib.jl

A Tensor Network Library (TeNLib.jl) built on top of [ITensors.jl](https://github.com/ITensor/ITensors.jl) for quantum many-body problems.

| **Build Status** | **Documentation** |
|:----------------:|:-----------------:|
| [![Build Status](https://github.com/titaschanda/TeNLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/titaschanda/TeNLib.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![Build Status](https://github.com/titaschanda/TeNLib.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/titaschanda/TeNLib.jl/actions/workflows/documentation.yml?query=branch%3Amain) |

The source code for TeNLib.jl can be found on [GitHub](https://github.com/titaschanda/TeNLib.jl)

The documentation for TeNLib.jl can be found [**here**](https://titaschanda.github.io/TeNLib.jl/dev/).

This library requires Julia 1.7+.

## Overview

TeNLib.jl features widely-used Tensor Network (TN) codes, designed with a **multi-layered abstraction**
to cater to diverse user needs. The library provides users with varying levels of control over their computations.
Currently, TeNLib.jl presents an array of functionalities for:
* *(a)* **Finite-size Matrix-Product States (MPS)**: Different variants of Density Matrix Renormalization Group (DMRG) and Time Dependent Variational Principle (TDVP) (including subspace expansion) methods.
* *(b)* **Tree Tensor Network (TTN)**: Variational search for the ground state and first few excited states.


## Installation

Currently, TeNLib.jl is not registered on Julia General Registry. To install the library (along with ITensors.jl), you can use the following steps:

```
$ julia

julia> ]

pkg> add ITensors

pkg> add https://github.com/titaschanda/TeNLib.jl
```

## Found an issue or bug?

> "Beware of bugs in the above code; I have only proved it correct, not tried it."
>    -- Donald Knuth


If you find bugs or a mistakes of any kind, please let us know by adding an issue to the
[GitHub issue tracker](https://github.com/titaschanda/TeNLib.jl/issues).
You are also welcome to submit a [pull request](https://github.com/titaschanda/TeNLib.jl/pulls).


## Future functionality?

Here is a list for future additions in the decreasing order of priority. Any help / suggestion is welcome.
* **Augmented Tree Tensor Network (aTTN)** for variational ground state search for 2D problems.
* **Infinite DMRG (iDMRG)** and/or **Variational Uniform Matrix Product States (VUMPS)** to tackle 1D / quasi-1D problems directly at the thermodynamic limit.
* **Projected Entangled Pair States (PEPS)** for 2D problems.
* Real-time evolution method using PEPS and TTN.

Also, please feel free to ask about a new feature by adding a new request to the
[GitHub issue tracker](https://github.com/titaschanda/TeNLib.jl/issues) labelled
`feature request`. Note that submitting a pull request, providing the needed changes to
introduced your requested feature, will speed up the process.



## Example: A simple DMRG code

The following code is for a simple DMRG run at **the highest level of abstraction without any additional control**.

```
using ITensors
using TeNLib

let
    N = 32
    sites = siteinds("S=1/2",N)
    os = OpSum()
    
    for j=1:N-1
        os += 1, "Sz", j,"Sz", j+1
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
    end
    
    H = MPO(os,sites)
    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    params = DMRGParams(;nsweeps = [5, 5], maxdim = [20, 50],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 3)

    # dmrg2 for two-site DMRG
    en, psi = dmrg2(psi0, H, params)
end
```

## Example: A simple TDVP code

The following code is for a simple TDVP run at **the highest level of abstraction without any additional control**.

```
using ITensors
using TeNLib

let
    N = 32
    sites = siteinds("S=1/2",N)
    os = OpSum()
    
    for j=1:N-1
        os += 1, "Sz", j,"Sz", j+1
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
    end
    
    H = MPO(os,sites)
    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    tau = -0.01im    
    engine = TDVPEngine(psi0, H)
    for ii = 1:100

    	# `nsite = "dynamic"` for dynamical selection between
	# single- and two-site variants at different bonds
        tdvpsweep!(engine, tau,
                   nsite = "dynamic";
                   maxdim = 200,
                   cutoff = 1E-12,
                   extendat = 5)

	psi = getpsi(engine)
	# DO STUFF
    end
end
```

## Example: A simple TTN ground-state optimzation code

The following code is for a simple TTN ground-state optimzation run at **the highest level of abstraction without any additional control**.
Here we use `OpStrings` and `CouplingModel` instead of `OpSum` and `MPO`.

```
using ITensors
using TeNLib

let
    N = 32
    sites = siteinds("S=1/2",N)
    os = OpStrings()
    
    for j=1:N-1
        os += 1, "Sz" => j,"Sz" => j+1
        os += 0.5, "S+" => j, "S-" => j+1
        os += 0.5, "S-"=> j, "S+" => j+1
    end
    
    H = CouplingModel(os,sites)
    psi0 = TTN(sites, 64, QN("Sz", 0))

    sweeppath = default_sweeppath(psi0)
    
    params = OptimizeParamsTTN(; maxdim = [64, 128], nsweeps = [5, 10], 
                               cutoff = 1e-14, noise = 1e-2, noisedecay = 5, 
                               disable_noise_after = 5)
			       
    en, psi = optimize(psi0, H, params, sweeppath)
end
```    