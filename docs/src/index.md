# TeNLib

A Tensor Network Library (TeNLib) built on top of [ITensors.jl](https://github.com/ITensor/ITensors.jl) for quantum many-body problems.

| **Build Status** | **Documentation** |
|:----------------:|:-----------------:|
| [![Build Status](https://github.com/titaschanda/TeNLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/titaschanda/TeNLib.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![Build Status](https://github.com/titaschanda/TeNLib.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/titaschanda/TeNLib.jl/actions/workflows/documentation.yml?query=branch%3Amain) |

The source code for TeNLib can be found on [GitHub](https://github.com/titaschanda/TeNLib.jl)

The documentation for TeNLib can be found [**here**](https://titaschanda.github.io/TeNLib.jl/dev/).

## Overview

Currently, TeNLib contains codes for
* *(a)* Finite-size Matrix-Product States (MPS): Different varaints of DMRG and TDVP (including subspace expansion).
* *(b)* Tree Tensor Network (TTN): Variational search for the ground state and first few excited states.