var documenterSearchIndex = {"docs":
[{"location":"dmrg.html#DMRG-Functions","page":"DMRG","title":"DMRG Functions","text":"","category":"section"},{"location":"dmrg.html","page":"DMRG","title":"DMRG","text":"dmrg!\ndmrg2\ndmrg1","category":"page"},{"location":"dmrg.html#TeNLib.dmrg!","page":"DMRG","title":"TeNLib.dmrg!","text":"dmrg!(sysenv::StateEnvs, params::DMRGParams, nsite::Int; kwargs...)::SweepData\n\nPerforms DMRG.\n\nArguments:\n\nsysenv::StateEnvs.\nparams::DMRGParams.\nnsite::Int of the environment. Either 1 or 2 for one-site or two-site update respectively.\n\nNamed arguments and their default values:\n\nnormalize::Bool = true: Whether to normalize after update.\nsvd_alg::String = \"divide_and_conquer\".\nweight::Float64 = -1.0: Weight for the excited state DMRG. Must be set to greater than 0 for the excited state DMRG.\noutputlevel::Int = 1. If 0 prints no information, for 1 outputs after every fullsweep, if 2 prints at every update step.\n\nConvergence criteria:\n\nenergyErrGoal: DMRG (at a particluar stage) stops when energy difference between two consecutive sweeps falls below this threshold and DMRG moves to the next stage. noise must be below Float64_threshold() to trigger this early stopping.\nentropyErrGoal: DMRG (at a particluar stage) stops when mid-chain entropy difference between two consecutive sweeps falls below this threshold and DMRG moves to the next stage. noise must be below Float64_threshold() to trigger this early stopping.\n\nWhen both energyErrGoal and entropyErrGoal are given, both conditions must be satisfied to trigger this early stopping.\n\nNamed arguments for solver and their default values:\n\nSee documentation of KrylovKit.jl.\n\nishermitian::Bool = true\nsolver_tol::Float64 = 1E-14.\nsolver_krylovdim::Int = 5.\nsolver_maxiter::Int = 2.\nsolver_outputlevel::Int = 0.: See verbosity in KrylovKit.jl.\nsolver_eager::Bool = false.\nsolver_check_convergence::Bool = false.\n\nReturn values:\n\nSweepData\n\n\n\n\n\n","category":"function"},{"location":"dmrg.html#TeNLib.dmrg2","page":"DMRG","title":"TeNLib.dmrg2","text":"dmrg2(psi0::MPS, H::MPO, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg2(psi0::MPS, H::CouplingModel, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg2(psi0::MPS, Hs::Vector{MPO}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg2(psi0::MPS, H::MPO, Ms::Vector{MPS}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg2(psi0::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg2(psi0::MPS, H::CouplingModel, Ms::Vector{MPS}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\n\nPerforms two-site DMRG.\n\nArguments:\n\npsi0::MPS: Initial MPS.\nH::MPO / H::CouplingModel / Hs:Vector{MPO} / H::MPO and Ms::Vector{MPS} / H::CouplingModel and Ms::Vector{MPS}.\nparams::DMRGParams.\nnsite::Int of the environment. Either 1 or 2 for one-site or two-site update respectively.\n\nNamed arguments and their default values:\n\nnormalize::Bool = true: Whether to normalize after update.\nsvd_alg::String = \"divide_and_conquer\".\nweight::Float64 = -1.0: Weight for the excited state DMRG. Must be set to greater than 0 for the excited state DMRG.\noutputlevel::Int = 1. If 0 prints no information, for 1 outputs after every fullsweep, if 2 prints at every update step.\n\nConvergence criteria:\n\nenergyErrGoal: DMRG (at a particluar stage) stops when energy difference between two consecutive sweeps falls below this threshold and DMRG moves to the next stage. noise must be below Float64_threshold() to trigger this early stopping.\nentropyErrGoal: DMRG (at a particluar stage) stops when mid-chain entropy difference between two consecutive sweeps falls below this threshold and DMRG moves to the next stage. noise must be below Float64_threshold() to trigger this early stopping.\n\nWhen both energyErrGoal and entropyErrGoal are given, both conditions must be satisfied to trigger this early stopping.\n\nNamed arguments for eig_solver and their default values:\n\nSee documentation of KrylovKit.jl.\n\nishermitian::Bool = true\nsolver_tol::Float64 = 1E-14.\nsolver_krylovdim::Int = 5.\nsolver_maxiter::Int = 2.\nsolver_outputlevel::Int = 0: See verbosity in KrylovKit.jl.\nsolver_eager::Bool = false.\nsolver_check_convergence::Bool = false.\n\nReturn values:\n\n::Float64: Energy.\n::MPS: The state psi.\n\n\n\n\n\n","category":"function"},{"location":"dmrg.html#TeNLib.dmrg1","page":"DMRG","title":"TeNLib.dmrg1","text":"dmrg1(psi0::MPS, H::MPO, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg1(psi0::MPS, H::CouplingModel, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg1(psi0::MPS, Hs::Vector{MPO}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg1(psi0::MPS, H::MPO, Ms::Vector{MPS}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg1(psi0::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\ndmrg1(psi0::MPS, H::CouplingModel, Ms::Vector{MPS}, params::DMRGParams; kwargs...)::Tuple{Float64, MPS}\n\nPerforms single-site DMRG. All other details are same as in dmrg2.\n\n\n\n\n\n","category":"function"},{"location":"index.html#TeNLib","page":"Introduction","title":"TeNLib","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"A Tensor Network Library (TeNLib) built on top of ITensors.jl for quantum many-body problems.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Build Status Documentation\n(Image: Build Status) (Image: Build Status)","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The source code for TeNLib can be found on GitHub","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The documentation for TeNLib can be found here.","category":"page"},{"location":"index.html#Overview","page":"Introduction","title":"Overview","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"TeNLib features popular Tensor Network (TN) codes with multi-layered abstraction, that provides varyling level of control to the user. Currently, TeNLib contains codes for","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"(a) Finite-size Matrix-Product States (MPS): Different varaints of DMRG and TDVP (including subspace expansion).\n(b) Tree Tensor Network (TTN): Variational search for the ground state and first few excited states.","category":"page"},{"location":"index.html#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Currently, TeNLib.jl is not registered on Julia General Registry. To install the library (along with ITensors.jl), you can use the following steps:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"$ julia\n\njulia> ]\n\npkg> add ITensors\n\npkg> add https://github.com/titaschanda/TeNLib.jl","category":"page"},{"location":"index.html#Future-functionality?","page":"Introduction","title":"Future functionality?","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Here is a list for future additions in the decreasing order of priority. Any help / suggestions will be helpful.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Augmented Tree Tensor Network (aTTN) for variational ground state search for 2D problems.\nInfinite DMRG (iDMRG) and/or Variational Uniform Matrix Product States (VUMPS) to tackle 1D / quasi-1D problems directly at the thermodynamic limit.\nProjected Entangled Pair States (PEPS) for 2D problems.\nReal-time evolution method using PEPS and TTN.","category":"page"},{"location":"index.html#Example:-A-simple-DMRG-code","page":"Introduction","title":"Example: A simple DMRG code","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The following code is for a simple DMRG run at the highest level of abstraction without any additional control.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"using ITensors\nusing TeNLib\n\nlet\n    N = 32\n    sites = siteinds(\"S=1/2\",N; conserve_qns = qn)\n    os = OpStrings()\n    \n    for j=1:N-1\n        os += 1, \"Sz\" => j,\"Sz\" => j+1\n        os += 0.5, \"S+\" =>j, \"S-\" => j+1\n        os += 0.5, \"S-\"=>j, \"S+\" => j+1\n    end\n    \n    H = CouplingModel(os,sites)\n    states = [isodd(n) ? \"Up\" : \"Dn\" for n in 1:N]\n    psi0 = MPS(sites, states)\n\n    params = DMRGParams(;nsweeps = [5, 5], maxdim = [20, 50],\n                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,\n                        disable_noise_after = 3)\n\n    # dmrg2 for two-site DMRG\n    en, psi = dmrg2(psi0, H, params)\nend","category":"page"}]
}
