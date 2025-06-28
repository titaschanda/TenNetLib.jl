push!(LOAD_PATH,"../src/")

using ITensors, ITensorMPS, DataStructures, KrylovKit
using Documenter, TenNetLib

############################################################################

sitename = "TenNetLib.jl"

settings = Dict(
    :pages =>
        [
            "Introduction" => "index.md"
            "Base functionalities" => [
                "Controlling global behaviors" => "base/global.md",
                "Dealing with Fermions" => "base/fermions.md",
                "OpStrings" => "base/opstrings.md",
                "CouplingModel" => "base/couplingmodel.md",
                "Solvers" => "base/solver.md",
                "The Graph object" => "base/graph.md",
                "Miscellaneous typedefs and functions" => "base/misc.md"
                ]
            "MPS based methods" => [
                "StateEnvs" => "mps/state_envs.md",
                "Perform local updates" => "mps/update_site.md",
                "Sweeping through the MPS" => "mps/sweep.md",
                "Performing DMRG" => "mps/dmrg.md",
                "Example: DMRG" => "mps/example_dmrg.md",
                "Performing TDVP" => "mps/tdvp.md",
                "Example: TDVP" => "mps/example_tdvp.md",
                "Measuring the MPS" => "mps/measure.md"
            ]
            "TTN based methods" => [
                "Tree Tensor Networks" => "ttn/ttn.md",
                "The TTN object" => "ttn/ttn_struct.md",
                "Generating TTN" => "ttn/ttn_gen.md",
                "StateEnvsTTN" => "ttn/state_envs_ttn.md",
                "Perform local updates" => "ttn/update_site_ttn.md",
                "Sweeping through the TTN" => "ttn/sweep_ttn.md",
                "Optimizing TTN" => "ttn/optimize_ttn.md",
                "Example: Optimizing TTN" => "ttn/example_optimize.md",
                "Measuring the TTN" => "ttn/measure_ttn.md"
            ]
        ],
    :format => Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        description = "A Tensor Network Library (TenNetLib.jl) built on top of ITensors.jl and ITensorMPS.jl for quantum many-body problems.",
        assets=["assets/favicon.ico"]),    
    :doctest => true,
    :checkdocs => :none,
)

                
############################################################################

makedocs(;sitename=sitename,
         authors = "Titas Chanda",
         settings...)

############################################################################

if get(ENV, "GITHUB_EVENT_NAME", nothing) == "workflow_dispatch"
  ENV["GITHUB_EVENT_NAME"] = "push"
end

deploydocs(;
  repo="github.com/titaschanda/TenNetLib.jl.git",
  devbranch="main",
  push_preview=true,
  deploy_config=Documenter.GitHubActions(),
)
