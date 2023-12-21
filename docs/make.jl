push!(LOAD_PATH,"../src/")

using ITensors, DataStructures
using Documenter, TeNLib

############################################################################

sitename = "TeNLib.jl"

settings = Dict(
    :pages =>
        [
            "Introduction" => "index.md"
            "Base functionalities" => [
                "Dealing with Fermions" => "base/fermions.md",
                "OpStrings" => "base/opstrings.md"
                ]
            "MPS based methods" => [
                "StateEnvs" => "mps/state_envs.md",
                "Perform local updates" => "mps/update_site.md",
                "Sweeping through the MPS" => "mps/sweep.md",
                "Performing DMRG" => "mps/dmrg.md",
                "Example: DMRG" => "mps/example_dmrg.md",
                "Performing TDVP" => "mps/tdvp.md"
                ]
        ],
    :format => Documenter.HTML(; assets=["assets/favicon.ico"], prettyurls=false),
    :doctest => true,
    :checkdocs => :none,
)

                
############################################################################

makedocs(; sitename=sitename, settings...)

############################################################################

if get(ENV, "GITHUB_EVENT_NAME", nothing) == "workflow_dispatch"
  ENV["GITHUB_EVENT_NAME"] = "push"
end

deploydocs(;
  repo="github.com/titaschanda/TeNLib.jl.git",
  devbranch="main",
  push_preview=true,
  deploy_config=Documenter.GitHubActions(),
)
