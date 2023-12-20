push!(LOAD_PATH,"../src/")

using Documenter, TeNLib

############################################################################

sitename = "TeNLib.jl"

settings = Dict(
    :pages =>
        [
            "Introduction" => "index.md"
            "DMRG" => "dmrg.md"
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
