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
