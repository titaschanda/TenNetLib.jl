using TeNLib
using Test

@testset "MPS" begin
#    @time include("test_MPS_DMRG.jl")
#    @time include("test_MPS_TDVP.jl")
end


@testset "TTN" begin
    @time include("test_TTN.jl")
end
