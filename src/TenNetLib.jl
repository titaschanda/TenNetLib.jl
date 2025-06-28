
module TenNetLib

##*********** MAIN Imports ************************

using Printf
using ITensors
using ITensorMPS
using LinearAlgebra

import DataStructures:
    Deque


import ITensors:
    QNBlocks,    
    scalartype,
    OneITensor


import KrylovKit:
    eigsolve,
    exponentiate


import ITensorMPS:
    AbstractProjMPO,
    orthocenter,
    set_nsite!,
    position!,
    setleftlim!,
    setrightlim!,
    nsite,
    site_range,
    truncate,
    rproj,
    lproj



##*********** FILES *****************************

include("base/global_variables.jl")
include("base/typedef.jl")
include("base/helper_internal_funcs.jl")
include("base/fermions.jl")
include("base/solver.jl")
include("base/graph.jl")
include("base/opstrings.jl")
include("base/opstringsmpo.jl")
include("base/couplingmodel.jl")


include("mps/projmps2.jl")
include("mps/projmpo_mps2.jl")
include("mps/projmposum2.jl")
include("mps/projmposum_mps.jl")
include("mps/projcouplingmodel.jl")
include("mps/projcouplingmodel_mps.jl")
include("mps/state_envs.jl")
include("mps/update_site.jl")
include("mps/sweep.jl")
include("mps/dmrg.jl")
include("mps/tdvp.jl")
include("mps/measure.jl")

include("ttn/typedef.jl")
include("ttn/ttn.jl")
include("ttn/helper_internal_funcs.jl")
include("ttn/ttn_generators.jl")
include("ttn/linktensors.jl")
include("ttn/linkproj.jl")
include("ttn/measure_ttn.jl")
include("ttn/environment.jl")
include("ttn/state_envs_ttn.jl")
include("ttn/update_site_ttn.jl")
include("ttn/sweep_ttn.jl")
include("ttn/optimize_ttn.jl")


##**************************************************

include("export.jl")

ITensors.enable_contraction_sequence_optimization()

##**************************************************

end


    


