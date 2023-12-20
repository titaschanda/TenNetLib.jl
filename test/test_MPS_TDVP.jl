using ITensors
using TeNLib


function coupling_model(;qn = true)
    
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

    return sites, H, psi0
end

function mpo(;qn = true)
    
    N = 32
    sites = siteinds("S=1/2",N; conserve_qns = qn)
    os = OpSum()
    
    for j=1:N-1
        os += 1, "Sz", j, "Sz", j+1
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
    end
    
    H = MPO(os,sites)
    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    return sites, H, psi0
end


function tdvp(sites, H, psi0)
    
    engine = TDVPEngine(psi0, H)

    for ii = 1:10
        tdvpsweep!(engine, -0.05im,
                   nsite = "dynamic";
                   maxdim = 200,
                   cutoff = 1E-12,
                   extendat = 5)
    end
end


let
    sites, H, psi0 = coupling_model(; qn = true)
    tdvp(sites, H, psi0)


    sites, H, psi0 = mpo(; qn = true)
    tdvp(sites, H, psi0)
end
