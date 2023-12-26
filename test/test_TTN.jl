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
    psi0 = default_randomTTN(sites, 12, QN("Sz", 0))
    
    return sites, H, psi0
end


function do_ttn_optimize(sites, H, psi0)
    
    sweeppath = default_sweeppath(psi0)
    
    params = OptimizeParamsTTN(; maxdim = [128], nsweeps = [5], 
                               cutoff = 1e-14, noise = 1e-2, noisedecay = 5, 
                               disable_noise_after = 2)

    en, psi = optimize(psi0, H, params, sweeppath; outputlevel=1)

    return en, psi
end

function do_ttn_optimize_ex(sites, H, psi0, psi_gr)

    sweeppath = default_sweeppath(psi0)
    
    params = OptimizeParamsTTN(; maxdim = [128], nsweeps = [10], 
                               cutoff = 1e-14, noise = 1e-2, noisedecay = 5, 
                               disable_noise_after = 6)

    en, psi = optimize(psi0, H, [psi_gr], params, sweeppath; weight = 10, outputlevel=1)

    return en, psi
end


let
    sites, H, psi0 = coupling_model()
    en, psi_gr = do_ttn_optimize(sites, H, psi0)

    @show measure(psi_gr, "Sz")    
    @show measure(psi_gr, ["Sz" => 1, "Sz" => 10])
    s1 = removeqns(sites[1])
    s10 = removeqns(sites[10])
    @show measure(psi_gr, [op("Sz", s1), op("Sz", s10)])
    @show measure(psi_gr, [op("Sx", s1), op("Sx", s10)])
                    
    en1, psi1 = do_ttn_optimize_ex(sites, H, psi0, psi_gr)
end
    
    
