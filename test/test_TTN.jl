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
    psi0 = TTN(sites, 64, QN("Sz", 0))
    
    return sites, H, psi0
end


function do_ttn_optimize(sites, H, psi0)
    
    sweeppath = default_sweeppath(psi0)
    
    params = OptimizeParamsTTN(; maxdim = [128], nsweeps = [10], 
                               cutoff = 1e-14, noise = 1e-2, noisedecay = 5, 
                               disable_noise_after = 6)

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

    @assert abs(en + 13.9973156) < 1E-7
    println("SUCCESS !!")

    en1, psi1 = do_ttn_optimize_ex(sites, H, psi0, psi_gr)

    @assert abs(en1 + 13.8795758) < 1E-7
    println("SUCCESS !!")
end
    
    
