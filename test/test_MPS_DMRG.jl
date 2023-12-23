
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


function dmrg_1(sites, H, psi0;
                nsite = 1)
   
    sysenv = StateEnvs(psi0, H)
    params = DMRGParams(;nsweeps = [5], maxdim = [20],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 2)      
    sw = dmrg!(sysenv, params, nsite)
    return sw.energy[end], sysenv.psi
end


function dmrg_2(sites, H, psi0;
                nsite = 1)
    
    params = DMRGParams(;nsweeps = [5], maxdim = [20],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 2)      

    en, psi = TeNLib.dmrg(psi0, H, params, nsite)

    return en, psi 
end


function dmrg_ex_1(sites, H, psi0, psi_gr;
                                  nsite = 1)

    sysenv = StateEnvs(psi0, H, [psi_gr]; weight = 10.0)
    
    params = DMRGParams(;nsweeps = [5], maxdim = [20],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 2)
    
    sw = dmrg!(sysenv, params, nsite)
    return sw.energy[end], sysenv.psi
end




function dmrg_ex_2(sites, H, psi0, psi_gr;
                                  nsite = 1)
    
    params = DMRGParams(;nsweeps = [5], maxdim = [20],
                        cutoff = 1e-14, noise = 1e-3, noisedecay = 2,
                        disable_noise_after = 2)      
    
    en, psi = TeNLib.dmrg(psi0, H, [psi_gr], params, nsite; weight = 10.0)
    
    return en, psi
end


let
    sites, H, psi0 = coupling_model(; qn = true)

    for nsite in [1, 2]
        @time en, psi_gr = dmrg_1(sites, H, psi0;
                                  nsite = nsite)
      
        @time en1, psi1 = dmrg_ex_1(sites, H, psi0, psi_gr;
                                    nsite = nsite)
    end

    for nsite in [1, 2]
        @time en, psi_gr = dmrg_2(sites, H, psi0;
                                  nsite = nsite)
        
        @time en1, psi1 = dmrg_ex_2(sites, H, psi0, psi_gr;
                              nsite = nsite)
    end

    
    sites, H, psi0 = mpo(; qn = true)

    for nsite in [1, 2]
        @time en, psi_gr = dmrg_1(sites, H, psi0;
                                  nsite = nsite)

        @time en1, psi1 = dmrg_ex_1(sites, H, psi0, psi_gr;
                                    nsite = nsite)
    end

    for nsite in [1, 2]
        @time en, psi_gr = dmrg_2(sites, H, psi0;
                                  nsite = nsite)

        @show measure(psi_gr, "Sz")
        @show measure(psi_gr, ["Sz" => 1, "Sz" => 10])
        
        @time en1, psi1 = dmrg_ex_2(sites, H, psi0, psi_gr;
                                    nsite = nsite)

    end

end
