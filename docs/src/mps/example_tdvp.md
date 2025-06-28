# Example: TDVP

## Time-evolution with time dependent variational principle (TDVP).

A straight-forward TDVP (with dynamical sweeps) can be performed as follows.
Note that TDVP computes `ψ' = exp(time_step * H) * ψ`. Therefore, for real-time dynamics
with step `dt`, `time_step` should be `-im * dt`.

```
using ITensors
using ITensorMPS
using TenNetLib

let
    N = 32
    sites = siteinds("S=1/2",N)
    os = OpSum()
    
    for j=1:N-1
        os += 1, "Sz", j, "Sz", j+1
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
    end
    
    H = MPO(os,sites)
    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    tau = -0.01im    
    engine = TDVPEngine(psi0, H)
    for ii = 1:100

    	# `nsite = "dynamic"` for dynamical selection between
	# single- and two-site variants at different bonds
        tdvpsweep!(engine, tau,
                   nsite = "dynamic";
                   maxdim = 200,
                   cutoff = 1E-12,
                   extendat = 5)

	# Errors in the last sweep and the total error till this point
	swerr = sweeperror(engine)
	totalerr = totalerror(engine)
	
	psi = getpsi(engine)
	
	# DO STUFF
    end
end
```

Here, we have used `OpSum` and `MPO`. Alternatively, standard `OpStrings` and `CouplingModel` can
be used.

Optionally, one can use pure single- or two-site updates. If the `TDVPEngine` is created with a
single `MPO` then Global Subspace Expansion can be performed.
```
for ii = 1:100

    # GSE at every 5th sweep.
    if ii % 5 == 1
        krylov_extend!(engine)
    end
    
    # `nsite = 1` for single-site update
    tdvpsweep!(engine, tau,
               nsite = 1;
               maxdim = 200,
               cutoff = 1E-12)

    # Errors in the last sweep and the total error till this point
    swerr = sweeperror(engine)
    totalerr = totalerror(engine)
	
    psi = getpsi(engine)
	
    # DO STUFF
end
```

## Time-evolution with time-dependent Hamiltonian

For time-dependent Hamiltonian, [`updateH!`](@ref updateH!(engine::TDVPEngine, H::T; recalcEnv::Bool = true) where T <: Union{MPO, Vector{MPO}, CouplingModel}) can be used.

```
using ITensors
using TenNetLib

function makeHt(sites, t)

    os = OpSum()    
    for j=1:N-1
        os += 1, "Sz", j, "Sz", j+1
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
    end
    for j=1:N
        os += t, "Sz", j
    end
    return MPO(os,sites)
end

let
    N = 32
    sites = siteinds("S=1/2",N)

    states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = MPS(sites, states)

    dt = 0.01
    tau = -im * dt

    H0 = makeHt(sites, 0.0)
    engine = TDVPEngine(psi0, H0)
    
    for ii = 1:100

    	# `nsite = "dynamic"` for dynamical selection between
	# single- and two-site variants at different bonds
        tdvpsweep!(engine, tau,
                   nsite = "dynamic";
                   maxdim = 200,
                   cutoff = 1E-12,
                   extendat = 5)

	# Errors in the last sweep and the total error till this point
	swerr = sweeperror(engine)
	totalerr = totalerror(engine)
	
	psi = getpsi(engine)
	
	# DO STUFF

	# update Hamiltonian for the next iteration
	Ht = makeHt(sites, ii*dt)
	updateH!(engine, Ht)
    end
end
```

**Note 1**: The above example for time-dependent Hamiltonian is very crude. In real situations,
the update in the Hamiltonian should be done with proper care.

**Note 2**: When  the `TDVPEngine` is created from a single `MPO` and `dt` is small,
one can use `recalcEnv = false` in [`updateH!`](@ref updateH!(engine::TDVPEngine, H::T; recalcEnv::Bool = true) where T <: Union{MPO, Vector{MPO}, CouplingModel}), so that environments from the
last step is reused.