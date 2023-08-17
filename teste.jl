using Pkg
using LinearAlgebra
using ITensors
using Plots

print("Comecei  ")


function EGSN(Lmax)
    x = zeros(0)
    L = zeros(0)

    for count in 1:(Lmax/5)
    append!(L, 5*count)
    end

    L = [floor(Int,x) for x in L]



    let
    for N in L
        #N = 5
        Npart = N
        t1 = 1.0
        U = 2.0
        mu = U/2

        sites = siteinds("Electron", N; conserve_qns=true)

        os = OpSum()
        for b in 1:(N - 1)
        os += -t1, "Cdagup", b, "Cup", b + 1
        os += -t1, "Cdagup", b + 1, "Cup", b
        os += -t1, "Cdagdn", b, "Cdn", b + 1
        os += -t1, "Cdagdn", b + 1, "Cdn", b

        end

        for i in 1:N
        os += U, "Nupdn", i
        end
        for i in 1:N
        os += -mu, "Nup", i
        os += -mu, "Ndn", i
        end

        H = MPO(os, sites)

        nsweeps = 6
        maxdim = [50, 100, 200, 400, 800, 800]
        cutoff = [1E-12]

        
        state = ["Emp" for n in 1:N]
        p = Npart
        for i in N:-1:1
        if p > i
            #println("Doubly occupying site $i")
            state[i] = "UpDn"
            p -= 2
        elseif p > 0
            #println("Singly occupying site $i")
            state[i] = (isodd(i) ? "Up" : "Dn")
            p -= 1
        end
        end
        
        # Initialize wavefunction to be bond 
        # dimension 10 random MPS with number
        # of particles the same as `state`
        psi0 = randomMPS(sites, state, 10)

        # Check total number of particles:
        @show flux(psi0)

        # Start DMRG calculation:
        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        println("\nGround State Energy = $energy")

        append!(x, energy/N)
    end
    end

    return L, x
end

#L, x = EGSN(20)

#print(L)

#plot(L, x, seriestype=:scatter)

let
    N = 5
    t = 1
    U = 2
    mu = U/2
    cutoff = 1E-8
    tau = 0.1
    ttotal = 5.0
  
    # Make an array of 'site' indices
    s = siteinds("Electron", N; conserve_qns=true)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    for t in 0.0:tau:ttotal
      Sz = expect(psi, "Sz"; sites=c)
      println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
  
    return
  end


print("Terminei")