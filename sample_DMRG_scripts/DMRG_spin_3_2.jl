using Random
using Distributions
using DelimitedFiles
using JSON
using ITensors

function ITensors.space(::SiteType"S=3/2";
    conserve_qns=false)
if conserve_qns
return [QN("Sz",3)=>1,QN("Sz",1)=>1,
QN("Sz",-1)=>1,QN("Sz",-3)=>1]
end
return 4
end

ITensors.op(::OpName"Sz",::SiteType"S=3/2") =
[+3/2   0    0    0
0  +1/2   0    0 
0    0  -1/2   0
0    0    0  -3/2]

ITensors.op(::OpName"S+",::SiteType"S=3/2") =
[0  √3  0  0
0   0  2  0
0   0  0 √3
0   0  0  0] 

ITensors.op(::OpName"S-",::SiteType"S=3/2") =
[0   0  0   0
√3  0  0   0
0   2  0   0
0   0  √3  0] 

ITensors.state(::StateName"Up", ::SiteType"S=3/2") = [1.0, 0, 0, 0]
ITensors.state(::StateName"HalfUp", ::SiteType"S=3/2") = [0.0, 1.0, 0, 0]
ITensors.state(::StateName"HalfDn", ::SiteType"S=3/2") = [0.0, 0.0, 1.0, 0]
ITensors.state(::StateName"Dn", ::SiteType"S=3/2") = [0.0, 0, 0, 1]

function H_MPO(sites, cs)
    N = length(sites)
    os = OpSum()
    for j=1:N-1
        s = -cs[j] * cs[j+1]
        os += s  , "Sz",j,"Sz",j+1
        os += s/2, "S+",j,"S-",j+1
        os += s/2, "S-",j,"S+",j+1
    end
    H = MPO(os,sites)
end

function entropies!(psi)
    N = length(psi)
    Ss = zeros(N)
    for b in 2:N
        orthogonalize!(psi, b)
        U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
        SvN = 0.0
        for n=1:ITensors.dim(S, 1)
            p = S[n,n]^2
            SvN -= p * log(p)
        end
        Ss[b] = SvN
        #@show b, SvN
    end
    return Ss
end

mutable struct SizeObserver <: AbstractObserver
    energy_tol::Float64
    last_energy::Float64
    min_sweeps::Int64
    SizeObserver(energy_tol=0.0, min_sweeps=0) = new(energy_tol,1000.0,min_sweeps)
 end

 function ITensors.checkdone!(o::SizeObserver;kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    if sw > o.min_sweeps && abs(energy-o.last_energy)/abs(energy) < o.energy_tol
      println("Stopping DMRG after sweep $sw")
      return true
    end
    # Otherwise, update last_energy and keep going
    o.last_energy = energy
    return false
  end

function ITensors.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
    GC.gc()
    if bond==1 && half_sweep==2
        psi_size =  Base.format_bytes(Base.summarysize(psi))
        PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
        println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
    end
end

let
    pars_json = read("input.json", String)
    parameters = JSON.parse(pars_json)


    chis = [parameters["chi1"], parameters["chi2"], parameters["chi3"]]
    nsweeps = 120
    obs = SizeObserver(1e-11, 30)
    #obs = DMRGObserver(energy_tol=1e-11, minsweeps=10)
    maxdims = Vector{Vector{Int64}}(undef,3)
    maxdims[1] = [10,10,20,20,50,50,100,100,100,200,200,200, chis[1]]
    maxdims[2] = [chis[2]]
    maxdims[3] = [chis[3]]
    cutoff = parameters["cutoff"]

    cs = vec(readdlm("configs.txt"))
    N = length(cs)

    sites = siteinds("S=3/2",N; conserve_qns=true)

    H = H_MPO(sites, cs)

    state = [cs[n]==1 ? "Up" : "Dn" for n=1:N]
    psi = MPS(sites,state)

    for chi_ind in 1:3

        E0, psi = dmrg(H,psi; nsweeps=nsweeps, maxdim=maxdims[chi_ind], cutoff=cutoff, noise=1e-9,
                    observer=deepcopy(obs))
        E0, psi = dmrg(H,psi; nsweeps=10, maxdim=chis[chi_ind], cutoff=cutoff, noise=0,
		    observer=deepcopy(obs))

        var_E0 = inner(H,psi,H,psi) - E0^2

        err_E0 = (var_E0>0 ? sqrt(var_E0) : 0)
        
        ####################################
        ####   measure GS properties    ####
        ####################################

        Szs = expect(psi,"Sz")
        open("Z_exp_chi_$(chis[chi_ind]).txt","a") do f
            writedlm(f, (cs.*Szs)')
        end

        #SvN = entropies!(psi)
        #open("entropies_chi_$(chis[chi_ind]).txt","a") do f
        #    writedlm(f, SvN')
        #end

        zz_corr = (cs * cs') .* correlation_matrix(psi,"Sz","Sz")
        xy_corr_no_c = 0.5 * real.(correlation_matrix(psi, "S+", "S-"))
        xy_corr_no_c += 0.5 * real.(correlation_matrix(psi, "S-", "S+"))
        xy_corr = (cs * cs') .* xy_corr_no_c


        open("ZZ_corr_chi_$(chis[chi_ind]).txt","w") do f
            writedlm(f, zz_corr)
        end
        open("XY_corr_chi_$(chis[chi_ind]).txt","w") do f
            writedlm(f, xy_corr)
        end
        open("XY_corr_no_c_chi_$(chis[chi_ind]).txt","w") do f
            writedlm(f, xy_corr_no_c)
        end
        open("energies_chi_$(chis[chi_ind]).txt","w") do f
            write(f, "$E0\t$err_E0\t$(maxlinkdim(psi))\n")
        end
    end

    state = [cs[n]==1 ? "Up" : "Dn" for n=1:N]
    #flip one spin to increase the spin
    I = findfirst(x->x==-sign(sum(cs)),cs)
    state[I] = (cs[I]==-1 ? "HalfDn" : "HalfUp")
    
    psi = MPS(sites,state)
    for chi_ind in 1:3

        ####################################
        #########   compute gap    #########
        ####################################

        E1, psi = dmrg(H,psi; nsweeps=nsweeps, maxdim=maxdims[chi_ind], cutoff=cutoff, noise=1e-9,
                    observer=deepcopy(obs))
        E1, psi = dmrg(H,psi; nsweeps=10, maxdim=[chis[chi_ind]], cutoff=cutoff, noise=0,
		    observer=deepcopy(obs))

        var_E1 = inner(H,psi,H,psi) - E1^2
        err_E1 = (var_E1>0 ? sqrt(var_E1) : 0)

        open("energies_chi_$(chis[chi_ind]).txt","a") do f
            write(f, "$E1\t$err_E1\t$(maxlinkdim(psi))\n")
        end

    end

end