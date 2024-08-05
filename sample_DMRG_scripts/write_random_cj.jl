using Random
using Distributions
using DelimitedFiles
using JSON


let
    pars_json = read("input.json", String)
    parameters = JSON.parse(pars_json)

    N = parameters["N"]
    rng = MersenneTwister()


    cs = zeros(N)
    while true
        cs = rand(rng, Bernoulli(0.5), N)
        cs = 2cs .- 1
        if sum(cs)!=0 && abs(sum(cs))<N
            break
        end
    end

    open("configs.txt","w") do f
        writedlm(f, cs')
    end


end
