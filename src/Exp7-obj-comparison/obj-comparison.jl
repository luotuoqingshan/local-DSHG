using Distributed
@everywhere begin 
    include("../header.jl")
    include("bulkeval.jl")
end

for dataset in [
    "walmart-trips",
    "trivago-clicks",
    "threads-ask-ubuntu",
    "threads-math-sx",
    ]
    # preprocess and load data
    preprocess_graph(dataset, true)
    data = load_graph(dataset, true)
    H = data["H"]
    Ht = sparse(H')

    # show some basic stats
    @show dataset, size(H)
    stats(H, Ht)

    # generate 100 Rs with |R| = 200 expanded from 10 seeds
    ntrial = 100 
    ϵs = Vector(0.0:0.1:1.2)  
    seedsize = 10 
    Rsize = 200 

    # generate R using RW-geo with stopping probability 0.3
    prob = 0.3
    genType = "RW-geo"
    Rdataname = "Rs-$seedsize-$Rsize-ntrial-$ntrial-$genType"
    datapath = homedir()*"/local-DHSG/data/"*dataset*"/"*Rdataname*".mat"

    # generate Rs
    #Rs = bulk_generate_R(H, Ht, seedsize, Rsize, ntrial, genType; dataset=dataset, prob=prob)

    # evaluate running time
    penaltyTypes = ["fracvol", "vol", "UCE", "WCE"]
    Volpenalizeddensities_results = []
    FractionalVolpenalizeddensities_results = []
    for penaltyType in penaltyTypes
        #res = bulkeval_ADHSG(H, Ht, Rs, ϵs, penaltyType; dataset=dataset, filename="$penaltyType-Rs-$seedsize-$Rsize-ntrial-$ntrial-$genType")
        res = matread(homedir()*"/local-DHSG/results/"*dataset*"/$penaltyType-Rs-$seedsize-$Rsize-ntrial-$ntrial-$genType.mat")
        push!(Volpenalizeddensities_results, res["Volpenalizeddensities"])
        push!(FractionalVolpenalizeddensities_results, res["FractionalVolpenalizeddensities"])
    end

    for i = [2, 3, 4, 1] 
        @. FractionalVolpenalizeddensities_results[i] /= FractionalVolpenalizeddensities_results[1]
    end

    for i = [1, 3, 4, 2] 
        @. Volpenalizeddensities_results[i] /= Volpenalizeddensities_results[2]
    end
    # plot
    penaltyNames = ["FracVol", "Vol", "UCE", "WCE"]
    cols = [:red, :blue, :green, :purple]
    savepath = homedir()*"/local-DHSG/figs/"*dataset*"/"*dataset*"-fractionalvolpenalized-density-comparison.pdf"
    plot_objective_ratio(ϵs, FractionalVolpenalizeddensities_results, penaltyNames, cols, savepath;)
    savepath = homedir()*"/local-DHSG/figs/"*dataset*"/"*dataset*"-volpenalized-density-comparison.pdf"
    plot_objective_ratio(ϵs, Volpenalizeddensities_results, penaltyNames, cols, savepath;)
end

