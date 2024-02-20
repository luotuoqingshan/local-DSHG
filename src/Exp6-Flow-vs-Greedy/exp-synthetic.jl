using Distributed
@everywhere begin 
    include("../header.jl")
    include("gen-synthetic.jl")
    include("bulkeval-synthetic.jl")
end

n = 1000
ncluster = 30 
m2 = 50000
hep1 = 0.8
hep2 = 0.8
seedratio = 0.05
Rratio = 1.5
ntrials = 10 
genType = "RW"

# pick of epsilon is based on some non-exhausting observations
# i.e. we observe on a small set of hypergraphs and decide 
# which epsilon is best for each method
ϵs = [1.0, 0.3]
penaltyTypes = ["fracvol", "vol"]

for m1 = 5000:5000:70000

    # generate the hypergraph and save it
    cluster_vertices, cluster_edges, vlabel = gen_cluster(n, m2, ncluster)
    H_edges = hypergraph_SBM(n, ncluster, vlabel, m1, cluster_edges, hep1, hep2) 
    clusters = Vector{Vector{Int64}}()
    for i = 1:ncluster
        C = findall(x -> vlabel[x] == i, 1:n)
        push!(clusters, C)
    end

    H = elist2inc(H_edges, n)
    Ht = sparse(H')
    stats(H, Ht)
    dataset = "HSBM-$n-$ncluster-$m1-$m2-$hep1-$hep2"
    matwrite(homedir()*"/local-DHSG/data/HSBM/"*dataset*".mat",
        Dict(
            "H" => H,
            "clusters" => clusters,
        )
    )

    # generate the Rs and save
    Rs_list = Vector{Vector{Int64}}[]
    for i = 1:ncluster
        C = clusters[i]
        Csz = length(C)
        Rs = Vector{Int64}[]
        for j = 1:ntrials
            seed = sample(C, Int64(round(seedratio*Csz)), replace=false)
            R = generate_R(H, Ht, seed, Int64(round(Rratio*Csz)), genType)
            push!(Rs, R)
        end
        push!(Rs_list, Rs)
    end

    matwrite(homedir()*"/local-DHSG/data/HSBM/"*dataset*"-Rs-$seedratio-$Rratio-$genType.mat",
        Dict(
            "Rs_list" => Rs_list,
        )
    )

    for i = eachindex(ϵs) 
        ϵ = ϵs[i]
        penaltyType = penaltyTypes[i]
        for j = 1:ncluster
            res = bulkeval_ADHSG(H, Ht, ϵs[i], clusters[j], Rs_list[j], penaltyTypes[i])
            matwrite(homedir()*"/local-DHSG/results/HSBM/"*dataset*"-$penaltyType-ϵ-$ϵ-cluster-$j.mat", res)
        end
    end
end

# Plot
ϵs = [1.0, 0.3, 1.0, 0.3]
penaltyTypes = ["fracvol", "vol", "fracvol", "vol"]
penaltyLabels = ["FracVol", "Vol", "FracVol(Greedy)", "Vol(Greedy)"]
colors = [:red, :blue, :green, :orange]
type_prefix = ["flow_", "flow_", "greedy_", "greedy_"]
fig = Figure(figure_padding=1.5, resolution=(530,265))
axis = Axis(fig[1,1], 
  xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false,
  xlabelpadding=10, ylabelpadding=10, ylabel="F1 Score",
  topspinevisible=false, rightspinevisible=false) 

ratio = Vector(0.1:0.1:1.4)
savepath=homedir()*"/local-DHSG/results/HSBM/"
for i = 1:length(colors)
    result = zeros(Float64, length(ratio), ncluster * ntrials)
    penaltyType = penaltyTypes[i]
    penaltyLabel = penaltyLabels[i]
    ϵ = ϵs[i]
    for j = 1:14 
        m1 = j * 5000
        dataset = "HSBM-$n-$ncluster-$m1-$m2-$hep1-$hep2"
        for clusterid = 1:ncluster
            res = matread(savepath*dataset*"-$penaltyType-ϵ-$ϵ-cluster-$clusterid.mat")
            f1score = res[type_prefix[i]*"F1scores"]
            result[j, (clusterid-1)*ntrials+1:clusterid * ntrials] .= f1score
        end
    end
    plot_ADHSG!(axis, ratio, result, penaltyLabel, colors[i]; yscale=identity)
end
axis.xticks = 0.0:0.25:1.5
axis.yticks = 0.0:0.2:0.8
axislegend(axis, position=(0.02, 0.7), merge = true, unique = true)
save(homedir()*"/local-DHSG/figs/HSBM/flow_vs_greedy_f1score-comparison.pdf", fig, pt_per_unit=1)