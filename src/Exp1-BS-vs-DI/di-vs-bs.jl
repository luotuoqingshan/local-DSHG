include("../header.jl")

flowtol = 1e-8

for dataset in [
    "walmart-trips",
    "trivago-clicks",
    "threads-math-sx",
    "threads-ask-ubuntu",
    "amazon-reviews",
]
    preprocess_graph(dataset, true) 
    data = load_graph(dataset, true)
    H = data["H"]
    Ht = sparse(H')
    @show dataset
    stats(H, Ht)
    @show size(H), mean(H_order(Ht))

    di_res = DHSG_flow_density_improvement(H, Ht; flowtol=flowtol)
    @show di_res["total_dt"], di_res["niter"], di_res["optval"]
    bs_res = DHSG_flow_binary_search(H, Ht; flowtol=flowtol)
    @show bs_res["total_dt"], bs_res["niter"], bs_res["optval"] 
end
