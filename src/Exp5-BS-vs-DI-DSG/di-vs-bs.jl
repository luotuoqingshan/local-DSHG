include("../header.jl")

flowtol = 1e-8

for dataset in [
    "ca-AstroPh",
    "ca-HepPh",
    "email-Enron",
    "com-amazon",
    "com-youtube",
]
    preprocess_graph(dataset) 
    data = load_graph(dataset)
    A = data["A"]
    @show dataset
    @show size(A, 1), div(nnz(A), 2)

    di_res = DSG_flow_density_improvement(A; flowtol=flowtol)
    @show di_res["total_dt"], di_res["niter"], di_res["optval"]
    bs_res = DSG_flow_binary_search(A; flowtol=flowtol)
    @show bs_res["total_dt"], bs_res["niter"], bs_res["optval"] 
end
