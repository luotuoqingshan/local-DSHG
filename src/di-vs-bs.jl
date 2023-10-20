include("readdata.jl")
include("DHSG.jl")
include("../include/Helper_Functions.jl")

dataset = "walmart-trips"

data = load_graph(dataset, true)

H = data["H"] 
Ht = sparse(H')
@show size(H)
@show mean(H_order(Ht))

di_res = DHSG_flow_density_improvement(H, Ht; flowtol=1e-8)
@show di_res["total_dt"], di_res["niter"]
bs_res = DHSG_flow_binary_search(H, Ht; flowtol=1e-8)
@show bs_res["total_dt"], bs_res["niter"] 


dataset = "trivago-clicks"

data = load_graph(dataset, true)

H = data["H"] 
Ht = sparse(H')

@show size(H)
@show mean(H_order(Ht))

#di_res = DHSG_flow_density_improvement(H, Ht; flowtol=1e-8)
#@show di_res["total_dt"], di_res["niter"]
#bs_res = DHSG_flow_binary_search(H, Ht; flowtol=1e-8)
#@show bs_res["total_dt"], bs_res["niter"] 

dataset = "amazon-reviews"

data = load_graph(dataset, true)

H = data["H"] 
Ht = sparse(H')

@show size(H)
@show mean(H_order(Ht))

#di_res = DHSG_flow_density_improvement(H, Ht; flowtol=1e-8)
#@show di_res["total_dt"], di_res["niter"]
#bs_res = DHSG_flow_binary_search(H, Ht; flowtol=1e-8)
#@show bs_res["total_dt"], bs_res["niter"] 

dataset = "threads-math-sx"
data = load_graph(dataset)
H = data["H"]
Ht = sparse(H')
stats(H, Ht)
@show size(H)
@show mean(H_order(Ht))
di_res = DHSG_flow_density_improvement(H, Ht; flowtol=1e-8)
@show di_res["total_dt"], di_res["niter"]
bs_res = DHSG_flow_binary_search(H, Ht; flowtol=1e-8)
@show bs_res["total_dt"], bs_res["niter"] 

dataset = "threads-ask-ubuntu"
data = load_graph(dataset)
H = data["H"]
Ht = sparse(H')
stats(H, Ht)
@show size(H)
@show mean(H_order(Ht))
di_res = DHSG_flow_density_improvement(H, Ht; flowtol=1e-8)
@show di_res["total_dt"], di_res["niter"]
bs_res = DHSG_flow_binary_search(H, Ht; flowtol=1e-8)
@show bs_res["total_dt"], bs_res["niter"] 
