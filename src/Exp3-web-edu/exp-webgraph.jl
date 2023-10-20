include("../header.jl")


"""
    print_domain_names(S, id2domain)

Print the domain names of the vertices in S.
"""
function print_domain_names(S, id2domain)
    domain_names = String[]
    for v in S
        push!(domain_names, id2domain[v])
    end
    @show domain_names
end


"""
    extract_region(domain, regions)

Extract the region of the domain name.
"""
function extract_region(domain::String, regions::Vector{String})
    domain_chunks = split(domain, ".")
    if domain_chunks[1] == "edu"
        return "us"
    elseif domain_chunks[2] == "edu" || domain_chunks[2] == "ac"
        if domain_chunks[1] in regions
            return domain_chunks[1]
        elseif domain_chunks[1] == "uk"
            return "gb"
        else
            return "other"
        end
    else
        error("not an edu domain")
    end
end


"""
    get_region_counts(S, id2domain, regions)

Get the counts of the regions of the vertices in S.
"""
function get_region_counts(
    S::Vector{Int64}, 
    id2domain::Dict{Int64, String},
    regions::Vector{String}
)
    region_counts = Dict(region_name => 0 for region_name in regions)
    for v in S
        domain = id2domain[v]
        region = extract_region(domain, regions) 
        if region == "other"
            continue
        end
        region_counts[region] += 1
    end
    return region_counts
end


"""
    get_region_vertices(region, domain_ids, regions)

Get the vertices in the region.
"""
function get_region_vertices(
    region::String,
    domain_ids,
    regions::Vector{String},
)
    R = Int64[]
    for str in keys(domain_ids)
        if extract_region(str, regions) == region
            push!(R, domain_ids[str])
        end
    end
    return R
end


data_path = homedir()*"/local-DHSG/data/webgraph/"
res_path = homedir()*"/local-DHSG/results/webgraph/"
# load processed data
filepath = data_path*"webgraph-edu-ac.json"
data = JSON.parsefile(filepath) 

I = data["I"]
J = data["J"]
chost2domain = data["chost2domain"]
domain_ids = data["domain_ids"]
id2domain = Dict(value => key for (key, value) in domain_ids)

H::SparseMatrixCSC{Float64, Int64} = sparse(I, J, 1, length(chost2domain), length(domain_ids))
Ht = sparse(H')

# do some preprocessing on the hypergraph
m, n = size(H)
elist = inc2elist(Ht)  
elist = unique(elist)
H = elist2inc(elist, n)
H = deduplicate(H)
Ht = sparse(H')
H, Ht = inc_del_selfloops(H, Ht)
@show size(H), nnz(H)
stats(H, Ht)

# load some country and region info
region = CSV.File(homedir()*"/local-DHSG/src/Exp3-web-edu/all.csv") |> DataFrame;
regions = lowercase.(region[:, "alpha-2"])

R_cn = get_region_vertices("cn", domain_ids, regions)
R_uk = get_region_vertices("gb", domain_ids, regions)
matwrite(data_path*"webgraph-edu-ac-R.mat", 
        Dict("R_cn" => R_cn,
             "R_uk" => R_uk,))


penaltyType = "fracvol"

cn_answers = Vector{Int64}[]
uk_answers = Vector{Int64}[]
epsvals = Vector(0.0:0.1:1.2) 
for (i, eps) = enumerate(epsvals)
    res_cn = ADHSG_flow_density_improvement(H, Ht, R_cn, eps, penaltyType;flowtol=1e-8, type="global")
    S_cn = res_cn["optsol"]
    push!(cn_answers, S_cn)
    res_uk = ADHSG_flow_density_improvement(H, Ht, R_uk, eps, penaltyType;flowtol=1e-8, type="global")
    S_uk = res_uk["optsol"]
    push!(uk_answers, S_uk)
end

uk_largest_eps = argmax(map(x->length(x), uk_answers))
cn_largest_eps = argmax(map(x->length(x), cn_answers))

uk_largest_S = uk_answers[uk_largest_eps]
cn_largest_S = cn_answers[cn_largest_eps]

@show edensity(H, H_order(Ht), uk_largest_S)
@show edensity(H, H_order(Ht), cn_largest_S)

intersect_cn_uk = intersect(cn_largest_S, uk_largest_S)
@show edensity(H, H_order(Ht), intersect_cn_uk)

uk_region_counts = get_region_counts(uk_largest_S, id2domain, regions)
cn_region_counts = get_region_counts(cn_largest_S, id2domain, regions)
inter_region_counts = get_region_counts(intersect_cn_uk, id2domain, regions)
total_region_counts = get_region_counts(Vector(1:length(id2domain)), id2domain, regions)

matwrite(res_path*"webgraph-edu-ac.mat",
    Dict("uk_region_counts" => uk_region_counts,
         "cn_region_counts" => cn_region_counts,
         "inter_region_counts" => inter_region_counts,
         "total_region_counts" => total_region_counts,))
