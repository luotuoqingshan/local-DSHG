include("../header.jl")

edges_file_list = readlines(homedir()*"/local-DHSG/src/Exp3-web-edu/cc-main-2020-21-oct-nov-jan-host-edges.paths")
vertices_file_list = readlines(homedir()*"/local-DHSG/src/Exp3-web-edu/cc-main-2020-21-oct-nov-jan-host-vertices.paths")

edges_file_list = map(x->split(x,'/')[end], edges_file_list)
vertices_file_list = map(x->split(x,'/')[end], vertices_file_list)

# as we only focus on edu domain, we need to compress the host id

# domain string => domain id
domain_ids = Dict{String, Int}()
# host id => compressed host id 
host2chost = Dict{Int, Int}()
# compressed host id => domain id
chost2domain = Dict{Int, Int}()

I = Int[]
J = Int[]

function update_domains_from_file!(
    domain_ids, 
    host2chost,
    chost2domain,
    I, 
    J,
    filename,
)
    open(GzipDecompressorStream, filename) do stream
        for line in eachline(stream)
            host_id, rev_host = split(line, "\t")
            host_id = parse(Int, host_id)
            host_id += 1
            rev_host_chunks = split(rev_host, ".")
            if rev_host_chunks[1] == "edu"
                domain = rev_host_chunks[1]*"."*rev_host_chunks[2]
            elseif rev_host_chunks[2] == "ac" || rev_host_chunks[2] == "edu"
                if length(rev_host_chunks) < 3
                    continue
                end
                domain = rev_host_chunks[1]*"."*rev_host_chunks[2]*"."*rev_host_chunks[3]
            else
                continue
            end
            host2chost[host_id] = length(host2chost) + 1
            if !haskey(domain_ids, domain)
                domain_ids[domain] = length(domain_ids) + 1 
            end
            chost_id = host2chost[host_id]
            chost2domain[chost_id] = domain_ids[domain]
            push!(I, chost_id)
            push!(J, chost2domain[chost_id])
        end
    end
end

for i = 1:12
    dt = @elapsed begin
        update_domains_from_file!(
            domain_ids, 
            host2chost,
            chost2domain,
            I,
            J,
            homedir()*"/local-DHSG/data/webgraph/"*vertices_file_list[i],
        )
    end
    @show (i, dt)
end


function build_graph!(
    I,
    J,
    host2chost,
    chost2domain,
    filename,
)  
    open(GzipDecompressorStream, filename) do stream
        for line in eachline(stream)
            u, v = split(line, "\t")
            u = parse(Int, u)
            v = parse(Int, v)
            u += 1
            v += 1
            if !haskey(host2chost, u) || !haskey(host2chost, v)
                continue
            end
            chost_u = host2chost[u]
            chost_v = host2chost[v]
            push!(I, chost_u)
            push!(J, chost2domain[chost_v])
        end
    end
end

for i = 1:24
    build_dt = @elapsed begin
        build_graph!(
            I,
            J,
            host2chost,
            chost2domain,
            homedir()*"/local-DHSG/data/webgraph/"*edges_file_list[i],
        )
    end
    @show (i, build_dt)
end

H_dt = @elapsed begin
    H = sparse(I, J, 1, length(host2chost), length(domain_ids))
    @show size(H), nnz(H)
end
@show H_dt

transform_dt = @elapsed begin
    data = Dict("I" => I, "J" => J, "host2chost"=>host2chost, "chost2domain" => chost2domain, "domain_ids" => domain_ids)
    stringdata = JSON.json(data)
end
@show transform_dt

save_dt = @elapsed begin
    open(homedir()*"/local-DHSG/data/webgraph/webgraph-edu-ac.json", "w") do f
        write(f, stringdata)
    end
end
@show save_dt


