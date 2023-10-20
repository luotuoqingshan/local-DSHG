"""
    deduplicate(H)

Remove duplicate nodes in one hyperedge.
"""
function deduplicate(
    H::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer}
    # deal with the corner case that one vertex 
    # appears in one hyperedge multiple times 
    I, J, _ = findnz(H)
    return sparse(I, J, ones(Tf, length(I)), size(H)...)
end


"""
    adj_del_selfloops(A)

Delete all self-loops in an adjacency matrix A.  
"""
function adj_del_selfloops(
    A::SparseMatrixCSC{Tf, Ti}
) where {Ti <: Integer, Tf}
    A[diagind(A)] .= zero(Tf) 
    return dropzeros(A)
end


"""
    inc_del_selfloops(H, Ht)

Delete all self-loops in an incidence matrix H. 
"""
function inc_del_selfloops(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Ti <: Integer, Tf} 
    m, _ = size(H)
    order = H_order(Ht)
    eids = findall(x->order[x]>1, 1:m)
    return H[eids, :], Ht[:, eids] 
end


"""
    adj2elist(A, [undirected=true])

Convert an adjacency matrix of an (un)directed graph to a list of edges.
"""
function adj2elist(
    A::SparseMatrixCSC{Tf, Ti},
    undirected=true,
) where {Tf, Ti <: Integer}
    n = size(A, 1)
    rowval = A.rowval    
    colptr = A.colptr
    nzval = A.nzval
    Hyperedges = Vector{Vector{Ti}}()
    for u = 1:n
        for nzi in colptr[u]:colptr[u+1]-1
            v = rowval[nzi]
            edge = Vector{Ti}()
            if u <= v  
                push!(edge, u)
                push!(edge, v)
                push!(Hyperedges, edge)
            end
            if (u > v && !undirected)
                push!(edge, u)
                push!(edge, v)
                push!(Hyperedges, edge)
            end
        end
    end
    return Hyperedges
end


"""
    elist2inc(Hyperedges, n)

Convert a list of hyperedges into an incidence matrix of a hypergraph with n nodes.
"""
function elist2inc(
    Hyperedges::Vector{Vector{Ti}},
    n::Ti,
) where {Ti <: Integer}
    eid = 0
    I = Ti[]
    J = Ti[]
    V = Float64[] 
    for i = eachindex(Hyperedges)
        eid += 1
        for v in Hyperedges[i]
            push!(I, eid)
            push!(J, v)
            push!(V, 1.0)
        end
    end
    return sparse(I, J, V, length(Hyperedges), n)
end


"""
    inc2adj(H, Ht)

Turn an incidence matrix into an adjacency matrix.
Only works for normal graphs without self-loops.
"""
function inc2adj(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer} 
    m, n = size(H)
    order = H_order(Ht) 
    # some sanity check
    # make sure this is a normal graph
    @assert maximum(order) == 2
    @assert minimum(order) == 2
    colptr = Ht.colptr
    rowval = Ht.rowval
    I = Ti[]; J = Ti[]
    for e = 1:m
        u = colptr[e]; v = colptr[e]+1
        u = rowval[u]; v = rowval[v]
        push!(I, u); push!(J, v)
        push!(I, v); push!(J, u)
    end
    return sparse(I, J, ones(Float64, length(I)), n, n)
end


"""
    inc2elist(Ht)

Turn the transpose of one incidence matrix into a list of hyperedges.
"""
function inc2elist(
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer}
    n, m = size(Ht)
    hyperedges = Vector{Vector{Ti}}()
    for e = 1:m
        st = Ht.colptr[e]
        ed = Ht.colptr[e+1]-1
        push!(hyperedges, Ht.rowval[st:ed])
    end
    return hyperedges
end


"""
    H_largest_component(H, Ht)

Compute the largest component of the hypergraph. Here we assume the hypergraph is undirected.
""" 
function H_largest_component(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Ti <: Integer, Tf} 
    m, n = size(H)
    # connected component id of each vertex
    cc_id = zeros(Ti, n) 
    # connected component count
    cc_cnt = 0
    for i = 1:n
        if cc_id[i] != 0 # already visited
            continue
        end
        cc_cnt += 1
        cc_id[i] = cc_cnt
        q = Queue{Ti}()
        enqueue!(q, i)
        while !isempty(q)
            u = dequeue!(q)
            for j = H.colptr[u]:(H.colptr[u+1]-1)
                e = H.rowval[j]
                for k = Ht.colptr[e]:(Ht.colptr[e+1]-1)
                    v = Ht.rowval[k]
                    if cc_id[v] == 0
                        cc_id[v] = cc_cnt
                        enqueue!(q, v)
                    end
                end
            end
        end
    end
    cc_sz = zeros(Ti, cc_cnt)
    for i = 1:n
        cc_sz[cc_id[i]] += 1
    end
    # id of the largest connected component
    lcc_id = argmax(cc_sz) 
    # vertices in the largest connected component
    lcc = findall(x->cc_id[x]==lcc_id, 1:n) 
    _, lcc_eid = etouchS(H, lcc)
    return H[lcc_eid, lcc], Ht[lcc, lcc_eid], lcc
end


"""
    inc_sortbysize(H, Ht)

Sort rows(hyperedges) of an incidence matrix according to their sizes.
"""
function inc_sortbysize(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Ti <: Integer, Tf}
    order = H_order(Ht)  
    e_ord = sortperm(order)
    return H[e_ord, :], Ht[:, e_ord]
end


"""
    precision(S, Target)

Computing the precision of S with Target as the ground truth.
"""
function precision(
    S::Vector{Ti},
    Target::Vector{Ti},
)where {Ti}
    return length(intersect(S, Target)) / length(S)
end


"""
    recall(S, Target)

Computing the recall of S with Target as the ground truth.
"""
function recall(
    S::Vector{Ti},
    Target::Vector{Ti},
) where {Ti}
    return length(intersect(S, Target)) / length(Target)
end


"""
    F1score(S, Target)

Computing the F1 score of S with Target as the ground truth.
"""
function F1score(
    S::Vector{Ti},
    Target::Vector{Ti},
) where {Ti}
    return 2 * length(intersect(S, Target)) / (length(S) + length(Target))
end


"""
    etouch(H, S)

Count and list those edges incident to a set S.
"""
function etouchS(
    H::SparseMatrixCSC{Tf, Ti},
    S::Vector{Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    HS = H[:, S] 
    rp_S = HS.rowval 
    Sedges = unique(rp_S) 
    return length(Sedges), Sedges
end


"""
    einS(H, order, S)

Given the incidence matrix of a hypergraph and hyperedge orders, and a vertex set S,
compute the number of hyperedges fully contained in S and list those hyperedges.
"""
function einS(
    H::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    S::Vector{Ti},
) where {Tf, Ti <: Integer}
    HS = H[:, S]
    rp_HS = HS.rowval 
    etouchS_list = unique(rp_HS) 
    @inbounds for i in rp_HS
        order[i] -= 1 
    end
    einS_list = findall(x->x==0, order[etouchS_list])
    einS_list = etouchS_list[einS_list]
    einS_count = length(einS_list)
    @inbounds for i in rp_HS
        order[i] += 1
    end
    return einS_count, einS_list 
end


"""
    edensity(H, order, p, S)

Given a set S, the incidence matrix of a hypergraph,
order of hyperedges, and the vertex penalty vector, 
compute the edge density (e[S] - p(S)) / |S|.
"""
function edensity(
    H::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    p::Vector{Tf},
    S::Vector{Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    @assert length(S) > 0
    eS, _ = einS(H, order, S)
    pS = sum(p[S])
    return (eS - pS) / length(S) 
end


"""
    edensity(A, p, S)

Given a set S and the adjacency matrix of a graph, and
the vertex penalty vector p, compute the edge density (e[S] - p(S)) / |S|.
"""
function edensity(
    A::SparseMatrixCSC{Tf, Ti},
    p::Vector{Tf},
    S::Vector{Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    @assert length(S) > 0
    eS = sum(A[S, S]) / 2.0
    pS = sum(p[S])
    return (eS - pS) / length(S)
end


"""
    edensity(H, order, S)

Given a set S, and the incidence for a hypergraph and its 
edge orders, compute the edge density e[S]] / |S|.
"""
function edensity(
    H::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    S::Vector{Ti},
)where {Tf <: AbstractFloat, Ti <: Integer}
    @assert length(S) > 0
    eS, _ = einS(H, order, S)
    return eS / length(S) 
end


"""
    edensity(A, S)

Given a set S and the adjacency matrix for a graph,
compute the edge density e[S]] / |S|.
"""
function edensity(
    A::SparseMatrixCSC{Tf, Ti},
    S::Vector{Ti},
)where {Tf <: AbstractFloat, Ti <: Integer}
    @assert length(S) > 0
    eS = sum(A[S, S]) / 2.0 
    return eS / length(S)
end


"""
    conductance(H, Ht, S)

Given a set S, the incidence matrix of a hypergraph, and its 
tranpose, compute the conductance of S, which is
    e(S, \bar{S}) / min(vol(S), vol(\bar{S})),
and the cut is all-or-nothing cut and volume is the sum of degrees.
"""
function conductance(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    S::Vector{Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    order = H_order(Ht) 
    deg = H_deg(H) 
    volH = sum(order)
    eS, _ = einS(H, order, S)  
    eneighhood, _ = etouchS(H, S)
    cut = eneighhood - eS 
    volS = sum(deg[S])
    return cut / min(volH - volS, volS)
end


"""
    expansion(H, Ht, S)

Given a set S, the incidence matrix of a hypergraph, and its 
tranpose, compute the expansion of S, which is
    e(S, \bar{S}) / min(|S|,  |\bar{S}|),
and the cut is all-or-nothing cut.
"""
function expansion(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    S::Vector{Ti},
) where{Tf, Ti <: Integer}
    order = H_order(Ht) 
    deg = H_deg(H) 
    m, n = size(H)
    eS, _ = einS(H, order, S)  
    eneighhood, _ = etouchS(H, S)
    cut = eneighhood - eS 
    return cut / min(n - length(S), length(S))
end


"""
    βcut(H, order, p, S, β)

Given the incidence matrix of a hypergraph, and its
hyperedge orders, and the vertex penalty vector p,
and a vertex set S, and a scalar β, 
compute β |S| - e[S] + p(S).
"""
function βcut(
    H::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    p::Vector{Tf},
    S::Vector{Ti},
    β::Tf,
) where {Tf <: AbstractFloat, Ti <: Integer}
    eS, _ = einS(H, order, S)
    return β * length(S) - eS + sum(p[S])
end


"""
    βcut(H, order, S, β)

Given the incidence matrix of a hypergraph, and its
hyperedge orders, and a vertex set S, and a scalar β, 
compute β |S| - e[S].
"""
function βcut(
    H::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    S::Vector{Ti},
    β::Tf,
) where {Tf <: AbstractFloat, Ti <: Integer}
    eS, _ = einS(H, order, S)
    return β * length(S) - eS
end


"""
    DHSG_expansion_inc(H, Ht, order, Inf_cap)

Construct the flow network 
For each hyperedge e, if |e| > 2, we add a new node v_e and
build one directed edge from v to v_e for each v in e with capacity 1 / |e|,
  and one directed edge from v_e to v with capacity Inf_cap.
If |e| = 2, we build one directed edge from v to u with capacity 1 / 2,
                 and one directed edge from u to v with capacity 1 / 2.
"""
function DHSG_expansion_inc(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    Inf_cap::Tf,
) where {Tf <: AbstractFloat, Ti <: Integer}
    (m, n) = size(H)
    BigEdges = length(findall(x->x>2, order))
    N = n + BigEdges 

    Hyperedges = inc2elist(Ht)

    ## Build the adjacency matrix
    ap = n + 1   # "auxiliary node pointer", points to next "free" aux node
    U = Vector{Ti}()
    V = Vector{Ti}()
    vals = Vector{Tf}()

    for e in eachindex(Hyperedges)
        edge = Hyperedges[e]
        nv = length(edge)
        # We remove self-loops in preprocessing
        @assert nv > 1 
        if nv == 2
            i = edge[1]; j = edge[2]
            push!(U, i); push!(V, j); push!(vals, 1 / nv)
            push!(U, j); push!(V, i); push!(vals, 1 / nv)
        else 
            for i = edge 
                push!(U, i); push!(V, ap); push!(vals, 1 / nv)
                push!(U, ap); push!(V, i); push!(vals, Inf_cap) 
            end
            ap += 1
        end
    end
    A = sparse(U, V, vals, N, N)
    return A
end


"""
    clique_expansion(H, Ht; [type="UCE"])

This function expands one hypergraph into a normal graph
via expanding each hyperedge into a clique. 
If type is weighted, then the clique is weighted by 1 / |e|.
Else, the clique is unweighted.
"""
function clique_expansion(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti};
    type::String="UCE", # weighted or unweighted
) where {Tf, Ti <: Integer}
    @assert type in ["UCE", "WCE"] println("unsupported type of clique expansion.")
    m, n = size(H)     
    I = Ti[]; J = Ti[]; V = Float64[];
    for e in 1:m
        e_vertices = Ti[]
        for nzv in Ht.colptr[e]:Ht.colptr[e+1]-1
            push!(e_vertices, Ht.rowval[nzv])
        end
        e_size = length(e_vertices) 
        for nzu = 1:e_size
            u = e_vertices[nzu]
            for nzv = nzu+1:e_size
                v = e_vertices[nzv]
                push!(I, u); push!(J, v)
                push!(I, v); push!(J, u)
                if type == "WCE"
                    push!(V, 1.0 / e_size)
                    push!(V, 1.0 / e_size)
                else
                    push!(V, 1.0)
                    push!(V, 1.0)
                end
            end
        end
    end
    return sparse(I, J, V, n, n)
end


"""
    A_deg(A)

Weighted degree of each vertex in A 
"""
function A_deg(
    A::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer}
    return sum(A, dims=1)[1, :]
end


"""
    H_deg(H)

Degree of each vertex in H
"""
function H_deg(
    H::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer}
    return Ti.(sum(H, dims=1)[1, :])
end


"""
    H_order(Ht)

Size of each hyperedge in H
"""
function H_order(Ht::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer}
    return Ti.(sum(Ht, dims=1)[1, :])
end


"""
    H_fracdeg(H, Ht)

Fractional degree of each vertex in H,
i.e. the sum of 1 / |e| over all hyperedges e incident to v
"""
function H_fracdeg(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Tf, Ti <: Integer}
    m, n = size(H)
    # 1 / |e|
    order_rev = ones(Tf, m) 
    order = H_order(Ht)
    order_rev ./= order 
    frac_deg = zeros(Tf, n)
    for i in axes(H, 2) 
        for j = H.colptr[i]:(H.colptr[i+1]-1)
            e = H.rowval[j]
            frac_deg[i] += order_rev[e]
        end
    end
    return frac_deg
end


"""
    get_immediate_neighbors(A, R)

Get those vertices that are within 1-hop neighborhood of R,
excluding vertices from R.
"""
function get_immediate_neighbors(
    A::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
) where {Tf, Ti <: Integer} 
    Rn = unique(A[:, R].rowval)
    Rn = setdiff(Rn, R)
    return Rn
end


"""
    neighborhood(A, R)

Get those vertices that are within 1-hop neighborhood of R,
including vertices from R.
"""
function neighborhood(
    A::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
) where {Tf, Ti <: Integer}
    Rn = unique(A[:, R].rowval)
    Rn = union(Rn, R)
    return Rn
end


"""
    get_immediate_neighbors(H, Ht, R)

Get those vertices that are within 1-hop neighborhood of R,
excluding vertices from R.
"""
function get_immediate_neighbors(
    H::SparseMatrixCSC{Tf,Ti},
    Ht::SparseMatrixCSC{Tf,Ti},
    R::Vector{Ti}
) where {Tf, Ti <: Integer}

    Hr = H[:,R]
    rp_r = Hr.rowval
    R_edges = unique(rp_r)

    He = Ht[:,R_edges]
    rp_e = He.rowval
    Rneighbs = unique(rp_e)
    Rn = setdiff(Rneighbs,R)

    return Rn

end


"""
    neighborhood(H, Ht, R)

Get those vertices that are within 1-hop neighborhood of R,
including vertices from R.
"""
function neighborhood(
    H::SparseMatrixCSC{Tf,Ti},
    Ht::SparseMatrixCSC{Tf,Ti},
    R::Vector{Ti}
) where {Tf, Ti <: Integer}
    Hr = H[:,R]
    rp_r = Hr.rowval
    R_edges = unique(rp_r)

    He = Ht[:,R_edges]
    rp_e = He.rowval
    Rn = unique(rp_e)

    return Rn
end


"""
    volume_penalty(H, R, n, ϵ)

Compute the vertex penalty function p(v) = ϵ vol(v ∩ \bar{R})/2
"""
function volume_penalty(
    H::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    n::Ti,
    ϵ::Tf, 
) where {Tf <: AbstractFloat, Ti <: Integer}
    inR = ones(n)
    inR[R] .= 0
    deg = H_deg(H)
    return @. deg * inR * ϵ / 2  
end


"""
    volume_penalty(A, R, ϵ)

Compute the vertex penalty function p(v) = ϵ vol(v ∩ \bar{R})/2
"""
function volume_penalty(
    A::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    ϵ::Tf,
) where {Tf <: AbstractFloat, Ti <: Integer}
    n = size(A, 1)
    inR = ones(n)
    inR[R] .= 0
    deg = A_deg(A)
    return @. deg * inR * ϵ / 2
end


"""
    frac_volume_penalty(H, Ht, R, n, ϵ)

Compute the vertex penalty function p(v) = ϵ \fracvol(v ∩ \bar{R})
"""
function frac_volume_penalty(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    n::Ti,
    ϵ::Tf,
) where {Tf <: AbstractFloat, Ti <: Integer}
    inR = ones(n)
    inR[R] .= 0
    fracdeg = H_fracdeg(H, Ht)
    return @. fracdeg * inR * ϵ
end


"""
    subgraph(H, Ht, R)

Extract the subgraph induced by R. 
Include those vertices in R and edges fully contained in R.
"""
function subgraph(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    order = H_order(Ht)
    _, einR_list = einS(H, order, R)
    return H[einR_list, R]
end


"""
    subgraph_var(H, R)

Extract the subgraph induced by R.
Include those vertices in R and edges containing vertices from R.
Shrink the hyperedges to only contain vertices in R.
"""
function subgraph_var(
    H::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    _, elist = etouchS(H, R)
    return H[elist, R]
end


"""
    bfs(H, Ht, s)

Perform a bfs starting from s.
"""
function bfs(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    s::Ti, 
) where {Tf <: AbstractFloat, Ti <: Integer}
    _, n = size(H)
    q = Queue{Ti}()
    dis = fill(Tf(Inf), n)
    dis[s] = 0
    enqueue!(q, s)
    while !isempty(q)
        u = dequeue!(q)
        for nzi in H.colptr[u]:H.colptr[u+1]-1 
            e = H.rowval[nzi]
            for nzj in Ht.colptr[e]:Ht.colptr[e+1]-1
                v = Ht.rowval[nzj]
                if dis[v] == Inf
                    dis[v] = dis[u] + 1
                    enqueue!(q, v)
                end
            end
        end
    end
    return dis
end


"""
    diameter(H, Ht)

Compute the diameter of a hypergraph, i.e. the longest shortest path.
We assume the hypergraph is connected.
"""
function diameter(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    dis = bfs(H, Ht, 1)
    s = argmax(dis) 
    dis = bfs(H, Ht, s)
    return maximum(dis)
end


"""
    diameter(H, Ht)

Show some basic statistics for a graph.
"""
function stats(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Tf <: AbstractFloat, Ti <: Integer}
    deg = H_deg(H)
    @printf("Vertex degree: mean = %.1lf, max = %d, min = %d\n",
        mean(deg), maximum(deg), minimum(deg))
    
    fracdeg = H_fracdeg(H, Ht)
    @printf("Fractional vertex degree: mean = %.1lf, max = %.1lf, min = %.1lf\n",
        mean(fracdeg), maximum(fracdeg), minimum(fracdeg))

    order = H_order(Ht)
    @printf("Edge order: mean = %.1lf, max = %d, min = %d\n",
        mean(order), maximum(order), minimum(order))
    
    Delta = diameter(H, Ht)
    @printf("Diameter = %d\n", Delta)
end


"""
    randomwalk(H, Ht, seed, k)

Perform a simple length-k random walk starting from seed.
"""
function randomwalk(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    seed::Ti,
    k::Ti;
) where {Tf <: AbstractFloat, Ti <: Integer} 
    RW = Ti[]
    push!(RW, seed)
    for i = 1:k
        nze = rand(H.colptr[seed]:H.colptr[seed+1]-1)
        e = H.rowval[nze]
        nzv = rand(Ht.colptr[e]:Ht.colptr[e+1]-1)
        v = Ht.rowval[nzv]
        seed = v
        push!(RW, v)
    end
    return RW
end


"""
    generate_R(H, Ht, R, k, [type="random", prob=0.5])

Expand an anchor set of size k from R. 
If type is RW, then we keep randomly selecting a seed from R 
and perform length-2 random walks until visited k vertices.
If type is RW-geo, then we keep randomly selecting a seed from R
and perform a random walk with geometrically distributed length
until visited k vertices.
"""
function generate_R(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    k::Ti,
    type::String="random";
    prob::Tf = 0.5,
)where {Tf <: AbstractFloat, Ti <: Integer}
    @assert length(R) <= k
    if type == "RW"
        inR = Dict{Ti, Bool}()
        Rn = Ti[]
        for v in R
            inR[v] = true
            push!(Rn, v)
        end
        while length(Rn) < k
            seed = rand(R)
            RW = randomwalk(H, Ht, seed, 2)
            for v in RW
                if length(Rn) >= k
                    break
                end
                if !haskey(inR, v)
                    inR[v] = true
                    push!(Rn, v)
                end
            end
        end
        return Rn
    elseif type == "RW-geo"
        inR = Dict{Ti, Bool}()
        Rn = Ti[]
        p = map(i->((1 - prob)^(i-1) * prob), 1:40) 
        for v in R
            inR[v] = true
            push!(Rn, v)
        end
        while length(Rn) < k
            seed = rand(R)
            len = wsample(1:40, p)
            RW = randomwalk(H, Ht, seed, len) 
            for v in RW
                if length(Rn) >= k
                    break
                end
                if !haskey(inR, v)
                    inR[v] = true
                    push!(Rn, v)
                end
            end
        end
        return Rn
    else
        error("not supported generator.")
    end
end