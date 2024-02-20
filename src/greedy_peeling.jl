"""
Naive implementation of Greedy++ using BBST
original C++ code: 
https://www.dropbox.com/s/jzouo9fjoytyqg3/code-greedy%2B%2B.zip?dl=0&file_subpath=%2Fcode-greedy%2B%2B
"""
function greedy_peeling(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    p::Vector{Tf},
) where {Ti <: Integer, Tf <: AbstractFloat} 
    total_dt = @elapsed begin
        # making it faster to access vertices inside one hyperedge
        m, n = size(H) # number of edges, vertices

        # compute the degrees of vertices
        deg = H_deg(H)

        # objective (e[S] - ϵ * Vol(S ∩ \bar{R})) / |S| 
        total_wts = m - sum(p) 
        best_ans = total_wts / n 
        best_size = n
        best_S = Vector(1:n)

        peeling_ord = zeros(Ti, n)

        # indicate  
        # whether one node has been peeled off or not
        # whether one edge has already not been fully contained
        exists_v = ones(Bool, n) 
        contained_e = ones(Bool, m)

        # keep vertices sorted by marginal reward 
        pq = PriorityQueue{Ti, Tf}() 
        curr_deg = zeros(n) 

        total_wts = m - sum(p)
        curr_wts = total_wts

        for u in 1:n 
            curr_deg[u] = deg[u] - p[u] 
            exists_v[u] = true
            enqueue!(pq, u=>curr_deg[u])
        end
        for e in 1:m
            contained_e[e] = true  
        end
        for i in 1:n
            # delete the vertex with the smallest marginal reward
            u = dequeue!(pq)
            exists_v[u] = false
            peeling_ord[i] = u
            curr_wts -= curr_deg[u]
            for nzi in H.colptr[u]:(H.colptr[u+1]-1)
                e = H.rowval[nzi]
                # if this edge has already been not fully contained
                # then we skip
                if !contained_e[e]
                    continue
                end
                contained_e[e] = false 
                for nzj in Ht.colptr[e]:(Ht.colptr[e+1]-1) 
                    v = Ht.rowval[nzj]
                    # we only process those vertices that
                    # haven't been peeled off
                    if !exists_v[v]
                        continue
                    end
                    # we stick to this simple implementation first
                    # the number of delete/push is bounded by vol(H)
                    curr_deg[v] -= Ht.nzval[nzj] 
                    pq[v] = curr_deg[v]
                end
            end
            if i < n
                curr_density = curr_wts / (n - i)
                if curr_density > best_ans
                    best_ans = curr_density
                    best_size = n - i
                end
            end
        end
        best_S = peeling_ord[n - best_size + 1:n]
    end
    return Dict(
        "optval" => best_ans, 
        "optsol" => best_S,
        "peeling_ord" => peeling_ord,
        "total_dt" => total_dt, 
    )
end


