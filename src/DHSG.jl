"""
    DHSG_flow_binary_search(H, Ht; [flowtol=zero(Tf)])

Given the incidence matrix of a hypergraph and its transpose,
compute the densest subhypergraph using binary search
on the answer, and then solve each subproblem using flow. 

# Arguments
- H: the incidence matrix of a hypergraph
- Ht: the transpose of H
- flowtol: the tolerance for the flow algorithm
"""
function DHSG_flow_binary_search(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti};
    flowtol::Tf = zero(Tf),
) where {Tf <: AbstractFloat, Ti <: Integer}
    total_dt = @elapsed begin
        m, n = size(H)

        # compute edge order, vertex fractional degree  
        order = H_order(Ht)
        fracdeg = H_fracdeg(H, Ht)

        # initialize binary search lower/upper bound
        # and optimal answer
        lb = m / n # when S = V, density is m / n 
        ub = maximum(fracdeg) # when each vertex has the same degree.
        optS = Vector(1:n)

        Inf_cap = 2.0 * m * n # infinity capacity

        # build the auxiliary graph for flow method
        A = DHSG_expansion_inc(H, Ht, order, Inf_cap)
        iter_dts = Tf[]
        niter = 0
        # termination condition
        # stop when the interval is short enough
        while ub - lb > 1 / n / (n - 1)  
            niter += 1
            iter_dt = @elapsed begin
                β = (ub + lb) / 2
                S, _ = DHSG_flow_step(A, Vector(1:n), fracdeg, β; flowtol=flowtol)
                # check if β|S| - e[S] < 0
                if βcut(H, order, S, β) < 0  
                    lb = β 
                    optS = S
                else
                    ub = β
                end
            end
            push!(iter_dts, iter_dt)
        end
    end
    return Dict(
        "optval" => edensity(H, order, zeros(Tf, n), optS),
        "optsol" => optS,
        "total_dt" => total_dt,
        "dt_per_iter" => mean(iter_dts),
        "niter" => niter,
    )
end


"""
    DHSG_flow_density_improvement(H, Ht; [flowtol=zero(Tf)])

Given the incidence matrix of a hypergraph and its transpose,
compute the densest subhypergraph using our density improvement
framework, and then solve each subproblem using flow. 

# Arguments
- H: the incidence matrix of a hypergraph
- Ht: the transpose of H
- flowtol: the tolerance for the flow algorithm
"""
function DHSG_flow_density_improvement(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti};
    flowtol::Tf = zero(Tf),
) where{Tf <: AbstractFloat, Ti <: Integer}
    total_dt = @elapsed begin
        # precompute β
        m, n = size(H)

        # compute edge order, vertex fractional degree
        order = H_order(Ht)
        fracdeg = H_fracdeg(H, Ht)

        Inf_cap = 2.0 * m * n

        A = DHSG_expansion_inc(H, Ht, order, Inf_cap)
        still_improving = true
        β = m / n
        S = Vector{Ti}(1:n)
        iter_dts = Tf[] 
        niter = 0
        # stop when no improvement can be made
        while still_improving
            niter += 1
            iter_dt = @elapsed begin
                still_improving = false
                S_new, _ = DHSG_flow_step(A, Vector{Ti}(1:n), fracdeg, β; flowtol=flowtol)
                # check if it's possible to find a set with higher density
                if length(S_new) > 0 && edensity(H, order, S_new) > β 
                    S = S_new
                    β = edensity(H, order, S_new)
                    still_improving = true
                end
            end
            push!(iter_dts, iter_dt)
        end
    end
    return Dict(
        "optval" => edensity(H, order, S),
        "optsol" => S,
        "total_dt" => total_dt,
        "dt_per_iter" => mean(iter_dts),
        "niter" => niter,
    )
end


"""
    DHSG_flow_step(A, V, fracdeg, β; [flowtol=zero(Tf)])

Solve the decision problem that if there is a set S with 
    e[S] > β|S| using flow method.

# Arguments
- A: the auxiliary graph
- V: the set of vertices from the original hypergraph
- fracdeg: the fractional degree of each vertex in V
- β: parameter for the decision problem
- flowtol: the tolerance for the flow algorithm
"""
function DHSG_flow_step(
    A::SparseMatrixCSC{Tf, Ti},
    V::Vector{Ti},
    fracdeg::Vector{Tf},
    β::Tf;
    flowtol::Tf = zero(Tf),
)where {Tf <: AbstractFloat, Ti <: Integer}
    N = size(A, 1)
    sVec = zeros(N)
    tVec = zeros(N)

    # set up terminal edges
    sVec[V] .= fracdeg  
    tVec[V] .= β

    # pre-route the flow
    minVec = min.(sVec, tVec) 
    sVec .-= minVec
    tVec .-= minVec

    # remove those zero edges
    sVec = sparse(sVec)
    tVec = sparse(tVec)
    dropzeros!(sVec)
    dropzeros!(tVec)
    C = [spzeros(1, 1) sVec' spzeros(1, 1);
         spzeros(N, 1) A tVec;
         spzeros(1, 1) spzeros(1, N) spzeros(1, 1);]
    F = hlpp.maxflow(C, flowtol)
    # notice that our hlpp only computes a preflow,
    # so it's hard to implement the source_nodes_min
    Src = hlpp.source_nodes(F)[2:end] .- 1
    S = intersect(Src, V)
    # add the pre-routed flow back
    return S, F.cutvalue + sum(minVec)
end
