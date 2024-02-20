"""
    DSG_flow_binary_search(A; [flowtol=zero(Tf)])

Given the adjacency matrix of a graph, 
compute the densest subgraph using binary search 
on the answer, and then solve each subproblem using flow.

# Arguments
- A: the incidence matrix of a graph 
- flowtol: the tolerance for the flow algorithm
"""
function DSG_flow_binary_search(
    A::SparseMatrixCSC{Tf, Ti};
    flowtol::Tf = zero(Tf),
) where {Tf <: AbstractFloat, Ti <: Integer}
    total_dt = @elapsed begin
        n = size(A, 1)
        deg = A_deg(A)
        m = div(sum(A), 2)

        # initialize binary search lower/upper bound
        # and optimal answer
        lb = m / n # when S = V, density is m / n 
        ub = maximum(deg) # when each vertex has the same degree.
        optS = Vector(1:n)

        # build the auxiliary graph for flow method
        iter_dts = Tf[]
        niter = 0
        # termination condition
        # stop when the interval is short enough
        while ub - lb > 1 / n / (n - 1)  
            niter += 1
            iter_dt = @elapsed begin
                β = (ub + lb) / 2
                S, _ = DSG_flow_step(A, Vector(1:n), deg, β; flowtol=flowtol)
                # check if β|S| - e[S] < 0
                if  length(S) > 0 && edensity(A, S) > β  
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
        "optval" => edensity(A, optS),
        "optsol" => optS,
        "total_dt" => total_dt,
        "dt_per_iter" => mean(iter_dts),
        "niter" => niter,
    )
end


"""
    DSG_flow_binary_search(A; [flowtol=zero(Tf)])

Given the adjacency matrix of a graph, 
compute the densest subgraph using our density improvement
framework, and then solve each subproblem using flow.

# Arguments
- A: the adjacency matrix of a graph
- flowtol: the tolerance for the flow algorithm
"""
function DSG_flow_density_improvement(
    A::SparseMatrixCSC{Tf, Ti};
    flowtol::Tf = zero(Tf),
) where{Tf <: AbstractFloat, Ti <: Integer}
    total_dt = @elapsed begin
        # precompute β
        n = size(A, 1)
        deg = A_deg(A)
        still_improving = true
        β = sum(A) / 2 / n
        # initialize the solution as V
        S = Vector{Ti}(1:n)
        iter_dts = Tf[] 
        niter = 0
        # stop when no improvement can be made
        while still_improving
            niter += 1
            iter_dt = @elapsed begin
                still_improving = false
                S_new, _ = DSG_flow_step(A, Vector{Ti}(1:n), deg, β; flowtol=flowtol)
                # check if it's possible to find a set with higher density
                if length(S_new) > 0 && edensity(A, S_new) > β 
                    S = S_new
                    β = edensity(A, S_new)
                    still_improving = true
                end
            end
            push!(iter_dts, iter_dt)
        end
    end
    return Dict(
        "optval" => edensity(A, S),
        "optsol" => S,
        "total_dt" => total_dt,
        "dt_per_iter" => mean(iter_dts),
        "niter" => niter,
    )
end


"""
    DSG_flow_step(A, V, deg, β; [flowtol=zero(Tf)])

Solve the decision problem that if there is a set S with 
    e[S] > β|S| using flow method.

# Arguments
- A: the auxiliary graph
- V: the set of vertices from the original hypergraph
- deg: the degree of each vertex in V
- β: parameter for the decision problem
- flowtol: the tolerance for the flow algorithm
"""
function DSG_flow_step(
    A::SparseMatrixCSC{Tf, Ti},
    V::Vector{Ti},
    deg::Vector{Tf},
    β::Tf;
    flowtol::Tf = zero(Tf),
)where {Tf <: AbstractFloat, Ti <: Integer}
    N = size(A, 1)
    sVec = zeros(N)
    tVec = zeros(N)

    # set up the terminal edges
    @. sVec[V] = deg / 2 
    tVec[V] .= β
    minVec = min.(sVec, tVec) 
    # pre-route the flow
    sVec .-= minVec
    tVec .-= minVec

    # remove those zero edges
    sVec = sparse(sVec)
    tVec = sparse(tVec)
    dropzeros!(sVec)
    dropzeros!(tVec)
    C = [spzeros(1, 1) sVec' spzeros(1, 1);
         spzeros(N, 1) (A./2) tVec;
         spzeros(1, 1) spzeros(1, N) spzeros(1, 1);]
    F = hlpp.maxflow(C, flowtol)
    # notice that hlpp only computes a preflow,
    # so it's hard to implement the source_nodes_min
    Src = hlpp.source_nodes(F)[2:end] .- 1
    S = intersect(Src, V)
    # add the pre-routed flow back
    return S, F.cutvalue + sum(minVec)
end
