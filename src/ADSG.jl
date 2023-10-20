"""
    ADSG_flow_density_improvement(A, R, ϵ; [flowtol=zero(Tf), type="global"])

Compute the set S maximizing the anchored densest subgraph objective
max_{S} (e[S] - ϵ Vol(S cap R)/2) / |S|.

# Arguments
- A: the adjacency matrix of a graph
- R: the set of anchor vertices
- ϵ: the volume penalty parameter
- flowtol: the tolerance for the flow algorithm
- type: "global" indicates using the global flow method, 
    and "local" indicates using the local flow method.
"""
function ADSG_flow_density_improvement(
    A::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    ϵ::Tf;
    flowtol=zero(Tf),
    type::String="global",
) where{Tf <: AbstractFloat, Ti <: Integer}
    total_dt = @elapsed begin
        # precompute beta
        n = size(A, 1)
        deg = A_deg(A)

        p = volume_penalty(A, R, ϵ) 
        p = Tf.(p)
        
        # initialize the solution as the densest subgraph 
        # of the graph induced by R
        A_R = A[R, R] 
        A_R_res = DSG_flow_density_improvement(A_R; flowtol=flowtol)

        optS = R[A_R_res["optsol"]] 
        @assert length(optS) > 0 
        @assert edensity(A, optS) == A_R_res["optval"] println("Make sure that R does not contain duplicate vertices.")

        niter = 0
        iter_dts = Vector{Tf}()
        still_improving = true
        β = max(0, A_R_res["optval"])

        if type == "local"
            # set up S and its proper neighborhood
            # S are those vertices with terminal edges from source s. 
            if ϵ < 1.0
                S = findall(x->deg[x]/2 - p[x] - β > flowtol, 1:n) 
                Sn = setdiff(neighborhood(A, S), S)
            else
                S = R
                Sn = setdiff(neighborhood(A, S), S)
            end
        end
        # stop when no improvement can be made
        while still_improving
            niter += 1
            iter_dt = @elapsed begin
                still_improving = false
                if type == "global"
                    S_new, _ = ADSG_flow_step_global(
                        A, Vector(1:n), p, deg, β; flowtol=flowtol,    
                    )
                else 
                    S_new = ADSG_flow_step_local(
                        A, deg, p, S, Sn, β; flowtol=flowtol,
                    )
                end
                if length(S_new) > 0 && edensity(A, p, S_new) > β 
                    β = edensity(A, p, S_new)
                    optS = S_new
                    still_improving = true
                end
            end
            push!(iter_dts, iter_dt)
        end
    end
    return Dict(
        "optval" => β,
        "optsol" => optS,
        "total_dt" => total_dt, 
        "dt_per_iter" => mean(iter_dts),
        "niter" => niter,
    )
end


"""
    ADSG_flow_density_improvement(H, Ht, R, ϵ; [flowtol=zero(Tf), type="global", expansion="WCE"])

Compute the set S maximizing the anchored densest subgraph objective 
where the graph is the clique expansion of the hypergraph given by incidence matrices.

# Arguments
- H: the incidence matrix of a hypergraph
- expansion: type of clique expansion
"""
function ADSG_flow_density_improvement(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    ϵ::Tf;
    flowtol=zero(Tf),
    type::String="global",
    expansion::String="WCE",
) where{Tf <: AbstractFloat, Ti <: Integer}
    preprocess_dt = @elapsed begin
        A = clique_expansion(H, Ht; type = expansion)
    end
    di_res = ADSG_flow_density_improvement(
        A, R, ϵ; flowtol=flowtol, type=type)
    di_res["total_dt"] += preprocess_dt
    return di_res
end


"""
    ADSG_flow_step_global(A, V, p, deg, β; [flowtol=zero(Tf)])

Solve the decision problem that if there is a set S with
    e[S] - ϵ Vol(S cap R)/2 > β|S| using flow method.

# Arguments
- A: the adjacency matrix of a graph
- V: the set of vertices from the original hypergraph
- p: the penalty vector
- deg: the degree vector
- β: parameter for the decision problem
- flowtol: the tolerance for the flow algorithm
"""
function ADSG_flow_step_global(
    A::SparseMatrixCSC{Tf, Ti},
    V::Vector{Ti},
    p::Vector{Tf},
    deg::Vector{Tf},
    β::Tf;
    flowtol::Tf = zero(Tf),
) where {Tf <: AbstractFloat, Ti <: Integer}
    N = size(A, 1)
    sVec = zeros(N)
    tVec = zeros(N)

    # set up the terminal edges
    @. sVec[V] = deg / 2  
    @. tVec[V] = p + β
    minVec = min.(sVec, tVec) 
    # pre-route the flow
    sVec .-= minVec
    tVec .-= minVec

    sVec = sparse(sVec)
    tVec = sparse(tVec)
    dropzeros!(sVec)
    dropzeros!(tVec)
    C = [spzeros(1, 1) sVec' spzeros(1, 1);
         spzeros(N, 1) (A./2) tVec;
         spzeros(1, 1) spzeros(1, N) spzeros(1, 1);]
    F = hlpp.maxflow(C, flowtol)

    Src = hlpp.source_nodes(F)[2:end].-1
    S = intersect(Src, V)
    return S, F.cutvalue + sum(minVec) 
end


"""
    ADSG_flow_step_local(A, deg, p, R, Rn, β; [flowtol=zero(Tf)])

Solve the decision problem that if there is a set S with
e[S] - ϵ Vol(S cap R)/2 > β|S| using local flow method,
which starts from a local  graph and gradually grow 
based on the flow solution.

# Arguments
- A: the adjacency matrix of a graph
- deg: the degree vector
- p: the penalty vector
- R: the set of vertices explored, detailed definition is in the paper
- Rn: the proper neighborhood of R
- β: parameter for the decision problem 
- flowtol: the tolerance for the flow algorithm
"""
function ADSG_flow_step_local(
    A::SparseMatrixCSC{Tf, Ti},
    deg::Vector{Tf},
    p::Vector{Tf},
    R::Vector{Ti},
    Rn::Vector{Ti},
    β::Tf;
    flowtol::Tf = 0.0,
) where {Tf <: AbstractFloat, Ti <: Integer}
    # id mapping from local graph to global graph
    Local2Global = [R; Rn]

    Lsize = length(Local2Global)
    # C is the set of explored nodes
    # I is the set of unexplored nodes and they are amoung the neighbors of C
    # local and global denote the local and global indices
    C_global = R
    I_global = Rn

    C_local = collect(1:length(R))
    I_local = collect(length(R)+1:Lsize)
    R_local = collect(1:length(R))

    # the local graph 
    A_L = A[Local2Global, Local2Global]
    n_L = length(Local2Global)

    # solve the flow problem on the local graph
    S_local, _ = ADSG_flow_step_global(
        A_L, Vector{Ti}(1:n_L), p[Local2Global],
        deg[Local2Global], β; flowtol=flowtol)
    
    E_local = intersect(S_local, I_local) # new explored
    E_global = Local2Global[E_local]

    while length(E_local) > 0

        # Update which nodes are explored and which are unexplored but adjacent to C
        C_local = [C_local; E_local]
        C_global = Local2Global[C_local]

        # Take these away from I_local
        I_local = setdiff(I_local,E_local)

        # This is better
        Nbs_of_E = get_immediate_neighbors(A,E_global)
        Lnew = setdiff(Nbs_of_E,Local2Global)
        numNew = length(Lnew)
        # Update the set of indices in L
        Local2Global = [Local2Global; Lnew]

        # Store local indices for new nodes added to L
        Lnew_local = collect((Lsize+1):(Lsize+numNew))
        Lsize = length(Local2Global)

        # These are going to be "unexplored" nodes
        I_local = [I_local; Lnew_local]
        I_global = Local2Global[I_local]

        # Expand into a directed graph
        A_L = A[Local2Global, Local2Global] 
        n_L = length(Local2Global)  # number of non-auxiliary nodes in A_L

        # Find the mincut on the local graph
        S_local, _ = ADSG_flow_step_global(
            A_L, Vector{Ti}(1:n_L), p[Local2Global],
            deg[Local2Global], β; flowtol=flowtol)
        # Find nodes to "expand" around:
        #   any nodes in the cut set tha are "unexplored" still
        E_local = intersect(S_local,I_local)
        E_global = Local2Global[E_local]
    end

    return Local2Global[S_local]
end

