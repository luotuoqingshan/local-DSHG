"""
    ADHSG_flow_density_improvement(H, Ht, R, ϵ, [penalty="vol"]; [type="global", flowtol=zero(Tf)])

Compute the set S maximizing the anchored densest subhypergraph objective(Problem 5, 6 in the paper).

# Arguments
- H: the incidence matrix of a hypergraph of shapse |E|x|V|
- Ht: the transpose of H
- R: the set of anchor vertices
- ϵ: the parameter for the penalty
- penalty: the type of penalty, either "vol" or "fracvol"
- type: the type of algorithm, either "global" or "local"
- flowtol: the tolerance for the flow algorithm
"""
function ADHSG_flow_density_improvement(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    R::Vector{Ti},
    ϵ::Tf,
    penalty::String="vol";
    type::String="global", # "global" or "local"
    flowtol=zero(Tf),
) where{Tf <: AbstractFloat, Ti <: Integer}
    total_dt = @elapsed begin
        # precompute beta
        m, n = size(H)
        order = H_order(Ht)
        fracdeg = H_fracdeg(H, Ht)
        Inf_cap = 2.0 * m * n

        # compute the penalty vector
        if penalty == "vol"
            p = volume_penalty(H, R, n, ϵ)
        elseif penalty == "fracvol"
            p = frac_volume_penalty(H, Ht, R, n, ϵ)
        else
            error("Penalty type not supported.")
        end
        p = Tf.(p)
        
        inR = zeros(n)
        inR[R] .= Inf_cap

        # initialize the optimal solution
        H_R = subgraph(H, Ht, R)
        H_R_res = DHSG_flow_density_improvement(H_R, sparse(H_R'); flowtol=flowtol)
        optS = R[H_R_res["optsol"]] 
        @assert length(optS) > 0 
        @assert edensity(H, order, optS) == H_R_res["optval"] println("Make sure R does not contain duplicate vertices.") 

        niter = 0
        iter_dts = Vector{Tf}()
        still_improving = true
        β = max(0, edensity(H, order, optS))

        if type == "global"
            A = DHSG_expansion_inc(H, Ht, order, Inf_cap)
        elseif type == "local"
            # set up S and its proper neighborhood
            # S are those vertices with terminal edges from source s. 
            if ϵ < 1.0
                S = findall(x->fracdeg[x] - p[x] - β > flowtol, 1:n) 
                Sn = setdiff(neighborhood(H, Ht, S), S)
            else
                S = R
                Sn = setdiff(neighborhood(H, Ht, S), S)
            end
        else
            error("type should be either global or local")
        end
        while still_improving
            niter += 1
            iter_dt = @elapsed begin
                still_improving = false
                if type == "global"
                    S_new, _ = ADHSG_flow_step_global(
                        A, Vector(1:n), p, fracdeg,  β; flowtol=flowtol,   
                    )
                else 
                    S_new = ADHSG_flow_step_local(
                        H, Ht, order, fracdeg, p, S, Sn, β, Inf_cap; flowtol=flowtol,
                    )
                end
                if length(S_new) > 0 && edensity(H, order, p, S_new) > β 
                    β = edensity(H, order, p, S_new)
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
    ADHSG_flow_step_global(A, V, p, fracdeg, β; flowtol=zero(Tf))

Solve the decision problem of Problem 5,6 in the paper.

# Arguments
- A: the auxiliary graph
- V: the set of vertices from the original hypergraph
- p: the penalty vector
- fracdeg: the fractional degree vector
- β: parameter for the decision problem
- flowtol: the tolerance for the flow algorithm
"""
function ADHSG_flow_step_global(
    A::SparseMatrixCSC{Tf, Ti},
    V::Vector{Ti},
    p::Vector{Tf},
    fracdeg::Vector{Tf},
    β::Tf;
    flowtol::Tf = zero(Tf),
) where {Tf <: AbstractFloat, Ti <: Integer}
    N = size(A, 1)
    sVec = zeros(N)
    tVec = zeros(N)

    # set up the terminal edges
    sVec[V] .= fracdeg  
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
         spzeros(N, 1) A tVec;
         spzeros(1, 1) spzeros(1, N) spzeros(1, 1);]
    F = hlpp.maxflow(C, flowtol)

    Src = hlpp.source_nodes(F)[2:end].-1
    S = intersect(Src, V)
    return S, F.cutvalue + sum(minVec) 
end


"""
    ADHSG_flow_step_local(H, Ht, order, fracdeg, p, R, Rn, β, Inf_cap; [flowtol=zero(Tf)])

Solve the decision problem of Problem 5,6 in the paper using the local flow method.

# Arguments
- H: the incidence matrix of a graph
- Ht: the transpose of H
- order: the order of the edges
- fracdeg: the fractional degree vector
- p: the penalty vector
- R: the set of vertices explored, detailed definition is in the paper
- Rn: the proper neighborhood of R
- β: parameter for the decision problem 
- Inf_cap: the capacity of edges which will never be saturated 
- flowtol: the tolerance for the flow algorithm
"""
function ADHSG_flow_step_local(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    order::Vector{Ti},
    fracdeg::Vector{Tf},
    p::Vector{Tf},
    R::Vector{Ti},
    Rn::Vector{Ti},
    β::Tf,
    Inf_cap::Tf;
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
    Hc = H[:, C_global]
    rp_c = Hc.rowval
    L_edges = unique(rp_c)

    # the local hypergraph
    HL = H[L_edges, Local2Global]
    order_L = order[L_edges]
    A_L = DHSG_expansion_inc(HL, sparse(HL'), order_L, Inf_cap)
    n_L = length(Local2Global)

    # solve the flow problem on the local hypergraph
    S_local, _ = ADHSG_flow_step_global(
        A_L, Vector{Ti}(1:n_L), p[Local2Global],
        fracdeg[Local2Global], β; flowtol=flowtol)
    
    E_local = intersect(S_local, I_local) # new explored
    E_global = Local2Global[E_local]

    while length(E_local) > 0

        # Update which nodes are explored and which are unexplored but adjacent to C
        C_local = [C_local; E_local]
        C_global = Local2Global[C_local]

        # Take these away from I_local
        I_local = setdiff(I_local,E_local)

        # This is better
        Nbs_of_E = get_immediate_neighbors(H,Ht,E_global)
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

        # Now we have a new set of explored and unexplored edges,
        # we do the same thing over again to find a localize min-cut
        Hc = H[:,C_global]
        rp_c = Hc.rowval
        L_edges = unique(rp_c)

        # Binary indicence matrix for the local hypergraph (without terminal edges)
        HL = H[L_edges,Local2Global]
        order_L = order[L_edges]

        # Expand into a directed graph
        A_L = DHSG_expansion_inc(HL, sparse(HL'), order_L,Inf_cap)
        n_L = length(Local2Global)  # number of non-auxiliary nodes in A_L

        # Find the first mincut, which can be done by calling HLC_Step
        # with localized objects
        S_local, _ = ADHSG_flow_step_global(
            A_L, Vector{Ti}(1:n_L), p[Local2Global],
            fracdeg[Local2Global], β; flowtol=flowtol)
        # Find nodes to "expand" around:
        #   any nodes in the cut set tha are "unexplored" still
        E_local = intersect(S_local,I_local)
        E_global = Local2Global[E_local]
    end

    return Local2Global[S_local]
end