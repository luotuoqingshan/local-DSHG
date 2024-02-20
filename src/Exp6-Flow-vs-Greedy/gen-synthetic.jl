"""
  _draw_with_duplicate_check(e, rng, elements, duplicates)

draw one element from elements using random generator rng, and add it to e.
If duplicates is true, then we allow duplicates, 
otherwise we redraw.
"""
function _draw_with_duplicate_check(e, rng, elements, duplicates)
  while true
    val = rand(rng, elements) # draw a random elements
    if duplicates == true 
      return val 
    else 
      duplicate = false 
      for curval in e
        if val == curval # found a duplicate 
          duplicate = true 
          break 
        end
      end
      if duplicate
        continue # continue the outer loop
      else
        return val
      end 
    end
  end 
end 


"""
  random_hyperedge!(e, elements, max_size, p; rng, duplicates)

Generate a random hyperedge where the hyperedge keeps growing with probability p.
If hyperedge size is already max_size, then we stop growing.
Vertices are drawn from elements.
"""
function random_hyperedge!(e, elements, max_size, p; rng, duplicates)
  push!(e, rand(rng, elements))
  push!(e, _draw_with_duplicate_check(e, rng, elements, duplicates))
  for i in 3:max_size
    rv = rand(rng) 
    if rv <= p 
      # then the hyperedge keeps growing... 
      push!(e, _draw_with_duplicate_check(e, rng, elements, duplicates))
    else
      # we stop the hyperedge... 
      break
    end 
  end 
  return e
end 
"""
    random_hyperedge(elements, max_size, p; [rng=GLBOAL_RNG, duplicates=false])

Generate a random hyperedge where the hyperedge keeps growing with probability p. 
(Seeting p=0 will generate a random edge.) The hyperedge will grow up to max_size. 
Elements of the hyperedge are chosen by `rand(rng, elements)`.

We do not allow duplicate elements in the hyperedge if `duplicates` is false. 
**This may cause problems for sets of small elements and cause the algorithm to run forever.**
""" 
function random_hyperedge(elements, max_size::Int, p::Real; rng = Random.GLOBAL_RNG, duplicates::Bool=false)
  e = Vector{eltype(elements)}()
  return random_hyperedge!(e, elements, max_size, p; rng, duplicates)
end 

"""
    random_hypergraph(n, nedges, hep; [rng = Random.GLOBAL_RNG, duplicates=false, max_size])
    random_hypergraph(elements, nedges, hep; [rng = Random.GLOBAL_RNG, duplicates=false, max_size])

Generate a random hypergraph where elements are chosen from 1:n. 
- the default value of max_size is based on twice the expected length of the hyperedge continuation probability.
  So if this is 0.4 (i.e. short hyperedges), then max_size = 2*ceil(Int, 0.4/0.6 ) = 4. 
- duplicates refers to duplicate values within a hyperedge. This function may produce 
  duplicate hyperedges and does not check for uniqueness. If you need that, call 
  `unique!` on the final list of hyperedges. If you need unique sorted hyperedges, then
  use `unique!(sort!.(random_hypergraph()))`
""" 
function random_hypergraph(elements, nedges::Int, hep::Real; rng = Random.GLOBAL_RNG, duplicates::Bool=false, 
  max_size = 2+2*ceil(Int, hep/(1-hep)))

  if duplicates == false && max_size >= length(elements) 
    @warn "This may run forever because max_size is longer than the total list of elements but we don't allow duplicates"
  end 
  edges = Vector{Vector{Int}}() 
  for i=1:nedges
    push!(edges, random_hyperedge(elements, max_size, hep))
  end 
  return edges
end 

random_hypergraph(n::Int, nedges::Int, hep::Real; kwargs...) = random_hypergraph(1:n, nedges, hep; kwargs... )


"""
  planted_partition(n, m, nedges1, nedges2, hep1, hep2; [rng=Random.GLOBAL_RNG, kwargs...])

Generate a planted partition model with n vertices with a planted partition on the first
m of them (region: 1:m).

We generate nedges1 in the full region and nedges2 in the planted region. 
The edges in the full region have hyperedge length probability hep1 and the 
edges in the planted region have length probablitiy hep2 

See random_hypergraph for other kwargs... 
""" 
function planted_partition(n, m, nedges1, nedges2, hep1, hep2; kwargs...)
  m < n || throw(ArgumentError("m must be strictly less than n"))
  edges1 = random_hypergraph(1:n, nedges1, hep1; kwargs...)
  edges2 = random_hypergraph(1:m, nedges2, hep2; kwargs...)
  edges = append!(edges1, edges2) 
  return unique!(sort!.(edges))
end 


"""
  hypergraph_SBM(n, cluster_vertices, vlabel, nedges1, cluster_edges, hep1, hep2; kwargs...)

Generate a hypergraph using a simplified SBM model.
The graph has n vertices, nedges1 + nedges2 edges, ncluster clusters.
nedges1 edges are edges between different clusters and 
nedges2 edges are edges within the same cluster.

# Arguments
n: number of vertices
ncluster: number of clusters
vlabel: cluster label for each vertex
nedges1: number of edges between clusters
cluster_edges: number of edges within each cluster
hep1: hyperedge extension probability for edges between clusters
hep2: hyperedge extension probability for edges within clusters
"""
function hypergraph_SBM(n, ncluster, vlabel, nedges1, cluster_edges, hep1, hep2; kwargs...) 
  @assert ncluster <= n "number of clusters must not exceed number of vertices"

  # generate edges between clusters
  edges = random_hypergraph(1:n, nedges1, hep1; kwargs...)
  for i = 1:ncluster
    # generate edges within each cluster
    # vlabel is the cluster label for each vertex
    # cluster_edges mark the number of edges within each cluster
    C = findall(x->vlabel[x]==i, 1:n)
    edgesi = random_hypergraph(C, cluster_edges[i], hep2; kwargs...)
    append!(edges, edgesi)
  end
  return unique!(sort!.(edges))
end


"""
  gen_cluster(n, m, ncluster)

Generating the cluster sizes and edge sizes for the hypergraph SBM. 
Splitting n vertices and m hyperedges into ncluster clusters.
Currently the edge number is proportional to the vertex number, which is a bit unrealistic.
"""
function gen_cluster(n, m, ncluster) 
  cluster_vertices = zeros(Int64, ncluster)
  cluster_edges = zeros(Int64, ncluster)
  vlabel = rand(1:ncluster, n)
  elabel = rand(1:ncluster, m)
  for i = 1:ncluster
    cluster_vertices[i] = length(findall(x->vlabel[x]==i, 1:n))
    cluster_edges[i] = length(findall(x->elabel[x]==i, 1:m))
  end
  return cluster_vertices, cluster_edges, vlabel
end

