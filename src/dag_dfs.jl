struct DFSNode
  parent::Vector{DFSNode}
  visited_set_indices::Set{Int64}
  final_node_idx::Int64
  key::Tuple{Set{Int64}, Int64}
end

function DFSNode(parent::Vector{DFSNode}, visited_set_indices::Set{Int64}, final_node_idx::Int64)
  return DFSNode(parent, visited_set_indices, final_node_idx, (visited_set_indices, final_node_idx))
end

function dag_dfs(dist::AbstractArray{Int64, 2}, sets::Vector{Vector{Int64}}, membership::Vector{Int64}, inf_val::Int64, stop_time::Float64)
  bt = time_ns()
  closed_list = Set{Tuple{Set{Int64}, Int64}}()
  dfs_stack = [DFSNode(Vector{DFSNode}(), Set(membership[1]), 1)]
  set_indices = Set(1:length(sets))

  # First dimension is node index, second dimension is set index. inf_cost_to_sets[i, j] = true if the cost from node i to all nodes of set j is infinity
  inf_cost_to_sets = stack([mapslices(minimum, dist[:, set], dims=2) for set=sets], dims=2) .== inf_val
  before = [setdiff(Set(findall(inf_cost_to_sets[node_idx, :])), membership[node_idx]) for node_idx=1:size(dist, 1)]

  num_outgoing_edges = sum(dist .!= inf_val, dims=2)

  solved = false
  cost = inf_val
  while length(dfs_stack) != 0
    if time() > stop_time
      println("Timeout during initial tour generation")
      break
    end

    pop = pop!(dfs_stack)
    if pop.key in closed_list
      continue
    end

    push!(closed_list, pop.key)

    if length(pop.visited_set_indices) == length(sets)
      tour = [pop.final_node_idx]
      node_tmp = pop.parent
      while length(node_tmp) != 0
        node = node_tmp[1]
        push!(tour, node.final_node_idx)
        node_tmp = node.parent
      end
      reverse!(tour)
      at = time_ns()
      println("Found initial tour after ", (at - bt)/1.0e9, " s")
      return tour
    else
      unvisited_set_indices = setdiff(set_indices, pop.visited_set_indices)
      neighbors_mask = zeros(Bool, size(dist, 1))
      for set_idx=unvisited_set_indices
        neighbors_mask[sets[set_idx]] .= true
      end
      neighbors_mask = neighbors_mask .& (dist[pop.final_node_idx, :] .!= inf_val)
      neighbors = findall(neighbors_mask)
    end

    if length(neighbors) == 0
      continue
    end

    sort_idx = sortperm(num_outgoing_edges[neighbors])
    for node_idx=neighbors[sort_idx]
      next_unvisited_set_indices = setdiff(unvisited_set_indices, membership[node_idx])
      # Dumas test 2
      if length(intersect(next_unvisited_set_indices, before[node_idx])) != 0
        continue
      end
      neighbor_node = DFSNode([pop], union(pop.visited_set_indices, membership[node_idx]), node_idx)
      if neighbor_node.key in closed_list
        continue
      end
      push!(dfs_stack, neighbor_node)
    end
  end

  at = time_ns()
  println("Failed to generate initial tour after ", (at - bt)/1.0e9, " s")
  return Vector{Int64}()
end
