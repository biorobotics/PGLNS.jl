using DataStructures
using FastPriorityQueues

include("compact_astar_insertion.jl")

mutable struct VDNodeDFS
  tour_idx::Int64
  parent::Vector{VDNodeDFS}
  visited_removed_sets::Vector{Bool}
  final_node_idx::Int64
  key::Tuple{Int64, Vector{Bool}, Int64}
  g_val::Int64
  h_val::Int64
  f_val::Int64
  unvisited_removed_sets::Vector{Int64}
  num_nonremoved_visited::Int64
  next_nonremoved_idx::Int64
  next_cand_successor_node::Int64
end

function VDNodeDFS(tour_idx::Int64, parent::Vector{VDNodeDFS}, visited_removed_sets::Vector{Bool}, final_node_idx::Int64, g_val::Int64, h_val::Int64, removed_set_indices_per_tour_idx::Vector{Vector{Int64}}, partial_tour::Vector{Int64}, membership::Vector{Int64}, sets_to_insert::Vector{Int64}, vd_info::VDInfo, sets::Vector{Vector{Int64}})
  if tour_idx == length(sets)
    unvisited_removed_sets = Vector{Int64}()
  else
    unvisited_removed_sets = [removed_set_idx for removed_set_idx in removed_set_indices_per_tour_idx[tour_idx + 1] if !visited_removed_sets[removed_set_idx]]
  end

  # Handle next non-removed set
  num_nonremoved_visited = tour_idx - sum(visited_removed_sets)
  next_nonremoved_idx = num_nonremoved_visited + 1
  if next_nonremoved_idx <= length(partial_tour)
    push!(unvisited_removed_sets, -1)
  end

  next_cand_successor_node = 1

  key = (tour_idx, visited_removed_sets, final_node_idx)
  return VDNodeDFS(tour_idx, parent, visited_removed_sets, final_node_idx, key, g_val, h_val, g_val + h_val, unvisited_removed_sets, num_nonremoved_visited, next_nonremoved_idx, next_cand_successor_node) 
end

# Note: this whole function assumes open tsp, doesn't account for cost of returning to depot.
# Not a fundamental limitation, I just didn't implement it to handle closed TSP
function dfs_insertion!(sets_to_insert::Vector{Int64}, dist::AbstractArray{Int64, 2}, sets::Vector{Vector{Int64}}, membership::Vector{Int64}, inf_val::Int64, stop_time::Float64, vd_info::VDInfo, partial_tour::Vector{Int64}, ub::Int64)
  # When the next unvisted index in partial_tour is tour_idx, the h value is dist[node_idx, partial_tour[tour_idx]] + h_vals[tour_idx]
  h_vals = Vector{Int64}(undef, length(partial_tour))
  h_vals[length(partial_tour)] = 0
  for tour_idx=length(partial_tour)-1:-1:1
    h_vals[tour_idx] = dist[partial_tour[tour_idx], partial_tour[tour_idx + 1]] + h_vals[tour_idx + 1]
  end

  removed_set_indices_per_tour_idx = compute_removed_set_indices_per_tour_idx(partial_tour, sets_to_insert, vd_info)

  open_list = Vector{VDNodeDFS}()

  # Generate root node and push to open
  tour_idx = 1
  visited_removed_sets = zeros(Bool, length(sets_to_insert))
  node_idx = 1
  key = (tour_idx, visited_removed_sets, node_idx)
  g_val = 0
  h_val = 0
  root_node = VDNodeDFS(tour_idx, Vector{VDNodeDFS}(), visited_removed_sets, node_idx, g_val, h_val, removed_set_indices_per_tour_idx, partial_tour, membership, sets_to_insert, vd_info, sets)
  push!(open_list, root_node)

  # vd_info.num_nodes_pushed_to_open += 1

  seen_nodes = Dict{Tuple{Int64, Vector{Bool}, Int64}, VDNodeDFS}()
  seen_nodes[root_node.key] = root_node

  bt = time_ns()
  removed_sets_mask = zeros(Bool, length(sets))
  removed_sets_mask[sets_to_insert] .= true
  removed_set_indices = cumsum(removed_sets_mask)
  cand_successor_nodes = cat(partial_tour, sets[sets_to_insert]..., dims=1)
  shuffle!(cand_successor_nodes)
  at = time_ns()
  vd_info.prealloc_time += (at - bt)/1e9

  goal_node_vec = Vector{VDNodeDFS}()
  while length(open_list) != 0 && time() < stop_time
    # bt = time_ns()

    top = open_list[end]

    if top.next_cand_successor_node == length(cand_successor_nodes) + 1
      pop!(open_list)
      continue
    end

    # at = time_ns()
    # vd_info.open_pop_time += (at - bt)/1e9

    node_idx = cand_successor_nodes[top.next_cand_successor_node]
    top.next_cand_successor_node += 1
    if top.next_nonremoved_idx != length(partial_tour) + 1 && node_idx == partial_tour[top.next_nonremoved_idx]
      removed_set_idx = -1
    elseif removed_sets_mask[membership[node_idx]]
      set_idx = membership[node_idx]
      removed_set_idx = removed_set_indices[membership[node_idx]]
      if top.visited_removed_sets[removed_set_idx] || (top.next_nonremoved_idx != length(partial_tour) + 1 && vd_info.after[partial_tour[top.next_nonremoved_idx], set_idx])
        continue
      end
    else
      continue
    end

    # bt = time_ns()
    h_val = top.num_nonremoved_visited == length(partial_tour) ? 0 : dist[node_idx, partial_tour[top.next_nonremoved_idx]] + h_vals[top.next_nonremoved_idx]
    if top.g_val + dist[top.final_node_idx, node_idx] + h_val >= ub
      # at = time_ns()
      # vd_info.inf_and_prune_check_time += (at - bt)/1e9
      continue
    end

    if dist[top.final_node_idx, node_idx] == inf_val
      # at = time_ns()
      # vd_info.inf_and_prune_check_time += (at - bt)/1e9
      continue
    end

    # Check if unvisited removed customer is in BEFORE[node_idx]. If so, prune
    prune = false
    for removed_set_idx2 in top.unvisited_removed_sets
      if removed_set_idx2 != -1 && vd_info.before[node_idx, sets_to_insert[removed_set_idx2]]
        prune = true
        break
      end
    end
    if prune
      # at = time_ns()
      # vd_info.inf_and_prune_check_time += (at - bt)/1e9
      continue
    end

    # Check if unvisited nonremoved node is unreachable from node_idx. If so, prune
    if top.next_nonremoved_idx <= length(partial_tour) && node_idx != partial_tour[top.next_nonremoved_idx] && dist[node_idx, partial_tour[top.next_nonremoved_idx]] == inf_val
      # at = time_ns()
      # vd_info.inf_and_prune_check_time += (at - bt)/1e9
      continue
    end

    # at = time_ns()
    # vd_info.inf_and_prune_check_time += (at - bt)/1e9

    tour_idx = top.tour_idx + 1
    visited_removed_sets = copy(top.visited_removed_sets)
    if removed_set_idx != -1
      visited_removed_sets[removed_set_idx] = true
    end
    neighbor_key = (tour_idx, visited_removed_sets, node_idx)
    g_val = top.g_val + dist[top.final_node_idx, node_idx]

    # bt = time_ns()
    if haskey(seen_nodes, neighbor_key) && seen_nodes[neighbor_key].g_val <= g_val
      # at = time_ns()
      # vd_info.seen_key_time += (at - bt)/1e9
      continue
    end
    # at = time_ns()
    # vd_info.seen_key_time += (at - bt)/1e9

    # bt = time_ns()
    neighbor_node = VDNodeDFS(tour_idx, [top], visited_removed_sets, node_idx, g_val, h_val, removed_set_indices_per_tour_idx, partial_tour, membership, sets_to_insert, vd_info, sets)
    # at = time_ns()
    # vd_info.succ_gen_time += (at - bt)/1e9

    # bt = time_ns()
    seen_nodes[neighbor_node.key] = neighbor_node
    # at = time_ns()
    # vd_info.seen_update_time += (at - bt)/1e9

    # bt = time_ns()
    push!(open_list, neighbor_node)
    # at = time_ns()
    # vd_info.open_push_time += (at - bt)/1e9
    # vd_info.num_nodes_pushed_to_open += 1
    # if neighbor_node.g_val + h_val == pop.f_val
    #   vd_info.num_nodes_pushed_to_open_with_same_f_val_as_parent += 1
    # end

    # bt = time_ns()
    if neighbor_node.tour_idx == length(sets)
      if length(goal_node_vec) == 0
        push!(goal_node_vec, neighbor_node)
      elseif neighbor_node.f_val < goal_node_vec[1].f_val
        goal_node_vec[1] = neighbor_node
        ub = neighbor_node.f_val
      end
    end
    # at = time_ns()
    # vd_info.goal_check_time += (at - bt)/1e9
  end
  # println("max tour idx: ", max_tour_idx)

  # No solution
  if length(goal_node_vec) == 0
    return Vector{Int64}()
  end

  goal_node = goal_node_vec[1]

  # println("dfs goal cost = ", goal_node.f_val)
  # throw("intentional error")

  tour = [goal_node.final_node_idx]
  node_tmp = goal_node.parent
  while length(node_tmp) != 0
    node = node_tmp[1]
    push!(tour, node.final_node_idx)
    node_tmp = node.parent
  end
  reverse!(tour)

  return tour
end
