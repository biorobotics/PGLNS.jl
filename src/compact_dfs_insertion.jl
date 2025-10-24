using DataStructures
using FastPriorityQueues

include("compact_astar_insertion.jl")

struct VDNodeDFS
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
  cand_successors::Vector{Tuple{Int64,Int64}} # node_idx, and index in sets_to_insert it corresponds to
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
  next_nonremoved_set_idx = -1
  if next_nonremoved_idx <= length(partial_tour)
    push!(unvisited_removed_sets, -1)
    next_nonremoved_set_idx = membership[partial_tour[next_nonremoved_idx]]
  end

  cand_successors = Vector{Tuple{Int64, Int64}}()

  for removed_set_idx in unvisited_removed_sets
    # bt = time_ns()
    set_idx = removed_set_idx == -1 ? next_nonremoved_set_idx : sets_to_insert[removed_set_idx]

    if next_nonremoved_set_idx != -1 && vd_info.after[partial_tour[next_nonremoved_idx], set_idx]
      # at = time_ns()
      # vd_info.inf_and_prune_check_time += (at - bt)/1e9
      continue
    end
    # at = time_ns()
    # vd_info.inf_and_prune_check_time += (at - bt)/1e9

    this_set = removed_set_idx == -1 ? [partial_tour[next_nonremoved_idx]] : sets[set_idx]
    for node_idx in this_set
      # TODO: seems like I should at least do the dist pruning here
      push!(cand_successors, (node_idx, removed_set_idx))
    end
  end

  shuffle!(cand_successors)

  key = (tour_idx, visited_removed_sets, final_node_idx)
  return VDNodeDFS(tour_idx, parent, visited_removed_sets, final_node_idx, key, g_val, h_val, g_val + h_val, unvisited_removed_sets, num_nonremoved_visited, next_nonremoved_idx, cand_successors)
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

  goal_node_vec = Vector{VDNodeDFS}()
  while length(open_list) != 0 && time() < stop_time
    # bt = time_ns()

    top = open_list[end]

    if length(top.cand_successors) == 0
      pop!(open_list)
      continue
    end

    # at = time_ns()
    # vd_info.open_pop_time += (at - bt)/1e9

    (node_idx, removed_set_idx) = pop!(top.cand_successors)

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
