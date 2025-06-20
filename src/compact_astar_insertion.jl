using DataStructures
using FastPriorityQueues

# VD is visitation DAG
struct VDNode
  tour_idx::Int64
  parent::Vector{VDNode}
  visited_removed_sets::Vector{Bool}
  final_node_idx::Int64
  key::Tuple{Int64, Vector{Bool}, Int64}
  g_val::Int64
  h_val::Int64
  f_val::Int64
end

function VDNode(tour_idx::Int64, parent::Vector{VDNode}, visited_removed_sets::Vector{Bool}, final_node_idx::Int64, g_val::Int64, h_val::Int64)
  return VDNode(tour_idx, parent, visited_removed_sets, final_node_idx, (tour_idx, visited_removed_sets, final_node_idx), g_val, h_val, g_val + h_val)
end

struct VDInfo
  before::Array{Bool, 2}
  before_set_to_set::Array{Bool, 2}
  before_time::Float64
  before_update_time::Float64
  open_pop_time::Float64
  closed_check_time::Float64
  closed_push_time::Float64
  inf_and_prune_check_time::Float64
  succ_gen_time::Float64
  succ_closed_time::Float64
  seen_key_time::Float64
  seen_update_time::Float64
  open_push_time::Float64
  goal_check_time::Float64
end

function VDInfo(dist::AbstractArray{Int64, 2}, sets::Vector{Vector{Int64}}, membership::Vector{Int64}, inf_val::Int64)
  vd_info = VDInfo(ones(Bool, length(membership), length(sets)), zeros(Bool, length(sets), length(sets)), zeros(12)...)
  if length(sets) == 0
    return vd_info
  end

  # bt = time_ns()
  for set_idx=1:length(sets)
    for node_idx=1:length(membership)
      for node_idx2 in sets[set_idx]
        if dist[node_idx, node_idx2] != inf_val
          vd_info.before[node_idx, set_idx] = false
          break
        end
      end
    end
  end

  for node_idx = 2:length(membership)
    vd_info.before[node_idx, membership[node_idx]] = false
  end

  for set_idx1=1:length(sets)
    for set_idx2=1:length(sets)
      if set_idx1 == set_idx2
        continue
      end
      vd_info.before_set_to_set[set_idx1, set_idx2] = all(vd_info.before[sets[set_idx1], set_idx2])
    end
  end

  # at = time_ns()
  # vd_info.before_time += (at - bt)/1e9

  return vd_info
end

# Note: this whole function assumes open tsp, doesn't account for cost of returning to depot.
# Not a fundamental limitation, I just didn't implement it to handle closed TSP
function dp_insertion!(sets_to_insert::Vector{Int64}, dist::AbstractArray{Int64, 2}, sets::Vector{Vector{Int64}}, membership::Vector{Int64}, inf_val::Int64, stop_time::Float64, vd_info::VDInfo, partial_tour::Vector{Int64}, ub::Int64)
  # When the next unvisted index in partial_tour is tour_idx, the h value is dist[node_idx, partial_tour[tour_idx]] + h_vals[tour_idx]
  h_vals = Vector{Int64}(undef, length(partial_tour))
  h_vals[length(partial_tour)] = 0
  for tour_idx=length(partial_tour)-1:-1:1
    h_vals[tour_idx] = dist[partial_tour[tour_idx], partial_tour[tour_idx + 1]] + h_vals[tour_idx + 1]
  end

  removed_set_indices_per_tour_idx = Vector{Vector{Int64}}(undef, length(sets))
  for tour_idx=1:length(sets)
    removed_set_indices_per_tour_idx[tour_idx] = Vector{Int64}()
  end

  for (removed_set_idx, set_idx) in enumerate(sets_to_insert)
    Pmin = sum(vd_info.before_set_to_set[set_idx, :]) + 1
    Pmax = length(sets) - sum(vd_info.before_set_to_set[:, set_idx])
    for tour_idx=Pmin:Pmax
      push!(removed_set_indices_per_tour_idx[tour_idx], removed_set_idx)
    end
  end

  closed_list = Set{Tuple{Int64, Vector{Bool}, Int64}}()

  open_list = HeapPriorityQueue{VDNode, Int64}()
  root_node = VDNode(1, Vector{VDNode}(), zeros(Bool, length(sets_to_insert)), 1, 0, 0)
  enqueue!(open_list, root_node, 0)

  seen_nodes = Dict{Tuple{Int64, Vector{Bool}, Int64}, VDNode}()
  seen_nodes[root_node.key] = root_node

  goal_node = VDNode(0, Vector{VDNode}(), zeros(Bool, 1), 1, typemax(Int64), 0)
  while length(open_list) != 0 && goal_node.f_val > peek(open_list).first.f_val
    if time() >= stop_time
      println("Timeout during A*")
      return Vector{Int64}()
    end

    # bt = time_ns()

    pop = dequeue!(open_list)

    # at = time_ns()
    # vd_info.open_pop_time += (at - bt)/1e9

    # bt = time_ns()
    if pop.key in closed_list
      # at = time_ns()
      # vd_info.closed_check_time += (at - bt)/1e9
      continue
    end

    # at = time_ns()
    # vd_info.closed_check_time += (at - bt)/1e9


    # bt = time_ns()

    push!(closed_list, pop.key)

    # at = time_ns()
    # vd_info.closed_push_time += (at - bt)/1e9

    unvisited_removed_sets = [removed_set_idx for removed_set_idx in removed_set_indices_per_tour_idx[pop.tour_idx + 1] if !pop.visited_removed_sets[removed_set_idx]]

    # Handle next non-removed set
    num_nonremoved_visited = pop.tour_idx - sum(pop.visited_removed_sets)
    next_nonremoved_idx = num_nonremoved_visited + 1
    next_nonremoved_set_idx = -1
    if next_nonremoved_idx <= length(partial_tour)
      push!(unvisited_removed_sets, -1)
      next_nonremoved_set_idx = membership[partial_tour[next_nonremoved_idx]]
    end

    for removed_set_idx in unvisited_removed_sets
      # bt = time_ns()
      set_idx = removed_set_idx == -1 ? next_nonremoved_set_idx : sets_to_insert[removed_set_idx]

      if next_nonremoved_set_idx != -1 && vd_info.before_set_to_set[set_idx, next_nonremoved_set_idx]
        continue
      end

      this_set = removed_set_idx == -1 ? [partial_tour[next_nonremoved_idx]] : sets[set_idx]
      for node_idx in this_set
        h_val = num_nonremoved_visited == length(partial_tour) ? 0 : dist[node_idx, partial_tour[next_nonremoved_idx]] + h_vals[next_nonremoved_idx]
        if pop.g_val + dist[pop.final_node_idx, node_idx] + h_val >= ub
          continue
        end

        if dist[pop.final_node_idx, node_idx] == inf_val
          continue
        end

        # Check if unvisited removed customer is in BEFORE[node_idx]. If so, prune
        prune = false
        for removed_set_idx2 in unvisited_removed_sets
          if removed_set_idx2 != -1 && vd_info.before[node_idx, sets_to_insert[removed_set_idx2]]
            prune = true
            break
          end
        end
        if prune
          # vd_info.inf_and_prune_check_time += (at - bt)/1e9
          continue
        end

        # Check if unvisited nonremoved node is unreachable from node_idx. If so, prune
        if next_nonremoved_idx <= length(partial_tour) && node_idx != partial_tour[next_nonremoved_idx] && dist[node_idx, partial_tour[next_nonremoved_idx]] == inf_val
          # vd_info.inf_and_prune_check_time += (at - bt)/1e9
          continue
        end

        # at = time_ns()
        # vd_info.inf_and_prune_check_time += (at - bt)/1e9

        # bt = time_ns()
        neighbor_node = VDNode(pop.tour_idx + 1, [pop], copy(pop.visited_removed_sets), node_idx, pop.g_val + dist[pop.final_node_idx, node_idx], h_val)
        if removed_set_idx != -1
          neighbor_node.visited_removed_sets[removed_set_idx] = true
          neighbor_node.key[2][removed_set_idx] = true
        end

        if neighbor_node.key in closed_list
          continue
        end

        if haskey(seen_nodes,  neighbor_node.key) && seen_nodes[neighbor_node.key].g_val <= neighbor_node.g_val
          continue
        end

        seen_nodes[neighbor_node.key] = neighbor_node

        enqueue!(open_list, neighbor_node, neighbor_node.f_val)

        if neighbor_node.f_val < goal_node.f_val && neighbor_node.tour_idx == length(sets)
          goal_node = neighbor_node
        end
      end
    end
  end
  # println("max tour idx: ", max_tour_idx)

  # No solution
  if length(open_list) == 0
    # println("No solution found by A*")
    if goal_node.f_val != typemax(Int64)
      throw("f_value of goal node should be typemax(Int64)")
    end
    return Vector{Int64}()
  end

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
