using DataStructures
using FastPriorityQueues

include("compact_dp_insertion.jl")

# VD is visitation DAG
struct VDNodeAstar
  tour_idx::Int64
  parent::Vector{VDNodeAstar}
  visited_removed_sets::Vector{Bool}
  final_node_idx::Int64
  key::Tuple{Int64, Vector{Bool}, Int64}
  g_val::Int64
  h_val::Int64
  f_val::Int64
end

function VDNodeAstar(tour_idx::Int64, parent::Vector{VDNodeAstar}, visited_removed_sets::Vector{Bool}, final_node_idx::Int64, g_val::Int64, h_val::Int64)
  # return VDNodeAstar(tour_idx, parent, visited_removed_sets, final_node_idx, (tour_idx, visited_removed_sets, final_node_idx), g_val, h_val, g_val + 10*h_val)
  return VDNodeAstar(tour_idx, parent, visited_removed_sets, final_node_idx, (tour_idx, visited_removed_sets, final_node_idx), g_val, h_val, g_val + h_val)
end

# Note: this whole function assumes open tsp, doesn't account for cost of returning to depot.
# Not a fundamental limitation, I just didn't implement it to handle closed TSP
function astar_insertion!(sets_to_insert::Vector{Int64}, dist::AbstractArray{Int64, 2}, sets::Vector{Vector{Int64}}, membership::Vector{Int64}, inf_val::Int64, stop_time::Float64, vd_info::VDInfo, partial_tour::Vector{Int64}, ub::Int64)
  # When the next unvisted index in partial_tour is tour_idx, the h value is dist[node_idx, partial_tour[tour_idx]] + h_vals[tour_idx]
  h_vals = Vector{Int64}(undef, length(partial_tour))
  h_vals[length(partial_tour)] = 0
  for tour_idx=length(partial_tour)-1:-1:1
    h_vals[tour_idx] = dist[partial_tour[tour_idx], partial_tour[tour_idx + 1]] + h_vals[tour_idx + 1]
  end

  # num_generated_nodes = 0
  # num_expanded_nodes = 0

  #=
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
  =#

  removed_set_indices_per_tour_idx = compute_removed_set_indices_per_tour_idx(partial_tour, sets_to_insert, vd_info)

  closed_list = Set{Tuple{Int64, Vector{Bool}, Int64}}()

  open_list = HeapPriorityQueue{VDNodeAstar, Int64}()
  root_node = VDNodeAstar(1, Vector{VDNodeAstar}(), zeros(Bool, length(sets_to_insert)), 1, 0, 0)
  enqueue!(open_list, root_node, 0)
  # num_generated_nodes += 1
  # vd_info.num_nodes_pushed_to_open += 1

  # numbers_of_successors = Vector{Int64}()

  seen_nodes = Dict{Tuple{Int64, Vector{Bool}, Int64}, VDNodeAstar}()
  seen_nodes[root_node.key] = root_node

  #=
  known_feas_tour = [1, 223, 393, 240, 155, 2, 36, 138, 172, 121, 19, 104, 189, 70, 87, 53, 257, 274, 291, 206, 376, 512, 631, 427, 529, 546, 308, 699, 359, 342, 563, 325, 614, 682, 410, 495, 478, 580, 444, 461, 767, 665, 869, 818, 1073, 597, 886, 648, 716, 937, 1022, 784, 1005, 750, 903, 1039, 733, 852, 1124, 835, 801, 971, 1141, 988, 1090, 954, 1294, 920, 1362, 1107, 1175, 1158, 1056, 1282, 1345, 1498, 1549, 1192, 1311, 1209, 1464, 1260, 1396, 1243, 1226, 1413, 1481, 1328, 1430, 1447, 1685, 1395, 1583, 1668, 1532, 1566, 1634, 1736, 1651, 1617, 1515, 1806, 1872, 1753, 1838, 1702, 1600, 1906, 1855, 1719, 2093, 1770, 1796, 1823, 1889, 1974, 1991, 2008, 2161, 2127, 1940, 2399, 2144, 1957, 1923, 2229, 2450, 2042, 2059, 2025, 2518, 2076, 2484, 2178, 2467, 2297, 2110, 2705, 2246, 2263, 2195, 2382, 2212, 2348, 2365, 2314, 2773, 2331, 2535, 2824, 2280, 2433, 2586, 2756, 2552, 2501, 2722, 2671, 2858, 2790, 2416, 2688, 2603, 2620, 3062, 2569, 2977, 2739, 2926, 2960, 2654, 2637, 3096, 2841, 2815, 2920, 3045, 3079, 2943, 3130, 2875, 3198, 2892, 3011, 3113, 3028, 3283, 3215, 2994, 3334, 3300, 3249, 3147, 3164, 3232, 3266, 3181, 3317, 3351, 3385, 3368]
  known_feas_vd_node_seq = [root_node]
  for (tour_idx, node_idx) in enumerate(known_feas_tour)
    if tour_idx == 1
      continue
    end
    visited_removed_sets = copy(known_feas_vd_node_seq[end].visited_removed_sets)
    removed_set_idx = findfirst(==(membership[node_idx]), sets_to_insert)
    if removed_set_idx != nothing
      visited_removed_sets[removed_set_idx] = true
      if !(removed_set_idx in removed_set_indices_per_tour_idx[tour_idx])
        set_idx = sets_to_insert[removed_set_idx]
        println("tour_idx = ", tour_idx)
        println("removed_set_idx = ", removed_set_idx)
        println("set_idx = ", set_idx)
        println("number of removed sets = ", length(sets_to_insert))

        for node_idx2 in known_feas_tour[tour_idx + 1:end]
          if vd_info.after[node_idx2, set_idx]
            throw("Error: after[node_idx2, set_idx] but set_idx comes before")
          end

          if vd_info.before_set_to_set[set_idx, membership[node_idx2]]
            throw("Error: before_set_to_set[set_idx, membership[node_idx2]] but set_idx comes before")
          end
        end

        for node_idx2 in known_feas_tour[1:tour_idx-1]
          if vd_info.before[node_idx2, set_idx]
            throw("Error: before[node_idx2, set_idx] but set_idx comes after")
          end

          if vd_info.before_set_to_set[membership[node_idx2], set_idx]
            throw("Error: before_set_to_set[membership[node_idx2], set_idx] but set_idx comes after")
          end
        end
        throw("Error: cannot reconstruct known feas tour")
      end
    end
    vd_node = VDNodeAstar(tour_idx, [known_feas_vd_node_seq[end]], visited_removed_sets, node_idx, 0, 0)
  end
  =#

  # num_cand_successors_per_expansion = Vector{Int64}()

  goal_node = VDNodeAstar(0, Vector{VDNodeAstar}(), zeros(Bool, 1), 1, typemax(Int64), 0)
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

    # num_expanded_nodes += 1

    # at = time_ns()
    # vd_info.closed_check_time += (at - bt)/1e9


    # bt = time_ns()

    push!(closed_list, pop.key)

    # push!(numbers_of_successors, 0)

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

    cand_successors = Vector{Tuple{Int64,Int64}}()
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
      cand_successors = cat(cand_successors, [(node_idx, removed_set_idx) for node_idx in this_set], dims=1)
    end

    # push!(num_cand_successors_per_expansion, length(cand_successors))

    num_threads = Threads.nthreads()
    neighbor_nodes_per_thread = Vector{Vector{VDNodeAstar}}(undef, num_threads)
    for thread_idx=1:num_threads
      neighbor_nodes_per_thread[thread_idx] = Vector{VDNodeAstar}()
    end
    # bt = time_ns()
    Threads.@threads :static for (node_idx, removed_set_idx) in cand_successors
      # bt = time_ns()
      h_val = num_nonremoved_visited == length(partial_tour) ? 0 : dist[node_idx, partial_tour[next_nonremoved_idx]] + h_vals[next_nonremoved_idx]
      if pop.g_val + dist[pop.final_node_idx, node_idx] + h_val >= ub
        # at = time_ns()
        # vd_info.inf_and_prune_check_time += (at - bt)/1e9
        continue
      end

      if dist[pop.final_node_idx, node_idx] == inf_val
        # at = time_ns()
        # vd_info.inf_and_prune_check_time += (at - bt)/1e9
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
        # at = time_ns()
        # vd_info.inf_and_prune_check_time += (at - bt)/1e9
        continue
      end

      # Check if unvisited nonremoved node is unreachable from node_idx. If so, prune
      if next_nonremoved_idx <= length(partial_tour) && node_idx != partial_tour[next_nonremoved_idx] && dist[node_idx, partial_tour[next_nonremoved_idx]] == inf_val
        # at = time_ns()
        # vd_info.inf_and_prune_check_time += (at - bt)/1e9
        continue
      end

      # at = time_ns()
      # vd_info.inf_and_prune_check_time += (at - bt)/1e9

      # bt = time_ns()
      neighbor_node = VDNodeAstar(pop.tour_idx + 1, [pop], copy(pop.visited_removed_sets), node_idx, pop.g_val + dist[pop.final_node_idx, node_idx], h_val)
      if removed_set_idx != -1
        neighbor_node.visited_removed_sets[removed_set_idx] = true
        neighbor_node.key[2][removed_set_idx] = true
      end
      # at = time_ns()
      # vd_info.succ_gen_time += (at - bt)/1e9

      # bt = time_ns()
      if neighbor_node.key in closed_list
        # at = time_ns()
        # vd_info.succ_closed_time += (at - bt)/1e9
        continue
      end
      # at = time_ns()
      # vd_info.succ_closed_time += (at - bt)/1e9

      # bt = time_ns()
      if haskey(seen_nodes,  neighbor_node.key) && seen_nodes[neighbor_node.key].g_val <= neighbor_node.g_val
        # at = time_ns()
        # vd_info.seen_key_time += (at - bt)/1e9
        continue
      end
      # at = time_ns()
      # vd_info.seen_key_time += (at - bt)/1e9
      push!(neighbor_nodes_per_thread[Threads.threadid()], neighbor_node)
    end
    # at = time_ns()
    # This now also includes succ_gen_time, succ_closed_time, and seen_key_time
    # vd_info.inf_and_prune_check_time += (at - bt)/1e9

    for neighbor_node in cat(neighbor_nodes_per_thread..., dims=1)
      # bt = time_ns()
      seen_nodes[neighbor_node.key] = neighbor_node
      # at = time_ns()
      # vd_info.seen_update_time += (at - bt)/1e9

      # bt = time_ns()
      enqueue!(open_list, neighbor_node, neighbor_node.f_val)
      # numbers_of_successors[end] += 1
      # num_generated_nodes += 1
      # at = time_ns()
      # vd_info.open_push_time += (at - bt)/1e9
      # vd_info.num_nodes_pushed_to_open += 1
      # if neighbor_node.g_val + h_val == pop.f_val
      #   vd_info.num_nodes_pushed_to_open_with_same_f_val_as_parent += 1
      # end

      # bt = time_ns()
      if neighbor_node.f_val < goal_node.f_val && neighbor_node.tour_idx == length(sets)
        goal_node = neighbor_node
      end
      # at = time_ns()
      # vd_info.goal_check_time += (at - bt)/1e9
    end
  end
  # println("max tour idx: ", max_tour_idx)
  # println("Number of expanded nodes = ", num_expanded_nodes, ", number of generated nodes = ", num_generated_nodes, ", ratio = ", num_expanded_nodes/num_generated_nodes)
  # println("min/median/max number of successors = ", minimum(numbers_of_successors), " ", median(numbers_of_successors), " ", maximum(numbers_of_successors))

  # println("min/median/max number of candidate successors per expansion = ", minimum(num_cand_successors_per_expansion), " ", median(num_cand_successors_per_expansion), " ", maximum(num_cand_successors_per_expansion))

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
