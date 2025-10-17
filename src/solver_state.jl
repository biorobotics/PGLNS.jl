
# Is it worth including the phase and temperature in here?
mutable struct SolverState
  inf_val::Int64
  num_vertices::Int64
  num_sets::Int64
  sets::Vector{Vector{Int64}}
  dist::AbstractArray{Int64,2}
  evaluated_edge_mat::AbstractArray{Bool,2}
  stop_at_first_improvement_with_unevaluated_edges::Bool
  membership::Vector{Int64}
  instance_read_time::Float64
  cost_mat_read_time::Float64
  max_threads::Int64
  powers::Dict{String,Any}
  update_powers::Bool
  pin_cores::Vector{Int64}

  setdist_time::Float64
  sets_deepcopy_time::Float64

  param
  count
  setdist::Distsv
  sets_unshuffled::Vector{Vector{Int64}}
end

function initialize_solver_state(args, inf_val::Int64, given_initial_tours::AbstractArray{Int64,1}, dist::AbstractArray{Int64,2}, evaluated_edge_mat::AbstractArray{Bool,2}, stop_at_first_improvement_with_unevaluated_edges::Bool, max_threads::Int64, pin_cores::Vector{Int64}=Vector{Int64}())
  problem_instance, optional_args = parse_cmd(args)
  problem_instance = String(problem_instance)

  read_start_time = time_ns()
  num_vertices, num_sets, sets, tmp_dist, membership = read_file(problem_instance, size(dist, 1) == 0)
  read_end_time = time_ns()
  instance_read_time = (read_end_time - read_start_time)/1.0e9
  if get(optional_args, :verbose, 0) == 3
    println("Reading GTSPLIB file took ", instance_read_time, " s")
  end

  cost_mat_read_time = 0.

  powers = Dict{String, Any}()
  update_powers = false

  if size(dist, 1) == 0
    dist = tmp_dist
  end

  if length(sets[membership[1]]) != 1
    for set_idx=1:length(sets)
      sets[set_idx] .+= 1
    end
    pushfirst!(sets, [1])
    membership = vcat(1, membership .+ 1)
    dist = hcat(zeros(Int64, length(membership)), vcat(0, dist))
    num_sets += 1
    num_vertices += 1
  end

  args = optional_args

  # Random.seed!(1234)

	param = parameter_settings(num_vertices, num_sets, sets, problem_instance, args)
  if length(given_initial_tours) != 0
    @assert(length(given_initial_tours)%num_sets == 0)
    param[:cold_trials] = div(length(given_initial_tours), num_sets)
  end

	init_time = time()

	count = Dict(:latest_improvement => 1,
	  			 :first_improvement => false,
	 		     :warm_trial => 0,
	  		     :cold_trial => 1,
				 :total_iter => 0,
         :print_time => init_time)

  bt = time_ns()
	# compute set distances which will be helpful
	setdist = set_vertex_dist(dist, num_sets, membership)
  at = time_ns()
  setdist_time = (at - bt)/1e9

  if length(powers) == 0
    powers = initialize_powers(param)
  elseif update_powers
    power_update!(powers, param)
  end

  bt = time_ns()
  sets_unshuffled = deepcopy(sets)
  at = time_ns()
  sets_deepcopy_time = (at - bt)/1e9
  # sets_unshuffled = sets # Need to use this to match GLNS

  solver_state = SolverState(inf_val,
                             num_vertices,
                             num_sets,
                             sets,
                             dist,
                             evaluated_edge_mat, 
                             stop_at_first_improvement_with_unevaluated_edges,
                             membership,
                             instance_read_time,
                             cost_mat_read_time,
                             max_threads,
                             powers,
                             update_powers,
                             pin_cores,
                             setdist_time,
                             sets_deepcopy_time,
                             param,
                             count,
                             setdist,
                             sets_unshuffled)
end

function solve_with_state!(solver_state::SolverState, new_time_limit::Float64, updated_edges::AbstractArray{Int64,2}, given_initial_tours::AbstractArray{Int64,1})
  start_time_for_tour_history = time_ns()

  Random.seed!(1234)

  inf_val = solver_state.inf_val
  num_vertices = solver_state.num_vertices
  num_sets = solver_state.num_sets
  sets = solver_state.sets
  # sets = deepcopy(solver_state.sets_unshuffled) # If we want to reproduce the behavior of the version that doesn't use a solver state, we need to use this line instead of the above. From the profiling I did on 10/16/2025, deepcopying the sets isn't expensive, but I didn't see a point in doing it
  dist = solver_state.dist
  evaluated_edge_mat = solver_state.evaluated_edge_mat
  stop_at_first_improvement_with_unevaluated_edges = solver_state.stop_at_first_improvement_with_unevaluated_edges
  membership = solver_state.membership
  max_threads = solver_state.max_threads
  powers = solver_state.powers
  update_powers = solver_state.update_powers
  pin_cores = solver_state.pin_cores
  param = solver_state.param
  count = solver_state.count
  setdist = solver_state.setdist
  sets_unshuffled = solver_state.sets_unshuffled

  nthreads = min(Threads.nthreads(), max_threads)
  if nthreads != 1
    if length(pin_cores) != 0
      pinthreads(pin_cores)
    else
      pinthreads(:cores)
    end
  end

  param[:max_time] = new_time_limit

	#####################################################
	init_time = time()

  # Note: there may be some value in not resetting all these counters after each edge evaluation step. Haven't explored
	count[:latest_improvement] = 1
	count[:first_improvement] = false
	count[:warm_trial] = 0
	count[:cold_trial] = 1
	count[:total_iter] = 0
  count[:print_time] = init_time

	lowest = Tour(Int64[], typemax(Int64))

	start_time = time_ns()
  start_proc_time = CPUtime_us()

  # Assumes that edge costs only increase
  bt = time_ns()
  for edge_idx in 1:size(updated_edges, 1)
    node_idx1 = updated_edges[edge_idx, 1]
    node_idx2 = updated_edges[edge_idx, 2]
    prev_cost = updated_edges[edge_idx, 3]
    if dist[node_idx1, node_idx2] < prev_cost
      println("Edge: ", node_idx1, " ", node_idx2)
      println("Original cost: ", prev_cost)
      println("Updated cost: ", dist[node_idx1, node_idx2])
      throw("Error: updated edge cost is smaller than original cost")
    end

    set1 = membership[node_idx1]
    set2 = membership[node_idx2]

    if setdist.set_vert[set1, node_idx2] == prev_cost
      setdist.set_vert[set1, node_idx2] = typemax(Int64)
      for i in sets[set1]
        if dist[i, node_idx2] < setdist.set_vert[set1, node_idx2]
          setdist.set_vert[set1, node_idx2] = dist[i, node_idx2]
        end
      end
    end

    if setdist.min_sv[set1, node_idx2] == prev_cost
      setdist.min_sv[set1, node_idx2] = typemax(Int64)
      for i in sets[set1]
        if dist[i, node_idx2] < setdist.min_sv[set1, node_idx2]
          setdist.min_sv[set1, node_idx2] = dist[i, node_idx2]
        end
        if dist[node_idx2, i] < setdist.min_sv[set1, node_idx2]
          setdist.min_sv[set1, node_idx2] = dist[node_idx2, i]
        end
      end
    end

    if setdist.vert_set[node_idx2, set1] == prev_cost
      setdist.vert_set[node_idx2, set1] = typemax(Int64)
      for i in sets[set1]
        if dist[node_idx2, i] < setdist.vert_set[node_idx2, set1]
          setdist.vert_set[node_idx2, set1] = dist[node_idx2, i]
        end
      end
    end

    if setdist.set_vert[set2, node_idx1] == prev_cost
      setdist.set_vert[set2, node_idx1] = typemax(Int64)
      for j in sets[set2]
        if dist[j, node_idx1] < setdist.set_vert[set2, node_idx1]
          setdist.set_vert[set2, node_idx1] = dist[j, node_idx1]
        end
      end
    end

    if setdist.min_sv[set2, node_idx1] == prev_cost
      setdist.min_sv[set2, node_idx1] = typemax(Int64)
      for j in sets[set2]
        if dist[j, node_idx1] < setdist.min_sv[set2, node_idx1]
          setdist.min_sv[set2, node_idx1] = dist[j, node_idx1]
        end
        if dist[node_idx1, j] < setdist.min_sv[set2, node_idx1]
          setdist.min_sv[set2, node_idx1] = dist[node_idx1, j]
        end
      end
    end

    if setdist.vert_set[node_idx1, set2] == prev_cost
      setdist.vert_set[node_idx1, set2] = typemax(Int64)
      for j in sets[set2]
        if dist[node_idx1, j] < setdist.vert_set[node_idx1, set2]
          setdist.vert_set[node_idx1, set2] = dist[node_idx1, j]
        end
      end
    end
  end
  at = time_ns()
  update_setdist_time = at - bt
  #=
  setdist_from_scratch = set_vertex_dist(dist, num_sets, membership)
  if !all(setdist_from_scratch.set_vert .== setdist.set_vert)
    throw("Error: setdist from scratch does not match incremental setdist for set_vert")
  end

  if !all(setdist_from_scratch.vert_set .== setdist.vert_set)
    throw("Error: setdist from scratch does not match incremental setdist for vert_set")
  end

  if !all(setdist_from_scratch.min_sv .== setdist.min_sv)
    throw("Error: setdist from scratch does not match incremental setdist for min_sv")
  end
  =#

  
  tour_history = Array{Tuple{Float64, Array{Int64,1}, Int64},1}()
  num_trials_feasible = 0
  num_trials = 0
  # time_per_trial = Vector{Float64}()
  num_trials_lock = ReentrantLock()

  set_locks = [ReentrantLock() for set=sets]

  # rngs = [Future.randjump(Random.default_rng(), thread_idx*big(10)^20) for thread_idx=1:nthreads]

  powers_lock = ReentrantLock()
  best_lock = ReentrantLock()
  current_lock = ReentrantLock()

  phase_lock = ReentrantLock()

  count_lock = ReentrantLock()

  iter_count_lock = ReentrantLock()

  temperature_lock = ReentrantLock()

  lock_times = zeros(nthreads)

  stop_upon_budget = param[:budget] != typemin(Int64)

  time_spent_waiting_for_termination = 0.

	while true
    if count[:cold_trial] > param[:cold_trials] && !stop_upon_budget
      break
    end

    if length(given_initial_tours) != 0
      start_idx = (count[:cold_trial] - 1)*num_sets + 1
      end_idx = count[:cold_trial]*num_sets
      initial_tour = given_initial_tours[start_idx:end_idx]
    else
      initial_tour = Vector{Int64}()
    end
    best = initial_tour!(lowest, dist, sets, setdist, count[:cold_trial], param, num_sets, membership, initial_tour, inf_val, init_time + param[:max_time])
    if length(best.tour) == 0
      return inf_val, powers
    end
    timer = (time_ns() - start_time)/1.0e9
		# print_cold_trial(count, param, best)
		phase = :early

    if count[:cold_trial] > 1 && update_powers
      power_update!(powers, param)
    end

    unevaluated_edge_in_this_cold_trial = false

    while count[:warm_trial] <= param[:warm_trials]
      best_update_time = time_ns()
      iter_count = 1
      current = tour_copy(best)
      temperature = 1.442 * param[:accept_percentage] * best.cost
      # accept a solution with 50% higher cost with 0.05% change after num_iterations.
      cooling_rate = ((0.0005 * lowest.cost)/(param[:accept_percentage] *
                  current.cost))^(1/param[:num_iterations])

      if count[:warm_trial] > 0	  # if warm restart, then use lower temperature
        temperature *= cooling_rate^(param[:num_iterations]/2)
        phase = :late
      end
      this_phase = phase
      thread_broke = false
      @threads for thread_idx=1:nthreads
      # for thread_idx=1:1 # Need to use this to match GLNS
        this_num_trials_feasible = 0
        this_num_trials = 0
        # this_time_per_trial = Vector{Float64}()
        budget_met = false
        unevaluated_edge = false
        while true
          bt_trial_time = time()
          do_break = false
          bt = time()
          lock(count_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            unevaluated_edge_in_this_cold_trial = unevaluated_edge
            if budget_met || unevaluated_edge || thread_broke || 
               (count[:latest_improvement] > (count[:first_improvement] ?
                                              param[:latest_improvement] : param[:first_improvement]))
              thread_broke = true
              do_break = true
            end
          finally
            unlock(count_lock)
          end

          if do_break
            break
          end

          this_iter_count = 0
          bt = time()
          lock(iter_count_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            this_iter_count = iter_count
          finally
            unlock(iter_count_lock)
          end

          bt = time()
          lock(phase_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            if this_iter_count > param[:num_iterations]/2 && phase == :early
              phase = :mid  # move to mid phase after half iterations
            end
            this_phase = phase
          finally
            unlock(phase_lock)
          end
          trial = remove_insert(current, dist, membership, setdist, sets, sets_unshuffled, powers, param, this_phase, powers_lock, current_lock, set_locks, update_powers, lock_times, thread_idx)

          trial_infeasible = dist[trial.tour[end], trial.tour[1]] == inf_val
          @inbounds for i in 1:length(trial.tour)-1
            if trial_infeasible
              break
            end
            trial_infeasible = dist[trial.tour[i], trial.tour[i+1]] == inf_val
          end
          if ~trial_infeasible
            this_num_trials_feasible += 1
          end
          this_num_trials += 1

          # decide whether or not to accept trial
          this_temperature = 0.
          bt = time()
          lock(temperature_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            this_temperature = temperature
          finally
            unlock(temperature_lock)
          end

          accept = false
          if param[:mode] == "slow"
            bt = time()
            lock(current_lock)
            at = time()
            lock_times[thread_idx] += at - bt
            try
              accept = accepttrial_noparam(trial.cost, current.cost, param[:prob_accept]) || accepttrial(trial.cost, current.cost, this_temperature)
            finally
              unlock(current_lock)
            end

            if accept
              opt_cycle!(trial, dist, sets_unshuffled, membership, param, setdist, "full")
              bt = time()
              lock(current_lock)
              at = time()
              lock_times[thread_idx] += at - bt
              try
                current = tour_copy(trial)
              finally
                unlock(current_lock)
              end
            else
              bt = time()
              lock(current_lock)
              at = time()
              lock_times[thread_idx] += at - bt
              try
                trial = tour_copy(current) # I don't remember if there was a point to doing this
              finally
                unlock(current_lock)
              end
            end
          else
            bt = time()
            lock(current_lock)
            at = time()
            lock_times[thread_idx] += at - bt
            try
              if accepttrial_noparam(trial.cost, current.cost, param[:prob_accept]) ||
                 accepttrial(trial.cost, current.cost, this_temperature)
                accept = true
                current = tour_copy(trial)
              else
                trial = tour_copy(current) # I don't remember if there was a point to doing this
              end
            finally
              unlock(current_lock)
            end
          end

          updated_best = false
          bt = time()
          lock(best_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            if accept && trial.cost < best.cost
              best_update_time = time_ns()
              updated_best = true
              best = tour_copy(trial)
              timer = (time_ns() - start_time)/1.0e9
              if param[:print_output] == 3
                println("Thread ", thread_idx, " found new best tour after ", timer, " s with cost ", best.cost, " (before opt cycle)")
              end

              if stop_at_first_improvement_with_unevaluated_edges
                for (node_idx1, node_idx2) in zip(best.tour[1:end-1], best.tour[2:end])
                  if !evaluated_edge_mat[node_idx1, node_idx2]
                    unevaluated_edge = true
                    break
                  end
                end

                node_idx1 = best.tour[end]
                node_idx2 = best.tour[1]
                if !evaluated_edge_mat[node_idx1, node_idx2]
                  unevaluated_edge = true
                end
              end
            end

            budget_met = best.cost <= param[:budget]
          finally
            unlock(best_lock)
          end

          if budget_met
            continue
          end

          if unevaluated_edge
            continue
          end

          bt = time()
          lock(count_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            if updated_best
              count[:latest_improvement] = 1
              count[:first_improvement] = true
              if count[:cold_trial] > 1 && count[:warm_trial] > 1
                count[:warm_trial] = 1
              end
            else
              count[:latest_improvement] += 1
            end
          finally
            unlock(count_lock)
          end

          if updated_best
            opt_cycle!(trial, dist, sets_unshuffled, membership, param, setdist, "full")

            bt = time()
            lock(best_lock)
            at = time()
            lock_times[thread_idx] += at - bt
            try
              # Uncomment || nthreads = 1 to match GLNS
              if trial.cost < best.cost # || nthreads == 1
                best_update_time = time_ns()
                best = tour_copy(trial)
                # print_best(count, param, best, lowest, init_time)
                timer = (time_ns() - start_time)/1.0e9
                if param[:print_output] == 3
                  println("Thread ", thread_idx, " found new best tour after ", timer, " s with cost ", best.cost)
                end

                if param[:output_file] != "None"
                  push!(tour_history, (round((time_ns() - start_time_for_tour_history)/1.0e9, digits=3), best.tour, best.cost))
                end

                if stop_at_first_improvement_with_unevaluated_edges
                  for (node_idx1, node_idx2) in zip(best.tour[1:end-1], best.tour[2:end])
                    if !evaluated_edge_mat[node_idx1, node_idx2]
                      unevaluated_edge = true
                      break
                    end
                  end

                  node_idx1 = best.tour[end]
                  node_idx2 = best.tour[1]
                  if !evaluated_edge_mat[node_idx1, node_idx2]
                    unevaluated_edge = true
                  end
                end
              end
              budget_met = best.cost <= param[:budget]
            finally
              unlock(best_lock)
            end

            if budget_met
              continue
            end

            if unevaluated_edge
              continue
            end

            bt = time()
            lock(current_lock)
            at = time()
            lock_times[thread_idx] += at - bt
            try
              if trial.cost < current.cost
                current = tour_copy(trial)
              end
            finally
              unlock(current_lock)
            end
          end

          bt = time()
          lock(temperature_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            temperature *= cooling_rate  # cool the temperature
          finally
            unlock(temperature_lock)
          end

          bt = time()
          lock(iter_count_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            iter_count += 1
            count[:total_iter] += 1
          finally
            unlock(iter_count_lock)
          end

          at_trial_time = time()
          trial_time = at_trial_time - bt_trial_time
          # push!(this_time_per_trial, trial_time)

          if time() - init_time > param[:max_time]
            break
          end
        end
        bt = time()
        lock(num_trials_lock)
        at = time()
        lock_times[thread_idx] += at - bt
        try
          num_trials += this_num_trials
          num_trials_feasible += this_num_trials_feasible
          # append!(time_per_trial, this_time_per_trial)
        finally
          unlock(num_trials_lock)
        end
      end
      print_warm_trial(count, param, best, iter_count)
      count[:warm_trial] += 1
      count[:latest_improvement] = 1
      count[:first_improvement] = false
      param[:budget_met] = best.cost <= param[:budget]
      if param[:budget_met]
        break
      end

      if unevaluated_edge_in_this_cold_trial
        break
      end

      time_spent_waiting_for_termination += (time_ns() - best_update_time)/1e9
      if time() - init_time > param[:max_time]
        param[:timeout] = true
        break
      end
    end
		lowest.cost > best.cost && (lowest = best)
		count[:warm_trial] = 0
		count[:cold_trial] += 1
    if param[:budget_met]
      break
    end
    if unevaluated_edge_in_this_cold_trial
      break
    end
    if time() - init_time > param[:max_time]
      break
    end
	end
	timer = (time_ns() - start_time)/1.0e9
  if param[:output_file] != "None"
    push!(tour_history, (round((time_ns() - start_time_for_tour_history)/1.0e9, digits=3), lowest.tour, lowest.cost))
  end

  proc_timer = (CPUtime_us() - start_proc_time)/1e6

  @assert(lowest.cost == tour_cost(lowest.tour, dist))
  @assert(length(lowest.tour) == num_sets)
  if param[:print_output] == 3
    println(lock_times)
  end

  return lowest.tour, timer, proc_timer, num_trials_feasible, num_trials, param[:timeout], lock_times, time_spent_waiting_for_termination, tour_history, update_setdist_time
end
