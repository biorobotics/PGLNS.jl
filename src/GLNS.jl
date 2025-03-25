# Copyright 2017 Stephen L. Smith and Frank Imeson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
module GLNS
export solver
using Random
using Sockets
using Printf
using NPZ
using Polyester: @batch
using Base.Threads
using ThreadPinning
import Future
include("utilities.jl")
include("parse_print.jl")
include("tour_optimizations.jl")
include("adaptive_powers.jl")
include("insertion_deletion.jl")
include("parameter_defaults.jl")

function solver(problem_instance::String, given_initial_tours::AbstractArray{Int64,1}, start_time_for_tour_history::UInt64, inf_val::Int64, num_vertices::Int64, num_sets::Int64, sets::Vector{Vector{Int64}}, dist::AbstractArray{Int64,2}, membership::Vector{Int64}, instance_read_time::Float64, cost_mat_read_time::Float64, max_threads::Int64, powers::Dict{String,Any}=Dict{String,Any}(), update_powers::Bool=true, pin_cores::Vector{Int64}=Vector{Int64}(); args...)
  Random.seed!(1234)

  nthreads = min(Threads.nthreads(), max_threads)
  if nthreads != 1
    if length(pin_cores) != 0
      pinthreads(pin_cores)
    else
      pinthreads(:cores)
    end
  end

	param = parameter_settings(num_vertices, num_sets, sets, problem_instance, args)
  if length(given_initial_tours) != 0
    @assert(length(given_initial_tours)%num_sets == 0)
    param[:cold_trials] = div(length(given_initial_tours), num_sets)
  end
	#####################################################
	init_time = time()

	count = Dict(:latest_improvement => 1,
	  			 :first_improvement => false,
	 		     :warm_trial => 0,
	  		     :cold_trial => 1,
				 :total_iter => 0,
				 :print_time => init_time)
	lowest = Tour(Int64[], typemax(Int64))

	start_time = time_ns()
	# compute set distances which will be helpful
	setdist = set_vertex_dist(dist, num_sets, membership)
  if length(powers) == 0
    powers = initialize_powers(param)
  elseif update_powers
    power_update!(powers, param)
  end
  
  sets_unshuffled = deepcopy(sets)
  # sets_unshuffled = sets # Need to use this to match GLNS

  tour_history = Array{Tuple{Float64, Array{Int64,1}, Int64},1}()
  num_trials_feasible = 0
  num_trials = 0
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

	while true
    if count[:cold_trial] > param[:cold_trials] && lowest.cost <= param[:budget]
      break
    end

    if length(given_initial_tours) != 0
      start_idx = (count[:cold_trial] - 1)*num_sets + 1
      end_idx = count[:cold_trial]*num_sets
      initial_tour = given_initial_tours[start_idx:end_idx]
    else
      initial_tour = given_initial_tours
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

    while count[:warm_trial] <= param[:warm_trials]
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
        while true
          do_break = false
          bt = time()
          lock(count_lock)
          at = time()
          lock_times[thread_idx] += at - bt
          try
            if thread_broke || 
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
              updated_best = true
              best = tour_copy(trial)
              timer = (time_ns() - start_time)/1.0e9
              if param[:print_output] == 3
                println("Thread ", thread_idx, " found new best tour after ", timer, " s with cost ", best.cost, " (before opt cycle)")
              end
            end
          finally
            unlock(best_lock)
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
                best = tour_copy(trial)
                # print_best(count, param, best, lowest, init_time)
                timer = (time_ns() - start_time)/1.0e9
                if param[:print_output] == 3
                  println("Thread ", thread_idx, " found new best tour after ", timer, " s with cost ", best.cost)
                end

                if param[:output_file] != "None"
                  push!(tour_history, (round((time_ns() - start_time_for_tour_history)/1.0e9, digits=3), best.tour, best.cost))
                end
              end
            finally
              unlock(best_lock)
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
        finally
          unlock(num_trials_lock)
        end
      end
      print_warm_trial(count, param, best, iter_count)
      count[:warm_trial] += 1
      count[:latest_improvement] = 1
      count[:first_improvement] = false

      if time() - init_time > param[:max_time]
        param[:timeout] = true
        break
      end
    end
		lowest.cost > best.cost && (lowest = best)
		count[:warm_trial] = 0
		count[:cold_trial] += 1
    if time() - init_time > param[:max_time]
      break
    end
	end
	timer = (time_ns() - start_time)/1.0e9
  if param[:output_file] != "None"
    push!(tour_history, (round((time_ns() - start_time_for_tour_history)/1.0e9, digits=3), lowest.tour, lowest.cost))
  end

  print_summary(lowest, timer, membership, param, tour_history, cost_mat_read_time, instance_read_time, num_trials_feasible, num_trials, param[:timeout], lock_times)

  @assert(lowest.cost == tour_cost(lowest.tour, dist))
  @assert(length(lowest.tour) == num_sets)
  if param[:print_output] == 3
    println(lock_times)
  end
  return lowest.cost, powers
end

function parse_cmd(ARGS)
	if isempty(ARGS)
		println("no input instance given")
		exit(0)
	end
	if ARGS[1] == "-help" || ARGS[1] == "--help"
		println("Usage:  GTSPcmd.jl [filename] [optional flags]\n")
		println("Optional flags (vales are give in square brackets) :\n")
		println("-mode=[default, fast, slow]      (default is default)")
		println("-max_time=[Int]                  (default set by mode)")
		println("-trials=[Int]                    (default set by mode)")
		println("-restarts=[Int]                  (default set by mode)")
		println("-noise=[None, Both, Subset, Add] (default is Both)")
		println("-num_iterations=[Float]          (default set by mode. Number multiplied by # of sets)")
		println("-verbose=[0, 1, 2, 3]            (default is 3. 0 is no output, 3 is most.)")
		println("-output=[filename]               (default is None)")
		println("-epsilon=[Float in [0,1]]        (default is 0.5)")
		println("-reopt=[Float in [0,1]]          (default is 1.0)")
		println("-budget=[Int]                    (default has no budget)")
		println("-socket_port=[Int]               (default is 65432)")
		println("-new_socket_each_instance=[filename]    (default is 0)")
		exit(0)
	end
	int_flags = ["-max_time", "-trials", "-restarts", "-verbose", "-budget", "-socket_port", "-new_socket_each_instance"]
	float_flags = ["-epsilon", "-reopt", "-num_iterations", "-latest_improvement", "-first_improvement"]
	string_flags = ["-mode", "-output", "-noise", "-devel"]
	filename = ""
	optional_args = Dict{Symbol, Any}()
	for arg in ARGS
		temp = split(arg, "=")
		if length(temp) == 1 && filename == ""
			filename = temp[1]
		elseif length(temp) == 2
			flag = temp[1]
			value = temp[2]
			if flag in int_flags
				key = Symbol(flag[2:end])
				optional_args[key] = parse(Int64, value)
			elseif flag in float_flags
				key = Symbol(flag[2:end])
				optional_args[key] = parse(Float64, value)
			elseif flag in string_flags
				key = Symbol(flag[2:end])
				optional_args[key] = value
			else
				println("WARNING: skipping unknown flag ", flag, " in command line arguments")
			end
		else
			error("argument ", arg, " not in proper format")
		end
	end
	return filename, optional_args
end

function main(args::Vector{String}, max_time::Float64, inf_val::Int64, given_initial_tours::AbstractArray{Int64,1}, dist::AbstractArray{Int64,2}, max_threads::Int64, pin_cores::Vector{Int64}=Vector{Int64}())
  start_time_for_tour_history = time_ns()
  problem_instance, optional_args = parse_cmd(args)
  problem_instance = String(problem_instance)

  if max_time >= 0
    optional_args[Symbol("max_time")] = max_time
  end

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

  timing_result = @timed GLNS.solver(problem_instance, given_initial_tours, start_time_for_tour_history, inf_val, num_vertices, num_sets, sets, dist, membership, instance_read_time, cost_mat_read_time, max_threads, powers, update_powers, pin_cores; optional_args...)
  if get(optional_args, :verbose, 0) == 3
    println("Compile time: ", timing_result.compile_time)
  end
  return timing_result.value
end

end
