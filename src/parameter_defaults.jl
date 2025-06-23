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

"""
Three different parameter settings for the GTSP solver
default -- moderate runtime.  ~ 10-20 seconds for 100 sets
fast -- very quick, lower quality solution
slow -- slower, but potentially higher quality solution
"""	
function parameter_settings(num_vertices, num_sets, sets, problem_instance, args)	
	args = Dict(args)
	mode = get(args, :mode, "default")	
	############ default setting #####################
	if mode == "default"
		num_iterations = get(args, :num_iterations, 60) * num_sets
    latest_improvement = get(args, :latest_improvement, num_iterations/2)
    first_improvement = latest_improvement/2
		param = Dict(
		:cold_trials => get(args, :trials, 5),
		:warm_trials => get(args, :restarts, 3),
		:max_time => get(args, :max_time, 360),
		:init_tour => get(args, :init_tour, "rand"),

		:prob_reopt => get(args, :reopt, 1.0),
		:accept_percentage => 0.05,
		:prob_accept => 10.0/num_iterations,

		:num_iterations => num_iterations,
		:latest_improvement => latest_improvement,
		:first_improvement => first_improvement,
		:max_removals => min(100, max(round(Int64, 0.3*num_sets), 1)),
		:insertions => ["randpdf", "cheapest"],
		)
			
	################## very_fast  ##########################
	elseif mode == "fast"
    num_iterations_multiplier = get(args, :num_iterations, 60)
		num_iterations = num_iterations_multiplier * num_sets
    latest_improvement_multiplier = get(args, :latest_improvement, num_iterations_multiplier/4)
    latest_improvement = latest_improvement_multiplier * num_sets
    first_improvement_multiplier = get(args, :first_improvement, num_iterations_multiplier/6)
    first_improvement = first_improvement_multiplier * num_sets
    max_removal_fraction = get(args, :max_removal_fraction, 0.1)
    max_removals_cap = get(args, :max_removals_cap, 20)
		param = Dict(
		:cold_trials => get(args, :trials, 3),
		:warm_trials => get(args, :restarts, 2),
		:max_time => get(args, :max_time, 300),
		:init_tour => get(args, :init_tour, "insertion"),

		:prob_reopt => get(args, :reopt, 0.2),
		:accept_percentage => 0.05,
		:prob_accept => 10.0/num_iterations,

		:num_iterations => num_iterations,
		:latest_improvement => latest_improvement,
		:first_improvement => first_improvement,
		:max_removals => min(max_removals_cap, max(round(Int64, max_removal_fraction*num_sets), 1)),
		:insertions => ["randpdf"],
		)
	
	################## attempt slow search  ##########################
	elseif mode == "slow"
		num_iterations = get(args, :num_iterations, 150) * num_sets
    latest_improvement = get(args, :latest_improvement, num_iterations/3)
    first_improvement = latest_improvement/2
		param = Dict(
		:cold_trials => get(args, :trials, 10),
		:warm_trials => get(args, :restarts, 5),
		:max_time => get(args, :max_time, 1200),
		:init_tour => get(args, :init_tour, "rand"),

		:prob_reopt => get(args, :reopt, 1.0),
		:accept_percentage => 0.05,
		:prob_accept => 10.0/num_iterations,

		:num_iterations => num_iterations,
		:latest_improvement => latest_improvement,
		:first_improvement => first_improvement,
		:max_removals => max(round(Int64, 0.4*num_sets), 1),
		:insertions => ["randpdf", "cheapest"],
		)
		
	else
		error("mode not recognized.  Use default, fast, or slow")
	end

  # println("Latest improvement = ", param[:latest_improvement], ", first improvement = ", param[:first_improvement])
	
	# insertion algs
	if get(args, :insertion_algs, "default") == "cheapest"
		param[:insertions] = ["randpdf", "cheapest"]
	elseif get(args, :insertion_algs, "default") == "randpdf"
		param[:insertions] = ["randpdf"]
	end
	
	# insertion powers
	if get(args, :insertion_algs, "default") == "classic"
		param[:insertion_powers] = [-10.0, 0.0, 10.0]
	else
		param[:insertion_powers] = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 10.0]
	end
	
	# removal algs and powers
	if get(args, :removal_algs, "default") == "classic"
		param[:removals] = ["worst"]
		param[:removal_powers] = [0.0, 10.0]
	else
		param[:removals] = ["distance", "worst", "segment"]
		param[:removal_powers] = [-0.5, 0.0, 0.5, 1.0, 10.0]
	end
	
	# parameters that are common for all modes
	param[:mode] = mode
	param[:problem_instance] = split(problem_instance, "/")[end]
	param[:num_sets] = num_sets
	param[:num_vertices] = num_vertices
	param[:output_file] = get(args, :output, "None")
	param[:print_output] = get(args, :verbose, 0)
	param[:epsilon] = get(args, :epsilon, 0.5)
	param[:noise] = get(args, :noise, "Add")
	param[:adaptive_iter] = 1
	param[:print_time] = 5
	param[:budget] = get(args, :budget, typemin(Int64))
	param[:socket_port] = get(args, :socket_port, 65432)
	param[:lazy_edge_eval] = get(args, :lazy_edge_eval, 1)
	param[:new_socket_each_instance] = get(args, :new_socket_each_instance, 0)
	param[:timeout] = false
	param[:budget_met] = false
	param[:min_set] = min_set(sets)
	param[:min_removals] = (param[:max_removals] > 1 ? 2 : 1)
	param[:search_order] = get(args, :search_order, "astar")
	print_params(param)
	
	return param
end
