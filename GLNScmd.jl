#!/usr/bin/env julia

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

using Sockets
import Pkg
# Pkg.activate(expanduser("~/GLNS_lazy_edge_eval.jl"))
# using GLNS
using Printf
using Dates

using NPZ
include("src/utilities.jl")
include("src/parse_print.jl")
include("src/tour_optimizations.jl")
include("src/adaptive_powers.jl")
include("src/insertion_deletion.jl")
include("src/parse_print.jl")

"""
Optional Flags -- values are given in square brackets []  
	-max_time=[Int]				 (default set by mode)
	-trials=[Int]				 (default set by mode)
	-restarts=[Int]              (default set by mode)
	-mode=[default, fast, slow]  (default is default)
	-verbose=[0, 1, 2, 3]        (default is 3.  0 is no output, 2 is most verbose)
	-output=[filename]           (default is output.tour)
	-epsilon=[Float in [0,1]]	 (default is 0.5)
	-reopt=[Float in [0,1]]      (default is set by mode)
	-socket_port=[Int>1023]      (default is 65432)
	-lazy_edge_eval=[Int in [0,1]]      (default is 1)
"""
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
		println("-num_iterations=[Int]            (default set by mode. Number multiplied by # of sets)")
		println("-verbose=[0, 1, 2, 3]            (default is 3. 0 is no output, 3 is most.)")
		println("-output=[filename]               (default is None)")
		println("-epsilon=[Float in [0,1]]        (default is 0.5)")
		println("-reopt=[Float in [0,1]]          (default is 1.0)")
		println("-budget=[Int]                    (default has no budget)")
		println("-socket_port=[Int]               (default is 65432)")
		println("-lazy_edge_eval=[Int]            (default is 1)")
		println("-new_socket_each_instance=[filename]    (default is 0)")
		exit(0)
	end
	int_flags = ["-max_time", "-trials", "-restarts", "-verbose", "-budget", "-num_iterations", "-socket_port", "-lazy_edge_eval", "-new_socket_each_instance"]
	float_flags = ["-epsilon", "-reopt"]
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

function main()
  problem_instance, optional_args = parse_cmd(ARGS)
  problem_instance = String(problem_instance)

	output_file = get(optional_args, :output, "None")
  if output_file != "None"
    f = open(output_file, "w")
    write(f, "\n")
    close(f)
  end

  if haskey(optional_args, Symbol("socket_port"))
    PORT = optional_args[Symbol("socket_port")]
  else
    PORT = 65432
  end

  @printf("Server attempting to listen on port %d\n", PORT)
  server = TCPSocket()
  try
    server = listen(PORT)
  catch e
    @printf("Server on port %d failed to listen\n", PORT)
    exit()
  end
  @printf("Server listening on port %d\n", PORT)

  client_socket = accept(server)

  try
    iter_count = 0
    while true
      if iter_count != 0 && haskey(optional_args, Symbol("new_socket_each_instance")) && optional_args[Symbol("new_socket_each_instance")] == 1
        client_socket = accept(server)
      end
      msg = readline(client_socket)
      start_time_for_tour_history = time_ns()
      if msg == "terminate"
        @printf("Server on port %d received termination signal", PORT)
        break
      end
      if length(msg) == 0
        iter_count += 1
        continue # Assume a client just closed its connection
      end
      if !isfile(problem_instance)
        println("the problem instance  ", problem_instance, " does not exist")
        break
      end
      msg_split = split(msg, " ")
      optional_args[Symbol("max_time")] = parse(Float64, msg_split[1])
      inf_val = parse(Int64, msg_split[2])
      given_initial_tours = Vector{Int64}()
      for node_idx_str in msg_split[3:end]
        push!(given_initial_tours, parse(Int64, node_idx_str))
      end

      # Get already evaluated edges
      evaluated_edges = Vector{Tuple{Int64, Int64}}()
      open_tsp = false
      if optional_args[Symbol("lazy_edge_eval")] == 1
        msg = readline(client_socket)
        if msg == "terminate\n"
          @printf("Server on port %d received termination signal", PORT)
          break
        end
        if length(msg) == 0
          iter_count += 1
          continue # Assume a client just closed its connection
        end
        msg_split = split(msg, " ")
        for edge_str in msg_split
          if edge_str == "o"
            open_tsp = true
          else
            node_strs = split(edge_str, "-")
            push!(evaluated_edges, (parse(Int64, node_strs[1]), parse(Int64, node_strs[2])))
          end
        end
      end

      # do_perf = true
      do_perf = false
      perf_file = ""
      if do_perf
        msg = readline(client_socket)
        if msg == "terminate\n"
          @printf("Server on port %d received termination signal", PORT)
          break
        end
        if length(msg) == 0
          iter_count += 1
          continue # Assume a client just closed its connection
        end
        perf_file = msg
      end

      read_start_time = time_ns()
      num_vertices, num_sets, sets, membership = read_file(problem_instance)
      read_end_time = time_ns()
      instance_read_time = (read_end_time - read_start_time)/1.0e9
      @printf("Reading GTSPLIB file took %f s\n", instance_read_time)

      # Read cost matrix from npy file
      read_start_time = time_ns()
      npyfile = first(problem_instance, length(problem_instance) - length(".gtsp")) * ".npy"
      dist = npzread(npyfile)
      read_end_time = time_ns()
      cost_mat_read_time = (read_end_time - read_start_time)/1.0e9
      @printf("Reading cost mat file took %f s\n", cost_mat_read_time)

      GLNS.solver(problem_instance, client_socket, given_initial_tours, start_time_for_tour_history, inf_val, evaluated_edges, open_tsp, num_vertices, num_sets, sets, dist, membership, instance_read_time, cost_mat_read_time, do_perf, perf_file; optional_args...)
      #=
      timing_result = @timed GLNS.solver(problem_instance, client_socket, given_initial_tours, start_time_for_tour_history, inf_val, evaluated_edges, open_tsp, num_vertices, num_sets, sets, dist, membership, instance_read_time, cost_mat_read_time, do_perf, perf_file; optional_args...)
      println(timing_result)
      if timing_result.compile_time >= 1e-3
        filename = "err_time.txt"
        if !isfile(filename)
          f = open(filename, "w")
          write(f, Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"), "\n")
          write(f, string(timing_result), "\n")
          close(f)
        end
      end
      =#
      # @assert(timing_result.compile_time <= 1e-3)
      write(client_socket, "solved\n")
      iter_count += 1
    end
  finally
    close(server)
    @printf("Closed server on port %d\n", PORT)
  end
end

main()
