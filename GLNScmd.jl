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
# Pkg.activate(expanduser("."))
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

function main()
  problem_instance, optional_args = GLNS.parse_cmd(ARGS)
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

  println("Server attempting to listen on port ", PORT)
  server = TCPSocket()
  try
    server = listen(PORT)
  catch e
    println("Server on port ", PORT, " failed to listen")
    exit()
  end
  println("Server listening on port ", PORT)

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
        println("Server on port ", PORT, " received termination signal")
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

      read_start_time = time_ns()
      num_vertices, num_sets, sets, _, membership = read_file(problem_instance, false)
      read_end_time = time_ns()
      instance_read_time = (read_end_time - read_start_time)/1.0e9
      println("Reading GTSPLIB file took ", instance_read_time, " s")

      # Read cost matrix from npy file
      read_start_time = time_ns()
      npyfile = first(problem_instance, length(problem_instance) - length(".gtsp")) * ".npy"
      dist = npzread(npyfile)
      read_end_time = time_ns()
      cost_mat_read_time = (read_end_time - read_start_time)/1.0e9
      println("Reading cost mat file took ", cost_mat_read_time, " s")

      GLNS.solver(problem_instance, given_initial_tours, start_time_for_tour_history, inf_val, num_vertices, num_sets, sets, dist, membership, instance_read_time, cost_mat_read_time, 10; optional_args...)
      write(client_socket, "solved\n")
      iter_count += 1
    end
  finally
    close(server)
    println("Closed server on port ", PORT)
  end
end

main()
