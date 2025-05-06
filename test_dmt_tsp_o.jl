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
Pkg.activate(expanduser("~/PGLNS_dp_repair.jl"))
using GLNS
using Printf
using NPZ
using PythonCall
using Profile, FileIO

function main()
  # instance_folder = "debug"
  instance_folder = "dmt_tsp_o"
  # instance_folder = ""

  i = 0
  ARGS = [expanduser("~/GLKH-1.1/GTSPLIB/"*instance_folder*"/custom"*string(i)*".gtsp"), "-output=custom.tour", "-socket_port=65432", "-lazy_edge_eval=0", "-new_socket_each_instance=0", "-verbose=3", "-mode=fast", "-num_iterations=8", "-latest_improvement=2", "-first_improvement=1.33", "-max_removal_fraction=0.1", "-max_removals_cap=20"]

  problem_instance = ARGS[1]
  npyfile = first(problem_instance, length(problem_instance) - length(".gtsp")) * "_initial_tour.npy"
  given_initial_tours = npzread(npyfile)[1:end-1] .+ 1

  npyfile = first(problem_instance, length(problem_instance) - length(".gtsp")) * ".npy"
  dist = npzread(npyfile)
  inf_val = maximum(dist)

  GLNS.main(ARGS, 100., inf_val, PyArray{Int64, 1, true, true, Int64}(given_initial_tours), PyArray{Int64, 2, true, true, Int64}(dist), 8)
  GLNS.main(ARGS, 100., inf_val, PyArray{Int64, 1, true, true, Int64}(given_initial_tours), PyArray{Int64, 2, true, true, Int64}(dist), 8)
end

main()
