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
Pkg.activate(expanduser("."))
using GLNS
using Printf
using NPZ
include("src/utilities.jl")
include("src/parse_print.jl")
include("src/tour_optimizations.jl")
include("src/adaptive_powers.jl")
include("src/insertion_deletion.jl")
include("src/parse_print.jl")

function main()
  # Should be in the /home/cobra/GLKH-1.1/GTSPLIB folder
  instance_folder = "debug"

  for mode=["slow", "default", "fast"]
    for i=1:2
      ARGS = ["/home/cobra/GLKH-1.1/GTSPLIB/"*instance_folder*"/custom"*string(i)*".gtsp", "-output=custom.tour", "-socket_port=65432", "-new_socket_each_instance=0", "-verbose=3", "-mode="*mode]

      problem_instance, optional_args = GLNS.parse_cmd(ARGS)
      problem_instance = String(problem_instance)

      read_start_time = time_ns()
      num_vertices, num_sets, sets, _, membership = read_file(problem_instance, false)
      read_end_time = time_ns()
      instance_read_time = (read_end_time - read_start_time)/1.0e9

      read_start_time = time_ns()
      npyfile = first(problem_instance, length(problem_instance) - length(".gtsp")) * ".npy"
      dist = npzread(npyfile)
      read_end_time = time_ns()
      cost_mat_read_time = (read_end_time - read_start_time)/1.0e9

      given_initial_tours = [   0,    1,   35,   86,   18,   69,  120,  137,  103,   52,  154,  171,  222,  205,  188,  273,  239,  256,  307,  290,  341,  324,  392,  358,  375,  426,  460,  409,  443,  494,  477,  528,  562,  511,  596,  545,  579,  613,  681,  647,  630,  664,  715,  698,  749,  732,  766,  834,  817,  783,  800,  851,  885,  868,  936,  902,  919,  970,  987,  953, 1004, 1021, 1055, 1089, 1072, 1038, 1106, 1123, 1140, 1191, 1208, 1157, 1174, 1242, 1225, 1276, 1310, 1344, 1293, 1259, 1327, 1378, 1412, 1429, 1361, 1395, 1463, 1480, 1497, 1531, 1548, 1446, 1565, 1514, 1582, 1616, 1599, 1650, 1633, 1684, 1718, 1735, 1752, 1667, 1769, 1701, 1837, 1854, 1786, 1803, 1820, 1888, 1905, 1871, 1922, 1939, 1956, 1990, 2007, 1973, 2024, 2041, 2058, 2075, 2126, 2092, 2160, 2109, 2143, 2177, 2194, 2228, 2211, 2245, 2279, 2296, 2313, 2262, 2347, 2364, 2330, 2381, 2415, 2432, 2398, 2483, 2466, 2500, 2449, 2517, 2551, 2585, 2602, 2534, 2636, 2568, 2653, 2619, 2670, 2738, 2704, 2687, 2789, 2772, 2721, 2806, 2755, 2823, 2840, 2857, 2908, 2925, 2891, 2874, 3010, 2942, 2959, 2993, 2976, 3061, 3027, 3044, 3095, 3112, 3146, 3078, 3163, 3180, 3197, 3129, 3231, 3248, 3265, 3214, 3282, 3299, 3316, 3333, 3367, 3350, 3384] .+ 1

      inf_val = maximum(dist)

      @time GLNS.solver(problem_instance, given_initial_tours, time_ns(), inf_val, num_vertices, num_sets, sets, dist, membership, instance_read_time, cost_mat_read_time, 10; optional_args...)
    end
  end
end

main()
