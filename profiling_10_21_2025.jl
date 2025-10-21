import Pkg
Pkg.activate(".")
using NPZ
using GLNS

load_path = expanduser("~/catkin_ws/src/mapf/data/10_21_2025/pglns_astar_profiling/")
ARGS = [load_path*"custom0.gtsp", "-output=custom.tour", "-verbose=3", "-mode=fast", "-num_iterations=4", "-search_order=astar"]

dist_file = load_path*"dist.npy"
dist = npzread(dist_file)

initial_tour_file = load_path*"initial_tour.npy"
given_initial_tours = npzread(initial_tour_file)

inf_val = maximum(dist)

@time GLNS.main(ARGS, 1., inf_val, given_initial_tours, dist, 1)
@time GLNS.main(ARGS, 10., inf_val, given_initial_tours, dist, 1)
