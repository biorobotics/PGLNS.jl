import Pkg
Pkg.activate(".")
using NPZ
using GLNS

# Need to specify trials=1 to avoid assertion failure, since if we don't pass initial tours, we use the default number of cold trials which is 5, and we don't support multiple cold trials yet
ARGS = ["/home/noopygbhat/GLKH-1.1/GTSPLIB/debug/custom0.gtsp", "-output=custom.tour", "-socket_port=65432", "-lazy_edge_eval=0", "-new_socket_each_instance=0", "-verbose=3", "-mode=fast", "-trials=1"]

npyfile = first(ARGS[1], length(ARGS[1]) - length(".gtsp")) * ".npy"
dist = npzread(npyfile)

given_initial_tours = Vector{Int64}()
@time GLNS.main(ARGS, 10., 298309430, given_initial_tours, dist)
@time GLNS.main(ARGS, 10., 298309430, given_initial_tours, dist)
