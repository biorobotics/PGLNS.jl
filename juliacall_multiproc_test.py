# run "JULIA_NUM_THREADS=10 PYTHON_JULIACALL_HANDLE_SIGNALS=yes python3 juliacall_multiproc_test.py" without quotes in a terminal

from multiprocessing import Pool, Process
import time
import os
import numpy as np
import psutil
import sys
sys.path.append(os.path.expanduser("~/catkin_ws/src/mapf/scripts/"))
from solve_gtsp_glkh import solve_gtsp_glns_lazy_edge_eval

def proc_fn(i):
  bt = time.perf_counter()
  from juliacall import Main as jl
  from juliacall import convert
  at = time.perf_counter()
  print(at - bt)

  # If we're not using a sysimage. Results in ~0.07% compilation time even after solving the dummy GTSP below.
  # Pkg activate and import run pretty fast, though
  jl.seval("import Pkg")
  jl.seval('Pkg.activate(expanduser("~/PGLNS.jl"))')
  jl.seval("import GLNS")

  # Despite using a sysimage, something is still getting compiled at runtime, and I'm not sure what. Solve this dummy GTSP to trigger remaining compilation
  gtsp_cost_mat = np.array([[0, 1], [1, 0]])
  solve_gtsp_glns_lazy_edge_eval(gtsp_cost_mat, [[0], [1]], file_suffix=str(0))
  ARGS = [os.path.expanduser("~/GLKH-1.1/GTSPLIB/custom" + str(0) + ".gtsp"), "-output=custom.tour", "-verbose=3", "-mode=fast", "-num_iterations=60"]
  problem_instance = ARGS[0]
  npyfile = problem_instance[:-len(".gtsp")] + ".npy"
  dist = np.load(npyfile).astype(int)
  ARGS_jl = convert(jl.Vector[jl.String], ARGS)
  jl.GLNS.main(ARGS_jl, 10., 298309430, np.array([1, 2]), np.asfortranarray(dist), 10)

  instance_folder = "debug"
  i = 0
  ARGS = [os.path.expanduser("~/GLKH-1.1/GTSPLIB/" + instance_folder + "/custom" + str(i) + ".gtsp"), "-output=custom.tour", "-verbose=3", "-mode=fast", "-num_iterations=60"]

  given_initial_tours = np.array([   0,    1,   35,   86,   18,   69,  120,  137,  103,   52,  154,  171,  222,  205,  188,  273,  239,  256,  307,  290,  341,  324,  392,  358,  375,  426,  460,  409,  443,  494,  477,  528,  562,  511,  596,  545,  579,  613,  681,  647,  630,  664,  715,  698,  749,  732,  766,  834,  817,  783,  800,  851,  885,  868,  936,  902,  919,  970,  987,  953, 1004, 1021, 1055, 1089, 1072, 1038, 1106, 1123, 1140, 1191, 1208, 1157, 1174, 1242, 1225, 1276, 1310, 1344, 1293, 1259, 1327, 1378, 1412, 1429, 1361, 1395, 1463, 1480, 1497, 1531, 1548, 1446, 1565, 1514, 1582, 1616, 1599, 1650, 1633, 1684, 1718, 1735, 1752, 1667, 1769, 1701, 1837, 1854, 1786, 1803, 1820, 1888, 1905, 1871, 1922, 1939, 1956, 1990, 2007, 1973, 2024, 2041, 2058, 2075, 2126, 2092, 2160, 2109, 2143, 2177, 2194, 2228, 2211, 2245, 2279, 2296, 2313, 2262, 2347, 2364, 2330, 2381, 2415, 2432, 2398, 2483, 2466, 2500, 2449, 2517, 2551, 2585, 2602, 2534, 2636, 2568, 2653, 2619, 2670, 2738, 2704, 2687, 2789, 2772, 2721, 2806, 2755, 2823, 2840, 2857, 2908, 2925, 2891, 2874, 3010, 2942, 2959, 2993, 2976, 3061, 3027, 3044, 3095, 3112, 3146, 3078, 3163, 3180, 3197, 3129, 3231, 3248, 3265, 3214, 3282, 3299, 3316, 3333, 3367, 3350, 3384]) + 1

  problem_instance = ARGS[0]
  npyfile = problem_instance[:-len(".gtsp")] + ".npy"
  dist = np.load(npyfile).astype(int)

  ARGS_jl = convert(jl.Vector[jl.String], ARGS)

  bt = time.perf_counter()
  jl.GLNS.main(ARGS_jl, 10., 298309430, given_initial_tours, np.asfortranarray(dist), 10)
  at = time.perf_counter()
  print(at - bt)

num_proc = 1 # psutil.cpu_count(logical=False)

procs = []
for proc_idx in range(num_proc):
  procs.append(Process(target=proc_fn, args=(proc_idx,)))
  procs[-1].start()

for proc_idx in range(num_proc):
  procs[proc_idx].join()
