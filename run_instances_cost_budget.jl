import Pkg
Pkg.activate(".")
using GLNS
using Base.Threads

BAF_folder = expanduser("~/gtsp_instances/BAF/INSTANCES/")
MOM_folder = expanduser("~/gtsp_instances/MOM/INSTANCES/")
gtsplib_a_folder = expanduser("~/gtsp_instances/GTSP_all/a/")
gtsplib_s_folder = expanduser("~/gtsp_instances/GTSP_all/s/")
gtsplib_plus_folder = expanduser("~/gtsp_instances/gtsp+_lib/")
large_folder = expanduser("~/gtsp_instances/large_lib/large_lib/")
sat_folder = expanduser("~/gtsp_instances/sat_lib/sat_lib/")

# folders = [BAF_folder, MOM_folder, gtsplib_a_folder, gtsplib_s_folder, gtsplib_plus_folder, large_folder, sat_folder]
folders = [BAF_folder]
given_initial_tours = Vector{Int64}()
num_thread_options = [1, Threads.nthreads()]
for folder=folders
  for filename=readdir(folder)
    if !occursin(".gtsp", filename)
      continue
    end
    cost_1thread = 0
    for num_threads in num_thread_options
      soln_dir = "solutions_"*string(num_threads)*"threads/"
      if !isdir(soln_dir)
        mkdir(soln_dir)
      end
      ARGS = [folder*filename, "-output="*soln_dir*filename*".tour", "-verbose=3"]
      if num_threads != 1
        ARGS = [folder*filename, "-output="*soln_dir*filename*".tour", "-verbose=3", "-budget="*string(cost_1thread)]
      end
      timing_result = @timed GLNS.main(ARGS, -1., given_initial_tours, false, num_threads)
      while timing_result.compile_time != 0
        timing_result = @timed GLNS.main(ARGS, -1., given_initial_tours, false, num_threads)
      end
      if num_threads == 1
        cost_1thread = timing_result.value[1]
      end
    end
  end
end
