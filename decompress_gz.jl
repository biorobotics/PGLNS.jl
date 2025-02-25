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

folders = [BAF_folder, MOM_folder, gtsplib_a_folder, gtsplib_s_folder, gtsplib_plus_folder, large_folder, sat_folder]
given_initial_tours = Vector{Int64}()
num_thread_options = [1, Threads.nthreads()]
for folder=folders
  for filename=readdir(folder)
    if !occursin(".gtsp", filename)
      continue
    end
    if occursin(".gz", filename)
      filepath = folder*filename
      println(filepath)
      cmd = `gzip -d $filepath`
      run(cmd)
    end
  end
end
