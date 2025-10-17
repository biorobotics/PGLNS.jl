using PackageCompiler
using Pkg
Pkg.activate(expanduser("."))
using GLNS

create_sysimage(["GLNS"]; sysimage_path="GLNS_solver_state.so",precompile_statements_file="precompilation_statements_solver_state.jl")
