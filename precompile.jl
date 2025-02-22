using PackageCompiler
using Pkg
Pkg.activate(expanduser("."))
using GLNS

# create_sysimage(["GLNS"]; sysimage_path="GLNS.so",precompile_execution_file="precompilation_script.jl")
create_sysimage(["GLNS"]; sysimage_path="GLNS.so",precompile_statements_file="precompilation_statements.jl")
