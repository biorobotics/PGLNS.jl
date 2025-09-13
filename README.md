# PGLNS

Before using this package, install dependencies by running
```cmd
julia
```
in a terminal, then running
```julia
julia> import Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

This is meant to be used via the main function in src/GLNS.jl. Build a precompiled sysimage as follows:
```cmd
julia precompile.jl
julia -JGLNS.so
```
