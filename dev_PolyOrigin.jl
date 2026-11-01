
using Pkg
Pkg.activate()
try 
    Pkg.rm("PolyOrigin")
catch err
    @warn "PolyOrigin not installed"
    @warn err
end
Pkg.develop(path = abspath(@__DIR__))
# @time Pkg.build("PolyOrigin")
Pkg.activate()
