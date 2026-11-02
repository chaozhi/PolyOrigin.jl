using Pkg

function dev_pkg(pkg::AbstractString; path::AbstractString = abspath(@__DIR__))
    try 
        Pkg.activate()        
        if haskey(Pkg.project().dependencies, pkg)
            @info "Removing package: $pkg..."
            Pkg.rm(pkg)
        else
            @info "Package '$pkg' is already removed or not found."
        end    
        Pkg.develop(Pkg.PackageSpec(name=pkg, path=path))
        Pkg.activate()
        Pkg.instantiate() 
        Pkg.gc()
        0
    catch err
        @error "Failed to develop package $pkg"
        rethrow(err)
        -1
    end
end

dev_pkg("PolyOrigin")
