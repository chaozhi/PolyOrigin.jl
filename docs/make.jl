
using Plots
using Documenter
try 
    using PolyOrigin
catch
    using Pkg
    Pkg.develop(PackageSpec(path = abspath(@__DIR__,"..")))
    # Pkg.instantiate()
end


println("PolyOrigin_dir = ",pkgdir(PolyOrigin))
isdir(joinpath(@__DIR__,"build")) && error("delete old build!")


makedocs(
    modules=[PolyOrigin],
    authors="Chaozhi Zheng",    
    sitename="PolyOrigin.jl",    
    # repo =  Documenter.Remotes.GitHub("chaozhi", "PolyOrigin.jl"),    
    repo = "https://github.com/chaozhi/PolyOrigin.jl",
    # remotes = nothing,
    format=Documenter.HTML(;
        # prettyurls=get(ENV, "CI", nothing) == "true",
        prettyurls=true, # false for local browsing        
        canonical="https://chaozhi.github.io/PolyOrigin.jl",      
        repolink = "https://github.com/chaozhi/PolyOrigin.jl",  
    ),
    pages = [
        "Home" => "index.md", "Manual" => Any[
            "0 Get started" => "man/0_get_started.md",
            "1 Prepare input" => "man/1_prepare_input.md",
            "2 Run polyOrigin" => "man/2_run_polyOrigin.md",
            "3 Understand output" => "man/3_understand_output.md",
            "4 Examples" => "man/4_examples.md",
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md",
        ],
        "Index" => "symbols.md",
    ]
)

deploydocs(;
    repo="github.com/chaozhi/PolyOrigin.jl.git",           
    devbranch = "main" 
)


# cd(@__DIR__)