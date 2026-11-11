
using Plots
using DataFrames # 
using Documenter
try 
    using PolyOrigin
catch
    using Pkg
    Pkg.activate()
    pkgpath = abspath(@__DIR__,"..")
    Pkg.develop(PackageSpec(path = pkgpath))    
end

println("PolyOrigin_dir = ",pkgdir(PolyOrigin))
isdir(joinpath(@__DIR__,"build")) && error("delete old build!")

makedocs(
    modules=[PolyOrigin],
    authors="Chaozhi Zheng",    
    sitename="PolyOrigin2",
    repo =  Documenter.Remotes.GitHub("chaozhi", "PolyOrigin.jl"),        
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        # prettyurls=false, # false for local manual
        canonical="https://chaozhi.github.io/PolyOrigin.jl",
        # collapselevel = 1, 
    ),
    checkdocs = :none, 
    # draft = true, 
    pages = [
        "Home" => "index.md", "Manual" => Any[
            "0 Get started" => "man/0_get_started.md",
            "1 Prepare input" => "man/1_prepare_input.md",
            "2 Run polyorigin in Julia" => "man/2_run_polyorigin.md",
            "2* Run polyorigin in cmd" => "man/2_run_polyorigin_cmd.md",
            "3 Understand output" => "man/3_understand_output.md",
            "4 Examples" => "man/4_examples.md",                
        ],
        # "Library" => Any[
        #     "Public" => "lib/public.md",
        # ],
        # "Index" => "symbols.md",
    ]
)

deploydocs(;
    repo="github.com/chaozhi/PolyOrigin.jl.git",
    devbranch = "v2dev",
)
