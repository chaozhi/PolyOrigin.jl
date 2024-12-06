
using PolyOrigin
using Plots
using Documenter
println("PolyOrigin_dir = ",pkgdir(PolyOrigin))
isdir(joinpath(pkgdir(PolyOrigin),"docs","build")) && error("delete old build!")

makedocs(
    modules=[PolyOrigin],
    authors="Chaozhi Zheng",
    repo="https://github.com/chaozhi/PolyOrigin.jl/blob/{commit}{path}#L{line}",    
    sitename="PolyOrigin.jl",
    format=Documenter.HTML(;
        # prettyurls=get(ENV, "CI", "false") == "true",
        prettyurls=true,
        canonical="https://chaozhi.github.io/PolyOrigin.jl",
    ),
    pages = [
        "Home" => "index.md", "Manual" => Any[
            "0 Get started" => "man/0_get_started.md",
            "1 Prepare input" => "man/1_prepare_input.md",
            "2 Run polyOrigin" => "man/2_run_polyOrigin.md",
            "3 Examples" => "man/3_examples.md",
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
    devbranch = "main",    
)
