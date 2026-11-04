function tryusing(
    pkgname::AbstractString;
    pkgurl::Union{Nothing,AbstractString} = nothing,
)
    try
        @eval using $(Symbol(pkgname))
    catch
        @eval using Pkg
        if isnothing(pkgurl)
            Pkg.add(pkgname)
        else
            Pkg.add(PackageSpec(url = pkgurl, rev = "master"))
        end
        @eval using $(Symbol(pkgname))
    end
end

tryusing("PolyOrigin", pkgurl = "https://github.com/chaozhi/PolyOrigin.jl")
tryusing("ArgParse")

function parse_commandline()
    s = ArgParseSettings(
        prog = "polyorigin",
        version = string(pkgversion(PolyOrigin)), 
        add_version = false
    )
    s.description = "Haplotype reconstruction in polypoid multiparental populations"
    workdir = pwd()
    @add_arg_table! s begin
        "--genofile", "-g"
        help = "filename for genotypic data file"
        arg_type = AbstractString
        required = true
        "--pedfile", "-p"
        help = "filename for pedigree info"
        arg_type = AbstractString
        required = true
        "--delimchar"
        help = "text delimiter"
        arg_type = AbstractChar
        default = ','
        "--missingstring"
        help = "string code for missing value"
        arg_type = AbstractString
        default = "NA"
        "--commentstring"
        help = "rows that begins with commentstring will be ignored"
        arg_type = AbstractString
        default = "#"
        "--isphysmap"
        help = "if true, input markermap is physical map (location in bp)"
        arg_type = Bool
        default = false
        "--recomrate"
        help = "recombination rate in unit of cM/Mbp"
        arg_type = Float64
        default = 1.0  
        "--chrpairing_phase"
        help = "chromosome pairing in parental phasing, with 22 being only
        bivalent formations and 44 being bi- and quadri-valent formations"
        arg_type = Int
        default = 22
        "--chrpairing"
        help = "chromosome pairing in offspring decoding, with 22 being only
        bivalent formations and 44 being bivalent and quadrivalent formations"
        arg_type = Int
        default = 44
        "--chrsubset"
        help = "subset of chromosomes, with nothing denoting all chromosomes,
        e.g, \"[2,10]\" denotes the second and tenth chromosomes"
        arg_type = String
        default = "nothing"
        "--snpthin"
        help = "subset of markers by taking every snpthin-th markers"
        arg_type = Int
        default = 1
        "--maxstuck"
        help = "the max number of consecutive iterations that are rejected
        in a phasing run"
        arg_type = Int
        default = 5
        "--maxiter"
        help = "the max number of iterations in a phasing run"
        arg_type = Int
        default = 30
        "--minrun"
        help = "if the min number of phasing runs that are at the same local maximimum or
        have the same parental phases reaches minrun, phasing algorithm will stop before reaching the maxrun."
        arg_type = Int
        default = 3
        "--maxrun"
        help = "the max number of phasing runs"
        arg_type = Int
        default = 3
        "--refhapfile"
        help = "reference haplotype file
        for setting absolute parental phases. It has the same format as the input genofile,
        except that parental genotypes are phased and offspring genotypes are ignored if they exist."
        arg_type = AbstractString
        default = "nothing"
        "--trueorigfile"
        help = "true ancestral origin file produced by PolyoriginSim for setting absoluting phasing of ancestral orgin states"
        arg_type = AbstractString
        default = "nothing"
        "--isrefine"
        help = "if true, refine marker map"
        arg_type = Bool
        default = true
        "--isrefineorder"
        help = "if true, refine marker mordering, valid only if isrefine=true"
        arg_type = Bool
        default = false
        "--isrefineerror"
        help = "if true, estimating genotypic error parameters"
        arg_type = Bool
        default = true
        "--maxwinsize"
        help = "max size of sliding windown in map refinning"
        arg_type = Int
        default = 50
        "--inittemperature"
        help = "initial temperature of simulated annealing in map refinning"
        arg_type = Float64
        default = 4.0
        "--coolingrate"
        help = "cooling rate of annealing temperature in map refinning"
        arg_type = Float64
        default = 0.5
        "--stripdis"
        help = "a chromosome end in map refinement is removed if it has a distance gap > stripdis
        (centiMorgan) and it contains less than 5% markers."
        arg_type = Float64
        default = 20.0
        "--skeletonsize"
        help = "the number of markers in the skeleton map that is used
        to reduce map length inflation by subsampling markers"
        arg_type = Int
        default = 100
        "--nworker"
        help = "number of parallel workers for computing among chromosomes"
        arg_type = Int
        default = 1
        "--isplot"
        help = "if true, plot haploprob"
        arg_type = Bool
        default = false
        "--outstem", "-o"
        help = "stem of output filenames"
        arg_type = AbstractString
        default = "outstem"
        "--workdir", "-w"
        help = "directory for reading and writing files"
        arg_type = AbstractString
        default = workdir
        "--verbose", "-v"
        help = "if true, print messages on console"
        arg_type = Bool
        default = true
        "--version", "-V"
        action = :show_version
        help = "show PolyOrigin's version and exit"
    end
    return parse_args(s, as_symbols = true)
end

function string2vec(strvec::String, t::DataType)
    str = strip(strvec)
    @assert str[1] == '[' && str[end] == ']'
    parse.(t, split(str[2:end-1], ","))
end

function main(args::Vector{String})
    parsed_args = parse_commandline()
    verbose = parsed_args[:verbose]
    if verbose
        println("Parsed arguments:")
        for (arg, val) in parsed_args
            if isnothing(val)
                println("  $arg  =>  nothing ")
            else
                println("  $arg  =>  $val")
                val == "nothing" && (parsed_args[arg] = nothing)
            end
        end
    end
    outstem = parsed_args[:outstem]
    logfile = string(outstem, ".log")
    genofile = parsed_args[:genofile]
    pedfile = parsed_args[:pedfile]
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedfile)
    push!(parsed_args, :logfile => logfile)
    snpthin = parsed_args[:snpthin]
    # assum the max number of markers in a linkage group < 10^6
    snpsubset= snpthin<=1 ? nothing : 1:snpthin:10^6
    delete!(parsed_args, :snpthin)
    push!(parsed_args, :snpsubset => snpsubset)
    nworker = parsed_args[:nworker]
    isparallel = nworker <= 1 ? false : true
    if isparallel
        tryusing("Distributed")
        addprocs(nworker) # add worker processes on local machine
        @info string("#parallel workers =", nworkers())
        @eval @everywhere using PolyOrigin
    end
    delete!(parsed_args, :nworker)
    push!(parsed_args, :isparallel => isparallel)
    a = parsed_args[:chrsubset]
    if !isnothing(a) && occursin("[", a) && occursin("]", a)
        parsed_args[:chrsubset] = string2vec(a, Int)
    end
    @time polyOrigin(genofile, pedfile; parsed_args...)
    workdir = parsed_args[:workdir]
    outfiles = filter(x -> occursin(outstem, x), readdir(workdir))
    verbose && println("output files: ", join(outfiles, ","))
    if isparallel
        rmprocs(workers()...)
    end
    return 0
end

main(ARGS)
