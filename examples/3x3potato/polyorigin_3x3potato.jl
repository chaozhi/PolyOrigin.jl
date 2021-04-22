using Distributed
addprocs(4) # add worker processes on local machine
@info string("nworkers=", nworkers())
@everywhere  using PolyOrigin

using PolyOrigin
cd(@__DIR__)
pwd()

# run polyOrigin
genofile = "TableS2_dose.csv"
pedfile = "TableS3_ped.csv"
outstem = "potato_output"
chrsubset = 1:4
@time polyancestry = polyOrigin(genofile, pedfile;
    isphysmap = true,
    recomrate = 1.25, # average recombination rate (cM/Mbp) for potato
    snpsubset = 1:5:10000, # any integer >= #marker
    chrsubset = chrsubset,
    refinemap = true,
    refineorder = false, # estimate distances but not ordering
    # isplot = true,
    isparallel = true,
    outstem,
)
outfiles = filter(x -> occursin(outstem, x), readdir())
println("outfiles=", outfiles)

# plot relative frequencies of valent configurations
# polyancestry = readPolyAncestry(outstem*"_polyancestry.csv")
valentfreq = calvalentfreq(polyancestry)
plotvalentfreq(valentfreq)

# comparing maps
# polygeno.markermap was obtained from input physmap (bp) with
# default recomrate = 1 cM/Mbp
polygeno = readPolyGeno(genofile, pedfile; isphysmap = true)
fig = plotMapComp(
    polygeno.markermap[chrsubset],
    polyancestry.markermap;
    xlabel = "Physical poistion (Mbp)",
    ylabel = "Estimated position (cM)"
)

# delete output files
# rm.(filter(x->occursin(outstem,x), readdir()))
