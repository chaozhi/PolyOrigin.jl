# the following three lines are only for experimenting
# using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__, "..",".."))

# load the pakcage and set work directory
using PolyOrigin
cd(dirname(@__FILE__))
pwd()

# run polyOrigin
dataid = "sim4x"
genofile = dataid*"_geno.csv"
pedfile = dataid*"_ped.csv"
outstem = dataid*"_output"
@time polyancestry = polyOrigin(genofile, pedfile;
    refinemap=true,
    refineorder=true,
    outstem
)
outfiles = filter(x -> occursin(outstem, x), readdir())
println(outfiles)

# plot relative frequencies of valent configurations
# polyancestry = readPolyAncestry(outstem*"_polyancestry.csv")
valentfreq = calvalentfreq(polyancestry)
plotvalentfreq(valentfreq)

# calculate the accuracy of parental phasing and ancestral inference
truefile = dataid*"_true.csv"
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

# compare maps
plotMapComp(truegeno.truemap,polyancestry.markermap;
    xlabel = "True map(cM)",
    ylabel = "Estimated map(cM)"
)

# plot conditional probabilities
plotCondprob(polyancestry; truegeno = truegeno, offspring = 1)
animCondprob(polyancestry;
    truegeno = truegeno,
    fps = 0.5,
    # outfile = "outstem_condprob.gif",
)

# delete output files
# rm.(filter(x->occursin(outstem,x), readdir()))
