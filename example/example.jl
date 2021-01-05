using Revise
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using PolyOrigin
cd(@__DIR__)
pwd()

# run polyorigin
genofile = "geno.csv"
pedfile = "ped.csv"
@time polyancestry = polyOrigin(genofile, pedfile;
    refinemap=true,
    refineorder=true,
    # isplot=true,
)

polygeno = readPolyGeno(genofile, pedfile)
PolyOrigin.plotdesign(polygeno)

# plot relative frequencies of valent configurations
polyancestry = readPolyAncestry("outstem_polyancestry.csv")
valentfreq = calvalentfreq(polyancestry)
plotvalentfreq(valentfreq)

# calculate the accuracy of parental phasing and ancestral inference
truefile = "true.csv"
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

# plot conditional probabilities
plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(polyancestry;
    truegeno = truegeno,
    fps = 0.5,
    # outfile = "outstem_condprob.gif",
)

# delete output files
rm.(filter(x->occursin("outstem",x), readdir()))
