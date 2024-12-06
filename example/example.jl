# using Pkg
# Pkg.develop(path=abspath(@__DIR__, ".."))

using Revise
using PolyOrigin
cd(@__DIR__)
pwd()

# run polyorigin
genofile = "geno.csv"
pedfile = "ped.csv"
outstem = "example_output"
@time polyancestry = polyOrigin(genofile, pedfile;            
    isplot = true, 
    nplot_subpop = 3, 
    outstem,
);

# calculate the accuracy of parental phasing and ancestral inference
polyancestry = readPolyAncestry(outstem*"_polyancestry.csv")
truefile = "true.csv"
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

# plot conditional probabilities
plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(polyancestry;
    truegeno = truegeno,
    fps = 0.5,
    outfile = string(outstem,"_condprob.gif"),
)

# delete output files
cd(@__DIR__)
outstem = "example_output"
rm(outstem*"_plots",recursive=true)
rm.(filter(x->occursin(outstem,x), readdir()))