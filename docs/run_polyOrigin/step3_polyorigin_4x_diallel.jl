using PolyOrigin
using Plots
cd(@__DIR__)
pwd()

# run polyorigin
genofile = "geno_disturbed.csv"
pedfile = "ped.csv"
outstem = "outstem"
@time polyancestry = polyOrigin(genofile, pedfile;
    refinemap=true,
    refineorder=true,
    isphysmap = false,
    # isplot=true,
)

truefile = "true.csv"
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)


plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(
    polyancestry,
    truegeno = truegeno,
    fps = 0.5,
    outfile=string(outstem,"_condprob.gif"),
)

polygeno = readPolyGeno(genofile, pedfile)
fig = plotMapComp(
    truegeno.truemap,
    polygeno.markermap,
    xlabel = "True poistion (cM)",
    ylabel = "Input position (cM)",
)
fig2 = plotMapComp(
    truegeno.truemap,
    polyancestry.markermap,
    xlabel = "True poistion (cM)",
    ylabel = "Estimated position (cM)",
)
plot(fig, fig2)

# delete output files
rm.(filter(x->occursin(outstem,x), readdir()))
