using PolyOrigin
using Plots
cd(@__DIR__)
pwd()

# run polyorigin
genofile = "geno_disturbed.csv"
pedfile = "ped.csv"
outstem = "outstem"
@time polyancestry = polyorigin(genofile, pedfile;
    isrefine = true, 
    isrefineorder=true,
    isphysmap = false,
    isplot=true,
)

truefile = "true2.csv"
# trueancestryfile = "true.csv"
# PolyOrigin.generate_newtruefile(trueancestryfile, pedfile; outfile = truefile)

truegeno = readtruegeno!(truefile, polyancestry)
acc = calaccuracy!(truegeno, polyancestry)
println(acc)


plotcondprob(polyancestry, truegeno = truegeno, offspring = 1)
animcondprob(
    polyancestry,
    truegeno = truegeno,
    fps = 0.5,
    outfile=string(outstem,"_condprob.gif"),
)

polygeno = formpolygeno(genofile, pedfile)
fig = PolyOrigin.plotmapcomp(
    truegeno.truemap,
    polygeno.markermap,
    xlabel = "True poistion (cM)",
    ylabel = "Input position (cM)",
)
fig2 = PolyOrigin.plotmapcomp(
    truegeno.truemap,
    polyancestry.markermap,
    xlabel = "True poistion (cM)",
    ylabel = "Estimated position (cM)",
)
plot(fig, fig2)
savefig("mapcomp.png")

# delete output files
cd(@__DIR__)
rm.(filter(x->occursin("outstem",x), readdir()))
