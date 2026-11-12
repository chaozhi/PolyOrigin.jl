using Distributed
nprocs() < 4 && addprocs(4-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using PolyOrigin

using Revise
using Distributions
using PolyOrigin
cd(@__DIR__)
pwd()


outstem = "outstem_4x"
truegenofile = outstem*"_truegeno.vcf"
trueorigfile = outstem*"_trueorig.vcf"
genofile = outstem*"_geno.vcf"   
pedfile = outstem*"_ped.csv"
phasedgeno, polyancestry, joinedancestry = polyorigin(genofile,pedfile;     
    refhapfile = truegenofile, 
    trueorigfile,        
    samplesize = 10,        
    outstem = outstem*"_out",
    outext = ".vcf"
);

0


@time acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
acc

# trueorig = readtruegeno!(trueorigfile,joinedancestry);
# plotcondprob(joinedancestry; offspring=2,trueorig)

# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^out",x), readdir()))


