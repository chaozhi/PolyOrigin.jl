using Distributed
nprocs() < 4 && addprocs(4-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using PolyOrigin

using Revise
using Distributions
using PolyOrigin

cd(@__DIR__)
pwd()

exportall(PolyOrigin)

# julia --project=. -e "using PolyOrigin; PolyOrigin.comonicon_install()"
# julia --project=. -e 'using PolyOrigin' -- --help

outstem = "outstem_4x"
truegenofile = outstem*"_truegeno.vcf"
trueorigfile = outstem*"_trueorig.vcf"
genofile = outstem*"_geno.vcf"   
pedfile = outstem*"_ped.csv"
phasedgeno, polyancestry, joinedancestry = polyorigin(genofile,pedfile;     
    refhapfile = truegenofile, 
    trueorigfile,        
    # maxrun = 3,     
    # chrpairing_phase = 22,
    # isrefine = false, 
    # isparentdoseerror = false, 
    # doseerror = ErrorSetting(name="doseerror",init=1e-2),
    # phaseerror_midparent = 0.01,
    samplesize = 10,        
    # ismarkerspecific = true, 
    # isparallel = false, 
    # byneighbor = true, 
    outstem = outstem*"_out",
    outext = ".vcf"
);

0


@time acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
acc


# exportall(PolyOrigin)
# trueorig = readtruegeno!(trueorigfile,joinedancestry);
# plotcondprob(joinedancestry; offspring=2,trueorig)

# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^out",x), readdir()))


