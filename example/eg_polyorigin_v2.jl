# using Distributed
# nprocs() < 4 && addprocs(4-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using PolyOrigin

using Revise
using Distributions
using PolyOrigin
cd(@__DIR__)
pwd()

outstem = "outstem_4x"
truegenofile = outstem*"_truegeno.vcf.gz"
trueorigfile = outstem*"_trueorig.vcf.gz"
genofile = outstem*"_geno.vcf.gz"   
pedfile = outstem*"_ped.csv"
polygeno = formpolygeno(genofile, pedfile); 
fexch = 0.1
disturbjump = round(Int, fexch*size(polygeno.markermap[1],1)/2)
@info "" disturbjump
# disturbstd =  10 
PolyOrigin.disturbmap!(polygeno; 
    disturbjump,
    disturbstd = 0
);
0

polygeno2 = deepcopy(polygeno);
polyancestry, joinedancestry = polyorigin!(polygeno2; 
    # maxrun = 6, 
    refhapfile = truegenofile, 
    trueorigfile,    
    # doseerror = ErrorSetting(name="doseerror",hardthresh=0.5), 
    # ismarkerspecific = true, 
    # isdelmarker = true, 
    samplesize = 100,          
    outstem = outstem*"",         
);
0
@time acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
acc


# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^out",x), readdir()))


