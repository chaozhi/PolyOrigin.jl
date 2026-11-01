# using Distributed
# nprocs() < 4 && addprocs(4-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using PolyOrigin

using Revise
using PolyOrigin
cd(@__DIR__)
pwd()
exportall(PolyOrigin)

using Revise
outstem = "outstem_2x"
truegenofile = outstem*"_truegeno.vcf.gz"
trueorigfile = outstem*"_trueorig.vcf.gz"
genofile = outstem*"_geno.vcf.gz"   
genofile = string(outstem, "_out_imputed.vcf")
pedfile = outstem*"_ped.csv"    
@time polygeno = formpolygeno(genofile, pedfile; 
    isphysmap = true,     
);

polygeno.offspringgeno[1]
polygeno.parentgeno[1]


@time polyphase!(polygeno; 
    refhapfile = truegenofile, 
    # maxrun=1,     
    # isparallel = false, 
    chrpairing_phase = 44,      
    outstem = outstem*"_phase",
    outext = ".vcf", 
);
0

@time polyrefine!(polygeno; 
    outstem = outstem*"_refine",
);


phasedgeno = deepcopy(polygeno);

@time polyancestry, joinedancestry = polyreconstruct!(phasedgeno; 
    trueorigfile,
    samplesize = 10,
    outstem = outstem*"_reconstruct",
    outext = ".vcf", 
);
0


@time acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
acc


# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^outstem",x), readdir()))

