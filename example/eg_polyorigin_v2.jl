# using Distributed
# nprocs() < 4 && addprocs(4-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using PolyOrigin

using PolyOrigin
cd(@__DIR__)
pwd()

@info string("Example datafiles: \n", PolyOrigin.get_datafiles())
genofile = PolyOrigin.get_datafile("potato4x_geno.vcf.gz")
pedfile = PolyOrigin.get_datafile("potato4x_ped.csv")
phasedgeno, polyancestry, joinedancestry = polyorigin(genofile,pedfile;     
    isphysmap = true, 
    recomrate = 1.25,  
    chrsubset = 1:1,    
    maxrun = 1, 
    isplot = true, 
    nplot_subpop = 1, 
    samplesize = 100, # default 1000 random draws for producing jointpolyancestry
)
    
PolyOrigin.plotcondprob(joinedancestry;                 
    individual = 1,        
)


# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^out",x), readdir()))


