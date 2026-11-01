using Distributed
nprocs() < 6 && addprocs(6-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using PolyOrigin

using Revise
using Distributions
using PolyOrigin
cd(@__DIR__)
pwd()

outstem = "outstem_2x"
truegenofile = outstem*"_truegeno.vcf.gz"
trueorigfile = outstem*"_trueorig.vcf.gz"
genofile = outstem*"_geno.vcf.gz"   
pedfile = outstem*"_ped.csv"
isinferrlsls = [[false,false,false,false],[true, false,false,false],[true,true,false,false],[true,true,true,false],[true,true,true,true]]
isinferls = isinferrlsls[5]
# refhapfile = "outstem_2x_M5_parentphased.vcf.gz"
# isfile(refhapfile)
polyancestry, joinedancestry = polyorigin(genofile, pedfile;
    isphysmap = true,    
    recomrate = 1, # cm/Mbp
    # refhapfile,     
    maxrun = 6, 
    # refhapfile = truegenofile, 
    # trueorigfile, 
    # isdelmarker = true, 
    doseerror = ErrorSetting(name="doseerror",isinfer=isinferls[1]),            
    allelicbias = ErrorSetting(name="allelicbias",isinfer=isinferls[2]),            
    overdispersion = ErrorSetting(name="overdispersion",isinfer=isinferls[3]),            
    baseerror = ErrorSetting(name="baseerror",isinfer=isinferls[4]),            
    samplesize = 100,          
    outstem = outstem*"",
);
0
@time acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
acc


# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^out",x), readdir()))

