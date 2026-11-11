using Distributed
nprocs() < 4 && addprocs(4-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using PolyOrigin

using Revise
using Distributions
using PolyOrigin
cd(@__DIR__)
pwd()




PolyOrigin.get_datadir()
readdir(PolyOrigin.get_datadir())

pedfile = PolyOrigin.get_datafile("potato2x_ped.csv")

root_dir = artifact"polyorigin2_datasets" 
file_path = joinpath(root_dir, "polyorigin2_datasets", filename)
isfile(file_path)

using CSV, DataFrames
df = CSV.read(file_path, DataFrame)
    
readdir(data_dir)
# Target only the file needed for this specific example
genofile = joinpath(data_dir, "example_A_data.csv") 
df = CSV.read(data_path, DataFrame)


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


