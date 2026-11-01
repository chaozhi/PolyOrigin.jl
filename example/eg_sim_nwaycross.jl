using Distributions
using PolyOriginSim
cd(@__DIR__)
pwd()

ploidy = 4
nfounder = 3
# famsize = [4,20]
famsize = 50

# simdata
outstem = string("outstem_", ploidy, "x")
fhaplofile = PolyOriginSim.simfhaplo(;
    ploidy,
    nchr = 1, 
    chrlenls = 100, 
    nmarkerls = 100, 
    founderls = string.("P",1:nfounder), 
    outstem, 
)
pedfile = PolyOriginSim.simnwaycross(;     
    ploidy,
    founderls = string.("P",1:nfounder),
    famsizels = famsize, 
    outstem,    
)

polyoriginsim(fhaplofile,pedfile; 
    isphysmap = true, 
    recomrate  = 1,     
    doseerror = 0.01,
    missfreq = 0.01,    
    isgbs = false, 
    gbsavgdepth = 10, 
    gbsallelicbias = LogNormal(0,0.3),
    gbsoverdispersion = Exponential(0.1),
    # gbsallelicbias = 10,     
    # gbsoverdispersion = 0.1, 
    outext = ".vcf",
    outstem
)


