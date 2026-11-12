using Distributions
using PolyOriginSim
cd(@__DIR__)
pwd()

ploidy = 6
nfounder = 3
# famsize = [4,20]
famsize = 20

# simdata
outstem = string("outstem_", ploidy, "x")
fhaplofile = PolyOriginSim.simfhaplo(;
    ploidy,
    nchr = 1, 
    chrlenls = 100, 
    nmarkerls = 200, 
    founderls = string.("P",1:nfounder), 
    outstem, 
)
pedfile = PolyOriginSim.simconnectf1(;     
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
    gbsavgdepth = 40, 
    gbsallelicbias = LogNormal(0,0.5),
    gbsoverdispersion = Exponential(0.3),
    # gbsallelicbias = 0.5,     
    # gbsoverdispersion = 0.1, 
    # outext = ".csv", 
    outstem
)
