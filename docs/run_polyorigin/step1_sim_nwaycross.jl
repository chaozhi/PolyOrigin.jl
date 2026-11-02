using Distributions
using PolyOriginSim
cd(@__DIR__)
pwd()

ploidy = 4
nfounder = 3
famsize = [20,1]

# simdata
outstem = string("sim_", ploidy, "x")
fhaplofile = PolyOriginSim.simfhaplo(;
    ploidy,
    nchr = 1, 
    chrlenls = 100, 
    nmarkerls = 200, 
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
    # gbsavgdepth = 20, 
    # gbsallelicbias = LogNormal(0,0.3),
    # gbsoverdispersion = Exponential(0.1),
    outstem
)


