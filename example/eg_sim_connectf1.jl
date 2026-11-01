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



# using Plots
# d = LogNormal(0,0.5)
# plot(x->pdf(d,x); xrange=(0,10))
# quantile(d, [0.005,0.025,0.25,0.5,0.75,0.975, 0.995])
# cdf(d,2) - cdf(d,0.5)
# cdf(d,10) - cdf(d,0.1)

# cdf(d,2.664) - cdf(d,0.375)
# cdf(d,3.625) - cdf(d,0.276)


