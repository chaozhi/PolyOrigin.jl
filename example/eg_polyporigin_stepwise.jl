# using Distributed
# nprocs() < 4 && addprocs(4-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using PolyOrigin

using PolyOrigin
cd(@__DIR__)
pwd()
exportall(PolyOrigin)

using Revise
outstem = "outstem_4x"
truegenofile = outstem*"_truegeno.vcf.gz"
trueorigfile = outstem*"_trueorig.vcf.gz"
genofile = outstem*"_geno.vcf.gz"   
pedfile = outstem*"_ped.csv"    
@time polygeno = formpolygeno(genofile, pedfile; 
    isphysmap = true,     
);
@time polyphase!(polygeno; 
    refhapfile = truegenofile, 
    maxrun=1,     
    isparallel = false, 
    chrpairing_phase = 44,      
    outstem = outstem*"_phase",
);
0


chrpairing = 44
chr = 1
isrefineerror = true
byneighbor = false
doseerror=ErrorSetting(name="doseerror")   
allelicbias=ErrorSetting(name="allelicbias")
overdispersion=ErrorSetting(name="overdispersion")
baseerror=ErrorSetting(name="baseerror")     
io = nothing
verbose = true
bindis = 1e-3
chrpairing

polyped = polygeno.polyped
maxploidy_off = maximum(polyped.offspringinfo[!,:ploidy])
priorspace = get_priorstatespace(polyped,chrpairing)



idlen=max([length(i[1,:linkagegroup]) for i=polygeno.markermap]...)
chrid=lpad(polygeno.markermap[chr][1,:linkagegroup],idlen)
chrparentgeno= getchrparentgeno(polygeno.parentgeno[chr])
chrmultiallele = polygeno.markermap[chr][!,:multiallele]
derivedgeno = get_derivedgeno(chrparentgeno,chrmultiallele, priorspace)
chroffgeno = get_chroffgeno(polygeno, chr)

chroffformat = polygeno.markermap[chr][!,:offspringformat]        
targetls = doseerror.isinfer ?  ["doseerror"] : String[] 
if !all(occursin.(r"^GT_",chroffformat))         
    baseerror.isinfer && push!(targetls,"baseerror")
    allelicbias.isinfer && push!(targetls,"allelicbias")
    overdispersion.isinfer && push!(targetls,"overdispersion")
end   
isrefineerror = isrefineerror && !isempty(targetls)
snporder=collect(1:size(chroffgeno,1))    
ismissdoseerror = all(ismissing.(polygeno.markermap[chr][!,:doseerror]))
ismissbaseerror = all(ismissing.(polygeno.markermap[chr][!,:baseerror]))
ismissbias = all(ismissing.(polygeno.markermap[chr][!,:allelicbias]))
ismissod = all(ismissing.(polygeno.markermap[chr][!,:overdispersion]))
doseerrorls = reset_errorls(polygeno.markermap[chr][!,:doseerror],doseerror.init)
baseerrorls = reset_errorls(polygeno.markermap[chr][!,:baseerror],baseerror.init)    
allelicbiasls = reset_errorls(polygeno.markermap[chr][!,:allelicbias],allelicbias.init)    
overdispersionls = reset_errorls(polygeno.markermap[chr][!,:overdispersion],overdispersion.init)        

using StatsBase
using Distributions

if ismissdoseerror || !in("doseerror",targetls)
    doseerrorprior = doseerror.prior    
else
    avg = max(0.01,mean(doseerrorls))
    doseerrorprior = Beta(1,round(1/avg - 1,digits=1))         
    msg = string("chr=",chrid,", reset doseerror.prior = ", doseerrorprior)
    printconsole(io,verbose,msg)        
end
if ismissbaseerror || !in("baseerror",targetls) 
    baseerrorprior = baseerror.prior    
else        
    avg = max(0.005,mean(baseerrorls))
    baseerrorprior = Beta(1,round(1/avg - 1,sigdigits=2))         
    msg = string("chr=",chrid,", reset baseerror.prior = ", baseerrorprior)
    printconsole(io,verbose,msg)        
end    
if ismissbias || !in("allelicbias",targetls)  
    biasprior = allelicbias.prior    
else
    priorlogstd = round(std(log.(rand(allelicbias.prior,10000))),digits = 2)                 
    biasstd = round(std(log.(allelicbiasls)),digits = 2)
    if priorlogstd < biasstd
        biasprior = LogNormal(0, max(biasstd,0.1))
        msg = string("chr=",chrid,", reset allelicbias.prior = ", biasprior)
        printconsole(io,verbose,msg)       
    else
        biasprior = allelicbias.prior                 
        msg = string("chr=",chrid,", keep allelicbias.prior = ", biasprior)
        printconsole(io,verbose,msg)       
    end
end
if ismissod || !in("overdispersion",targetls) 
    overdispersionprior = overdispersion.prior    
else        
    avg = max(0.05,round(mean(overdispersionls),digits=2))
    overdispersionprior = Exponential(avg) 
    msg = string("chr=",chrid,", reset overdispersion.prior = ", overdispersionprior)
    printconsole(io,verbose,msg)        
end
errorprior = (doseerror=doseerrorprior, baseerror = baseerrorprior, 
        allelicbias = biasprior, overdispersion = overdispersionprior)
# errorprior = (doseerror=doseerror.prior, baseerror = baseerror.prior, 
#     allelicbias = allelicbias.prior, overdispersion = overdispersion.prior)
errorsoftthresh = (doseerror = doseerror.softthresh, baseerror = baseerror.softthresh, 
    allelicbias = allelicbias.softthresh, overdispersion = overdispersion.softthresh)
errorhardthresh = (doseerror = doseerror.hardthresh, baseerror = baseerror.hardthresh, 
    allelicbias = allelicbias.hardthresh, overdispersion = overdispersion.hardthresh)    

phaseerror_midparent = 0.01    

maxvalent = max(digits(chrpairing)...)
priorprocess = get_priorprocess(polygeno,chr,maxvalent; bindis)
bvpair = init_bvpair(priorspace,polyped)    
nowlogl =refinebvpair!(bvpair, chrparentgeno,chroffgeno,chroffformat, chrmultiallele, 
    priorspace,priorprocess, polyped; 
    phaseerror_midparent,byneighbor=false, snporder, isrand = false,
    doseerrorls, baseerrorls, allelicbiasls, overdispersionls, 
)    




# using JLD2
# res = load("temp.jld2");

# argls = res["argls"]
# optls = res["optls"]
# chrparentgeno2 = res["chrparentgeno"]
# chrparentgeno == chrparentgeno2



# b = chrparentgeno .!= chrparentgeno2
# [polygeno.parentgeno[chr][b] chrparentgeno[b] chrparentgeno2[b]]

# refhapfile = truegenofile
# refparenthap= readparenthaplo(refhapfile,polyped.parentinfo,polygeno.markermap);
# homologmapls,newestparenthap = calabsolutehap(polygeno.parentgeno,refparenthap,polyped)    
# chrparentgeno3 = newestparenthap[1]
# homologmapls[1]

# nowlogl

# nowlogl3 =refinebvpair!(bvpair, chrparentgeno3,chroffgeno,chroffformat, chrmultiallele, 
#     priorspace,priorprocess, polyped; 
#     phaseerror_midparent,byneighbor=false, snporder, isrand = false,
#     doseerrorls, baseerrorls, allelicbiasls, overdispersionls, 
# )    


# priorspace = get_priorstatespace(polyped,44)

# priorspace
# priorspace["fam1"]

# priorspace["fam1"]["valentneighbor"]

# priorspace["fam1"]["valentkey"]

# priorspace["fam1"]["state"]



# @time polyrefine!(polygeno; 
#     outstem = outstem*"_refine",
# );

# @time polyancestry, joinedancestry = polyreconstruct!(polygeno; 
#     trueorigfile,
#     outstem = outstem*"_reconstruct",
#     outext = ".vcf", 
# );

# @time acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
# acc


# delete output files
# cd(@__DIR__)
# rm.(filter(x->!isdir(x) && occursin(r"^outstem",x), readdir()))

