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


doseerror = 0.01
baseerror = 0.001
chrpairing = 44
bindis = 0.005
chr = 1
phasedgeno = deepcopy(polygeno);



nmarker = size(phasedgeno.markermap[chr],1)
doseerrorls = reset_errorls(phasedgeno.markermap[chr][!,:doseerror],doseerror)
baseerrorls = reset_errorls(phasedgeno.markermap[chr][!,:baseerror],baseerror)    
allelicbiasls = reset_errorls(phasedgeno.markermap[chr][!,:allelicbias],1.0)    
overdispersionls = reset_errorls(phasedgeno.markermap[chr][!,:overdispersion],0.0)
doseerrorls = allequal(doseerrorls) ? doseerrorls[1] : doseerrorls
baseerrorls = allequal(baseerrorls) ? baseerrorls[1] : baseerrorls
allelicbiasls  = allequal(allelicbiasls) ? allelicbiasls[1] : allelicbiasls 
overdispersionls = allequal(overdispersionls) ? overdispersionls[1] : overdispersionls  

priorspace = get_priorstatespace(phasedgeno.polyped,chrpairing)
polyped = phasedgeno.polyped
# dict2group = getdict2group(priorspace,polyped)    
dict2group = nothing # NOT transforming diploprob into genoprob
maxvalent = max(digits(chrpairing)...)
priorprocess = get_priorprocess(phasedgeno,chr,maxvalent; bindis)

# idlen=max([length(i[1,:linkagegroup]) for i=phasedgeno.markermap]...)
chrid = phasedgeno.markermap[chr][1,:linkagegroup]        
chroffgeno = get_chroffgeno(phasedgeno,chr)      
chroffformat = phasedgeno.markermap[chr][!,:offspringformat]
chrmultiallele = phasedgeno.markermap[chr][!,:multiallele]    
chrparentgeno= getchrparentgeno(phasedgeno.parentgeno[chr])
derivedgeno = get_derivedgeno(chrparentgeno,chrmultiallele, priorspace)
targetls = ["doseerror"]
if !all(occursin.("GT", unique(chroffformat))) 
    append!(targetls, ["baseerror", "allelicbias","overdispersion"])
end

targetls

using DataFrames
# decoding for leaf-offspring    
itmax=20
res=Vector(undef,4)
correctiondf = DataFrame()    
isparallelind = false
jltempdir = mktempdir()
verbose = true
io = nothing

startt = time()        
msg = string("chr=", chrid, ", #markers=",nmarker)
msg *= get_msg_err(targetls, priorprocess; 
    doseerror = doseerrorls, baseerror = baseerrorls, 
    allelicbias = allelicbiasls, 
    overdispersion = overdispersionls)        
chrdecode = chrhmmdecode(chroffgeno, derivedgeno,priorspace, priorprocess, 
    polyped, chroffformat, chrmultiallele; dict2group, 
    doseerror = doseerrorls, baseerror = baseerrorls, 
    allelicbias = allelicbiasls, overdispersion = overdispersionls,
    ismidparent =false, isinclviterbi = true, isparallelind, jltempdir, verbose)
for k in 1:4
    res[k] = [i[k] for i=chrdecode]
end
isaccept = false
ncorrect = 0
msg_correct = ", iscorrectparent=false"
phasedgeno.offspringgeno[chr]
isimputefirst = true

phasedgeno.parentgeno[chr] = Matrix{Any}([ismissing(i) ? i : string.(i) for i in chrparentgeno])
phasedgeno.offspringgeno[chr] .= cal_condprob2doseprob(chrdecode, chrparentgeno, priorspace,priorprocess)                        
chroffformat .= "GP"    
if isimputefirst
    # baseerror/allelicbias/od are not requred for GP format 
    chrdecode = chrhmmdecode(phasedgeno.offspringgeno[chr], derivedgeno,priorspace, priorprocess, 
        polyped, chroffformat, chrmultiallele; dict2group, 
        doseerror=1e-3, ismidparent =false, isinclviterbi = true, isparallelind, jltempdir, verbose)                            
    for k in 1:4
        res[k] = [i[k] for i=chrdecode]
    end
end            
msgtused = string(", tused=",round(time() - startt, digits=1),"s" )
msg *= string(msgtused, msg_correct)            
printconsole(io,verbose,msg)            

push!(res, correctiondf)


# decoding for intermidiate parent
nmidparent = sum(.!polyped.parentinfo[!,:isfounder])
if nmidparent == 0
    append!(res, [nothing for _ in 1:4])
else
    startt = time()
    chrformat = repeat(["GT_phased"],size(chrparentgeno,1))        
    chrdecode = chrhmmdecode(chrparentgeno, derivedgeno,priorspace, priorprocess, 
            polyped, chrformat, chrmultiallele; 
            # baseerror/allelicbias/od are not requred for GT_phased
            # doseerror is used to for phasing error rate of GT_phased; too small phased error rates might result in NaN or Inf in posterior decoding                 
            doseerror = 0.05,                                 
            isparallelind = isparallelind && (nmidparent > nprocs() > 1), 
            ismidparent =true, isinclviterbi = true, jltempdir, verbose)    
    for k in 1:4
        push!(res, [isnothing(i) ? i : i[k] for i=chrdecode]) # parent  valentprob
    end        
    msgtused = string(", tused=",round(time() - startt, digits=1),"s" )
    msg = string("chr=", chrid, ", #markers=",nmarker, ", midparent decoding", msgtused)
    printconsole(io, verbose, msg)
end

