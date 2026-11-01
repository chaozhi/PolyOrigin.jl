using Distributed
nprocs() < 4 && addprocs(4-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using PolyOrigin

using Revise
using Distributions
using PolyOrigin
cd(@__DIR__)
pwd()

outstem = "outstem_4x"
truegenofile = outstem*"_truegeno.vcf.gz"
trueorigfile = outstem*"_trueorig.vcf.gz"
genofile = outstem*"_geno.vcf.gz"   
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

exportall(PolyOrigin)

polygeno = phasedgeno;
chr = 1
isvcf,isunphased_2dosage, missingstring = true, false, "NA"

target = "all"  # "all","founder","parent","offspring"
in(target, ["all","founder", "parent","offspring"]) || @error string("unknown target=",
    target, ", target must be in [all, founder,parent, offspring]")
if in(target, ["all","parent","founder"])
    if isnothing(polygeno.parentgeno)
        parentmtx = nothing
    else                        
        parentmtx = Matrix{Any}(polygeno.parentgeno[chr])
        if target == "founder"
            b = polygeno.polyped.parentinfo[!,:isfounder]
            parentmtx = parentmtx[:,b]
        end
        formatls = polygeno.markermap[chr][!,:parentformat]                        
        inverse_parsegenomtx!(parentmtx,formatls,isvcf,isunphased_2dosage, missingstring)
    end
else
    parentmtx = nothing
end
if in(target, ["all","offspring"])
    if isnothing(polygeno.offspringgeno)
        offmtx = nothing
        calledoffmtx = nothing
    else
        offmtx = Matrix{Any}( polygeno.offspringgeno[chr])
        formatvec = polygeno.markermap[chr][!,:offspringformat]
        inverse_parsegenomtx!(offmtx,formatvec,isvcf,isunphased_2dosage, missingstring)
    end
else
    offmtx = nothing
end



formatvec = polygeno.markermap[chr][!,:offspringformat]

offmtx = Matrix{Any}( polygeno.offspringgeno[chr])
calledoffgeno = biallele_gp2gt.(offmtx)





formatvec = polygeno.markermap[chr][!,:offspringformat]
inverse_parsegenomtx!(offmtx,formatvec,isvcf,isunphased_2dosage, missingstring)



#  ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat", 
#   "multiallele", "alleles", "quality", "filter", "info",
#   "doseerror","baseerror","allelicbias","overdispersion"]     
mapmtx = Matrix{Any}(polygeno.markermap[chr]);   
submapmtx = view(mapmtx,:, vcat([3],_col_1sterr:_col_1stsample-1))    
submapmtx .= [ismissing(i) ? i : round(i,digits=8) for i in submapmtx]
    
# mapmtx =  tovcfgenomtx!(mapmtx, parentmtx, offmtx)

offmtx
parentmtx

isnooff =false
offspringmtx = offmtx;

snp = 1
offspringformat = sort(unique([i[1:2] for i in split(mapmtx[snp,7],"&")]))
length(offspringformat) > 1 && @warn string("multiple offspringformat=",offspringformat) maxlog=10    

parentformat = sort(unique([i[1:2] for i in split(mapmtx[snp,6],"&")]))
genoformat = isnooff ? parentformat : sort(unique(vcat(parentformat,offspringformat)))        
length(genoformat) == 2 && genoformat[2] == "GT" && reverse!(genoformat) 

# reset genodata
nformat = length(genoformat)


formatindex = findfirst(==(only(offspringformat)),genoformat)

offgenols = view(offspringmtx, snp,:)
for j in eachindex(offgenols)
    offgenols[j] = join([i==formatindex ? offgenols[j] : "." for i in 1:nformat],":") 
end


