

using PolyOrigin
cd(@__DIR__)
pwd()
exportall(PolyOrigin)


outstem = "outstem_4x"
truegenofile = outstem*"_truegeno.vcf.gz"
trueorigfile = outstem*"_trueorig.vcf.gz"
genofile = outstem*"_geno.vcf.gz"   
pedfile = outstem*"_ped.csv"    
@time polygeno = formpolygeno(genofile, pedfile; 
    isphysmap = true,     
);

chr = 1
chrid = polygeno.markermap[chr][1,:linkagegroup]
chrparentgeno = polygeno.parentgeno[chr] # keep phase genotypes in format of vector{Int8}  instead of vector{String}    
chrparentformat = polygeno.markermap[chr][!,:parentformat]         
chroffgeno = get_chroffgeno(polygeno,chr)      
chrmultiallele = polygeno.markermap[chr][!,:multiallele]
chroffformat = polygeno.markermap[chr][!,:offspringformat]        



chrpairing_phase = 44
bindis = 1e-3

polyped = polygeno.polyped    
priorspace = get_priorstatespace(polyped,chrpairing_phase)    
maxvalent = isnothing(chrpairing_phase) ? 4 : max(digits(chrpairing_phase)...)    
priorprocess = get_priorprocess(polygeno,chr,maxvalent; bindis)       

phaploset, phaploweight = init_phaploprior(chrparentgeno, chrparentformat, chrmultiallele,polyped;
    doseerror = 0.01,
    baseerror = 0.001,
    allelicbias = 1, 
    overdispersion = 0,
)
phaploindex= init_phaploindex(phaploweight)    
fixedparentls = findall([unique(length.(i)) == [1] for i in phaploset])    
popidls = collect(skipmissing([issubset(priorspace[i]["parentindex"],fixedparentls) ? missing : i for i in keys(priorspace)]))
bvpair = init_bvpair(priorspace,polyped)
siblogl = init_siblogl(bvpair)    



chrparentgeno = polygeno.parentgeno[chr] 
chrparentformat = polygeno.markermap[chr][!,:parentformat]         

chrphaplo = PolyOrigin.get_phaplo(phaploindex,phaploset)
chrparent_hiddenunphased = unphased_genocode_2_index.(join.(sort.(chrphaplo),"/"))
logl_parent = cal_logl_parent(chrparentgeno,chrparent_hiddenunphased,
    chrparentformat,chrmultiallele,priorprocess,polyped;
    # doseerror=doseerror2, baseerror=baseerror2, allelicbias=allelicbias2, overdispersion=overdispersion2       
    doseerror = 0.01,
    baseerror = 0.001,
    allelicbias = 1, 
    overdispersion = 0,
)
