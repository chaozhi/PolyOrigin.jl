
"""
    calaccuracy!(truegenofile,trueorigfile, polyest,workdir=pwd(),verbose=true)

calculate inference accuracies according to the true values in `truegenofile` and `trueorigfile` . 
`joinedancestry` is modified with absolute phases and consistent genoprob.

# Positional arguments

`truegenofile::AbstractString`: a file for true genotypes

`trueorigfile::AbstractString`: a file for true ancestral origins 

`joinedancestry::JoinedPolyAncestry`: estimated by `polyorigin`. 

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for writing outfile.

`verbose::Bool=true`: true if print messages on console.

"""
function calaccuracy!(truegenofile::AbstractString, trueorigfile::AbstractString,
    joinedancestry::JoinedPolyAncestry; 
    workdir::AbstractString=pwd(),
    verbose::Bool=false)
    truegeno = readtruegeno(truegenofile,joinedancestry.polyped.ped; workdir,verbose)        
    alignchromosome!(truegeno,joinedancestry) 
    res_marker = calmarkersum(truegeno,joinedancestry)    
    alignmarker!(truegeno,joinedancestry) # afterward, same sets of markers    
    res_pacc = calparentacc(truegeno,joinedancestry)
    res_offacc = caloffspring_phaseacc(truegeno,joinedancestry)    
    trueorig = readtruegeno!(trueorigfile,joinedancestry; workdir,verbose)        
    res_orig = calviterbiacc(trueorig,joinedancestry)
    OrderedDict(vcat(res_marker, res_pacc,res_offacc, res_orig))
end

function calmarkersum(truegeno::PolyGeno, polyest::JoinedPolyAncestry)
    truesnps = reduce(vcat, [i[!,:marker] for i in truegeno.markermap])
    estsnps =  reduce(vcat,[i[!,:marker] for i in polyest.markermap])
    # nmarkertrue: #markers in truegeno
    # nmarkerdel: #markers in truegeno but not in magicest
    nmarkertrue = length(truesnps)
    nmarkerdel = length(setdiff(truesnps, estsnps))
    ["markersum" => (nmarkertrue=nmarkertrue, nmarkerdel=nmarkerdel)]
end


function readtruegeno!(truegenofile::AbstractString, polyest::Union{PolyGeno,PolyAncestry,JoinedPolyAncestry};
    workdir=pwd(), 
    verbose::Bool=false)    
    truegeno = readtruegeno(truegenofile,polyest.polyped.ped; workdir,verbose)        
    alignchromosome!(truegeno,polyest) 
    alignmarker!(truegeno,polyest) # afterward, same sets of markers
    truegeno
end

function readtruegeno(truegenofile::AbstractString, ped::Pedigree; 
    workdir=pwd(), 
    verbose::Bool=false)
    tempid = tempname(workdir)
    temppedfile = string(tempid,"_ped.csv")
    try     
        PolyOrigin.saveped(temppedfile,ped)
        PolyOrigin.check_infiles(truegenofile,temppedfile; workdir, verbose)            
        truegeno = PolyOrigin.formpolygeno(truegenofile,temppedfile; delmultiallelic=false, formatpriority=["GT"]);
        for chr in eachindex(truegeno.markermap)
            truegeno.parentgeno[chr] .= [parse.(Int,i) for i in truegeno.parentgeno[chr]]
            truegeno.offspringgeno[chr] .= [parse.(Int,i) for i in truegeno.offspringgeno[chr]]
        end 
        truegeno
    finally
        rm(temppedfile, force=true)
    end    
end


function alignchromosome!(truegeno::PolyGeno,polyest::Union{PolyGeno,PolyAncestry,JoinedPolyAncestry})
    truechrls = [i[1,:linkagegroup] for i in truegeno.markermap]
    estchrls =  [i[1,:linkagegroup] for i in polyest.markermap]
    truechrls == estchrls && return nothing
    if issubset(estchrls,truechrls)
        chrsubset = [findfirst(==(i),truechrls) for i in estchrls]
        getsubpolygeno!(truegeno; chrsubset)
    else
        error("estchrls = ", estchrls, " not subset of truechrls = ", truechrls)
    end
    nothing
end



function alignmarker!(truegeno::PolyGeno, polyest::Union{PolyGeno,PolyAncestry,JoinedPolyAncestry})
    for chr in eachindex(truegeno.markermap)                
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = polyest.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            rule = Dict(snpidls .=> 1:length(snpidls))
            ii = [get(rule, i, missing) for i in estsnpidls]
            ii2 = collect(skipmissing(ii))                
            truegeno.markermap[chr] = truegeno.markermap[chr][ii2,:]
            truegeno.parentgeno[chr] = truegeno.parentgeno[chr][ii2,:]
            truegeno.offspringgeno[chr] = truegeno.offspringgeno[chr][ii2,:]
            nonmissls = .!ismissing.(ii)
            if !all(nonmissls)                    
                polyest.markermap[chr] = polyest.markermap[chr][nonmissls, :]
                polyest.parentgeno[chr] = polyest.parentgeno[chr][nonmissls, :]
                if isa(polyest,PolyGeno)
                    polyest.offspringgeno[chr] = polyest.offspringgeno[chr][nonmissls, :]
                else
                    # viterbipath[chr][:,offspring]: most probably path for offspring at chr.
                    polyest.viterbipath[chr] = polyest.viterbipath[chr][nonmissls, :]
                    polyest.parentviterbipath[chr] = polyest.parentviterbipath[chr][nonmissls, :]
                    for condprob in [polyest.diploprob,polyest.genoprob,polyest.haploprob]
                        if !isnothing(condprob)
                            chrprob = condprob[chr]
                            for o in eachindex(chrprob)
                                chrprob[o] = chrprob[o][nonmissls,:]
                            end   
                        end
                    end                       
                end
            end            
        end       
    end 
    truegeno, polyest
end


function calparentacc(truegeno::PolyGeno,polyest::Union{PolyGeno, JoinedPolyAncestry})    
    df = PolyOrigin.calphaseerr(truegeno.parentgeno,polyest)
    isfls = polyest.polyped.parentinfo[!,:isfounder]
    nf = length(isfls)
    accls = df[1:nf,:allchr]
    keyls = keys(first(accls))
    acc_f = sum([[values(i)...] for i in accls[isfls]])
    res = ["founder_phasing" => (;zip(keyls,acc_f)...)]
    if !all(isfls)
        acc_midp = sum([[values(i)...] for i in accls[.!isfls]])
        push!(res, "midparent_phasing" => (;zip(keyls,acc_midp)...))
    end
    res
end

function calviterbiacc(truegeno::PolyGeno, polyest::JoinedPolyAncestry)
    acc = calviterbiacc(truegeno.offspringgeno, polyest.viterbipath; ismidparent=false)
    ismidpls = .!polyest.polyped.parentinfo[!,:isfounder]
    if any(ismidpls)
        truegenols  = [chrgeno[:,ismidpls] for chrgeno in truegeno.parentgeno]
        viterbigenols = [chrviterbi[:,ismidpls] for chrviterbi in polyest.parentviterbipath]
        parentacc = calviterbiacc(truegenols, viterbigenols; ismidparent=true)
        acc = vcat(parentacc,acc)
    end
    acc
end

function calviterbiacc(truegenols::AbstractVector,viterbigenols::AbstractVector;
    ismidparent::Bool=true)    
    accls = [begin 
        chrtruehap = truegenols[chr]
        chresthap= viterbigenols[chr] 
        accls = [calphaseerr(chrtruehap[:,p],chresthap[:,p]) for p in 1:size(chrtruehap,2)]
        keyls = keys(first(accls))
        valls = sum([values(j)...] for j in accls)
        (;zip(keyls,valls)...)
    end for chr in eachindex(truegenols,viterbigenols)]    
    resid = ismidparent ? "midparent_origin" : "offspring_origin"
    keyls = Symbol.(string.(keys(first(accls))))
    valls = sum([values(j)...] for j in accls)
    [resid => (;zip(keyls,valls)...)]
end

function caloffspring_phaseacc(truegeno::PolyGeno, polyest::JoinedPolyAncestry)
    truegenols = truegeno.offspringgeno
    estgenols = polyest.offspringgeno
    accls = [begin 
        chrtruehap = truegenols[chr]
        chresthap= estgenols[chr] 
        chresthap= [parse.(Int,i) for i in estgenols[chr]]
        accls = [calphaseerr(chrtruehap[:,p],chresthap[:,p]) for p in 1:size(chrtruehap,2)]    
        keyls = keys(first(accls))
        valls = sum([values(j)...] for j in accls)
        (;zip(keyls,valls)...)
    end for chr in eachindex(truegenols,estgenols)]    
    keyls = Symbol.(string.(keys(first(accls))))
    valls = sum([values(j)...] for j in accls)
    ["offspring_phasing" => (;zip(keyls,valls)...)]
end






