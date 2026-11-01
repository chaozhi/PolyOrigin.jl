
function parentcorrect!(phaploindex::AbstractVector,    
    phaploset::AbstractVector,
    phaploweight::AbstractVector,            
    phaploposterior::AbstractVector, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict, 
    polyped::PolyPed;         
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},    
    overdispersion::Union{Real,AbstractVector},
    snporder::Union{Nothing,AbstractVector}=nothing,            
    bvpair::Union{Nothing, NamedTuple}=nothing, 
    correctthreshold::Real=0.1, 
    callthreshold::Real=0.7,
    minerrdiff::Integer=3,
    alwayscorrect::Bool=false)
    isnothing(snporder) && (snporder = 1:size(chroffgeno,1))    
    chrphaplo = get_phaplo(phaploindex, phaploset)    
    # badchanges: ["markerindex", "parentindex","weight", "new_nerr","new_genotype", "old_nerr", "old_genotype"]    
    deltlogl, badchanges = parentcorrect(chrphaplo, chroffgeno, chroffformat, chrmultiallele, 
        priorspace, priorprocess, polyped; 
        phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, snporder,  bvpair, 
        callthreshold, minerrdiff, correctthreshold,alwayscorrect
    )          
    isempty(badchanges) && return (0, badchanges)
    logl_parentdata = sum(cal_logl_parentdata(phaploindex, phaploweight))
    # 
    newphaploindex = deepcopy(phaploindex)    
    newphaploset = deepcopy(phaploset)    
    newphaploweight = deepcopy(phaploweight)    
    newphaploposterior = deepcopy(phaploposterior)    
    for row in eachrow(badchanges) 
        marker, p, newg_str = row[[:markerindex,:parentindex,:new_genotype]]
        newg = parse.(Int,split(newg_str,"|"))
        hh = newphaploset[p][marker]
        ww = newphaploweight[p][marker]
        pp = newphaploposterior[p][marker]
        hindex = findfirst(==(newg),hh)        
        if isnothing(hindex)
            # println("parentcorrect: newg=",newg,", hh=",hh, ",ww=",ww)
            push!(hh,newg)
            push!(ww, maximum(ww))
            push!(pp, maximum(pp))
            ww .= 1/length(ww) 
            pp .= 1/length(pp) 
            newphaploindex[p][marker] = length(hh)
        else
            newphaploindex[p][marker] = hindex
            # ww .= 1/length(ww)                 
        end        
    end    
    deltlogprior = sum(cal_logl_parentdata(newphaploindex, newphaploweight)) - logl_parentdata
    # @info "" deltlogl deltlogprior badchanges    
    deltlogp = deltlogl + deltlogprior         
    if deltlogp > 0 
        phaploindex .= newphaploindex                
        phaploset .= newphaploset                
        phaploweight .= newphaploweight                
        phaploposterior .= newphaploposterior                
    end
    (deltlogp, badchanges)
end

function parentcorrect(chrphaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, chrmultiallele::AbstractVector,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict, 
    polyped::PolyPed;         
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},    
    overdispersion::Union{Real,AbstractVector},
    snporder::Union{Nothing,AbstractVector}=nothing,        
    bvpair::Union{Nothing, NamedTuple}=nothing,     
    correctthreshold::Real=0.1, 
    callthreshold::Real=0.7,
    minerrdiff::Integer=3,
    alwayscorrect::Bool=false)
    isnothing(snporder) && (snporder =1:size(chroffgeno,1))    
    derivedgeno = get_derivedgeno(chrphaplo,chrmultiallele, priorspace)    
    if isnothing(bvpair)        
        bvpairprop = get_prop_bvpair(priorspace,polyped)
    else        
        bvpairprop = (offspring = [[i] for i in bvpair.offspring], parent = [[i] for i in bvpair.parent])
    end
    # offspring decoding
    dict2group = getdict2group(priorspace, polyped)    
    chrdecode = chrhmmdecode(chroffgeno, derivedgeno,priorspace, priorprocess, 
        polyped, chroffformat, chrmultiallele; dict2group, bvpairprop,
        doseerror, baseerror, allelicbias, overdispersion, isinclviterbi = false, 
        ismidparent =false, isparallelind=false, jltempdir="",verbose=false)  # jltempdir is not used for non-parallel    
    # midparent decoding
    chrpformat = ["GT_phased" for _ in eachindex(chroffformat)]
    chrdecode_midparent = chrhmmdecode(chrphaplo, derivedgeno,priorspace, priorprocess, 
        polyped, chrpformat, chrmultiallele; dict2group, bvpairprop,
        doseerror, baseerror, allelicbias, overdispersion, isinclviterbi = false, 
        ismidparent = true, isparallelind=false, jltempdir="",verbose=false)  # jltempdir is not used for non-parallel    
    # parental correction
    offploidyls = polyped.offspringinfo[!,:ploidy]
    chrcalledgeno = get_chrcalledgeno(chroffgeno, chroffformat, offploidyls; 
        doseerror, baseerror, allelicbias, overdispersion, callthreshold) 
    changedf = parentcorrect_decode(chrdecode, chrphaplo, chrcalledgeno, chrmultiallele, 
        priorspace, priorprocess, polyped;
        snporder, correctthreshold, callthreshold, minerrdiff
    )        
    changedf = reduce(vcat,[df[[rand(1:size(df,1))],:] for df in groupby(changedf,:markerindex)])
    isempty(changedf) && return (0, changedf)        
    alwayscorrect && return (0, changedf)

    # calculate logl before parent correction 
    logl_decode_off= sum(logsumexp(i[1][:,2]) for i in chrdecode)    
    logl_decode_midparent= sum(isnothing(i) ? 0 : logsumexp(i[1][:,2]) for i in chrdecode_midparent)
    logl_old = logl_decode_off + logl_decode_midparent
   
    # calculate logl after parent correction
    newchrphaplo = copy.(chrphaplo)
    correct_chrphalo!(newchrphaplo, changedf)    
    newderivedgeno = get_derivedgeno(newchrphaplo,chrmultiallele, priorspace)    
    if isnothing(bvpair)        
        loglls = calmarglogl_bvpairprop(newchrphaplo, chroffgeno,chroffformat, chrmultiallele, newderivedgeno, 
            priorspace, priorprocess, bvpairprop, polyped; 
            phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, snporder)    
        logl = sum(sum(ismissing(j) ? 0 : logsumexp(j) for j in i) for i in loglls)            
    else
        loglls = calmarglogl_bvpair(newchrphaplo, chroffgeno,chroffformat, chrmultiallele, newderivedgeno, 
            priorspace, priorprocess, bvpair, polyped; 
            phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, snporder)       
        logl = sum_margloglls(loglls)
    end    
    deltlogl = logl - logl_old
    deltlogl, changedf
end    

function correct_chrphalo!(chrphaplo::AbstractMatrix, changedf::AbstractDataFrame)
    for row = eachrow(changedf) 
        marker, p, newg = row[[:markerindex,:parentindex,:new_genotype]]
        newg2 = parse.(Int8,split(newg,"|"))
        chrphaplo[marker, p] = newg2 # Do NOT use .=newg2, otherwise it would change pointing to  phaploset
    end
    nothing
end

function get_chrcalledgeno(chroffgeno::AbstractMatrix, 
    chroffformat::AbstractVector, offploidyls::AbstractVector;    
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},    
    overdispersion::Union{Real,AbstractVector},
    callthreshold::Real=0.7) 
    res = copy(chroffgeno)
    markerls = findall(chroffformat .== "AD")
    if !isempty(markerls)
        b = [isa(i,Real) for i in [doseerror, baseerror,allelicbias]]
        if all(b)
            for j in eachindex(offploidyls)
                res[markerls,j] .= genocall_AD.(res[markerls,j], offploidyls[j]; 
                    doseerror, baseerror, allelicbias,callthreshold)
            end
        else
            nmarker = size(chroffgeno,1)            
            doseerrorls = isa(doseerror,Real) ? doseerror .* ones(nmarker) : doseerror
            baseerrorls = isa(baseerror,Real) ? baseerror .* ones(nmarker) : baseerror
            allelicbiasls = isa(allelicbias,Real) ? allelicbias .* ones(nmarker) : allelicbias
            overdispersionls = isa(overdispersion,Real) ? overdispersion .* ones(nmarker) : overdispersion
            for j in eachindex(offploidyls)
                for m in markerls
                    res[m,j] = genocall_AD(res[m,j], offploidyls[j];                         
                        doseerror = doseerrorls[m], baseerror = baseerrorls[m], 
                        allelicbias = allelicbiasls[m], overdispersion = overdispersionls[m],callthreshold)
                end
            end
        end
    end
    res
end

function parentcorrect_decode(chrdecode::AbstractVector, 
    chrphaplo::AbstractMatrix,
    chrcalledgeno::AbstractMatrix,    
    chrmultiallele::AbstractVector, 
    priorspace::AbstractDict,
    priorprocess::AbstractDict, 
    polyped::PolyPed;             
    snporder::Union{Nothing,AbstractVector}=nothing,        
    correctthreshold::Real=0.1,         
    callthreshold::Real=0.7,
    minerrdiff::Integer=3)
    # changedf: ["markerindex", "parentindex","weight", "new_nerr","new_genotype", "old_nerr", "old_genotype"]
    isnothing(snporder) && (snporder =1:size(chrcalledgeno,1))    
    doseprob = cal_condprob2doseprob(chrdecode, chrphaplo, priorspace,priorprocess;  snporder)
    bestdose = cal_doseprob2calldose(doseprob; callthreshold)
    badmarkers = get_badmarkers(chrcalledgeno,bestdose, priorspace,polyped; correctthreshold)        
    changedf = get_badchanges!(chrdecode, badmarkers, 
        chrphaplo, chrcalledgeno, chrmultiallele, priorspace,priorprocess, polyped; 
        minerrdiff, callthreshold,snporder)                
    changedf
end


function unphasedcalling(unphasedprob::Union{Missing,AbstractVector}; callthreshold::Real=0.7)
    ismissing(unphasedprob) && return  missing
    maxp, index = findmax(unphasedprob)
    maxp >= callthreshold ? index : missing        
end

function cal_doseprob2calldose(doseprob::AbstractMatrix; callthreshold::Real=0.7)
    unphasedcalling.(doseprob; callthreshold)    
end


function get_badmarkers(chrcalledgeno::AbstractMatrix, 
    bestdose::AbstractMatrix, 
    priorspace::AbstractDict, 
    polyped; 
    correctthreshold::Real=0.1)
    dosemismatch = map((x,y)->ismissing(x) || ismissing(y) ? missing : x!=y,chrcalledgeno,bestdose)            
    # mismatchrates = [mean(skipmissing(i)) for i in eachrow(dosemismatch)]
    # mismatch_markers =  findall(mismatchrates .> 0.1)
    nparent = size(polyped.parentinfo,1)
    res=[]
    for popid = keys(priorspace)
        poppriorspace = priorspace[popid]        
        offls = poppriorspace["offspringindex"]
        popmismatch = dosemismatch[:,offls]
        popmismatch2=[collect(skipmissing(i)) for i= eachrow(popmismatch)]
        nlen=length.(popmismatch2)
        nerr = [isempty(i) ? 0 : sum(i) for i=popmismatch2]
        # nerr = sum.(popmismatch2)
        pos = findall(nerr .>= 1)
        snps = pos
        keepat!(nerr, pos)
        keepat!(nlen, pos)
        pp = zeros(Int, nparent)
        pp[poppriorspace["parentindex"]] .= 1
        # println("popid=",popid,", pls=",poppriorspace["parentindex"])
        for i in eachindex(snps, nerr, nlen)        
            push!(res,[snps[i], nerr[i]*pp,nlen[i]*pp])
        end
    end
    # ArgumentError: reducing over an empty collection is not allowed    
    if isempty(res)
        res
    else
        res2=permutedims(reduce(hcat, res)) 
        res3=sortsplit(res2,1)
        # return res[i]: snp_index, nerr, nobs
        res4 = []
        for i in res3
            sumerr = sum(i[:,2])
            sumobs = sum(i[:,3])
            errrates =sumerr ./ sumobs            
            errrates[isnan.(errrates)] .= -1.0
            if maximum(errrates) >= correctthreshold
                push!(res4, [i[1,1], errrates, sumerr,sumobs])
            end
        end        
        res4
    end    
end


function get_phaplo_missp(chrphaplo::AbstractMatrix, chrmultiallele::AbstractVector, 
    priorspace::AbstractDict, polyped::PolyPed,p::Integer,snp::Integer)
    # hhset, a list of all possible phased genotypes for parent p witn ploidy n
    # fhaplodict, a dict for each subpopulation, all haplotype of two-parents at chr/snp
    ploidy=polyped.parentinfo[p,:ploidy]
    nallele = chrmultiallele[snp]
    hhset = get_phased_codeset(ploidy, nallele)
    fhaplodict=Dict([begin    
        parents = priorspace[popid]["parentindex"]
        if p  in parents
            hhset2=repeat(Vector{Int8}.(permutedims(chrphaplo[snp,:])),length(hhset))
            hhset2[:,p] .= hhset            
            hhset3=hhset2[:,parents]
            hhset4=reduce(hcat,[vcat(i...) for i=eachrow(hhset3)])'
            popid=>hhset4
        else
            popid=>[]
        end
    end for popid = keys(priorspace)])
    hhset, fhaplodict
end

function get_estdose_missp(resdecode::AbstractVector, 
    chrphaplo::AbstractMatrix, 
    chrmultiallele::AbstractVector, 
    priorspace::AbstractDict, 
    priorprocess::AbstractDict, 
    polyped::PolyPed,
    p::Integer,snp::Integer; 
    callthreshold::Real,
    snporder::AbstractVector)
    #  return hhset,offls,estdose
    # hhset, a list of all H possible phased genotypes for parent p witn ploidy n
    # offls is a list of offspring indices whose parent is p
    # estdose matrix with size: H x length(offls)
    hhset, fhaplodict = get_phaplo_missp(chrphaplo, chrmultiallele, priorspace, polyped, p, snp)
    noff = size(polyped.offspringinfo,1)
    res=Vector{Any}(missing,noff)    
    markerincl = first(values(priorprocess)).markerincl # markerincl is the same for all popid
    markerindices = snporder[markerincl]
    t = findfirst(markerindices .== snp)
    for popid in keys(priorspace)
        fhaplo=fhaplodict[popid]
        isempty(fhaplo) && continue
        groupstates = priorspace[popid]["groupstate"]
        offls = priorspace[popid]["offspringindex"]
        ploidy = priorspace[popid]["offspringploidy"]        
        nunphased = get_nunphased(ploidy, chrmultiallele[snp])
        nhaplo = size(fhaplo,1)
        for off=offls            
            genoprob = last(resdecode[off])
            I,V = findnz(genoprob[t,:])
            prob = [zeros(nunphased) for _ in 1:nhaplo]
            for i in eachindex(I, V)
                d = unphased_genocode_2_index.(eachrow(fhaplo[:,groupstates[I[i]]]))
                for j in 1:nhaplo
                    prob[j][d[j]]  += V[i]
                end
            end
            res[off] = unphasedcalling.(prob; callthreshold)    
        end
    end    
    offls = findall(.!ismissing.(res))
    estdose = res[offls]
    estdose2 = collect(eachrow(reduce(hcat,estdose)))
    hhset,offls,estdose2
end

function get_badchanges!(resdecode::AbstractVector, badmarkers::AbstractVector,    
    chrphaplo::AbstractMatrix, 
    chrcalledgeno::AbstractMatrix, 
    chrmultiallele::AbstractVector, 
    priorspace::AbstractDict, 
    prorprocess::AbstractDict, 
    polyped::PolyPed; 
    minerrdiff::Integer=3, 
    callthreshold::Real=0.7,
    snporder::AbstractVector)    
    res = []
    for m= eachindex(badmarkers)
        snp, ww = badmarkers[m][1:2] # snp_index, mismatch_rate, nerr, nobs        
        any(ww .> 0) || continue        
        pp = findall(ww .≈ maximum(ww))  
        # pp = sortperm(ww,rev=true)
        # pp =  pp[ww[pp] .> 0]        
        while true
            subres=[]
            for p=pp
                hhset, offls, estdose = get_estdose_missp(resdecode, chrphaplo, chrmultiallele, 
                    priorspace, prorprocess, polyped,p,snp;  callthreshold,snporder)
                obs=chrcalledgeno[snp,offls]
                mismatch = [map((x,y)->ismissing(x) || ismissing(y) ? missing : x != y, obs, i) for i= estdose]
                nerr = [all(ismissing.(i)) ? 0 : sum(skipmissing(i)) for i in mismatch]                    
                newerr,newindex = findmin(nerr)
                oldg = chrphaplo[snp,p]
                oldindex = findfirst([i == oldg for i=hhset])
                olderr = nerr[oldindex]                
                if olderr - newerr >= minerrdiff
                    newg = hhset[newindex]                    
                    push!(subres,[snp,p,ww[p], newerr,join(newg,"|"),olderr,join(oldg,"|")])
                end
            end
            if isempty(subres)
                break
            else
                ndiff = [i[6]-i[4] for i=subres]
                ii = argmax(ndiff)
                push!(res,subres[ii])
                snp,p,newg=subres[ii][[1,2,5]]                
                setdiff!(pp, [p])
                isempty(pp) && break
            end
        end
    end
    pushfirst!(res,["markerindex", "parentindex","weight", "new_nerr","new_genotype", "old_nerr", "old_genotype"])
    res2=permutedims(reduce(hcat,res))    
    df=DataFrame(res2[2:end,:],res2[1,:])
    df
end


function getchangesum(changedf::DataFrame)
    # changedf: ["markerindex", "parentindex","weight", "new_nerr","new_genotype", "old_nerr", "old_genotype"]
    res=[]
    for m=1:size(changedf,1)
        oldphase, newphase, oldnerr,  newnerr = changedf[m,["old_genotype","new_genotype","old_nerr","new_nerr"]]
        newphase2 = parse.(Int,split(newphase,"|"))
        oldphase2 = parse.(Int,split(oldphase,"|"))
        nerrdiff = oldnerr - newnerr
        isdosediff = (sum(newphase2)-sum(oldphase2)) != 0
        nflip = sum(abs.(newphase2 - oldphase2))
        push!(res,[isdosediff,nflip,nerrdiff])
    end
    res2=reduce(hcat,res)'
    ndosediff,nalleleflip,nerrdiff =sum(res2,dims=1)
    nphasediff = size(changedf,1)
    (ndosediff=ndosediff,nphasediff=nphasediff,nalleleflip=nalleleflip,nerrdiff=nerrdiff)
end

function calchroffgenoprob(chrgenoprob::AbstractVector,chrphaplo::AbstractMatrix,
    offinfo::DataFrame,priorspace::AbstractDict)    
    chrfhaplo= get_phaplo(chrphaplo)
    offploidy=offinfo[!,:ploidy]
    offpop=offinfo[!,:family]
    nind = length(offpop)
    nsnp = size(chrgenoprob[1],1)
    doseprob= Matrix{Vector{Float64}}(undef,nind,nsnp)
    for popid = keys(priorspace)
        fhaplo = chrfhaplo[popid] .- 1
        states = priorspace[popid]["groupstate"]
        offls=findall(offpop .== popid)
        doseprob[offls,:] = [begin
            I,V = findnz(chrgenoprob[off][m,:])
            d=[sum(fhaplo[m,i]) for i=states[I]]
            if any(ismissing.(d))
                prob = ones(offploidy[off]+1)/(offploidy[off]+1)
            else
                prob = zeros(offploidy[off]+1)
                for i=eachindex(d)
                    prob[d[i]+1] += V[i]
                end
                prob
            end
        end for off=offls,m=1:nsnp]
    end
    permutedims(doseprob)
end
