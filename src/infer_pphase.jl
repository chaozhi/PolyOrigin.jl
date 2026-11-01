function update_phaplo_bvpair!(siblogl::NamedTuple,
    phaploindex::AbstractVector,    
    phaploset::AbstractVector,
    phaploweight::AbstractVector,        
    phaploposterior::AbstractVector, 
    priorspace::AbstractDict,
    priorprocess::AbstractDict,    
    bvpair::NamedTuple,    
    chrparentgeno::AbstractMatrix,
    chrparentformat::AbstractVector, 
    chrliketempfile::AbstractString,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    polyped::PolyPed;         
    phaseerror_midparent::Union{Nothing,Real},
    byparent::Bool,
    byneighbor::Bool,
    fixedparentls::AbstractVector, 
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},        
    overdispersion::Union{Real,AbstractVector},    
    israndls::AbstractVector,
    isreversechr::Bool,    
    alwaysaccept::Bool,     
    upbyblock::Bool, 
    upbyhalf::Bool)    
    cache_dataprob_persite(chroffgeno,chroffformat,chrmultiallele,priorspace; 
        doseerror, baseerror, allelicbias, overdispersion, chrliketempfile
    )    
    # inputlogl: discrete uniform prior distribution for alent formation
    inputlogl = sum_margloglls(siblogl)
    inputlogl += sum(cal_logl_parentdata(phaploindex,phaploweight))   
    oldlogl = inputlogl

    inputphaploindex = deepcopy(phaploindex)    
    newphaploindex = deepcopy(phaploindex)    
    newbvpair,newsiblogl = deepcopy(bvpair),deepcopy(siblogl)         
    
    pindexls = last(get_pindexls(polyped; fixedparentls, byparent))    
    nparent = size(polyped.parentinfo, 1)      
    if isempty(pindexls)
        @warn string("all parents are phased")
    else
        issetequal(1:nparent, vcat(reduce(vcat, pindexls),fixedparentls)) || @error string("unexpected pindexls=",pindexls)
    end
    reversechrls = isreversechr ? [true, false] : [false, true]                 
    for pindex = pindexls             
        isrand = all(israndls[pindex])         
        # For midparent parents update phaploposterior, The update does not depend on the current values of phaploposterior
        update_midparent_prior!(newphaploindex, phaploset,phaploweight, phaploposterior, pindex,  
            priorspace,priorprocess, chrparentgeno, chrparentformat, chrmultiallele, polyped;             
            doseerror, baseerror, allelicbias, overdispersion, 
            isreversechr, samplesize = 1000, 
            bvpair = newbvpair, #  bvpair = nothing results much larger phasing error of midparents
        )          
        # update pphase
        if upbyblock            
            # upbyhalf == upbyblock
            newpphase = rand_pphase_byhalf(newphaploindex,phaploset, phaploposterior, pindex, newbvpair, 
                priorspace, priorprocess,chrliketempfile, chroffformat, chrmultiallele; 
                phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, 
                reversechrls,isrand
            ) 
        else           
            newpphase = rand_pphase(newphaploindex,phaploset, phaploposterior, pindex, newbvpair, 
                priorspace, priorprocess,chrliketempfile, chroffformat, chrmultiallele; 
                phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, 
                isreversechr,isrand
            )  
        end
        for snp in eachindex(newpphase)
            pphase_snp = newpphase[snp]
            ismissing(pphase_snp) && continue # phases are missing at deleted marekrs or at non-updated markers
            for j in eachindex(pphase_snp, pindex)
                newphaploindex[pindex[j]][snp] = pphase_snp[j]
            end
        end
        if upbyblock
            rand_update_bvpair!(newsiblogl, newbvpair,newphaploindex,phaploset, priorspace, priorprocess,
                chroffgeno,chroffformat, chrmultiallele, polyped; 
                popidls = popfrompindex(pindex,priorspace), 
                parentindexls = pindex,
                phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion,  
                byneighbor, isrand = true, incldelmarker =false,
            )  
            newlogl = sum_margloglls(newsiblogl)
            newlogl += sum(cal_logl_parentdata(newphaploindex,phaploweight))            
            deltlogl = newlogl - oldlogl    
            # println("pindex=",pindex, ", newlogl=", newlogl, ", dellogl=",deltlogl)            
            if deltlogl > 0        
                phaploindex[pindex] .= deepcopy.(newphaploindex[pindex])                
                bvpair.offspring .= newbvpair.offspring 
                bvpair.parent .= newbvpair.parent 
                siblogl.offspring .= newsiblogl.offspring
                siblogl.parent .= newsiblogl.parent
                oldlogl = newlogl
            else
                newphaploindex[pindex] .= deepcopy.(phaploindex[pindex])                
                newbvpair.offspring .= bvpair.offspring 
                newbvpair.parent .= bvpair.parent 
                newsiblogl.offspring .= siblogl.offspring
                newsiblogl.parent .= siblogl.parent                
            end            
        end
    end          
    # update bvpair    
    if !upbyblock
        if isinf(inputlogl)
            popidls = keys(priorspace)
        else
            popidls = collect(skipmissing([issubset(priorspace[i]["parentindex"],fixedparentls) ? missing : i for i in keys(priorspace)]))
        end
        rand_update_bvpair!(newsiblogl, newbvpair,newphaploindex,phaploset, priorspace, priorprocess,
            chroffgeno,chroffformat, chrmultiallele, polyped; popidls, 
            phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion,  
            byneighbor, isrand = true, incldelmarker =false,
        )      
    end
    markerincl = first(values(priorprocess)).markerincl
    ndiff,ntot = cal_ndiff_phaploindex(inputphaploindex, newphaploindex, phaploset, markerincl,fixedparentls)       
    newlogl = sum_margloglls(newsiblogl)
    newlogl += sum(cal_logl_parentdata(newphaploindex,phaploweight))
    deltlogl = newlogl - inputlogl    
    # update phaplo
    isaccept = alwaysaccept || deltlogl > 0                      
    if isaccept        
        phaploindex .= deepcopy.(newphaploindex)                
        bvpair.offspring .= newbvpair.offspring 
        bvpair.parent .= newbvpair.parent 
        siblogl.offspring .= newsiblogl.offspring
        siblogl.parent .= newsiblogl.parent
    end
    ndiff,ntot, deltlogl,isaccept
end

function cache_dataprob_persite(chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,    
    priorspace::AbstractDict;
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},     
    overdispersion::Union{Real,AbstractVector},
    chrliketempfile::AbstractString)
    nmarker, noff = size(chroffgeno)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nmarker) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nmarker) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nmarker) : allelicbias    
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nmarker) : overdispersion
    errorhash = hash((doseerrorls,baseerrorls,allelicbiasls,overdispersionls))
    if isfile(chrliketempfile)
        olderrorhash = JLD2.load(chrliketempfile,"errorhash")
        errorhash == olderrorhash  && return nothing
    end
    jldopen(chrliketempfile,"w") do file        
        write(file,"errorhash",errorhash)        
        for m = 1:nmarker
            siteoffgeno = chroffgeno[m,:]
            site_offformat = chroffformat[m]
            site_nallele = chrmultiallele[m]
            doseerror_m = doseerrorls[m]
            baseerror_m = baseerrorls[m]
            allelicbias_m = allelicbiasls[m]
            overdispersion_m = overdispersionls[m]
            dataprobls = Vector{Vector{Float64}}(undef, noff)
            for popid in keys(priorspace)
                poppriorspace =  priorspace[popid]
                ploidy = poppriorspace["offspringploidy"]                  
                offls = poppriorspace["offspringindex"] # leaf_offspring        
                for o in offls
                    dataprobls[o] = get_dataprob(siteoffgeno[o]; 
                        doseerror = doseerror_m, baseerror = baseerror_m, 
                        allelicbias = allelicbias_m, overdispersion = overdispersion_m, 
                        format = site_offformat, nallele = site_nallele, ploidy,
                    ) 
                end
            end
            file[string(m)] = dataprobls
        end
    end
    nothing
end

function findgoodparents(polyped::PolyPed; progenyperhaplo::Integer=5, minshare::Real=0.45)
    popsizels = [sum(polyped.offspringinfo[!,:family] .== fam) for fam in polyped.designinfo[!,:family]]
    names(polyped.designinfo)[2:end] == polyped.parentinfo[!,:individual] || error("parentinfo and designinfo are not consistent")
    parentpopls = [findall(polyped.designinfo[!,c] .> 0) for c in 2:size(polyped.designinfo, 2)]
    sharels = map(x -> 1 - maximum(popsizels[x]) / sum(popsizels[x]), parentpopls)
    nprogenyls = polyped.parentinfo[!,:immediate_progeny]
    avgploidy = mean(polyped.parentinfo[!,:ploidy])
    goodparents = findall([sharels[i] >= minshare && nprogenyls[i] >= progenyperhaplo*avgploidy for i in eachindex(nprogenyls,sharels)])
    goodparents
end


function get_pindexls(polyped; byparent::Bool=false,fixedparentls=[], isrand::Bool=false,
    progenyperhaplo::Integer=5, minshare::Real=0.45)
    if byparent
        nparent = size(polyped.parentinfo, 1)        
        pindexls = [[i] for i in 1:nparent]
        goodparents = nothing
    else
        parent2off = get_parent2progeny(polyped)
        parentls = polyped.parentinfo[!,:individual]
        indls = collect(eachindex(parentls))   
        ppopls = unique([indls[Vector(i) .> 0] for i in  eachrow(polyped.designinfo[!,1 .+ indls])]) # 1st = familiy
        goodparents = PolyOrigin.findgoodparents(polyped; progenyperhaplo, minshare)                
        pindexls = Vector{Int}[]
        for ppop in ppopls
            setdiff!(ppop, goodparents)         
        end    
        unique!(ppopls)
        deleteat!(ppopls, isempty.(ppopls))        
        remove_parentchild_block!(ppopls,polyped)
        while !isempty(ppopls)
            if isrand
                ii = rand(eachindex(ppopls))
            else
                # wwls = [(minimum([length(parent2off[parentls[p]]) for p in pop])) for pop in ppopls]
                # ii = argmin(wwls)
                wwls = [length(reduce(union,[parent2off[parentls[p]] for p in pop])) for pop in ppopls]                
                probls = 1 ./ wwls 
                probls ./= sum(probls)
                ii = rand(Categorical(probls))
            end
            pindex = ppopls[ii]
            savepindex = copy(pindex)
            # println("pindex=",savepindex, ",ii=", ii, ", w=",wwls[ii], ",wls=",wwls)
            # @info "" [ppopls wwls]          
            pushfirst!(pindexls, savepindex)
            for ppop in ppopls
                setdiff!(ppop, savepindex)         
            end    
            unique!(ppopls)
            deleteat!(ppopls, isempty.(ppopls))                    
        end        
        isempty(goodparents) || prepend!(pindexls, [[i] for i in goodparents])    
        wwls = [length(reduce(union,[parent2off[parentls[p]] for p in pop])) - maximum(pop) for pop in pindexls]         
        oo = sortperm(wwls;rev=true)
        # oo = sortperm(maximum.(pindexls))
        pindexls = pindexls[oo]
    
    end        
    if !isempty(fixedparentls)
        for i in pindexls
            setdiff!(i,fixedparentls)
        end
        deleteat!(pindexls, isempty.(pindexls))
    end
    goodparents, pindexls
end

function remove_parentchild_block!(ppopls,polyped)
    ped = PolyOrigin.toindexped(polyped.ped)
    res = []
    for i in eachindex(ppopls)
        block = ppopls[i]
        length(block) != 2 && continue
        p1,p2 = block
        ped.member[p1] == p1 || error("inconsistent member indices")
        isparentchild = in(p2, [ped.mother[p1],ped.father[p1]]) || in(p1, [ped.mother[p2],ped.father[p2]])
        isparentchild && push!(res, i)
    end
    isempty(res) && return ppopls
    badls = ppopls[res]
    @info string("split parentchild blocks: ", badls) maxlog = 1
    setdiff!(ppopls, badls)
    d = setdiff(reduce(vcat,badls), reduce(vcat,ppopls))
    for i in d
        push!(ppopls,[i])
    end
    ppopls
end


# phaploweight is modified
function update_midparent_prior!(phaploindex::AbstractVector,
    phaploset::AbstractVector,
    phaploweight::AbstractVector,        
    phaploposterior::AbstractVector, 
    parentls::AbstractVector,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict,    
    chrparentgeno::AbstractMatrix,
    chrparentformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    polyped::PolyPed;         
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},    
    overdispersion::Union{Real,AbstractVector},    
    snporder::Union{Nothing, AbstractVector}= nothing,    
    bvpair::Union{Nothing, NamedTuple}=nothing,     
    isreversechr::Bool,
    samplesize::Integer)
    pmidparentls = parentls[.!polyped.parentinfo[parentls, :isfounder]]
    isempty(pmidparentls) && return nothing
    isnothing(snporder) && (snporder = collect(1:size(chrparentgeno,1)))         
    if isreversechr
        reverse!(snporder)        
        for (key, val) in priorprocess
             reverseprior!(val)
        end        
    end
    if isnothing(bvpair)        
        bvpairprop = get_prop_bvpair(priorspace,polyped)
    else        
        bvpairprop = (offspring = [[i] for i in bvpair.offspring], parent = [[i] for i in bvpair.parent])
    end
    popidls = unique(polyped.parentinfo[pmidparentls, :family])    
    phaplo = get_phaplo(phaploindex, phaploset)  
    derivedgeno = get_derivedgeno(phaplo,chrmultiallele, priorspace; popidls)  
    nmarker = length(chrmultiallele)
    chrmidparentformat = chrparentformat 
    for popid in popidls        
        popderivedgeno = derivedgeno[popid]            
        poppriorspace = priorspace[popid]        
        nstate = length(poppriorspace["state"])        
        dataprobls = [ones(nstate) for _ in 1:nmarker]                
        midparentls = intersect(pmidparentls,poppriorspace["midparentindex"])        
        for p = midparentls            
            pgeno = map((x,y)-> occursin("GT",x) && isa(y,AbstractString) ? unphased_genocode_2_index(y) : y, chrmidparentformat, chrparentgeno[:,p])                        
            ploidy = polyped.parentinfo[p,:ploidy]
            cal_ind_dataprobls!(dataprobls, pgeno, chrmidparentformat, chrmultiallele, popderivedgeno; 
                ploidy,doseerror, baseerror, allelicbias, overdispersion, 
            )                   
            bvlist = bvpairprop.parent[p]
            samplels =  priorsample_midparent(dataprobls, poppriorspace, priorprocess, bvlist; snporder, samplesize)                        
            update_phaploweight!(samplels, phaplo, p, phaploset, phaploweight, phaploposterior, poppriorspace, priorprocess,snporder)
        end    
    end    
    if isreversechr        
        reverse!(snporder)
        for (key, val) in priorprocess
             reverseprior!(val)
        end        
    end
    nothing
end


function update_phaploweight!(samplels::AbstractVector, phaplo::AbstractMatrix, 
    midparent::Integer, 
    phaploset::AbstractVector, 
    phaploweight::AbstractVector, 
    phaploposterior::AbstractVector, 
    poppriorspace::AbstractDict, priorprocess::AbstractDict,    
    snporder::AbstractVector)
    markerls = snporder[first(values(priorprocess)).markerincl]    
    statespace = poppriorspace["state"]
    parentindex = poppriorspace["parentindex"]
    parent_haplols = reduce.(vcat, eachrow(phaplo[markerls,parentindex]))    
    phasesetls = view(phaploset[midparent],markerls)
    likeweightls = view(phaploweight[midparent],markerls)
    posteriordistls = view(phaploposterior[midparent],markerls)
    # samplels[i] = [bv,logl, samplemtx ((samplesize+1) x #markers)]; first row of samplemtx results from viterbipath
    probls = [i[2] for i in samplels]  
    probls .= normalize(exp.(probls .- maximum(probls)),1)
    originls = last.(samplels)    
    for t in eachindex(parent_haplols, phasesetls, posteriordistls)               
        length(posteriordistls[t]) == 1 && continue   # homozygous genotypes    
        phaseset = phasesetls[t] 
        parent_haplo = parent_haplols[t]                
        likeweight = likeweightls[t]
        posteriordist = posteriordistls[t]
        # if length(phaseset) != length(posteriordist)
        #     @info "" phaseset posteriordist likeweight maxlog = 20
        # end
        posteriordist .= 0.0         
        for k in eachindex(probls, originls)        
            ls = [parent_haplo[statespace[i]] for i in originls[k][:,t]]
            # push!(samplegeno, unique(ls))
            freqls = [count(==(i), ls) for i in phaseset] ./ length(ls)
            if sum(freqls) > 0
                posteriordist .+= probls[k] .* freqls
            end             
        end             
        # update
        weightsum = sum(posteriordist) 
        if weightsum > 0                    
            posteriordist ./= weightsum
            nphase = length(posteriordist)
            # mix error: account for uncertainties in random sampling; isrelavent to phaseerror_midparent                        
            mix = 0.01 
            posteriordist .=  [p * (1-mix) + mix/nphase for p in posteriordist]
            posteriordist .*= likeweight            
            posteriordist ./= sum(posteriordist)            
        else
            posteriordist .= 1/length(posteriordist)    # kind of correction by changing prior      
        end
    end
    nothing
end


function priorsample_midparent(dataprobls::AbstractVector,
    poppriorspace::AbstractDict,
    priorprocess::AbstractDict,
    bvlist::AbstractVector; 
    snporder::AbstractVector, 
    samplesize::Integer)
    markerincl = first(values(priorprocess)).markerincl # markerincl is the same for all popid    
    dataprobls2 = view(dataprobls, snporder[markerincl])    
    statespace = poppriorspace["state"]
    condstates = poppriorspace["condstate"]
    keyls = poppriorspace["valentkey"]
    loglmax = -Inf
    res=[]
    bvkeys = keyls[bvlist]
    keyset = unique(bvkeys)
    mindiff = log(1e-4)
    for strkey =keyset
        subbvlist = bvlist[bvkeys .== strkey]
        _, startprob, tranprobseq = strkey2prior(priorprocess, strkey) # tranprobseq has accounted fro markerincl
        for bv=subbvlist
            # @info "" [bv, strkey]
            condstates_bv = condstates[bv]
            dataprobseq = [i[condstates_bv] for i = dataprobls2]
            fwprob,fwscale = forward(startprob,tranprobseq,dataprobseq)
            logl = calloglike(fwscale)            
            if logl - loglmax >= mindiff
                bwsample=backwardsample(tranprobseq,dataprobseq,fwprob; samplesize)
                bwsample2 = condstates_bv[bwsample]
                (;path) = viterbidecode(startprob,tranprobseq,dataprobseq);
                vpath2 = condstates_bv[path]            
                bwsample3 = alignsample_absphase(bwsample2, vpath2, statespace; isalign=true) 
                push!(res,[bv, logl, bwsample3])
            end
            loglmax = max(logl, loglmax)            
        end
    end    
    bool = findall([i[2]-loglmax < mindiff for i=res])
    deleteat!(res,bool)        
    res
end

# vpath3 = setabsphase_vpath(vpath2, statespace)    # it does not work when chrpairing_phase = 22           
# function setabsphase_vpath(vpath::AbstractVector, statespace::AbstractVector)
#     state2index = Dict(statespace .=> eachindex(statespace))        
#     statevpath = statespace[vpath]
#     oo = sortperm(sum(statevpath))
#     [state2index[i[oo]] for i in statevpath]
# end

function alignsample_absphase(samplemtx::AbstractMatrix, vpath::AbstractVector, 
    statespace::AbstractVector; isalign::Bool=true)
    state2index = Dict(statespace .=> eachindex(statespace))        
    statevpath = reduce(hcat, statespace[vpath])'
    ploidy = size(statevpath,2)        
    # alight samples with statevpath    
    samplevec = collect.(eachrow(samplemtx)) # eachrow is a sample path
    if isalign
        for s in eachindex(samplevec)
            samplestates = statespace[samplevec[s]]
            statesamplepath = reduce(hcat, samplestates)'
            perm = alignsample_vpath(statesamplepath, statevpath,ploidy)        
            if perm != 1:ploidy        
                samplevec[s] .= [state2index[i[perm]] for i in samplestates]             
            end
        end
        allunique(statespace) || @error "statespace not unique"    
    end
    pushfirst!(samplevec, vpath) 
    reduce(hcat,samplevec)'
end

function alignsample_vpath(statesamplepath::AbstractMatrix,statevpath::AbstractMatrix, ploidy::Integer)
    perm=[argmin([sum(i .!= col)  for col =eachcol(statesamplepath)]) for i=eachcol(statevpath)]
    if length(unique(perm)) < ploidy
        permls = collect(permutations(1:ploidy))
        dls = [sum(statesamplepath[:,i] .!= statevpath) for i in permls] 
        perm = permls[argmin(dls)]
    end
    perm
end

function rand_pphase_byhalf(inputphaploindex::AbstractVector,
    phaploset::AbstractVector,
    phaploposterior::AbstractVector,    
    pindex::AbstractVector,    
    bvpair::NamedTuple,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chrliketempfile::AbstractString,
    chroffformat::AbstractVector,  
    chrmultiallele::AbstractVector;    
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector},     
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector}, 
    reversechrls::AbstractVector, 
    isrand::Bool=true,
    snporder::Union{Nothing, AbstractVector}= nothing)
    nsnp = length(chroffformat)
    isnothing(snporder) && (snporder = collect(1:nsnp))            
    resphase = Vector{Any}(missing, nsnp)    
    phaploindex = deepcopy(inputphaploindex)
    ncut = div(nsnp,2)
    for reversechr =reversechrls        
        phaploindex2 = deepcopy(phaploindex)        
        # phaploset2 = deepcopy(phaploset)                
        # phaploposterior2 = deepcopy(phaploposterior)
        phaploset2 = [copy.(i) for i in phaploset] # faster than deepcopy
        phaploposterior2 = [copy.(i) for i in phaploposterior]                
        snpls = reversechr ? snporder[1:ncut] : snporder[ncut+1:nsnp]        
        for p = pindex
            for j in snpls
                ismissing(phaploindex[p][j]) && continue
                phaploset2[p][j]= phaploset[p][j][[phaploindex[p][j]]]
                phaploposterior2[p][j] = [1.0]
                phaploindex2[p][j] = 1
            end
        end    
        halfphase = rand_pphase(phaploindex2,phaploset2, phaploposterior2, pindex, bvpair, 
            priorspace, priorprocess,chrliketempfile, chroffformat, chrmultiallele; 
            phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, 
            isreversechr = !reversechr,isrand, snporder, 
        ) 
        changels = setdiff(snporder,snpls)
        resphase[changels] .= halfphase[changels]        
        for j in eachindex(pindex)
            for i in changels
                ismissing(halfphase[i]) && continue
                phaploindex[pindex[j]][i] = halfphase[i][j] 
            end 
        end
    end    
    resphase
end

function rand_pphase(phaploindex::AbstractVector,
    phaploset::AbstractVector,
    phaploposterior::AbstractVector,
    pindex::AbstractVector,    
    bvpair::NamedTuple,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chrliketempfile::AbstractString, 
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector},
    isreversechr::Bool, 
    isrand::Bool=true, 
    snporder::Union{Nothing, AbstractVector}= nothing)
    isnothing(snporder) && (snporder = collect(1:length(chroffformat)))        
    if isreversechr
        reverse!(snporder)        
        for (key, val) in priorprocess
             reverseprior!(val)
        end        
    end
    fwlogorigprob, fwphaseset, fwlogpphase,bvdict = polyforward(phaploindex,phaploset, phaploposterior, 
        pindex, bvpair, priorspace, priorprocess,chrliketempfile, chroffformat, chrmultiallele; 
        phaseerror_midparent,doseerror, baseerror, allelicbias,overdispersion, snporder, 
    )
    pphases = polybackward!(fwlogorigprob,fwlogpphase,bvdict,priorprocess; isrand,snporder)    # phases are missing at deleted marekrs 
    # @info string("fw=",tfw," bw=",tbw)
    res = map((x,y)->ismissing(x) || ismissing(y) ? missing : x[y],fwphaseset,pphases)    
    if isreversechr        
        reverse!(snporder)
        for (key, val) in priorprocess
             reverseprior!(val)
        end        
    end
    res
end


function polyforward(phaploindex::AbstractVector,
    phaploset::AbstractVector,
    phaploposterior::AbstractVector,    
    pindex::AbstractVector,    
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chrliketempfile::AbstractString, 
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},     
    overdispersion::Union{Real,AbstractVector},
    snporder::Union{Nothing, AbstractVector}= nothing)    
    nseq=length(chroffformat)
    isnothing(snporder) && (snporder = 1:nseq)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nseq) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nseq) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nseq) : allelicbias    
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nseq) : overdispersion
    # Float16 to save memoryuse for log-scale orginprob
    fwlogorigprob = Vector{Union{Missing, Vector{Vector{Vector{Float16}}}}}(missing,nseq) 
    fwlogpphase = Vector{Union{Missing,Vector}}(missing,nseq)    
    fwphaseset = Vector{Union{Missing,Vector}}(missing,nseq)
    bvdict = get_bvdict(pindex,bvpair,priorspace)    
    jldopen(chrliketempfile,"r") do file        
        tseq=findall(first(values(priorprocess)).markerincl)
        popidls = popfrompindex(pindex,priorspace)
        snp = snporder[tseq[1]]    
        sitedatalike = file[string(snp)]        
        sitephaplo, sitephaseset, sitephaseweight =get_sitephase(snp,pindex,phaploindex, phaploset, phaploposterior; minweight = 1e-8) 
        sitederived = get_site_derivedgeno(sitephaplo, priorspace; site_nallele = chrmultiallele[snp], popidls)
        dataprob = cache_sitedataprob_multiphase(sitephaplo, sitederived, sitedatalike, bvpair,priorspace; 
            phaseerror_midparent, doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp], popidls)
        startprob =  get_startprobls(priorprocess, bvdict)
        origprob = cal_forward_origprob!(startprob, dataprob)
        phaseprob = cal_forward_pphase!(origprob, sitephaseset,sitephaseweight) # fwlogorigprob is normalized         
        fwphaseset[snp] = sitephaseset
        fwlogpphase[snp] = log.(phaseprob)        
        fwlogorigprob[snp] =  [map(x->log.(x),h) for h in origprob] 
        # starttime = time()
        # tusedls = zeros(5)
        for kk in 2:length(tseq)
            tranprob = get_tranprobls(priorprocess, bvdict,tseq[kk-1])
            snp = snporder[tseq[kk]]     
            sitedatalike = file[string(snp)]     
            sitephaplo, sitephaseset, sitephaseweight =get_sitephase(snp,pindex,phaploindex, phaploset, phaploposterior; minweight = 1e-8)         
            sitederived = get_site_derivedgeno(sitephaplo, priorspace; site_nallele = chrmultiallele[snp], popidls)
            dataprob = cache_sitedataprob_multiphase(sitephaplo, sitederived, sitedatalike, bvpair,priorspace; 
                phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
                allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
                site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp], popidls)
            # snppre = snporder[tseq[kk-1]] 
            # phaseprob is give by the exp.(fwlogpphase[snppre])
            # origprob is given by exp  of fwlogorigprob[snppre]
            prob = sum(phaseprob .* origprob) 
            startprob = [any(isnothing.(tranprob[ind])) ? prob[ind] : kronvec(tranprob[ind][1]', tranprob[ind][2]',prob[ind]) for ind in eachindex(prob, tranprob)]                  
            origprob = cal_forward_origprob!(startprob, dataprob) # dataprob is modified into origprob 
            phaseprob = cal_forward_pphase!(origprob, sitephaseset, sitephaseweight)   
            fwphaseset[snp] = sitephaseset
            fwlogpphase[snp] = log.(phaseprob)            
            fwlogorigprob[snp] = [map(x->log.(x),h) for h in origprob]                    
        end        
        # tall = time()-starttime
        # println("tall=", round(tall,digits=3), ", tusedls/tall=", round.(tusedls/tall,digits=3))
    end    
    fwlogorigprob, fwphaseset, fwlogpphase,bvdict
end

function cal_forward_origprob!(startprob::AbstractVector, dataprob::AbstractVector)    
    for dataprob_h in dataprob
        for k in eachindex(startprob, dataprob_h)
            dataprob_h[k] .*= startprob[k] 
        end
    end
    dataprob
end

function cal_forward_pphase!(origprob::AbstractVector, sitephaseset::AbstractVector, siteweight::AbstractVector)    
    pphase = [sum(log.(normprobvec!.(h))) for h in origprob]  # modify origprob onsite
    pphase .+= log.(siteweight)        
    pphase .= exp.(pphase .- maximum(pphase))
    b = pphase .>= 1e-20
    keepat!(pphase,b)
    keepat!(origprob, b)
    keepat!(sitephaseset,b)
    pphase ./ sum(pphase)
end

function get_bvdict(pindex,bvpair::NamedTuple, priorspace::AbstractDict)
    popidls=popfrompindex(pindex,priorspace)
    get_bvdict(bvpair, priorspace; popidls)    
end

function get_bvdict(bvpair::NamedTuple, priorspace::AbstractDict; popidls=keys(priorspace))    
    # using ordereddict so that keys(bvdict) = popidls
    # the ordering will be the same as startprob and transtionprob
    OrderedDict([begin
        offls = priorspace[popid]["offspringindex"]
        pls = priorspace[popid]["midparentindex"]
        bvls = vcat(bvpair.offspring[offls],bvpair.parent[pls])               
        bvkeys = priorspace[popid]["valentkey"][bvls]
        popid => (offls,pls, bvkeys)
    end for popid = popidls])
end


function get_bvdict_exclparent(bvpair::NamedTuple, priorspace::AbstractDict; popidls=keys(priorspace))    
    # using ordereddict so that keys(bvdict) = popidls
    # the ordering will be the same as startprob and transtionprob
    OrderedDict([begin
        offls = priorspace[popid]["offspringindex"]
        bvls = bvpair.offspring[offls]
        bvkeys = priorspace[popid]["valentkey"][bvls]
        popid => (offls,[], bvkeys)
    end for popid = popidls])
end

# fwlogorigprob is modified
function polybackward!(fwlogorigprob::AbstractVector, fwlogpphase::AbstractVector,
    bvdict::AbstractDict,priorprocess::AbstractDict;
    snporder::Union{Nothing, AbstractVector}=nothing,
    isrand::Bool=true)
    nseq = length(fwlogpphase)
    isnothing(snporder) && (snporder = 1:nseq)        
    pphases = Vector{Union{Missing,Int}}(missing,nseq)
    tseq=findall(first(values(priorprocess)).markerincl)
    snp = snporder[tseq[end]]        
    probls = exp.(fwlogpphase[snp])    
    if isrand         
        pphases[snp] = rand(Categorical(probls)) 
        # change from Float16 to Float64 to increase accuracy
        logprobls = Vector{Vector{Float64}}(fwlogorigprob[snp][pphases[snp]]) 
        orig = [begin 
            pls = exp.(i)
            pls ./= sum(pls)
            rand(Categorical(pls)) 
        end for i in logprobls]
    else
        pphases[snp] = randargmax(probls) 
        logprobls = fwlogorigprob[snp][pphases[snp]]
        orig = [randargmax(i)   for i in logprobls]
    end
    fwlogorigprob[snp] = missing  # save memoryuse 
    for kk=length(tseq)-1:-1:1            
        tranprob= get_tranprobls(priorprocess, bvdict, tseq[kk]; calkron=true)
        snp = snporder[tseq[kk]]                        
        logorigprob = fwlogorigprob[snp]
        nstatels = length.(first(logorigprob))
        bwls = [if isnothing(tranprob[ind])
            pls = zeros(nstatels[ind])
            pls[orig[ind]] = 1.0
            pls
        else        
            tranprob[ind][:,orig[ind]] 
        end for ind in eachindex(tranprob, orig, nstatels)]
        nphase = length(logorigprob)
        ww = zeros(nphase)
        logpost = Vector{Vector{Vector{Float64}}}(undef, nphase)
        for f in 1:nphase
            fwls = logorigprob[f] 
            logpost[f] = [fwls[ind] .+ log.(bwls[ind]) for ind in eachindex(fwls, bwls)]         
            ls = logsumexp!.(logpost[f])
            ls[isinf.(ls)] .= -1000.0
            ww[f] = sum(ls)            
        end 
        ww .+= fwlogpphase[snp]
        logsumexp!(ww)            
        probls = exp.(ww)        
        if isrand 
            pphases[snp] = rand(Categorical(probls)) 
            orig .= [rand(Categorical(exp.(i))) for i in logpost[pphases[snp]]]
        else
            pphases[snp] = randargmax(probls) 
            orig .= [randargmax(i)  for i in logpost[pphases[snp]]]
        end
        fwlogorigprob[snp] = missing  # save memoryuse 
    end    
    pphases
end



function randargmax(probls::AbstractVector; atol=1e-5)
    m = maximum(probls)
    rand(findall(isapprox.(probls, m; atol)))
end