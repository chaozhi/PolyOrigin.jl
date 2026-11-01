
function rand_update_bvpair!(siblogl::NamedTuple,  bvpair::NamedTuple,    
    phaploindex::AbstractVector,
    phaploset::AbstractVector,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    polyped::PolyPed;    
    phaseerror_midparent::Union{Nothing,Real}, 
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},     
    overdispersion::Union{Real,AbstractVector},
    byneighbor::Bool=false,     
    isrand::Bool=true,
    popidls=keys(priorspace),   
    parentindexls::AbstractVector=[], 
    snporder::Union{Nothing,AbstractVector}=nothing,
    incldelmarker::Bool=true)     
    if byneighbor
        update_bvpair!(siblogl, bvpair,phaploindex,phaploset, priorspace, priorprocess,
            chroffgeno,chroffformat, chrmultiallele, polyped; 
            phaseerror_midparent,doseerror, baseerror, allelicbias, overdispersion, popidls, 
            parentindexls, snporder, incldelmarker
        )
    else            
        rand_bvpair!(siblogl, bvpair,phaploindex,phaploset, priorspace, priorprocess,
            chroffgeno,chroffformat, chrmultiallele, polyped; 
            phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, popidls, 
            isrand, parentindexls, snporder, incldelmarker
        )                
    end         
end

function update_bvpair!(siblogl::NamedTuple,  bvpair::NamedTuple,    
    phaploindex::AbstractVector,
    phaploset::AbstractVector,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    polyped::PolyPed;    
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},         
    overdispersion::Union{Real,AbstractVector},
    popidls=keys(priorspace),   
    parentindexls::AbstractVector=[], 
    snporder::Union{Nothing,AbstractVector}=nothing,
    incldelmarker::Bool=true)
    isnothing(snporder) && (snporder = 1:size(chroffgeno,1))
    phaplo = get_phaplo(phaploindex, phaploset)    
    derivedgeno = get_derivedgeno(phaplo,chrmultiallele, priorspace; popidls, parentindexls)    
    nmarker = size(chroffgeno,1)
    parentformatls = repeat(["GT_phased"], nmarker)
    for popid in popidls
        poppriorspace = priorspace[popid]                
        offls = poppriorspace["offspringindex"]
        nstate = length(poppriorspace["state"])
        dataprobls = [ones(nstate) for _ in 1:nmarker]        
        for off=offls 
            offgeno = chroffgeno[:,off]
            ploidy = polyped.offspringinfo[off,:ploidy]                         
            cal_ind_dataprobls!(dataprobls, offgeno,chroffformat, chrmultiallele, derivedgeno[popid]; 
                ploidy,doseerror, baseerror, allelicbias,overdispersion)        
            bvpair.offspring[off],siblogl.offspring[off] = randbvprop(bvpair.offspring[off],siblogl.offspring[off],
                dataprobls,priorspace,priorprocess,popid,snporder)        
        end
        pls  = poppriorspace["midparentindex"]
        for p = pls        
            pgeno = phased_genocode_2_index.(phaplo[:, p], chrmultiallele)            
            ploidy = polyped.parentinfo[p,:ploidy]
            cal_ind_dataprobls!(dataprobls, pgeno,parentformatls, chrmultiallele, derivedgeno[popid]; 
                ploidy,doseerror=isnothing(phaseerror_midparent) ? doseerror : phaseerror_midparent, 
                baseerror, allelicbias,overdispersion)                
            bvpair.parent[p],siblogl.parent[p] = randbvprop(bvpair.parent[p],siblogl.parent[p],
                dataprobls,priorspace,priorprocess,popid,snporder)        
        end
    end
    if incldelmarker
        add_logl_delmarker!(siblogl, bvpair, phaplo, priorspace, priorprocess, 
            chroffgeno,chroffformat, chrmultiallele, polyped; 
            phaseerror_midparent,doseerror, baseerror, allelicbias,overdispersion, snporder)
    end
    bvpair,siblogl
end

function rand_bvpair!(siblogl::NamedTuple,  bvpair::NamedTuple,    
    phaploindex::AbstractVector,
    phaploset::AbstractVector,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    polyped::PolyPed;        
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},     
    overdispersion::Union{Real,AbstractVector},
    isrand::Bool=true,    
    popidls=keys(priorspace),   
    parentindexls::AbstractVector=[], 
    snporder::Union{Nothing,AbstractVector}=nothing,
    incldelmarker::Bool=true)
    isnothing(snporder) && (snporder = 1:size(chroffgeno,1))
    phaplo = get_phaplo(phaploindex, phaploset)    
    derivedgeno = get_derivedgeno(phaplo,chrmultiallele, priorspace; popidls, parentindexls)
    bvpairprop = get_prop_bvpair(priorspace,polyped)
    loglls = calmarglogl_bvpairprop(phaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
        priorspace, priorprocess, bvpairprop, polyped;         
        phaseerror_midparent,doseerror, baseerror, allelicbias, overdispersion, 
        popidls, parentindexls, snporder)
    # offspring
    offls = findall(.!(ismissing.(loglls.offspring)))         
    for o in offls    
        probls = normalizelogl(loglls.offspring[o])
        index = isrand ? rand(Categorical(probls)) : randargmax(probls) 
        bvpair.offspring[o] = bvpairprop.offspring[o][index]
        siblogl.offspring[o] = maximum(loglls.offspring[o]) # loglls.offspring[o][index]
    end            
    # parent    
    pls = findall(.!(ismissing.(loglls.parent)))            
    for p in pls    
        probls = normalizelogl(loglls.parent[p])
        index = isrand ? rand(Categorical(probls)) : randargmax(probls) 
        bvpair.parent[p] = bvpairprop.parent[p][index]
        siblogl.parent[p] = maximum(loglls.parent[p]) #  loglls.parent[p][index]
        # println("parent=",p,",isrand=", isrand, ", bvpair.parent=", bvpair.parent[p], " out of ", bvpairprop.parent[p], 
        #     ", logl.parent=", siblogl.parent[p],",loglls.parent=",loglls.parent[p])
    end   
    if incldelmarker
        add_logl_delmarker!(siblogl, bvpair, phaplo, priorspace, priorprocess, 
            chroffgeno,chroffformat, chrmultiallele, polyped; 
            phaseerror_midparent,doseerror, baseerror, allelicbias,overdispersion, snporder)
    end
    siblogl, bvpair
end

function add_logl_delmarker!(siblogl::NamedTuple,  bvpair::NamedTuple,        
    phaplo::AbstractMatrix, 
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    polyped::PolyPed;    
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},         
    overdispersion::Union{Real,AbstractVector},
    snporder::Union{Nothing,AbstractVector}=nothing)
    isnothing(snporder) && (snporder = collect(1:size(chroffformat,1)))
    markerincl = first(values(priorprocess)).markerincl
    markerls = snporder[.!markerincl]
    if !isempty(markerls)
        derivedgeno = get_derivedgeno(phaplo,chrmultiallele, priorspace)
        singlelogl = cal_singlelogl(phaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
            priorspace, priorprocess, bvpair; phaseerror_midparent, doseerror, 
            baseerror, allelicbias, overdispersion, markerls, 
        ) 
        bvdict = get_bvdict(bvpair, priorspace)    
        isologlls = sum(singlelogl, dims=2)[:,1]
        isologl = logl_vec2tuple!(isologlls,bvdict, polyped)
        siblogl.offspring .+= isologl.offspring        
        siblogl.parent .+= isologl.parent
    end
    nothing
end

function logl_vec2tuple!(loglvec::AbstractVector,bvdict::AbstractDict, polyped::PolyPed)
    noff = size(polyped.offspringinfo, 1)
    np = size(polyped.parentinfo, 1)
    offls = Vector{Union{Missing, Float64}}(missing, noff)
    pls = Vector{Union{Missing, Float64}}(missing, np)
    for (_, val) in bvdict
        if !isempty(val[1]) 
            noff = length(val[1])
            offls[val[1]] .= loglvec[1:noff]
            deleteat!(loglvec, 1:noff)
        end
        if !isempty(val[2]) 
            np = length(val[2])
            pls[val[2]] .= loglvec[1:np]
            deleteat!(loglvec, 1:np)
        end    
    end
    (offspring=offls, parent=pls)
end


function normalizelogl(loglls_ind)
    # logllist is a list of logl for a list of bvpair (bivalent and multi-vavelent) for a given individual     
    # equally probable for each of 9 or 16 the chromosome pairing
    # prior(V_0|H)=1/9 for bvModel, or 1/16 for fullModel for tetraploid    
    ls = loglls_ind .- log(length(loglls_ind)) 
    logsumexp!(ls)
    exp.(ls)    
end



function init_bvpair(priorspace::AbstractDict,polyped)
    bvpairprop = get_prop_bvpair(priorspace,polyped)
    bvpair_off = rand.(bvpairprop.offspring)
    bvpair_parent = [ismissing(i) ? missing : rand(i) for i in bvpairprop.parent]
    (offspring=bvpair_off, parent=bvpair_parent)    
end


function init_siblogl(bvpair::NamedTuple)
    NamedTuple{keys(bvpair)}([[ismissing(j) ? missing : -Inf for j in i] for i in bvpair])    
end

function randbvprop(priorspace::AbstractDict, popid::AbstractString,loglgrid::AbstractMatrix)
    valents = priorspace[popid]["valent"]
    isnonmiss = .!ismissing.(valents)    
    linear = LinearIndices(valents)
    bool = ismissing.(loglgrid) .* isnonmiss
    res = linear[bool]
    length(res) == 0 ? missing : rand(res)
end

function randbvprop(bvinput::Integer,loglinput::Real,dataprob::AbstractVector,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    popid::AbstractString,
    snporder::Union{Nothing,AbstractVector}=nothing)
    loglgrid = Matrix{Union{AbstractFloat,Missing}}(missing,size(priorspace[popid]["valent"]))
    minrun = 3
    maxrun = 10
    loglhis = [-Inf]
    bvhis = [0]
    isnothing(snporder) && (snporder=1:length(dataprob))
    for run=1:maxrun
        bvnow = randbvprop(priorspace,popid,loglgrid)
        ismissing(bvnow) && break
        bv = uploglgrid!(loglgrid,bvnow,dataprob,priorspace,priorprocess,popid,snporder)
        logl = loglgrid[bv]
        push!(bvhis,bv)
        push!(loglhis,logl)
        nrepeat =sum(loglhis .== max(loglhis...))
        # println("run=",run, ",logl0=", round(loglinput,digits=2),
        #     ",bv0=", bvinput, ", bv=", bv,
        #     ", nnonmissgrid=",length(loglgrid)-sum(ismissing.(loglgrid)),
        #     ",loglhis=", round.(loglhis,digits=2),",nrepeat=",nrepeat)
        nrepeat == minrun && break
    end
    bv,logl=last(bvhis), last(loglhis)
    if logl < loglinput
        bv = uploglgrid!(loglgrid,bvinput,dataprob,priorspace,priorprocess,popid,snporder)
        logl = loglgrid[bv]
    end
    bv, logl
end



function get_prop_bvpair(priorspace::AbstractDict,polyped::PolyPed)
    bvpair_off = [begin
        popid = polyped.offspringinfo[off,:family]
        valents = priorspace[popid]["valent"]
        isnonmiss = findall(.!ismissing.(valents))
        findall([all(issorted.(i)) for i in valents[isnonmiss]]) # ordring of valent does not matter for unphased offspring data
        # LinearIndices(valents)[isnonmiss] 
    end for off=1:size(polyped.offspringinfo,1)]
    bvpair_parent = [if polyped.parentinfo[p, :isfounder]
        missing
    else
        popid = polyped.parentinfo[p,:family]
        valents = priorspace[popid]["valent"]
        isnonmiss = .!ismissing.(valents)
        LinearIndices(valents)[isnonmiss]  # ordring of valent does matter for unphased offspring data
    end for p=1:size(polyped.parentinfo,1)]
    (offspring=bvpair_off, parent=bvpair_parent)
end

function get_prop_bvpair(priorspace::AbstractDict,popid::AbstractString,
    loglgrid::AbstractMatrix)
    valents = priorspace[popid]["valent"]
    isnonmiss = .!ismissing.(valents)
    cart = CartesianIndices(valents)
    linear = LinearIndices(valents) .* isnonmiss
    # ntop = round(Int,length(loglgrid)^0.25)
    nb=priorspace[popid]["valentneighbor"]
    ntop=round(Int,mean(length.(nb[1])))
    isloglnonmiss = .!ismissing.(loglgrid)
    p = sortperm(loglgrid[isloglnonmiss],rev=true)
    (length(p) >= ntop) && (p = p[1:ntop])
    top = cart[isloglnonmiss][p]
    xtop = unique([i[1] for i=top])
    ytop = unique([i[2] for i=top])
    res=reshape(linear[xtop,ytop],:)
    res[res .> 0]
end