mutable struct JoinedPolyAncestry
    polyped::PolyPed    
    markermap::Vector{DataFrame}    
    parentgeno::Vector{Matrix}        
    offspringgeno::Vector{Matrix}
    joinedspace::DataFrame
    diploprob::Union{Nothing, Vector{Matrix}}
    genoprob::Union{Nothing, Vector{Matrix}}
    haploprob::Union{Nothing,Vector{Matrix}}     
    parentdiploprob::Union{Nothing, Vector{Matrix}}
    parentgenoprob::Union{Nothing, Vector{Matrix}}
    parenthaploprob::Union{Nothing,Vector{Matrix}}     
    viterbipath::Union{Nothing, Vector{Matrix}}    
    parentviterbipath::Union{Nothing, Vector{Matrix}}    
    function JoinedPolyAncestry(polyped, markermap,parentgeno,offspringgeno,joinedspace,        
        diploprob,genoprob, haploprob,parentdiploprob,parentgenoprob, parenthaploprob,
        viterbipath,parentviterbipath)        
        dims = [length(markermap),size(parentgeno,1)]
        if length(union(dims))!=1
            error("inconsistent number of linkage groups ",dims)
        end
        dims = [size.(markermap,1),size.(parentgeno,1)]
        if length(union(dims))!=1
            error("inconsistent number of markers within linkage groups ",dims)
        end
        dims = [[size(polyped.offspringinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of offspring ",dims)
        end
        dims = [union(size.(parentgeno,2)), [size(polyped.parentinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of parents ",dims)
        end
        new(polyped, markermap,parentgeno,offspringgeno,joinedspace,
            diploprob,genoprob, haploprob,parentdiploprob,parentgenoprob, parenthaploprob,
            viterbipath,parentviterbipath)
    end
end

#############################################################################################

function JoinedPolyAncestry(polyancestry::PolyAncestry; samplesize=1000,isparallel=true)
    empty!(lru_state_parent) # global constant        
    joined_parentviterbi = get_joined_parentviterbi(polyancestry)
    joined_viterbi = get_joined_viterbi(polyancestry,joined_parentviterbi)       
    parentgeno = polyancestry.parentgeno
    offspringgeno = polyancestry.offspringgeno
    joinedspace = get_joinedspace(polyancestry.polyped)
    parentdiploprob = get_joined_parentdiploprob(polyancestry; samplesize) # if no midparents, samplesize is not used!

    diploprob = get_joined_diploprob(polyancestry, parentdiploprob; samplesize,isparallel)
    parentgenoprob = tojoined_genoprob(parentdiploprob)
    genoprob = tojoined_genoprob(diploprob)
    parenthaploprob = tojoined_haploprob(parentgenoprob)
    haploprob = tojoined_haploprob(genoprob)
    JoinedPolyAncestry(polyancestry.polyped, polyancestry.markermap,parentgeno,offspringgeno,joinedspace,                
        diploprob,genoprob, haploprob,parentdiploprob,parentgenoprob, parenthaploprob,
        joined_viterbi,joined_parentviterbi
    )
end

function get_joinedspace(polyped::PolyPed)
    pinfo = polyped.parentinfo
    pls = pinfo[pinfo[!,:isfounder],1]
    fglls = get_fgl_founder(polyped)
    DataFrame("founder" =>pls, "homologs"=> join.(fglls, "|"))
end

function get_foundergeno(polyancestry::PolyAncestry)
    isfls = polyancestry.polyped.parentinfo[!,:isfounder]
    [copy.(i[:,isfls]) for i in polyancestry.parentgeno]
end


#############################################################################################
# transform parent viterbi


function get_joined_parentviterbi(polyancestry::PolyAncestry)
    polyped = polyancestry.polyped
    fgl_founder = get_fgl_founder(polyped)
    statespace_parent = get_statespace_parent(polyancestry)
    peddf_parent = get_peddf_parent(polyped)
    vpathls = polyancestry.parentviterbipath
    if isnothing(vpathls) 
        vpathls = [nothing for _ in eachindex(polyancestry.markermap)]
    end
    joined_parentviterbi = [begin 
        chrparentviterbi = vpathls[chr]        
        nmarker = size(polyancestry.markermap[chr],1)
        get_joined_chrparentviterbi(chrparentviterbi, fgl_founder, statespace_parent, peddf_parent,nmarker)
    end for chr in eachindex(vpathls)]
    joined_parentviterbi
end


const lru_state_parent = LRU{UInt, Vector}(maxsize=10^6)

function get_joined_chrparentviterbi(chrparentviterbi,fgl_founder, 
    statespace_parent,peddf_parent,nmarker)        
    if isnothing(chrparentviterbi)
        nparent = length(statespace_parent)    
        midindices = []
    else
        nmarker2, nparent = size(chrparentviterbi)    
        nmarker == nmarker2 || error("nmarker mismatch")
        midindices = findall(.!ismissing.(chrparentviterbi[1,:]))
    end
    res = Vector{Any}(undef, nparent)
    for i in eachindex(fgl_founder)
        res[i] = repeat([fgl_founder[i]], inner=nmarker)
    end    
    for p in midindices        
        res[p] = statespace_parent[p][chrparentviterbi[:,p]]
    end
    if length(res) == 1
        res2 = reshape(res[1],:, 1)
    else        
        res2 = reduce(hcat, res)
    end
    for t in 1:nmarker
        state_parent = view(res2, t, :)
        lru_key = hash(state_parent)
        new_state_parent = get!(lru_state_parent,lru_key) do 
            get_new_state_parent(state_parent, peddf_parent)                
        end
        state_parent .= new_state_parent
    end
    res2
end


function get_fgl_founder(polyped::PolyPed)
    founderindices = findall(polyped.parentinfo[!,:isfounder])
    founderindices == 1:length(founderindices) || @error string("unexpected founder indices: ", founderindices)
    ploidyls = polyped.parentinfo[founderindices,:ploidy]
    nallele = sum(ploidyls)
    inttype = nallele < typemax(Int16) ? Int16 : Int32
    aploidyls = accumulate(+,ploidyls)
    pushfirst!(aploidyls,0)
    fgl_founder = [Tuple(range(aploidyls[i]+1, aploidyls[i+1])) for i in 1:length(aploidyls)-1]
    [inttype.(i) for i in fgl_founder]    
end

function get_statespace_parent(polyancestry::PolyAncestry)
    pinfo = polyancestry.polyped.parentinfo
    statespace_parent = Vector{Any}(nothing, size(pinfo,1))
    midindices = findall(.!pinfo[!,:isfounder])
    for i in midindices    
        fam = pinfo[i,:family]
        statespace_parent[i] = Tuple.(polyancestry.statespace[fam]["state"])
    end 
    statespace_parent
end

function get_peddf_parent(polyped::PolyPed)
    pls =  polyped.parentinfo[!,:individual]
    subped = getsubped(polyped.ped, pls)
    subdf = ped2df(subped)
    dict = Dict(zip(pls, eachindex(pls)))
    push!(dict, "notapplicable"=>0)
    for j in 1:3
        subdf[!,j] .= [dict[i] for i in subdf[!,j]]
    end
    subdf[!,1] == 1:size(subdf,1) || @error string("unexpected peddf=",subdf)
    subdf
end

function get_new_state_parent(state_parent,peddf_parent)  
    new_state_parent = copy(state_parent)
    midindices = findall(peddf_parent[:,:father] .> 0)
    for i = midindices
        moth = peddf_parent[i,:mother]
        fath = peddf_parent[i,:father]
        state_ind,state_moth, state_fath = new_state_parent[[i,moth,fath]]
        ploidy_moth = peddf_parent[moth,:ploidy]
        halfploidy = div(ploidy_moth,2)
        orig_moth = state_moth[[state_ind[1:halfploidy]...]]
        if fath == moth            
            orig_fath = state_fath[[state_ind[halfploidy+1:end]...]]
        else
            ls = state_ind[halfploidy+1:end] .- ploidy_moth
            orig_fath = state_fath[[ls...]]
        end
        new_state_parent[i] = tuple(orig_moth...,orig_fath...) 
    end
    new_state_parent    
end

#############################################################################################
# transform offspring viterbi

function get_joined_viterbi(polyancestry::PolyAncestry,joined_parentviterbi::AbstractVector)
    statespace = polyancestry.statespace
    [begin 
        joined_chrparentviterbi = joined_parentviterbi[chr]
        chrviterbi = polyancestry.viterbipath[chr]
        get_joined_chrviterbi(chrviterbi,joined_chrparentviterbi, statespace)
    end for chr = eachindex(polyancestry.viterbipath, joined_parentviterbi)]
end

function get_joined_chrviterbi(chrviterbi,joined_chrparentviterbi, statespace)   
    isnothing(chrviterbi) && return nothing
    nind = size(chrviterbi,2)
    res = Vector{Any}(undef, nind)
    for popid = keys(statespace)
        statels = Tuple.(statespace[popid]["state"])
        pls = statespace[popid]["parentindex"]
        isself = allequal(pls)
        ploidy_moth = statespace[popid]["parentploidy"][1]
        joined_moth = joined_chrparentviterbi[:,pls[1]]
        joined_fath = joined_chrparentviterbi[:,pls[2]]
        offls = statespace[popid]["offspringindex"]
        for off in offls
            offpath = statels[chrviterbi[:,off]]
            res[off] = PolyOrigin.get_new_offpath(offpath,isself, ploidy_moth,joined_moth,joined_fath)
        end
    end
    if nind == 1
        reshape(res[1],:, 1)
    else
        reduce(hcat,res)
    end
end


function get_new_offpath(offpath,isself, ploidy_moth,joined_moth,joined_fath)    
    res = similar(offpath)
    halfploidy = div(ploidy_moth,2)
    for i in eachindex(res, offpath, joined_moth)            
        rg = [offpath[i][1:halfploidy]...]
        rg2 = [offpath[i][halfploidy+1:end]...]
        isself || (rg2 .-= ploidy_moth)        
        res[i] = (joined_moth[i][rg]...,joined_fath[i][rg2]...)
    end
    res
end

function tojoined_genoprob(joined_diploprob::AbstractVector)
    [diplo2geno_probpair.(i) for i in joined_diploprob]
end

function diplo2geno_probpair(diploprobpair::Pair)
    ii,vv = diploprobpair
    ss = [sort([i...]) for i in ii]
    ss2 = unique(ss)
    ww = [sum(vv[[s == i for i in ss]]) for s in ss2]
    o = sortperm(ww; rev=true)
    Tuple(Tuple.(ss2[o])) => Tuple(ww[o])
end

function tojoined_haploprob(joined_genoprob::AbstractVector)
    [geno2haplo_probpair.(i) for i in joined_genoprob]
end

function geno2haplo_probpair(genoprobpair::Pair)
    ii,vv = genoprobpair
    ii2 = [[i...] for i in ii]
    alleles = sort(unique(reduce(vcat, ii2)))
    prob = SparseVector(maximum(alleles), alleles, zeros(eltype(vv),length(alleles)))
    for i in eachindex(ii2)
        for j in ii2[i]
            prob[j] += vv[i]
        end
    end
    ii,vv = findnz(prob)
    o = sortperm(vv; rev=true)
    Tuple(ii[o]) => Tuple(vv[o])
end

#############################################################################################
# save/read joinedancestry

function savejoinedancestry(outfile::AbstractString,joinedancestry::JoinedPolyAncestry;    
    missingstring = "NA",
    workdir::AbstractString = pwd())
    ext = last(splitext(outfile))
    myopen = ext == ".gz" ? GZip.open : open    
    outputfile2 = getabsfile(workdir,outfile)
    myopen(outputfile2, "w") do io
        initial = "PolyOrigin"        
        write(io,initial,",pedigree","\n")
        savepolyped(io,joinedancestry.polyped; workdir,delim=',')
        # parentgeno
        append_dataframe(io, joinedancestry.joinedspace; dfname="joinedspace")    
        append_dataframe(io, phasedgeno2df(joinedancestry; isparent=true); dfname="parentgeno")    
        append_dataframe(io, phasedgeno2df(joinedancestry; isparent=false); dfname="offspringgeno")    
        # viterbipath                
        df = joined_viterbipath2df(joinedancestry; ismidparent=true)
        df2 = joined_viterbipath2df(joinedancestry; ismidparent=false)
        df =  isnothing(df) ? df2 : hcat(df, df2[!,6:end]) # col1-5: markermap
        append_dataframe(io, df; missingstring,dfname="viterbipath")        
        flush(io)        
        # condprob 
        for probtype in ["diploprob", "genoprob", "haploprob"]
            df = joined_condprob2df(joinedancestry; probtype)
            dfname = probtype == "diploprob" ? "phasedgenoprob" : (probtype == "genoprob" ? "unphasedgenoprob" : "haploprob")
            append_dataframe(io, df; missingstring,dfname)        
            flush(io)        
        end                
    end        
    outputfile2
end

function readjoinedancestry(ancestryfile::AbstractString;
    missingstring::AbstractString="NA",
    workdir::AbstractString=pwd())
    ancestryfile2 = getabsfile(workdir,ancestryfile)
    isfile(ancestryfile2) || @error(string(ancestryfile," does not exist in workdir = ",workdir))
    res = readmultitable(ancestryfile2; missingstring)
    haskey(res, "pedigree") || error(string("key genodata is not in ", keys(res)))  
    ped = Pedigree(res["pedigree"])        
    polyped = PolyPed(ped)
    markermap, parentgeno = parsephasedgeno(res["parentgeno"])    
    markermap2, offspringgeno = parsephasedgeno(res["offspringgeno"])    
    joinedspace = res["joinedspace"]    
    nallele = sum(length(split(i,"|")) for i in joinedspace[!,:homologs])
    inttype = nallele < typemax(Int16) ? Int16 : Int32
    if haskey(res, "viterbipath")         
        parentviterbipath, viterbipath = parseviterbipath(res["viterbipath"],polyped)
        joined_viterbi = isnothing(viterbipath) ? nothing : [parseviterbistate.(i,inttype) for i in viterbipath]
        joined_parentviterbi = isnothing(viterbipath) ? nothing : [parseviterbistate.(i,inttype) for i in parentviterbipath]
    else
        joined_viterbi, joined_parentviterbi = nothing, nothing
    end
    if haskey(res, "phasedgenoprob") || haskey(res, "diploprob") 
        probkey = haskey(res, "phasedgenoprob") ? "phasedgenoprob" : "diploprob"
        parentdiploprob, diploprob = parse_joinedcondprob(res[probkey], polyped; delimitor="|", inttype)
    else
        parentdiploprob, diploprob = nothing, nothing
    end
    if haskey(res, "unphasedgenoprob") || haskey(res, "genoprob") 
        probkey = haskey(res, "unphasedgenoprob") ? "unphasedgenoprob" : "genoprob"
        parentgenoprob, genoprob = parse_joinedcondprob(res[probkey], polyped; delimitor="/", inttype)
    else
        parentgenoprob, genoprob = nothing, nothing
    end
    if haskey(res, "haploprob")         
        parenthaploprob, haploprob = parse_joinedcondprob(res["haploprob"], polyped; delimitor=nothing, inttype)
    else
        parenthaploprob, haploprob = nothing, nothing
    end
    joinedancestry=JoinedPolyAncestry(polyped, markermap,parentgeno,offspringgeno,joinedspace,
        diploprob,genoprob, haploprob,parentdiploprob,parentgenoprob, parenthaploprob,        
        joined_viterbi,joined_parentviterbi)
    joinedancestry
end

function parseviterbistate(vp::AbstractString,inttype::DataType)
    Tuple(parse.(inttype,split(vp,"|")))
end

function phasedgeno2df(polyancestry::Union{PolyAncestry,JoinedPolyAncestry}; isparent::Bool=true)
    genols = isparent ? polyancestry.parentgeno : polyancestry.offspringgeno
    nch = length(genols)
    fhaplo = reduce(vcat, [stringjoin.(genols[ch],"|") for ch=1:nch])    
    vcatmap = Matrix(reduce(vcat, polyancestry.markermap))
    fhaplo2=hcat(vcatmap,fhaplo)
    
    indinfo = isparent ? polyancestry.polyped.parentinfo : polyancestry.polyped.offspringinfo
    DataFrame(fhaplo2,[propertynames(polyancestry.markermap[1]); Symbol.(indinfo[!,1])])
end

# function foundergeno2df(joinedancestry::JoinedPolyAncestry)
#     nch = length(joinedancestry.foundergeno)
#     fhaplo = reduce(vcat, [stringjoin.(joinedancestry.foundergeno[ch],"|") for ch=1:nch])    
#     vcatmap = Matrix(reduce(vcat, joinedancestry.markermap))
#     fhaplo2=hcat(vcatmap,fhaplo)    
#     pinfo = joinedancestry.polyped.parentinfo
#     fls = pinfo[pinfo[!,:isfounder],1]
#     DataFrame(fhaplo2,[propertynames(joinedancestry.markermap[1]); Symbol.(fls)])
# end

function joined_viterbipath2df(joinedancestry::JoinedPolyAncestry;
    ismidparent::Bool=false)    
    if ismidparent
        isnothing(joinedancestry.parentviterbipath) && return nothing
        vpath = join.(reduce(vcat,joinedancestry.parentviterbipath),"|")
        pinfo = joinedancestry.polyped.parentinfo
        indidls = pinfo[!, :individual]
    else
        isnothing(joinedancestry.viterbipath) && return nothing
        vpath = join.(reduce(vcat,joinedancestry.viterbipath),"|")
        indidls = joinedancestry.polyped.offspringinfo[!,:individual]        
    end
    vcatmap = reduce(vcat, joinedancestry.markermap)[!,1:5]
    vpathdf = DataFrame(vpath, Symbol.(indidls))
    hcat(vcatmap, vpathdf)
end


function joined_condprob2df(joinedancestry::JoinedPolyAncestry; probtype::AbstractString)
    vcatmap = reduce(vcat, joinedancestry.markermap)[!,1:5]
    if probtype == "diploprob"
        parentcondprob = joinedancestry.parentdiploprob
        condprob = joinedancestry.diploprob
        delimitor = "|"
    elseif probtype == "genoprob"
        parentcondprob = joinedancestry.parentgenoprob
        condprob = joinedancestry.genoprob
        delimitor = "/"
    elseif probtype == "haploprob"
        parentcondprob = joinedancestry.parenthaploprob
        condprob = joinedancestry.haploprob
        delimitor = ""
    else
        error("probtype must be diploprob, genoprob or haploprob")
    end
    #
    stringprob = probpair2str.(reduce(vcat,parentcondprob); delimitor)
    indidls = joinedancestry.polyped.parentinfo[!,:individual]
    parentdiploprobdf = DataFrame(stringprob, Symbol.(indidls))
    # 
    stringprob = probpair2str.(reduce(vcat,condprob); delimitor)
    indidls = joinedancestry.polyped.offspringinfo[!,:individual]
    diploprobdf = DataFrame(stringprob, Symbol.(indidls))
    hcat(vcatmap, parentdiploprobdf, diploprobdf)
end


function probpair2str(probpair; delimitor="|")
    ii, vv = probpair
    join([join(join.(ii,delimitor),"&"), join(vv,"&")],"=>")
end


function parse_probpairstr(probpairstr::AbstractString; delimitor="|", inttype=Int16)
    ii,vv = split(probpairstr,"=>")
    if isnothing(delimitor)
        ii2 = Tuple(parse.(inttype,split(ii, "&")))        
    else
        ii2 = Tuple([Tuple(parse.(inttype,split(i,delimitor))) for i in split(ii, "&")])
    end
    vv2 = Tuple(parse.(Float16,split(vv,"&")))
    ii2 => vv2
end

function parse_joinedcondprob(condprobdf::AbstractDataFrame, indidls::AbstractVector; 
    delimitor::Union{Nothing, AbstractString}, inttype::DataType)
    samplenames = names(condprobdf)[6:end] # col1-5: markermap
    issubset(indidls,samplenames) || return nothing
    dict = Dict(samplenames .=> eachindex(samplenames))    
    colls = vcat(1:5,5 .+ [dict[i] for i in indidls])
    condprobdf2 = condprobdf[!,colls]
    gdf = groupby(condprobdf2, :linkagegroup)
    [parse_probpairstr.(Matrix(df[!,6:end]); delimitor, inttype) for df in gdf]
end

function parse_joinedcondprob(condprobdf::AbstractDataFrame, polyped::PolyPed;
    delimitor::Union{Nothing, AbstractString}, inttype::DataType)
    indidls = polyped.parentinfo[!,:individual] 
    parentcondprob = parse_joinedcondprob(condprobdf,indidls; delimitor,inttype)
    indidls = polyped.offspringinfo[!,:individual] 
    condprob = parse_joinedcondprob(condprobdf,indidls; delimitor,inttype)
    parentcondprob, condprob
end

#############################################################################################

function get_joined_parentdiploprob(polyancestry::PolyAncestry; samplesize=1000)    
    polyped = polyancestry.polyped
    fgl_founder = get_fgl_founder(polyped)
    statespace_parent = get_statespace_parent(polyancestry)
    peddf_parent = get_peddf_parent(polyped)
    withmidparent = !isnothing(polyancestry.parentdiploprob)    
    # if no midparents, samplesize is not used!
    [begin 
        chrparentdiploprob = withmidparent ? polyancestry.parentdiploprob[chr] : nothing        
        nmarker = size(polyancestry.markermap[chr],1)
        get_joined_chrparentdiploprob(chrparentdiploprob,fgl_founder, 
            statespace_parent,peddf_parent,nmarker; samplesize)
    end for chr in eachindex(polyancestry.markermap)]
end

function get_joined_chrparentdiploprob(chrparentdiploprob,fgl_founder, 
    statespace_parent,peddf_parent,nmarker; samplesize=1000)        
    if isnothing(chrparentdiploprob)
        nparent = length(statespace_parent)    
        midindices = []
    else
        nparent = length(chrparentdiploprob)        
        midindices = findall(.!isnothing.(chrparentdiploprob))
    end    
    res = Matrix{Any}(undef, nmarker, nparent)
    for i in eachindex(fgl_founder)
        res[:,i] = repeat([Pair([fgl_founder[i]],ProbabilityWeights([1.0]))], inner=nmarker)
    end    
    ninf = 0 
    ntot = nmarker * nparent 
    for p in midindices        
        ss = statespace_parent[p]
        indprob = chrparentdiploprob[p]
        sparsels = sparse.(eachrow(indprob))
        res[:,p] .= [begin 
            ii,vv = findnz(sparsels[t])
            b = isnan.(vv)
            b2 = isinf.(vv)
            if any(b) || any(b2)
                ninf += 1
                vv0 = copy(vv)
                any(b) && (vv[b] .= 1/length(vv))            
                any(b2) && (vv[b2] .= 100)
                vv ./= sum(vv)
                @warn string("wegiths contain Inf or NaN:  oldweights=",vv0,",newweights=",vv) maxlog=1
            end
            Pair(ss[ii],ProbabilityWeights(vv))
        end for t in eachindex(sparsels)]
    end
    ninf > 0 && @warn string(ninf," weight vectors (out of ", ntot, ") containning Inf or NaN")
    if !isempty(midindices)
        for t in 1:nmarker
            postprob_parent = view(res,t,:)
            joinedspace_postprob_parent!(postprob_parent,peddf_parent, midindices; samplesize)
        end
    end
    res
end

function joinedspace_postprob_parent!(postprob_parent,peddf_parent, midindices; samplesize=1000)    
    samplels = [begin 
        state_parent = sample_state_parent(postprob_parent,midindices)
        lru_key = hash(state_parent)
        get!(lru_state_parent,lru_key) do 
            get_new_state_parent(state_parent, peddf_parent)                
        end
    end for _ in 1:samplesize]
    samplemtx = permutedims(reduce(hcat,samplels))
    for p = midindices
        freq = countmap(samplemtx[:,p])
        a = Tuple.(keys(freq))
        w = collect(values(freq))
        o = sortperm(w; rev=true)
        a2 = a[o]
        w2 = ProbabilityWeights(w[o] ./ sum(w))
        postprob_parent[p] = Pair(a2,w2)    
    end
    postprob_parent
end

function sample_state_parent(postprob_parent,midindices)
    res = copy(first.(first.(postprob_parent)))
    for p = midindices
        a,w = postprob_parent[p]
        res[p] = sample(a,w)
    end
    res
end


#############################################################################################


function get_joined_diploprob(polyancestry::PolyAncestry,joined_parentdiploprob::AbstractVector; 
    samplesize=1000,isindex::Bool=false,isparallel::Bool=true)    
    isnothing(polyancestry.diploprob) && return nothing
    statespace = polyancestry.statespace
    pinfo = polyancestry.polyped.parentinfo
    nallele = sum(pinfo[pinfo[!,:isfounder],:ploidy])
    [begin 
        joined_chrparentdiploprob = joined_parentdiploprob[chr]
        chrdiploprob = polyancestry.diploprob[chr]
        get_joined_chrdiploprob(chrdiploprob, joined_chrparentdiploprob,statespace,nallele; samplesize,isindex,isparallel)
    end for chr in eachindex(polyancestry.diploprob,joined_parentdiploprob)]
end

function get_joined_chrdiploprob(chrdiploprob, joined_chrparentdiploprob,statespace,nallele; 
    samplesize=1000,isindex::Bool=false,isparallel::Bool=true)
    noff = length(chrdiploprob)
    res = Vector{Any}(undef, noff)
    for popid = keys(statespace)        
        statels = Tuple.(statespace[popid]["state"])
        pls = statespace[popid]["parentindex"]
        isself = allequal(pls)
        ploidy_moth = statespace[popid]["parentploidy"][1]
        joined_moth = joined_chrparentdiploprob[:,pls[1]]
        joined_fath = joined_chrparentdiploprob[:,pls[2]]        
        isf1fam = length(unique(joined_moth)) == 1 && length(unique(joined_fath)) == 1
        offls = statespace[popid]["offspringindex"]
        funmap = isparallel && nworkers() > 1 ? pmap : map
        res[offls] .= funmap(x->get_joined_offdiploprob(x, statels, isf1fam, isself, ploidy_moth,
            joined_moth,joined_fath,nallele; samplesize,isindex),chrdiploprob[offls])
    end
    if noff == 1
        reshape(res[1],:, 1)
    else
        reduce(hcat,res)
    end
end


function get_joined_offdiploprob(offdiploprob, statels, isf1fam, isself, ploidy_moth,joined_moth,joined_fath,nallele; 
    samplesize=1000,isindex::Bool=false)     
    inttype = nallele < typemax(Int16) ? Int16 : Int32 
    offdiploprob2 = permutedims(offdiploprob)    
    if isf1fam
        # only one state if mother/father is a founder
        state_moth = only(first(only(unique(joined_moth))))
        state_fath = only(first(only(unique(joined_fath))))
        [begin             
            ii,vv = findnz(col)
            halfploidy = div(ploidy_moth,2)
            ii_state = [begin 
                s = statels[i]
                rg = [s[1:halfploidy]...]
                rg2 = [s[halfploidy+1:end]...]
                isself || (rg2 .-= ploidy_moth)   
                (state_moth[rg]...,state_fath[rg2]...)
            end for i in ii]
            oo = sortperm(vv;rev=true)
            ii_state2 = ii_state[oo]            
            g2 = isindex ? diploindex_2_index.(ii_state2,nallele) : Tuple(Tuple(inttype.(i))  for i in ii_state2)
            Pair(g2,vv[oo])
        end for col in eachcol(offdiploprob2)]
    else
        offprobls = [begin 
            ii,vv = findnz(col)
            Pair(ii, ProbabilityWeights(vv))
        end for col in eachcol(offdiploprob2)]
        resls = [begin 
            offpath = statels[sample_ind_path(offprobls)]
            joined_moth2 = sample_ind_path(joined_moth)
            joined_fath2 = sample_ind_path(joined_fath)
            get_new_offpath(offpath,isself, ploidy_moth,joined_moth2,joined_fath2)        
        end for i in 1:samplesize]
        resmtx = permutedims(reduce(hcat,resls))        
        [begin 
            freq = countmap(col)
            w = collect(values(freq))
            o = sortperm(w; rev=true)
            w2 = Tuple(Float16.(w[o] ./ sum(w)))
            g = collect(keys(freq))      
            g2 = isindex ? diploindex_2_index.(g[o],nallele) : Tuple(Tuple(inttype.(i))  for i in g[o])
            Pair(g2,w2)
        end for col in eachcol(resmtx)]
    end
end


function sample_ind_path(indprobls)
    [sample(i[1],i[2]) for i in indprobls]
end

#############################################################################################

function index_2_genoindex(index, ploidy, nallele; iszerobased::Bool=false)
    iszerobased || (index -= 1)
    genoindex = Int[]
    for _ in ploidy:-1:1    
        ls = [binomial(k + ploidy - 1, ploidy-1) for k in 0:nallele-1]    
        ls = accumulate(+,ls)
        # println("ploidy=", ploidy, ", ls=",ls, ",genoindex=",genoindex, ", index=",index)
        k = findfirst(x -> x > index, ls) - 1
        pushfirst!(genoindex, k)
        ploidy -= 1
        k >= 1 && (index -= ls[k])
    end
    index  == 0 || @error string("unexpected index=", index)
    iszerobased || (genoindex .+= 1)
    genoindex
end

function genoindex_2_index(genoindex::Union{Tuple,AbstractVector}; iszerobased::Bool=false)
    if iszerobased       
        # starting from 0
        sum(binomial(genoindex[i] + i -1, i) for i in eachindex(genoindex))
    else
        # starting from 1
        sum(binomial(genoindex[i] + i - 2, i) for i in eachindex(genoindex)) + 1
    end
end

function index_2_diploindex(index::Integer, ploidy::Integer, nallele::Integer;
    iszerobased::Bool=true)
    iszerobased || (index -= 1)
    diploindex = Int[]
    for i = ploidy:-1:1
        y = nallele^(i-1)
        k = div(index, y)
        pushfirst!(diploindex, k)        
        index -= k*y        
    end
    index == 0 || @error string("unexpected index=",index)
    iszerobased || (diploindex .+= 1)
    diploindex
end


function diploindex_2_index(diploindex::Union{Tuple,AbstractVector}, nallele::Integer;
    iszerobased::Bool=false)
    if iszerobased
        # starting from 0
        sum(nallele^(i-1) * diploindex[i] for i in eachindex(diploindex))
    else
        # starting from 1
        sum(nallele^(i-1) * (diploindex[i] -1) for i in eachindex(diploindex)) +1
    end
end

#############################################################################################