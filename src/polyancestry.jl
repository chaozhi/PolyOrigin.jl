"""
    PolyAncestry

mutable struct the results of haplotype reconstruction.

    PolyAncestry(polyped, markermap,parentgeno,offspringgeno, 
        statespace,valentprob,parentvalentprob, 
        diploprob,genoprob,haploprob,
        parentdiploprob, parentgenoprob, parenthaploprob,
        viterbipath, parentviterbipath, misc)

Constructor from the results of haplotype reconstruction.

# Fields

`polyped::PolyPed`: pedigree information

`markermap::Vector{DataFrame}`: marker map for each chromosome. `markermap[c]`
gives the dataframe of chromosome c with columns [:marker, :linkagegroup, :poscm].

`parentgeno::Vector{Matrix}`: genotypic data in parents. `parentgeno[c][m,i]`
gives the phased genotype of parent i at marker m of chromosome c.

`offspringgeno::Vector{Matrix}`: genotypic data in offspring. `offspringgeno[c][m,o]`
gives the phased genotype of offspring o at marker m of chromosome c.

`statespace::Dict`: specify valents and origin-genotypes for each population.
Each key is a population id, and its value is in turn a dict with keys:
"parent", "parentindex", "valent", "valentneighbor", and "groupstate".

`valentprob::Union{Nothing,Vector}`, posterior valent probability for each offspring in each chromosome.
valentprob[c][o] for offspring o in chromsome c is a matrix with three columns being
valentindex, loglike, and posterior probability. Here posterior probability is calcuated from loglike,
assuming a discrete uniform prior distribution of valent configurations.

`genoprob::Vector`: posterior origin-genotype probability for each offspring
in each chromosome. genoprob[c][o] for offspring o in chromsome c is a sparse matrix,
with element (m,s) being the probability of origin-genotype s at marker m.

`haploprob::Union{Nothing,Vector}`: posterior origin-haplotype probability for each offspring
in each chromosome. haploprob[c][o] for offspring o in chromsome c is a sparse matrix,
with element (m,s) being the probability of origin-haplotype s at marker m.

`misc::Dict{String, DataFrame}`: contains information such as (1) "deletion" dataframe
for markers that were removed from markermap, (2) "correction" dataframe for
parental error correction.

"""
mutable struct PolyAncestry
    polyped::PolyPed    
    markermap::Vector{DataFrame}    
    parentgeno::Vector{Matrix}        
    offspringgeno::Union{Nothing,Vector{Matrix}}        
    statespace::AbstractDict    
    valentprob::Vector{Vector{Matrix}} # valentprob[chr]: chr, ind, pop, valentindex, valent, loglike, valentprob    
    parentvalentprob::Union{Nothing,Vector{Vector}} 
    diploprob::Vector{Vector{SparseMatrixCSC}} # diplorob[chr][offspring]: a sparse matrix with dimension nsnp x nstate    
    genoprob::Union{Nothing, Vector}
    haploprob::Union{Nothing,Vector}     
    parentdiploprob::Union{Nothing, Vector}
    parentgenoprob::Union{Nothing, Vector}
    parenthaploprob::Union{Nothing,Vector}     
    viterbipath::Vector{Matrix}
    parentviterbipath::Union{Nothing, Vector}
    misc::Dict{String, DataFrame}
    function PolyAncestry(polyped, markermap,parentgeno,offspringgeno,statespace,        
        valentprob,parentvalentprob, 
        diploprob, genoprob,haploprob, 
        parentdiploprob, parentgenoprob, parenthaploprob,
        viterbipath, parentviterbipath, misc)        
        dims = [length(markermap),size(parentgeno,1),size(diploprob,1)]
        if length(union(dims))!=1
            error("inconsistent number of linkage groups ",dims)
        end
        dims = [size.(markermap,1),size.(parentgeno,1),size.(first.(diploprob),1)]
        if length(union(dims))!=1
            error("inconsistent number of markers within linkage groups ",dims)
        end
        dims = [union(length.(diploprob)),[size(polyped.offspringinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of offspring ",dims)
        end
        dims = [union(size.(parentgeno,2)), [size(polyped.parentinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of parents ",dims)
        end
        new(polyped, markermap,parentgeno,offspringgeno,statespace,
            valentprob,parentvalentprob, 
            diploprob, genoprob,haploprob,
            parentdiploprob, parentgenoprob, parenthaploprob, 
            viterbipath, parentviterbipath, misc)
    end
end


"""
    savepolyancestry(outfile,polyancestry,missingstring="NA",workdir=pwd())

saves polyancestry into outfile in CSV format.

# Positional arguments

`outfile::AbstractString`: file for saving polyancestry in the directory workdir.
The saved outfile consists of several tables: the first row of each table has two cells:
[polyorigin-PolyAncestry, nameoftable], and the rest rows denote a dataframe with column names.
The following is a list of table names and their descriptions.

* `pedigree`: pedigree information 

* `valentlist`: list of possible bi- or multi-valent formations for each population.
   The dataframe has columns [population, parentindex, parent, valentindex, valent].

* `valentprob`: posterior valent probability for each offspring in each chromosome.
   The dataframe has columns
   [chromosome, individual, population, valentindex, valent, loglike, valentprob].

* `parentgeno`: phased parental genotypes. The dataframe has columns
   [marker, chromosome, position, parent1, parent2, ...].

* `ancestralgenotype`: list of origin-genotypes for each population. The dataframe
   has columns [population, parentindex, parent, stateindex, state]

* `genoprob`: marginal posterior probabilities for origin-genotypes. For each offspring
   at each marker, the posterior probability vector is represented by a sparse vector
   in the form of I=>V: state index vector I and non-zero probability vector V such
   that I[k]=V[k] for k=1, ..., K. The elements of vector I or V are delimited by "|".

`polyancestry::PolyAncestry`: results of haplotype reconstruction returned
 from [`polyorigin`](@ref).

# Keyward arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for writing outfile.

"""
function savepolyancestry(outfile::AbstractString,polyancestry::PolyAncestry;
    missingstring::AbstractString="NA",
    workdir::AbstractString = pwd())
    ext = last(splitext(outfile))
    myopen = ext == ".gz" ? GZip.open : open    
    outputfile2 = getabsfile(workdir,outfile)
    myopen(outputfile2, "w") do io
        initial = "PolyOrigin"        
        write(io,initial,",pedigree","\n")
        savepolyped(io,polyancestry.polyped; workdir,delim=',')
        append_dataframe(io, polyancestry.polyped.offspringinfo; 
            missingstring,dfname = "offspringinfo")               

        # valent
        dfvalents, dfgroupstates, dfstates = statespace2df(polyancestry; missingstring)
        append_dataframe(io, dfvalents; missingstring,dfname="valentlist")             
        df = valentprob2df(polyancestry; ismidparent = true)
        df2 = valentprob2df(polyancestry; ismidparent = false)
        df = isnothing(df) ? df2 : vcat(df,df2)
        append_dataframe(io, df; missingstring,dfname="valentprob")                

        # parentgeno and states
        append_dataframe(io, phasedgeno2df(polyancestry; isparent=true); missingstring,dfname="parentgeno")  
        append_dataframe(io, phasedgeno2df(polyancestry; isparent=false); missingstring,dfname="offspringgeno")  
         
        append_dataframe(io, dfgroupstates; missingstring,dfname="unphasedstate")     
        append_dataframe(io, dfstates; missingstring,dfname="phasedstates")     
        
        # viterbipath
        df = viterbipath2df(polyancestry; ismidparent=true)
        df2 = viterbipath2df(polyancestry; ismidparent=false)
        df =  isnothing(df) ? df2 : hcat(df, df2[!,6:end]) # col1-5: markermap
        append_dataframe(io, df; missingstring,dfname="viterbipath")               

        # condprob 
        for probtype = ["diploprob","genoprob","haploprob"]
            df = condprob2df(polyancestry; ismidparent = true,  probtype)
            df2 = condprob2df(polyancestry; ismidparent = false,  probtype)
            df =  isnothing(df) ? df2 : hcat(df, df2[!,6:end]) # col1-5: markermap
            dfname = probtype == "diploprob" ? "phasedgenoprob" : (probtype == "genoprob" ? "unphasedgenoprob" : "haploprob")
            append_dataframe(io, df; missingstring,dfname)               
        end    
        # misc
        misc = polyancestry.misc
        if !isempty(misc)
            for (key,val) in misc
                key == "filecomment" && continue         
                append_dataframe(io, val; missingstring, dfname=string("misc|",key))                
            end
        end
        flush(io)
    end        
    outputfile2
end

function append_dataframe(io::IO, df::AbstractDataFrame; 
    initial="PolyOrigin", dfname="table",
    missingstring="NA", 
    delim=",")
    write(io,initial,delim,dfname,"\n")
    CSV.write(io,df; header=true,append=true,delim,missingstring)    
end



function statespace2df(polyancestry::PolyAncestry; missingstring::AbstractString)
    statespace=polyancestry.statespace
    popls = polyancestry.polyped.designinfo[!,:family]
    resvalents = []    
    resgroupstates = []
    resstates = []
    for pop = popls        
        chrpairing = statespace[pop]["chrpairing"]
        # valent
        valents0 = [[unique(j) for j=i] for i = statespace[pop]["valent"]]
        # gameteindex=[string(i,"|",j) for i=1:size(valents0,1), j=1:size(valents0,2)]        
        valents = vec(valents0)        
        v = [ismissing(k) ? missingstring : stringjoin([stringjoin(stringjoin.(j,":"),"-") for j=k],"|") for k=valents]
        valentindex = collect(eachindex(v))
        parent = stringjoin(statespace[pop]["parent"],"|")
        parentindex = stringjoin(statespace[pop]["parentindex"],"|")        
        resv = hcat(repeat([pop chrpairing parentindex parent],length(v)),valentindex, v)
        push!(resvalents,resv)
        # groupstate
        s = stringjoin.(statespace[pop]["groupstate"],"/")
        ns = length(s)
        ress = hcat(repeat([pop parentindex parent],ns),1:ns,s)
        push!(resgroupstates,ress)
        # state
        s = stringjoin.(statespace[pop]["state"],"|")
        ns = length(s)
        ress = hcat(repeat([pop parentindex parent],ns),1:ns,s)
        push!(resstates,ress)
    end
    dfvalents = DataFrame(reduce(vcat,resvalents),Symbol.(["family", "chrpairing", "parentindex", "parent", "valentindex", "valent"]))
    dfgroupstates = DataFrame(reduce(vcat,resgroupstates),Symbol.(["family", "parentindex", "parent", "stateindex", "unphasedstate"]))
    dfstates = DataFrame(reduce(vcat,resstates),Symbol.(["family","parentindex", "parent", "stateindex", "phasedstate"]))
    dfvalents,dfgroupstates, dfstates
end

function condprob2df(polyancestry::PolyAncestry; ismidparent::Bool=false,
    probtype::AbstractString="diploprob")
    if ismidparent
        pinfo = polyancestry.polyped.parentinfo
        indidls = pinfo[!, :individual]
        indls = findall(.!pinfo[!,:isfounder])
        isempty(indls) && return nothing                
        condprob = probtype == "diploprob" ? polyancestry.parentdiploprob : (probtype == "genoprob" ? polyancestry.parentgenoprob : polyancestry.parenthaploprob)
        isnothing(condprob) && return nothing
    else        
        condprob = probtype == "diploprob" ? polyancestry.diploprob : (probtype == "genoprob" ? polyancestry.genoprob : polyancestry.haploprob)
        indidls = polyancestry.polyped.offspringinfo[!,:individual]
        indls = eachindex(indidls)
    end    
    vcatmap = reduce(vcat, polyancestry.markermap)[!,1:5]
    nsnp,nmapcol = size(vcatmap)
    res= Matrix{Any}(missing,nsnp,length(indidls)+nmapcol)
    res[:,1:nmapcol] .= Matrix(vcatmap)
    nchr=length(condprob)
    for ind in eachindex(indidls)
        in(ind, indls) || continue
        indpr=[begin
            prob = condprob[ch][ind]
            [begin
                I, V = findnz(prob[m,:])
                join([join(I,"&"),join(V,"&")],"=>")
            end for m = 1:size(prob,1)]
        end for ch=1:nchr]
        res[:,ind+nmapcol] =  reduce(vcat, indpr)
    end
    rowid = [propertynames(vcatmap); Symbol.(indidls)]
    DataFrame(res,rowid)
end


function viterbipath2df(polyancestry::PolyAncestry;
    ismidparent::Bool=false)    
    if ismidparent
        pinfo = polyancestry.polyped.parentinfo
        indidls = pinfo[!, :individual]
        indls = findall(.!pinfo[!,:isfounder])
        isempty(indls) && return nothing        
        vpath = polyancestry.parentviterbipath
        isnothing(vpath) && return nothing
    else        
        vpath = polyancestry.viterbipath
        isnothing(vpath) && return nothing
        indidls = polyancestry.polyped.offspringinfo[!,:individual]
        indls = eachindex(indidls)
    end        
    vcatmap = reduce(vcat, polyancestry.markermap)[!,1:5]
    hcat(vcatmap, DataFrame(reduce(vcat,vpath),Symbol.(indidls)))
end

function valentprob2df(polyancestry::PolyAncestry;
    ismidparent::Bool=false)
    statespace = polyancestry.statespace
    if ismidparent
        pinfo = polyancestry.polyped.parentinfo
        indidls = pinfo[!, :individual]        
        indls = findall(.!pinfo[!,:isfounder])
        isempty(indls) && return nothing
        indpopls = pinfo[!,:family]
        valentprob = polyancestry.parentvalentprob
        isnothing(valentprob) && return nothing
    else
        valentprob = polyancestry.valentprob
        indidls = polyancestry.polyped.offspringinfo[!,:individual]
        indls = eachindex(indidls)
        indpopls = polyancestry.polyped.offspringinfo[!,:family]
    end
    markermap= polyancestry.markermap
    chrls = [markermap[i][1,2] for i=1:length(markermap)]
    rowid=["chromosome", "individual", "family", "valentindex", "valent", "loglike", "valentprob"]    
    res = []
    for chr= eachindex(chrls)
        for ind= eachindex(indidls)
            if in(ind, indls) 
                vv = valentprob[chr][ind]
                valentindex, logl, postprob = [stringjoin(vv[:,j],"&") for j=1:3]
                valent00 = statespace[indpopls[ind]]["valent"]
                valent0=  [[unique(j) for j=i] for i =valent00]
                valent = valent0[vv[:,1]]
                valent2 = [stringjoin([stringjoin([stringjoin(j,":") for j=k],"-") for k=i],"|") for i=valent]
                valent3=stringjoin(valent2,"&")
                resrow=[chrls[chr],indidls[ind],indpopls[ind],valentindex,valent3,logl,postprob]                
            else
                resrow=[chrls[chr],indidls[ind],indpopls[ind], missing, missing, missing, missing]                
            end            
            push!(res,resrow)
        end
    end
    res2= permutedims(reduce(hcat,res))
    DataFrame(res2,Symbol.(rowid))
end

function parse_saved_genoprob(saved_genoprob::AbstractDataFrame, polyped::PolyPed; 
    chrpairing, missingstring)
    nparent = size(polyped.parentinfo,1)
    noffspring = size(polyped.offspringinfo,1)
    nmidparent = sum(.!polyped.parentinfo[!,:isfounder])
    nmapcol = get_nmapcol(saved_genoprob)
    nprobcol = nmapcol+nparent+noffspring
    nmidparent > 0 && (nprobcol += nparent)
    if size(saved_genoprob,2) != nprobcol
        error(string(genoprobfile, " contains unexpected #column = ",size(saved_genoprob,2)))
    end    
    markermap = saved_genoprob[!, 1:nmapcol]
    parsemarkermap!(markermap)
    markermap2=[DataFrame(i) for i=groupby(markermap,:linkagegroup)]    

    # parentgeno
    dictcol = Dict(names(saved_genoprob) .=> 1:size(saved_genoprob,2))
    pls = polyped.parentinfo[!,:individual]
    pcols = [dictcol[i] for i in pls]
    snpindexlist=splitindex(markermap[!,:linkagegroup])
    parentgeno = [Matrix(saved_genoprob[i,pcols]) for i=snpindexlist]
    parentgeno2=parseinputgeno(parentgeno; missingstring)

    # condprob
    priorspace = get_priorstatespace(polyped,chrpairing)
    statespace = getstatespace(priorspace)
    nstaterule=Dict([key=>length(statespace[key]["groupstate"]) for key=keys(statespace)])
    # condprob_leafoffspring
    offls = polyped.offspringinfo[!,:individual]
    offcols = [dictcol[string(i, "_condprob")] for i in offls]
    if isempty(offcols)
        offgenoprob2 = nothing
    else
        offgenoprob = [Matrix(saved_genoprob[i,offcols]) for i=snpindexlist]    
        nstatels=[nstaterule[i] for i=polyped.offspringinfo[!,:family]]
        offgenoprob2=[begin
            mtx=offgenoprob[ch]
            [sparse(reduce(hcat,parseprobcell.(mtx[:,off],nstatels[off]))') for off=eachindex(nstatels)]
        end for ch = eachindex(offgenoprob)]
    end
    # condprob_midparentparent
    pcols = [dictcol[string(i, "_condprob")] for i in pls]
    if isempty(pcols)
        midparentgenoprob2 = nothing
    else
        parentgenoprob = [Matrix(saved_genoprob[i,pcols]) for i=snpindexlist]
        nstatels=[i[:isfounder] ? nothing : (haskey(nstaterule,i[:family]) ? nstaterule[i[:family]] : nothing) for i= eachrow(polyped.parentinfo)]
        midparentgenoprob2=[begin
            mtx=parentgenoprob[ch]
            [if isnothing(nstatels[ind])
                nothing
            else
                sparse(reduce(hcat,parseprobcell.(mtx[:,ind],nstatels[ind]))') 
            end for ind=eachindex(nstatels)]
        end for ch = eachindex(parentgenoprob)]
    end
    statespace, markermap2, parentgeno2, midparentgenoprob2, offgenoprob2
end



"""
    readpolyancestry(ancestryfile, missingstring="NA",workdir=pwd())

return a struct of type PolyAncestry from  ancestryfile in the
    directory workdir.

# Positional argument

`ancestryfile`: file storing polyancestry that is generated by [`savepolyancestry`](@ref)
It is one of output files by [`polyorigin`](@ref).

# Keyward arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for reading ancestryfile,

"""
function readpolyancestry(ancestryfile::AbstractString;
    missingstring::AbstractString="NA",
    workdir::AbstractString=pwd())
    ancestryfile2 = getabsfile(workdir,ancestryfile)
    isfile(ancestryfile2) || @error(string(ancestryfile," does not exist in workdir = ",workdir))
    res = readmultitable(ancestryfile2; missingstring)
    haskey(res, "pedigree") || error(string("key genodata is not in ", keys(res)))  
    ped = Pedigree(res["pedigree"])        
    polyped = PolyPed(ped) 
    if haskey(res, "offspringinfo")                   
        offinfo = res["offspringinfo"]
        offinfo[!,:family] .= string.(offinfo[!,:family])
        polyped = PolyPed(polyped.ped, polyped.designinfo,polyped.parentinfo,offinfo)
    end
    markermap, parentgeno = parsephasedgeno(res["parentgeno"])     
    if haskey(res, "offspringgeno")   
        markermap2, offspringgeno = parsephasedgeno(res["offspringgeno"])        
    else
        offspringgeno = nothing
    end
    valentdf = res["valentlist"]
    chrpairingls = unique(valentdf[!,"chrpairing"])
    length(chrpairingls) == 1 || @error "non-unique chrpairing!" valentdf
    priorspace = get_priorstatespace(polyped,chrpairingls[1])
    statespace = getstatespace(priorspace)
    if haskey(res, "valentprob") 
        parentvalentprob, valentprob = parsevalentprob(res["valentprob"],polyped) 
    else
        parentvalentprob, valentprob = nothing, nothing
    end    
    if haskey(res, "phasedgenoprob") || haskey(res, "diploprob") 
        probkey = haskey(res, "phasedgenoprob") ? "phasedgenoprob" : "diploprob"
        parentdiploprob, diploprob = parsecondprob(res[probkey],polyped, statespace; probtype="diploprob")
    else
        parentdiploprob, diploprob = nothing, nothing
    end
    if haskey(res, "unphasedgenoprob") || haskey(res, "genoprob") 
        probkey = haskey(res, "unphasedgenoprob") ? "unphasedgenoprob" : "genoprob"
        parentgenoprob, genoprob = parsecondprob(res[probkey],polyped, statespace; probtype="genoprob")        
    else
        parentgenoprob, genoprob = nothing, nothing
    end
    if haskey(res, "viterbipath") 
        parentviterbipath, viterbipath = parseviterbipath(res["viterbipath"],polyped)
    else
        parentviterbipath, viterbipath = nothing, nothing
    end
    misc = Dict{String, DataFrame}()
    for (key,val) in res
        if key[1:5]=="misc|"                         
            push!(misc, key[6:end]=>val)            
        end
    end
    polyancestry=PolyAncestry(polyped, markermap,parentgeno, offspringgeno, statespace,        
        valentprob,parentvalentprob,
        diploprob, genoprob,nothing, parentdiploprob, parentgenoprob,nothing,
        viterbipath, parentviterbipath, misc)    
    sethaploprob!(polyancestry)
    polyancestry
end


function parseviterbipath(vpathdf::AbstractDataFrame, indidls::AbstractVector)
    samplenames = names(vpathdf)[6:end] # col1-5: markermap
    issubset(indidls,samplenames) || return nothing
    dict = Dict(samplenames .=> eachindex(samplenames))    
    colls = vcat(1:5,5 .+ [dict[i] for i in indidls])
    vpathdf2 = vpathdf[!,colls]
    gdf = groupby(vpathdf2, :linkagegroup)
    [Matrix(df[!,6:end]) for df in gdf]
end

function parseviterbipath(vpathdf::AbstractDataFrame, polyped::PolyPed) 
    indidls = polyped.parentinfo[!,:individual] 
    parentviterbipath = parseviterbipath(vpathdf,indidls)
    indidls = polyped.offspringinfo[!,:individual] 
    viterbipath = parseviterbipath(vpathdf,indidls)
    parentviterbipath, viterbipath
end

function get_nmapcol(df::AbstractDataFrame)
    namels = names(df)
    nmapcol = findfirst(==("doublereduction"),namels)
    if isnothing(nmapcol ) 
        nmapcol =  _col_1stsample - 1
    else
        nmapcol ==  _col_1stsample || @error string("unexpected doublereduction col=",nmapcol)
    end        
    nmapcol
end


function parsephasedgeno(parentgeno)    
    nmapcol = get_nmapcol(parentgeno)
    markermap = parentgeno[:,1:nmapcol]
    for i=1:2
        markermap[!,i]=string.(markermap[!,i])
    end
    groupdf = groupby(parentgeno, :linkagegroup)
    parentgeno = [begin
        haplo = split.(Matrix(i[!,nmapcol+1:end]),"|")
        # [[x=="missing" ? missing : parse(Int,x) for x in h] for h in haplo]
        [mytryparse.(Int, h) for h in haplo]
    end for i = groupdf]
    markermap2=[DataFrame(i) for i=groupby(markermap,:linkagegroup)]
    markermap2, parentgeno
end

mytryparse(T, str) = something(tryparse(T, str), missing)

function parseprobcell(x::AbstractString,nstate::Integer)
    ls=split.(split(x,"=>"),"&")
    sparsevec(parse.(Int,ls[1]),parse.(Float64,ls[2]),nstate)
end


function parsecondprob(condprobdf::DataFrame,nstatels::AbstractVector)
    nmapcol = 5
    groupdf = groupby(condprobdf,2)    
    [begin
        mtx=Matrix(groupdf[ch][:,nmapcol+1:end])
        [if isnothing(nstatels[ind])
            nothing
        else
            sparse(reduce(hcat,parseprobcell.(mtx[:,ind],nstatels[ind]))') 
        end for ind = eachindex(nstatels)]
    end for ch=1:length(groupdf)]
end


function parsecondprob(condprobdf::DataFrame,polyped::PolyPed, 
    statespace::AbstractDict; probtype::AbstractString)
    if probtype == "diploprob" 
        statekind = "state" 
    elseif probtype == "genoprob" 
        statekind = "groupstate"
    else
        error(string("unknown probtype: ", probtype))
    end
    fam2nstate=Dict([key=>length(statespace[key][statekind]) for key=keys(statespace)])
    samplenames = names(condprobdf)[6:end] # cols1-5: markermap
    sampledict = Dict(samplenames .=> eachindex(samplenames))

    parentidls = polyped.parentinfo[!,:individual]
    if issubset(parentidls,samplenames)
        cols = vcat(1:5, 5 .+ [sampledict[i] for i in parentidls])
        df = condprobdf[!,cols]
        nstatels=[get(fam2nstate,i,nothing) for i = polyped.parentinfo[!,:family]]
        parentcondprob = parsecondprob(df, nstatels)
    else
        parentcondprob = nothing
    end

    offspringidls = polyped.offspringinfo[!,:individual]
    if issubset(offspringidls,samplenames)
        cols = vcat(1:5, 5 .+ [sampledict[i] for i in offspringidls])
        df = condprobdf[!,cols]
        nstatels=[get(fam2nstate,i,nothing) for i = polyped.offspringinfo[!,:family]]
        condprob = parsecondprob(df, nstatels)
    else
        condprob = nothing
    end
    parentcondprob, condprob
end

function parsevalentprob(valentprob::AbstractDataFrame, polyped::PolyPed)
    parentidls = polyped.parentinfo[!,:individual]
    offspringidls = polyped.offspringinfo[!,:individual]
    df = filter(row->in(row["individual"],parentidls), valentprob)
    parentvalentprob = parsevalentprob(df)
    df = filter(row->in(row["individual"],offspringidls), valentprob)
    valentprob = parsevalentprob(df)
    parentvalentprob, valentprob
end

function parsevalentprob(valentprob)
    groupdf = groupby(valentprob,1)
    [begin
        # col 4,6,7= valuentindex, loglike,valentprob
        mtx=Matrix{Any}(groupdf[ch][!,[4,6,7]])        
        for j=1:3
            mtx[:,j]=map(x->ismissing(x) ? missing : parse.(Float64, split(string(x),"&")),mtx[:,j])
        end
        [if any(ismissing.(i))
            nothing
        else
            vp=Matrix{Real}(reduce(hcat,i))
            vp[:,1]=Int.(vp[:,1])
            vp
        end for i=eachrow(mtx)]
    end for ch=1:length(groupdf)]
end

function repeatmultivalent(zygotevalent::AbstractVector)
    [reduce(vcat,[repeat([j],div(length(j), 2)) for j=i]) for i=zygotevalent]
end


function setgenoprob!(polyancestry::PolyAncestry)    
    dict2group = getdict2group(polyancestry.statespace,polyancestry.polyped)    
    if !isnothing(polyancestry.diploprob)
        polyancestry.genoprob = [begin
            [begin
                indpop = polyancestry.polyped.offspringinfo[ind,:family]
                tran2group = dict2group[indpop]                
                round.(chgeno[ind] * tran2group,digits=3)
            end for ind = eachindex(chgeno)]
        end for chgeno = polyancestry.diploprob]
    end
    if !isnothing(polyancestry.parentdiploprob)
        polyancestry.parentgenoprob = [begin             
            [if polyancestry.polyped.parentinfo[ind,:isfounder]
                nothing
            else 
                indpop = polyancestry.polyped.parentinfo[ind,:family]
                tran2group = dict2group[indpop]                
                round.(chgeno[ind] * tran2group,digits=3)
            end for ind=eachindex(chgeno)]
        end for chgeno = polyancestry.parentdiploprob]
    end
    nothing
end

function sethaploprob!(polyancestry::PolyAncestry)
    dict=Dict{String,SparseMatrixCSC}()
    if !isnothing(polyancestry.genoprob)
        polyancestry.haploprob = [begin
            [begin
                indpop = polyancestry.polyped.offspringinfo[ind,:family]
                indkey = string("offspring_",indpop)
                if !haskey(dict,indkey)
                    genotypes = polyancestry.statespace[indpop]["groupstate"]
                    g2h=geno2haplo(genotypes)
                    push!(dict,indkey =>sparse(g2h))
                end
                mtx = dict[indkey]
                chgeno[ind] * mtx
            end for ind = eachindex(chgeno)]
        end for chgeno = polyancestry.genoprob]
    end
    if !isnothing(polyancestry.parentgenoprob)
        polyancestry.parenthaploprob = [begin             
            [if polyancestry.polyped.parentinfo[ind,:isfounder]
                nothing
            else 
                indpop = polyancestry.polyped.parentinfo[ind,:family]
                indkey = string("parent_",indpop)
                if !haskey(dict,indkey)
                    genotypes = polyancestry.statespace[indpop]["groupstate"]
                    g2h=geno2haplo(genotypes)
                    push!(dict,indkey =>sparse(g2h))
                end
                mtx = dict[indkey]
                chgeno[ind] * mtx
            end for ind=eachindex(chgeno)]
        end for chgeno = polyancestry.parentgenoprob]
    end
    nothing
end

function geno2haplo(genotypes::AbstractVector)
    nrow =size(genotypes,1)
    ncol = maximum(reduce(vcat,genotypes))
    mtx = zeros(nrow,ncol)
    ploidy = length(genotypes[1])
    for i=1:nrow for j=1:ploidy
            mtx[i,genotypes[i][j]] += 1
        end
    end
    mtx ./= ploidy
    sparse(mtx)
end

function getsubpolyancestry(polyancestry::PolyAncestry;
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing)   
    isnothing(chrsubset) && return nothing
    polyped = deepcopy(polyancestry.polyped)
    nchr=length(polyancestry.markermap)    
    chrsubset2 = chrsubset[chrsubset.<= nchr]
    markermap = deepcopy(polyancestry.markermap[chrsubset2])
    parentgeno= deepcopy(polyancestry.parentgeno[chrsubset2])    
    statespace = deepcopy(polyancestry.statespace)
    valentprob = deepcopy(polyancestry.valentprob[chrsubset2])
    parentvalentprob = deepcopy(polyancestry.valentprob[chrsubset2])    
    # offspring
    if isnothing(polyancestry.diploprob) # diplorob[chr][offspring]: a sparse matrix with dimension nsnp x nstate    
        diploprob = nothing
    else
        diploprob = deepcopy(polyancestry.diploprob[chrsubset2])
    end
    if isnothing(polyancestry.genoprob)
        genoprob = nothing
    else
        genoprob = deepcopy(polyancestry.genoprob[chrsubset2])
    end
    if isnothing(polyancestry.haploprob)
        haploprob = nothing
    else
        haploprob = deepcopy(polyancestry.haploprob[chrsubset2])
    end
    # parent
    if isnothing(polyancestry.parentdiploprob) # diplorob[chr][offspring]: a sparse matrix with dimension nsnp x nstate    
        parentdiploprob = nothing
    else
        parentdiploprob = deepcopy(polyancestry.parentdiploprob[chrsubset2])
    end
    if isnothing(polyancestry.parentgenoprob)
        parentgenoprob = nothing
    else
        parentgenoprob = deepcopy(polyancestry.parentgenoprob[chrsubset2])
    end
    if isnothing(polyancestry.parenthaploprob)
        parenthaploprob = nothing
    else
        parenthaploprob = deepcopy(polyancestry.parenthaploprob[chrsubset2])
    end
    if isnothing(polyancestry.viterbipath)
        viterbipath = nothing
    else        
        viterbipath = deepcopy(polyancestry.viterbipath[chrsubset2])
    end
    if isnothing(polyancestry.parentviterbipath)
        parentviterbipath = nothing
    else        
        parentviterbipath = deepcopy(polyancestry.parentviterbipath[chrsubset2])
    end
    PolyAncestry(polyped, markermap,parentgeno,statespace,        
        valentprob,parentvalentprob, 
        diploprob, genoprob,haploprob, 
        parentdiploprob, parentgenoprob, parenthaploprob,
        viterbipath, parentviterbipath, deepcopy(polyancestry.misc))      
end