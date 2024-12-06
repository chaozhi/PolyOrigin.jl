

"""
    readTruegeno!(truefile,polyancestry, keyargs...)

returned truegeno of NamedTuple with truegeno.parentgeno storing phased parental genotypes,
and truegeno.offspringgeno storing true parental origins. `polyancestry` is modified
(1) extra markers in truefile are put into polyancestry.delmarker, (2)
absolute phase of polyancestry.parentgeno is set according to truegeno.parentgeno,
and (3) polyancestry.genoprob and polyancestry.haploprob are set to be consistent
with the absolute phase.

# Positional arguments

`truefile::AbstractString`: file storing true values, which has the same requirement
as the `genofile`. The parental genotypes must be given in the `phasedgeno`  format
with alleles being "1" and "2". The offspring origin-genotypes must be given in the
`phasedgeno`  format with allels being parental homologs.

`polyancestry::PolyAncestry`: for checking the consistency of marker IDs,
chromosome IDs, and individual IDs.

# Keyword arguments

`delimchar::AbstractChar=','`:  text delimiter.

`commentstring::AbstractString="#"`: rows that begins with commentstring will be ignored.

`workdir::AbstractString = pwd()`: directory for reading truevalue file.

"""
function readTruegeno!(truefile::AbstractString,polyancestry::PolyAncestry;
    delimchar::AbstractChar=',',
    missingstring::AbstractString="NA",
    commentstring::AbstractString="#",
    workdir::AbstractString=pwd())
    truefile2 =getabsfile(workdir,truefile)
    isfile(truefile2) || @error(string(truefile," does not exist in workdir = ",workdir))
    setAbsPhase!(truefile2,polyancestry, workdir=workdir,
        delim=delimchar,comment=commentstring,verbose=false)
    truedf=CSV.read(truefile2,DataFrame; delim=delimchar,comment=commentstring,missingstring=missingstring)
    for i=1:2
        truedf[!,i] =string.(strip.(string.(truedf[!,i])))
    end
    # marker indices
    dict=Dict(truedf[!,1] .=> 1:size(truedf,1))
    vcatmap = vcat(polyancestry.markermap...)
    iimarker=[get(dict,i,missing) for i=vcatmap[!,:marker]]
    bool=ismissing.(iimarker)
    any(bool) && @error string("markers that are not in truevalue file: ",truedf[iimarker[bool],1])
    # vcatmap[!,:marker] == truedf[iimarker,1]
    # check chromosome
    if strip.(vcatmap[!,:chromosome]) != truedf[iimarker,2]
        @error string("inconsistent chromosome between polyancestry and truefile")
    end
    # push extra markers of truedef into polyancestry.delmarker
    del=setdiff(1:size(truedf,1),iimarker)
    if !isempty(del)
        deldf = deepcopy(truedf[del,1:3])
        deldf = parsemarkermap!(deldf)
        b=[!in(i,polyancestry.delmarker[!,:marker]) for i=deldf[!,:marker]]
        polyancestry.delmarker =vcat(polyancestry.delmarker,deldf[b,:])
    end
    # get ind dict
    ind=string.(strip.(string.(names(truedf))))
    length(ind) < 5 && @error string("At least 5 coloumns are required")
    ind=ind[4:end]
    dict=Dict(ind .=> 1:size(ind,1))
    # parent indices
    jjparent=[get(dict,i,missing) for i=polyancestry.parentinfo[!,:individual]]
    bool=ismissing.(jjparent)
    any(bool) && @error string("offspring that are not in truevalue file: ",ind[jjparent[bool]])
    jjparent2=vcat(1:3,jjparent .+ 3)
    parentdf = truedf[iimarker,jjparent2]
    parentphase = [Matrix(i[!,4:end]) for i=groupby(parentdf,2)]
    parentphase2=[map(x->tryparse.(Int,replace(strip.(x),missingstring=>"-1")),split.(string.(i),"|")) for i=parentphase]
    parentphase3 = [replace.(i,-1=>missing) for i=parentphase2]
    # offsprign indices
    jjoff=[get(dict,i,missing) for i=polyancestry.offspringinfo[!,:individual]]
    bool=ismissing.(jjoff)
    any(bool) && @error string("offspring that are not in truevalue file: ",ind[jjoff[bool]])
    jjoff2=vcat(1:3,jjoff .+ 3)
    offdf = truedf[iimarker,jjoff2]
    offorig = [Matrix(i[!,4:end]) for i=groupby(offdf,2)]
    offorig2=[map(x->tryparse.(Int,replace(strip.(x),missingstring=>"-1")),split.(string.(i),"|")) for i=offorig]
    offorig3 = [replace.(i,-1=>missing) for i=offorig2]
    # markermap
    inputdf = rename(truedf[!,1:3],[:marker,:chromosome,:position])
    truemap=[DataFrame(i) for i=groupby(inputdf[!,1:3],2)]
    estmap=[DataFrame(i) for i=groupby(inputdf[iimarker,1:3],2)]
    (truemap=truemap,estmap=estmap,parentgeno=parentphase3,offspringgeno = offorig3)
end

function calparentacc(trueparentgeno::AbstractVector,estparentgeno::AbstractVector)
    if !all(size.(trueparentgeno) .== size.(estparentgeno))
        @error string("inconsistent size between true and estimated parentgeno")
    end
    d = trueparentgeno .- estparentgeno
    d2=[sum.(map(x->abs.(x),i)) for i=d]
    nalleleerr = sum([sum(skipmissing(i)) for i=d2])
    nphaseerr = sum([sum(skipmissing(i .> 0)) for i=d2])
    ndoseerr=sum(map((x,y)->sum(skipmissing(sum.(x) .!= sum.(y))),trueparentgeno,estparentgeno))
    nparentgeno = sum([sum(.!ismissing.(i)) for i=d2])
    (ndoseerr=ndoseerr,nphaseerr=nphaseerr,nalleleerr=nalleleerr,nparentgeno=nparentgeno)
end

function toindexancestry(trueoffgeno::AbstractVector, polyancestry::PolyAncestry)
    trueancestry = [sort.(i) for i=trueoffgeno]
    ancestry = [Matrix{Union{Int,Missing}}(missing,size(i)) for i=trueoffgeno]
    popls = polyancestry.designinfo[!,:population]
    nchr=length(trueancestry)
    for popid=popls
        gs = polyancestry.statespace[popid]["groupstate"]
        rule = Dict(gs .=> 1:length(gs))
        offls = findall(polyancestry.offspringinfo[!,:population] .== popid)
        for chr=1:nchr
            ancestry[chr][:,offls]=[get(rule,i, missing) for i=trueancestry[chr][:,offls]]
        end
    end
    # if any([ismissing.(i) for i=ancestry]...)
    #     wrong=unique(vcat([trueoffgeno[chr][ismissing.(ancestry[chr])] for chr=1:nchr]...))
    #     @info string("ancestry missing code: ", wrong)
    # end
    ancestry
end

function calancestryacc(trueoffgeno::AbstractVector, polyancestry::PolyAncestry)
    estgenoprob=polyancestry.genoprob
    trueoffancestry=toindexancestry(trueoffgeno,polyancestry)
    nchr = length(trueoffancestry)
    assignerr = [begin
        nsnp,noff = size(trueoffancestry[chr])
        acc = [begin
            trueanc = trueoffancestry[chr][:,off]
            prob = estgenoprob[chr][off]
            mls=findall(.!ismissing.(trueanc))
            isempty(mls) ? missing : mean([prob[m,trueanc[m]] for m=mls])
        end for off=1:noff]
        1-mean(skipmissing(acc))
    end for chr=1:nchr]
    # ancestrycall in polyreconstruct.jl
    doseprob = ancestrycall(estgenoprob,minprob=0.0)
    callerr=[begin
        nnonmiss=length(trueoffancestry[chr]) - sum(ismissing.(trueoffancestry[chr]))
        1.0-sum(doseprob[chr] .=== trueoffancestry[chr])/nnonmiss
    end for chr=1:nchr]
    (assignerr=round(mean(assignerr),sigdigits=6),callerr=round(mean(callerr),sigdigits=6))
end

function calabsaccuracy0(truegeno::NamedTuple,polyancestry::PolyAncestry)
    parentacc = calparentacc(truegeno.parentgeno,polyancestry.parentgeno)
    offacc=calancestryacc(truegeno.offspringgeno,polyancestry)
    ndel =size(polyancestry.delmarker,1)
    nrest = sum(size.(polyancestry.markermap,1))
    fdel = ndel/(ndel+nrest)
    fdel=round(fdel,sigdigits=6)
    merge(parentacc,offacc,(delfraction=fdel,))
end

"""
    calAccuracy!(truefile,polyancestry,workdir=pwd(),io=nothing,verbose=true)

calculate phasing and ancestral inference accuracies accoring to the true values in
`truefile`. The `polyancestry` is modified with absolute parental phases and
consistent genoprob.

# Positional arguments

`truefile::AbstractString`: contains true values for parental phases and
ancestral genotypes.

`polyancestry::PolyAncestry`: a struct returned by [`polyOrigin`](@ref).

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for writing outfile.

`io::Union{Nothing,IOStream}`: stream for writing log.

`verbose::Bool=true`: true if print messages on console.

"""
function calAccuracy!(truefile::AbstractString,polyancestry::PolyAncestry;
    workdir::AbstractString=pwd(),
    io::Union{Nothing,IOStream}=nothing,verbose::Bool=true)
    truegeno = readTruegeno!(truefile,polyancestry,workdir=workdir)
    calAccuracy!(truegeno,polyancestry)
end

"""
    calAccuracy!(truegeno,polyancestry,io=nothing,verbose=true)

calculate phasing and ancestral inference accuracies accoring to the true values in
`truegeno`. The `polyancestry` is modified with absolute parental phases and
consistent genoprob.

# Positional arguments

`truegeno::NamedTuple`: contains true values for parental phases and
ancestral genotypes, returned by [`readTruegeno!`](@ref).

`polyancestry::PolyAncestry`: a struct returned by [`polyOrigin`](@ref).

# Keyward arguments

`io::Union{Nothing,IOStream}`: stream for writing log.

`verbose::Bool=true`: true if print messages on console.

"""
function calAccuracy!(truegeno::NamedTuple, polyancestry::PolyAncestry;
    io::Union{Nothing,IOStream}=nothing,verbose::Bool=true)
    setabsphase0!(truegeno.parentgeno,polyancestry,io=io,verbose=verbose)
    calabsaccuracy0(truegeno,polyancestry)
end
