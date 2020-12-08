# """
#     parentCorrect!(polygeno::PolyGeno, polyancestry::PolyAncestry,keyargs...)
#
# performs parental genotype correction and return a dataframe
# listing of error corrections.
#
# # Positional arguments
#
# `polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.
#
# `chr::Integer`: chromosome index
#
# `chrgenoprob::AbstractVector`: ancestral probabilities for each offspring
# at chr-th chromsome
#
# `priorspace::AbstractDict`: state space for each subpopulation.
#
# `polyancestry::PolyAncestry`: a struct that stores results of ancestral inference.
#
# # Keyword arguments
#
# `correctthreshold::AbstractFloat=0.15`: a candidate marker is selected for
# parental error correction if the fraction of offspring genotypic error >= correctthreshold.
#
# """
function parentCorrect!(phasedgeno::PolyGeno,chr::Integer,
    chrgenoprob::AbstractVector,priorspace::AbstractDict;
    correctthreshold::AbstractFloat=0.15)
    # correctthreshold: fraction of offspring genotypes being incosistent
    # with  estimated parental phased genotypes
    if correctthreshold<1.0
        callthreshold=0.5
        chrdoseprob = calchrdoseprob(chrgenoprob,phasedgeno.parentgeno[chr],
            phasedgeno.offspringinfo,phasedgeno.designinfo,priorspace)
        chrcalldose = getchrcalldose(chrdoseprob,callthreshold=callthreshold)
        chrbadsnp = getchrbadsnp(phasedgeno,chr,chrcalldose,correctthreshold)
        # ["marker", "chromosome","parent", "old_genotype", "new_genotype","old_nerr", "new_nerr"]
        changedf = getchrbadchange!(chrgenoprob,chrbadsnp,priorspace,
            phasedgeno,chr,callthreshold=callthreshold)
    else
        res=[["marker", "chromosome", "parent", "old_genotype", "new_genotype",
            "old_nerr", "new_nerr"]]
        res2=permutedims(hcat(res...))
        changedf = DataFrame(res2[2:end,:],Symbol.(res2[1,:]))
    end
    changedf
end

function offspringCall(phasedgeno::PolyGeno,chr::Integer,
    chrgenoprob::AbstractVector,priorspace::AbstractDict;
    callthreshold::AbstractFloat=0.5)
    chrdoseprob = calchrdoseprob(chrgenoprob,phasedgeno.parentgeno[chr],
        phasedgeno.offspringinfo,phasedgeno.designinfo,priorspace)
    chrcalldose = getchrcalldose(chrdoseprob,callthreshold=callthreshold)
    chrcalldose
end

function getchangesum(changedf::DataFrame)
    res=[]
    for m=1:size(changedf,1)
        oldphase, newphase, oldnerr,  newnerr = Vector(changedf[m,4:end])
        newphase2 = parse.(Int,split(newphase,"|"))
        oldphase2 = parse.(Int,split(oldphase,"|"))
        nerrdiff = oldnerr - newnerr
        isdosediff = (sum(newphase2)-sum(oldphase2)) != 0
        nflip = sum(abs.(newphase2 - oldphase2))
        push!(res,[isdosediff,nflip,nerrdiff])
    end
    res2=hcat(res...)'
    ndosediff,nalleleflip,nerrdiff =sum(res2,dims=1)
    nphasediff = size(changedf,1)
    (ndosediff=ndosediff,nphasediff=nphasediff,nalleleflip=nalleleflip,nerrdiff=nerrdiff)
end

function getchrfhaplo(desinfo::DataFrame,chrparentgeno::AbstractMatrix)
    Dict([begin
        popid=desinfo[i,1]
        hh=chrparentgeno[:,Vector(desinfo[i,2:end]) .> 0]
        popid=>hcat([vcat(i...) for i=eachrow(hh)]...)'
    end for i=1:size(desinfo,1)])
end

function calchrdoseprob(chrgenoprob::AbstractVector,chrparentgeno::AbstractMatrix,
    offinfo::DataFrame,desinfo::DataFrame,priorspace::AbstractDict)
    chrfhaplo=getchrfhaplo(desinfo,chrparentgeno)
    offploidy=offinfo[!,:ploidy]
    offpop=offinfo[!,:population]
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

function getchrcalldose(chrdoseprob::AbstractMatrix;callthreshold::AbstractFloat=0.5)
    [begin
        val,index=findmax(i)
        val> callthreshold ? index-1 : missing
    end for i=chrdoseprob]
end

function gbserrorscore(x::AbstractVector,y::Integer)
    # x[y+1] ≤ 0.05
    x[y+1] ≤ 0.05 ? 1-x[y+1] : 0.0
end

function getchrbadsnp(phasedgeno::PolyGeno, chr::Integer,
    chrcalldose::AbstractMatrix,
    correctthreshold::AbstractFloat)
    chroffgeno = phasedgeno.offspringgeno[chr]
    bestdose = chrcalldose
    kind = kindofgeno([chroffgeno])
    nch = length(chroffgeno)
    if kind=="dosage"
        dosemismatch = map((x,y)->ismissing(x) || ismissing(y) ? missing : x!=y,chroffgeno,bestdose)
    elseif kind=="probability"
        # replace 0.05 with 10*seqerr? but the dose probability may also be calculated from SNP array data
        dosemismatch= map((x,y)->ismissing(y) ? missing : gbserrorscore(x,y),chroffgeno,bestdose)
    else
        @error string("unknown offspring data type: ",kind)
    end
    popparent = Dict([i[1]=>Int.(i[2:end].>0) for i=eachrow(Matrix(phasedgeno.designinfo))])
    offpop=phasedgeno.offspringinfo[!,:population]
    res=[]
    for popid = unique(offpop)
        offls = offpop .== popid
        popmismatch = dosemismatch[:,offls]
        popmismatch2=[collect(skipmissing(i)) for i= eachrow(popmismatch)]
        nlen=length.(popmismatch2)
        nerr = [isempty(i) ? 0 : sum(i) for i=popmismatch2]
        # nerr = sum.(popmismatch2)
        pos = findall(nerr .>= 3)
        ww = nerr[pos] ./ nlen[pos]
        bool = ww .> correctthreshold
        snps = pos[bool]
        ww = ww[bool]
        pp= popparent[popid]
        for i=1:length(ww)
            push!(res,[chr,snps[i], ww[i]*pp])
        end
    end
    res2=permutedims(hcat(res...))
    if isempty(res2)
        Matrix(undef,0,0)
    else
        res3=sortsplit(res2,2)
        res4=[[i[1,1],i[1,2],sum(i[:,3])] for i=res3]
        # return matrix with cols: ch_index, snp_index, parent_weights
        permutedims(hcat(res4...))
    end
end

function getfhaplomissp(phasedgeno::PolyGeno,p::Integer,chr::Integer,snp::Integer)
    # hhset, a vector of length 2^n, includes all possible phased genotypes of parent p witn ploidy n
    # fhaplodict, a dict for each subpopulation, all haplotype of two-parents at chr/snp
    chrparentgeno=phasedgeno.parentgeno[chr]
    popparent = Dict([i[1]=>findall(i[2:end].>0) for i=eachrow(Matrix(phasedgeno.designinfo))])
    ploidyls=phasedgeno.parentinfo[!,:ploidy]
    hhset=getdoserule(ploidyls[p])[missing]
    fhaplodict=Dict([begin
        parents = popparent[popid]
        if p  in parents
            hhset2=repeat(permutedims(chrparentgeno[snp,:]),length(hhset))
            hhset2[:,p ] = hhset
            hhset3=hhset2[:,parents]
            hhset4=hcat([vcat(i...) for i=eachrow(hhset3)]...)'
            popid=>hhset4
        else
            popid=>[]
        end
    end for popid = keys(popparent)])
    hhset, fhaplodict
end

function getdosemissp(chrgenoprob::AbstractVector,priorspace::AbstractDict,
    phasedgeno::PolyGeno,p::Integer,ch::Integer,snp::Integer;
    callthreshold::AbstractFloat=0.5)
    # return hhset,offls,estdose
    # hhset, a vector of length 2^n, includes all possible phased genotypes of parent p witn ploidy n
    # offls is a list of offspring indices whose parent is p
    # estdose matrix with size: 2^n x length(offls)
    offploidy=phasedgeno.offspringinfo[!,:ploidy]
    offpop=phasedgeno.offspringinfo[!,:population]
    nind = length(offpop)
    hhset, fhaplodict = getfhaplomissp(phasedgeno,p,ch,snp)
    res=Vector{Any}(missing,nind)
    for popid=keys(fhaplodict)
        fhaplo=fhaplodict[popid]
        isempty(fhaplo) && continue
        fhaplo =fhaplo .- 1
        states = priorspace[popid]["groupstate"]
        offls=findall(offpop .== popid)
        for off=offls
            I,V = findnz(chrgenoprob[off][snp,:])
            d=[sum(fhaplo[:,i],dims=2) for i=states[I]]
            dd=hcat(d...)
            prob = zeros(size(dd,1),offploidy[off]+1)
            for i=1:size(dd,1),j=1:size(dd,2)
                prob[i,dd[i,j]+1] += V[j]
            end
            res[off]=[begin
                val,ind = findmax(i)
                val > callthreshold ? ind-1 : missing
            end for i=eachrow(prob)]
        end
    end
    offls = findall(.!ismissing.(res))
    estdose = hcat(res[offls]...)
    hhset,offls,estdose
end

function getchrbadchange!(chrgenoprob::AbstractVector,chrbadsnp::AbstractMatrix,
    priorspace::AbstractDict, phasedgeno::PolyGeno,chr::Integer;
    callthreshold::AbstractFloat=0.5)
    chrobsgeno = phasedgeno.offspringgeno[chr]
    kind = kindofgeno([chrobsgeno])
    chrmarkermap = phasedgeno.markermap[chr]
    res=[]
    # println("chrbadsnp=", DataFrame(chrbadsnp))
    for m=1:size(chrbadsnp,1)
        snp, ww = chrbadsnp[m,[2,3]]
        pp = sortperm(ww,rev=true)
        pp =  pp[ww[pp] .> 0]
        # to update at the ordering of decreasing number of errors, for a given snp
        while true
            subres=[]
            for p=pp
                hhset, offls, estdose = getdosemissp(chrgenoprob,priorspace,phasedgeno,p,chr,snp,
                    callthreshold=callthreshold)
                obs=chrobsgeno[snp,offls]
                if kind=="dosage"
                    mismatch=[map((x,y)->ismissing(x) || ismissing(y) ? missing : x!=y, obs, i) for i=eachrow(estdose)]
                elseif kind=="probability"
                    mismatch=[map((x,y)->ismissing(y) ? missing : gbserrorscore(x,y), obs, i) for i=eachrow(estdose)]
                else
                    @error string("unknown offspring data type: ",kind)
                end
                mismatch2=[collect(skipmissing(i)) for i=mismatch]
                nlen=length.(mismatch2)
                nerr=sum.(mismatch2)
                newerr,newindex = findmin(nerr)
                old = phasedgeno.parentgeno[chr][snp,p]
                oldindex = findfirst([i == old for i=hhset])
                olderr = nerr[oldindex]
                minerrdiff = 3
                if olderr - newerr >= minerrdiff
                    new = hhset[newindex]
                    chid = chrmarkermap[1,:chromosome]
                    snpid = chrmarkermap[snp,:marker]
                    pid = phasedgeno.parentinfo[p,:individual]
                    push!(subres,[snpid,chid,pid,chr,snp,p,newerr,new,olderr,old])
                end
            end
            if length(subres)==0
                break
            else
                ndiff = [i[9]-i[7] for i=subres]
                ii = argmax(ndiff)
                push!(res,subres[ii])
                snp,p,new=subres[ii][[5,6,8]]
                phasedgeno.parentgeno[chr][snp,p] = new
                pp=pp[pp .!= p]
                # println("snp=",snp, ", pp=", pp, ", p=",p,", resp=",subres[ii])
                isempty(pp) && break
            end
        end
    end
    pushfirst!(res,["marker", "chromosome", "parent",
        "chrindex", "markerindex", "parentindex",
        "new_nerr","new_genotype", "old_nerr", "old_genotype"])
    res2=permutedims(hcat(res...))
    for j=[8,10]
        res2[2:end,j]=stringjoin.(res2[2:end,j],"|")
    end
    df=DataFrame(res2[2:end,:],Symbol.(res2[1,:]))
    df[!,[1,2,3,10,8,9,7]]
end
