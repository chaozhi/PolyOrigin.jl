## functions for updating chromosome pairings, gives parental phases
# a rule mapping dosages to phased genotype
function getdoserule(ploidy::Integer)
    genols=[begin
                ls=collect(powerset(1:ploidy,d,d))
                geno = [ones(Int8,ploidy) for i=1:length(ls)]
                for i=1:length(ls)
                    geno[i][ls[i]] .=2
                end
                geno
            end for d=0:ploidy]
    push!(genols,vcat(genols...))
    dls = Vector{Union{Missing,Int8}}(missing,ploidy+2)
    dls[1:ploidy+1] .= 0:ploidy
    Dict(dls .=> genols)
end

function calfhaploset(polygeno::PolyGeno, chr::Integer)
    fdose=polygeno.parentgeno[chr]
    ploidyls=polygeno.parentinfo[!,:ploidy]
    fhaploset=[begin
        doserule=getdoserule(ploidyls[f])
        [get(doserule, d, missing) for d = fdose[:,f]]
    end for f=1:length(ploidyls)]
    fhaploweight = [[repeat([1.0/length(i)],length(i)) for i=hh] for hh=fhaploset]
    fhaploset, fhaploweight
end


function callogbinprob(ploidy::Integer,seqerr::Real)
    p=[[dot([1-seqerr,seqerr], i) for i=[[1-i/ploidy,i/ploidy],[i/ploidy,1-i/ploidy]]] for i=0:ploidy]
    # return matrix size: 2x(ploidy+1)
    log.(reduce(hcat,p))
end

function calreadlogl(reads::AbstractVector,logbinprob::AbstractMatrix)
    # reads is a vector of integer
    ls = reads' * logbinprob
    ls[1,:]
end

function calreadpostprob(reads::AbstractVector,logbinprob::AbstractMatrix)
    # reads is a vector of integer
    # assume a discrete uniform prior
    ls= calreadlogl(reads,logbinprob)
    w = logsumexp(ls)
    exp.(ls .- w)
end

function calreadpostprob(reads::AbstractVector,ploidy::Integer,seqerr::Real)
    # reads is a vector of vector
    logp = callogbinprob(ploidy,seqerr)
    [any(ismissing.(i)) ? ones(ploidy+1)/(ploidy+1) : calreadpostprob(i,logp) for i=reads]
end

function calreadpostprob(reads::AbstractMatrix,ploidys::AbstractVector,
    seqerr::Real;digits=6)
    res=[calreadpostprob(reads[:,j],ploidys[j],seqerr) for j=1:length(ploidys)]
    map(x->round.(x,digits=digits),reduce(hcat,res))
end

function calfhaploset(doseprob::AbstractMatrix, ploidyls::AbstractVector)
    ind = [findall(i .> 0) for i=doseprob]
    prob = map((x,y)->x[y], doseprob, ind)
    fhaploset = [begin
        doserule=getdoserule(ploidyls[f])
        [vcat([get(doserule, i-1, missing) for i=d]...) for d = ind[:,f]]
    end for f=1:length(ploidyls)]
    fhaploweight = [begin
        n=ploidyls[f]
        doserule=getdoserule(n)
        funcount = [length(get(doserule, i,[])) for i=0:n]
        ncount = [funcount[i] for i= ind[:,f]]
        probf = prob[:,f]
        [vcat(map((x,y) -> repeat([x],y),probf[i],ncount[i])...) for i=1:length(ncount)]
    end for f=1:length(ploidyls)]
    fhaploweight=[[LinearAlgebra.normalize(i,1) for i=j] for j= fhaploweight]
    fhaploset, fhaploweight
end

function calparenthaploset(polygeno::PolyGeno, chr::Integer)
    kind = kindofgeno(polygeno.parentgeno)
    if kind == "dosage"
        calfhaploset(polygeno,chr)
    elseif kind == "probability"
        ploidyls = polygeno.parentinfo[!,:ploidy]
        doseprob = polygeno.parentgeno[chr]
        calfhaploset(doseprob, ploidyls)
    else
        @error string("unknow kind of parental geno: ",kind)
    end
end

function getfhaplo(fhaploindex::AbstractVector,fhaploset::AbstractVector)
    hh0=[map((x,y)->ismissing(y) ? x[1] .* missing : x[y],fhaploset[i],fhaploindex[i]) for i=1:length(fhaploset)]
    hh1 = reduce(hcat,hh0)
    hh2=[vcat(i...) for i=eachrow(hh1)]
    # return a matrix with size rxc, r=# markers for a given linkage group
    # c = sum of ploidy over founders
    reduce(hcat,hh2)'
end

function getderiveddose(fhaplo::AbstractMatrix, priorspace::AbstractDict,
    polygeno::PolyGeno,popidls=keys(priorspace))
    design=polygeno.designinfo
    nn0 = polygeno.parentinfo[!,:ploidy]
    nn1=accumulate(+, nn0)
    pushfirst!(nn1,0)
    rgls=[nn1[i]+1:nn1[i+1] for i=1:length(nn1)-1]
    Dict([begin
            states = get(priorspace, popid, missing)["state"]
            pp0 = design[design[!,:population] .== popid,2:end]
            pp1 = Vector(pp0[1,:])
            pp2 =[i for i=1:length(pp1) for j=1:pp1[i]]
            pp=vcat(rgls[pp2]...)
            subfhaplo = fhaplo[:,pp] .- Int8(1)
            derived= reduce(hcat,[Matrix{Int8}(sum(subfhaplo[:,j],dims=2)) for j=states])
            # derived is a matrix with size rxc, r=#markers, c=#states
            # matrix element is an integer doseage
            popid => derived
    end for popid = popidls])
end

function caldataprob_prob(probls::AbstractVector,truedose::AbstractArray,
    ploidy::Integer,doseerr::Real)
    probls2 = probls[truedose .+ 1]
    # (1-doseerr)*P_i + doseerr/n * \sum_{j \neq i}P_j
    (1-doseerr-doseerr/ploidy) .* probls2 .+ sum(probls)*doseerr/ploidy
end

# for one offspring in one linkage group
# offdose is a list of read pairs (AbstractVector) or a list of dosages (Integer)
# popderiveddose = deriveddose[popid], #snps x #states
function caldataprob(offdose::AbstractVector,ploidy::Integer,
    popderiveddose::AbstractMatrix,doseerrls::AbstractVector)
    ddose = popderiveddose
    type = eltype(typeof(offdose))
    nsnp, nstate = size(ddose)
    if type <: AbstractVector
        # if derived dose is missing, all offspring doses are missing,thus like =1
        res = [ones(Float64,nstate) for i=1:nsnp]
        for i=1:nsnp
            ismiss = any(ismissing.(ddose[i,:])) || ismissing(doseerrls[i])
            if !ismiss
                res[i]=caldataprob_prob(offdose[i],ddose[i,:],ploidy,doseerrls[i])
            end
        end
        return res
    elseif type <: Union{Missing,Integer}
        res=[[begin
                d = ddose[i,j] - offdose[i]
                cond = ismissing(d) || ismissing(doseerrls[i])
                cond ? 1.0 : (d == 0 ? (1-doseerrls[i]) : (doseerrls[i]/ploidy))
            end for j=1:nstate] for i=1:nsnp]
        return Vector{Vector{Float64}}(res)
    else
        @error string("unknow offspring genodata type: ",type)
    end
end

function upoffmarglogl!(loglgrid::AbstractMatrix,
    dataprob::AbstractVector,popid::AbstractString,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    bvprop::AbstractVector,
    snporder::AbstractVector)
    condstates = priorspace[popid]["condstate"]
    keyls = priorspace[popid]["valentkey"]
    loglgrid[bvprop] = [if ismissing(loglgrid[bv])
        markerincl, startprob, tranprobseq = strkey2prior(priorprocess, keyls[bv])
        # dataprobseq = [i[condstates[bv]] for i = dataprob[markerincl]]
        tseq = findall(markerincl)
        snpseq = snporder[tseq]
        dataprobseq = [i[condstates[bv]] for i=dataprob[snpseq]]
        calloglike(startprob, tranprobseq, dataprobseq)
    else
        loglgrid[bv]
    end for bv=bvprop]
end

# for one offspring in one linkage group
function caloffmarglogl(dataprob::AbstractVector,popid::AbstractString,
    priorspace::AbstractDict,priorprocess::AbstractDict,bvprop::AbstractVector,
    snporder::AbstractVector)
    condstates = priorspace[popid]["condstate"]
    keyls = priorspace[popid]["valentkey"]
    [begin
        markerincl, startprob, tranprobseq = strkey2prior(priorprocess, keyls[bv])
        tseq = findall(markerincl)
        snpseq = snporder[tseq]
        dataprobseq = [i[condstates[bv]] for i=dataprob[snpseq]]
        calloglike(startprob, tranprobseq, dataprobseq)
    end for bv=bvprop]
end


function calmarglogl(doseerrls::AbstractVector,deriveddose::AbstractDict,
    chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,bvpairprop::AbstractVector;
    snporder::Union{Nothing,AbstractVector}=nothing)
    noff = size(chrdose,2)
    res = Vector{Union{Missing,Vector}}(missing,noff)
    isnothing(snporder) && (snporder=1:size(chrdose,1))
    for popid = keys(deriveddose)
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        for off=offls
            # println("off = ", off)
            offdose = chrdose[:,off]
            # popid = polygeno.offspringinfo[off,:population]
            ploidy = polygeno.offspringinfo[off,:ploidy]
            offdataprob = caldataprob(offdose,ploidy,deriveddose[popid],doseerrls)
            res[off] = caloffmarglogl(offdataprob, popid, priorspace,priorprocess,bvpairprop[off],snporder)
        end
    end
    res
end

function calmarglogl(bvpair::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    doseerr::Real,chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    fhaplo=getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    bvpairprop = [[i] for i=bvpair]
    logllist = calmarglogl(doseerr,deriveddose,chrdose,priorspace,
        priorprocess,polygeno,bvpairprop)
    sum(vcat(logllist...))
end

## functions for updating parental phases, given chromosome pairing

function getfhaplofindex(findex::Integer,fhaploindex::AbstractVector,fhaploset::AbstractVector)
    getfhaplofindex([findex],fhaploindex, fhaploset)
end

function getfhaplofindex(findex::AbstractVector,fhaploindex::AbstractVector,fhaploset::AbstractVector)
    hh0=[map((x,y)->ismissing(y) ? x[1] .* missing : x[y],fhaploset[i],fhaploindex[i]) for i=1:length(fhaploset)]
    hh = reduce(hcat,hh0)
    ff0=fhaploset[findex]
    ff = reduce(hcat,ff0)
    nsnp = size(hh,1)
    nfindex = length(findex)
    if nfindex == 1
        hhset0 = [[(hh[m,findex] = [i]; vcat(hh[m,:]...)) for i=ff[m,1]] for m=1:nsnp]
    elseif nfindex ==2
        hhset0 = [[(hh[m,findex] = [i,j]; vcat(hh[m,:]...)) for i=ff[m,1] for j=ff[m,2]] for m=1:nsnp]
    else
        @error(string("too many #parents = ", nfindex, " in the ", popindex, "-th subpopulation"))
    end
    # return res=a list of matrix with length = #markers; res[i] size: rxc
    # r =# possible founder haplotypes at marker i; c= sum of ploidy over founders
    [reduce(hcat,i)' for i = hhset0]
end


# list of all subpopulations whose parnets contain findex
function popfromfindex(findex::Integer,polygeno::PolyGeno)
     popfromfindex([findex],polygeno)
end

# list of all subpopulations whose parnets are contained in the vector findex
function popfromfindex(findex::AbstractVector,polygeno::PolyGeno)
    polygeno.designinfo[Bool.(sign.(sum(Matrix(polygeno.designinfo[:,findex .+ 1]),dims=2)[:,1])),:population]
end

# findex: Integer or AbstractVector
function getderiveddose(findex,fhaploindex::AbstractVector,
    fhaploset::AbstractVector, priorspace::AbstractDict,
    polygeno::PolyGeno,bvpair::AbstractVector)
    # fhaplo =a list of matrices with length = #markers for a given linkage group
    # fhaplo[i] size: rxc; r =# founder haplotypes at marker i
    # c= sum of ploidy over founders
    fhaplo=getfhaplofindex(findex,fhaploindex, fhaploset)
    popidls = popfromfindex(findex,polygeno)
    design=polygeno.designinfo
    nn0 = polygeno.parentinfo[!,:ploidy]
    nn1=accumulate(+, nn0)
    pushfirst!(nn1,0)
    rgls=[nn1[i]+1:nn1[i+1] for i=1:length(nn1)-1]
    Dict([begin
        # caculate subfhaplo
        pp0 = design[design[!,:population] .== popid,2:end]
        pp1 = Vector(pp0[1,:])
        pp2 =[i for i=1:length(pp1) for j=1:pp1[i]]
        pp=vcat(rgls[pp2]...)
        subfhaplo = [i[:,pp] .- 1 for i=fhaplo]
        # calculate states
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        condstate=priorspace[popid]["condstate"]
        states0 = priorspace[popid]["state"]
        # derived is a matrix with size rxc, r=#markers, c=#states
        # matrix element is a vector of doseages, correspoinng each founder haplotype
        statecodes = sort(unique(vcat(condstate[bvpair[offls]]...)))
        states=states0[statecodes]
        nrow =size(subfhaplo,1)
        ncol = length(statecodes)
        I =repeat(1:nrow,outer=ncol)
        J=repeat(statecodes,inner=nrow)
        V0=[Vector{Int8}(sum(subfhaplo[i][:,j],dims=2)[:,1]) for i=1:nrow,j=states]
        V=reshape(V0,:)
        derived=sparse(I,J,V,nrow,length(states0))
        popid => derived
    end for popid = popidls])
end

function calsitedataprob(snp::Integer,findex, doseerr::Real,
    bvpair::AbstractVector,deriveddose::AbstractDict,chrdose::AbstractMatrix,
    priorspace::AbstractDict,polygeno::PolyGeno)
    popidls = popfromfindex(findex,polygeno)
    # return dataprob for the offspring of parents[findex], at a site snp
    # dataprobls is a list of prob matrix for each offspring
    # i-th row of the matrix corresonds to i-th  founder haplotype
    # the ij element is the prob of that offspring for  founder haplotype i and state j
    [begin
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        ploidy=polygeno.offspringinfo[first(offls),:ploidy]
        offdose = chrdose[snp,offls]
        offspace=priorspace[popid]["condstate"][bvpair[offls]]
        type = eltype(typeof(offdose))
        derived = deriveddose[popid][snp,:]
        if type <: Union{Missing,Integer}
            [begin
                ddose=reduce(hcat,derived[offspace[i]])
                ismiss = any(ismissing.(ddose)) || ismissing(offdose[i])
                if ismiss
                    ones(Float64,size(ddose))
                else
                    b = ddose .== offdose[i]
                    prob = ((1-doseerr) .* b) .+ ((doseerr/ploidy) .* (1 .- b))
                    Matrix{Float64}(prob)
                end
            end for i=1:length(offls)]
        elseif type <: AbstractVector
            [begin
                ddose=reduce(hcat,derived[offspace[i]])
                ismiss = any(ismissing.(ddose))
                if ismiss
                    ones(Float64,size(ddose))
                else
                    prob = caldataprob_prob(offdose[i],ddose,ploidy,doseerr)
                    Matrix{Float64}(prob)
                end
            end for i=1:length(offls)]
        else
            @error string("unknow offspring genodata type: ",type)
        end
    end for popid = popidls]
end

function calsitedataprobls(snp::Integer,doseerr::Real,
    bvpair::AbstractVector,deriveddose::AbstractDict,chrdose::AbstractMatrix,
    priorspace::AbstractDict,phasedgeno::PolyGeno)
    findex=1:size(phasedgeno.parentinfo,1)
    dataprob = calsitedataprob(snp,findex,doseerr,bvpair,deriveddose,chrdose,
        priorspace,phasedgeno)
    popidls=popfromfindex(findex,phasedgeno)
    noff=size(chrdose,2)
    res=Vector{Vector{Float64}}(undef,noff)
    for i=1:length(popidls)
        offls = findall(phasedgeno.offspringinfo[!,:population] .== popidls[i])
        res[offls] = [j[1,:] for j=dataprob[i]]
    end
    res
end
