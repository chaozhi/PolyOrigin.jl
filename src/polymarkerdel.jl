
function calsinglelogl(chrdose::AbstractMatrix, deriveddose::AbstractDict,
    bvpair::AbstractVector,priorspace::AbstractDict,priorprocess::AbstractDict,
    doseerrls::AbstractVector, polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing)
    isnothing(snporder) && (snporder = 1:size(chrdose,1))
    popidls = polygeno.designinfo[!,:population]
    markerincl = first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    noff = size(polygeno.offspringinfo,1)
    singlelogl = Matrix{Union{Missing,Float64}}(missing,length(markerincl),noff)
    for popid = popidls
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        # ploidy is the same for individuals in a population
        ploidy = polygeno.offspringinfo[first(offls),:ploidy]
        condstates = priorspace[popid]["condstate"]
        keyls = priorspace[popid]["valentkey"]
        for off = offls
            bv=bvpair[off]
            startprob = only(getstartprobls(priorprocess,[keyls[bv]]))
            ddose = deriveddose[popid][:, condstates[bv]]
            dataprob = caldataprob(chrdose[:,off],ploidy,ddose,doseerrls)
            dataprob = dataprob[snporder[tseq]]
            singlelogl[tseq,off] = [log(dot(startprob,i)) for i=dataprob]
        end
    end
    singlelogl
end

function callogbackward(chrdose::AbstractMatrix, deriveddose::AbstractDict,
    bvpair::AbstractVector,priorspace::AbstractDict,priorprocess::AbstractDict,
    doseerrls::AbstractVector, polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing)
    popidls = polygeno.designinfo[!,:population]
    res=Vector{Matrix{Union{Float64,Missing}}}(undef,size(polygeno.offspringinfo,1))
    nsnp = size(chrdose,1)
    isnothing(snporder) && (snporder = 1:nsnp)
    for popid = popidls
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        condstates = priorspace[popid]["condstate"]
        keyls = priorspace[popid]["valentkey"]
        for off = offls
            bv = bvpair[off]
            markerincl, startprob, tranprobseq = strkey2prior(priorprocess, keyls[bv])
            ploidy = polygeno.offspringinfo[off,:ploidy]
            ddose = deriveddose[popid][:, condstates[bv]]
            dataprob = caldataprob(chrdose[:,off],ploidy,ddose,doseerrls)
            dataprob = dataprob[snporder]
            dataprobseq = dataprob[markerincl]
            bw0 = logbackward(tranprobseq, dataprobseq)
            bw=reduce(hcat,bw0)'
            resbw = Matrix{Union{Float64,Missing}}(missing,nsnp,size(bw,2))
            resbw[markerincl,:] = bw
            res[off]=resbw
        end
    end
    res
end

function calinitforward(tinit::Integer,dataprobls::AbstractVector,
    bvpair::AbstractVector, priorspace::AbstractDict,
    priorprocess::AbstractDict, phasedgeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing)
    nsnp = length(first(values(priorprocess)).markerid)
    snporder == nothing && (snporder=1:nsnp)
    bvkeyls = getbvkeyls(bvpair,priorspace,phasedgeno)
    startprobls = getstartprobls(priorprocess,bvkeyls)
    popidls = phasedgeno.designinfo[!,:population]
    noff = size(phasedgeno.offspringinfo,1)
    fwprob=Vector{Matrix{Union{Float64,Missing}}}(undef,noff)
    fwlogl=Vector{Vector{Union{Float64,Missing}}}(undef,noff)
    for popid = popidls
        offls = findall(phasedgeno.offspringinfo[!,:population] .== popid)
        for off = offls
            dataprob = dataprobls[off]
            startprob = startprobls[off]
            ls = startprob .* dataprob
            scale = sum(ls)
            prob = Matrix{Union{Float64,Missing}}(missing,nsnp,length(startprob))
            prob[tinit,:] = ls ./ scale
            fwprob[off] = prob
            logscale = Vector{Union{Float64,Missing}}(missing,nsnp)
            logscale[tinit] = log(scale)
            fwlogl[off] = logscale
        end
    end
    Dict("fwprob"=>fwprob,"fwlogl"=>fwlogl)
end

function calnextforward!(fwdict::AbstractDict,tnow::Integer,tnext::Integer,
    dataprobls::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict, priorprocess::AbstractDict,phasedgeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    tnowdis::Union{Nothing,Real}=nothing)
    fwprob=fwdict["fwprob"]
    fwlogl=fwdict["fwlogl"]
    nsnp=length(fwlogl[1])
    snporder == nothing && (snporder=1:nsnp)
    if tnowdis == nothing
        tranprobdict = Dict([strkey=>priorprocess[strkey].tranprobseq[tnow]
            for (strkey, pri) in priorprocess])
    else
        tranprobdict = Dict([strkey=>getgametetran(tnowdis,pri.nvalent)
            for (strkey, pri) in priorprocess])
    end
    popidls = phasedgeno.designinfo[!,:population]
    for popid = popidls
        offls = findall(phasedgeno.offspringinfo[!,:population] .== popid)
        keyls = priorspace[popid]["valentkey"]
        for off = offls
            fwoff = fwdict["fwprob"][off]
            dataprob = dataprobls[off]
            tranprob = [tranprobdict[i] for i= split(keyls[bvpair[off]],"|")]
            pnow =fwoff[tnow,:]
            ls= kronvec(tranprob[1]', tranprob[2]', pnow) .* dataprob
            scale = sum(ls)
            fwoff[tnext,:] = ls ./ scale
            fwoff[tnow+1:tnext-1,:] .= missing
            fwlogl[off][tnext] = fwlogl[off][tnow] + log(scale)
            fwlogl[off][tnow+1:tnext-1] .=  missing
        end
    end
    fwdict
end

function calindlogl(fwdict::AbstractDict,logbw::AbstractVector,t::Integer)
    fwprob=fwdict["fwprob"]
    fwlogl=fwdict["fwlogl"]
    [begin
        lp=logbw[i][t,:]
        lpmax=max(lp...)
        log(dot(exp.(lp .- lpmax), fwprob[i][t,:])) + lpmax + fwlogl[i][t]
    end for i=1:length(logbw)]
end

function calvuongts(logldiff::AbstractVector,dfdiff::Integer)
    n = length(logldiff)
    v = var(logldiff)
    vuong = sum(logldiff)
    # correction for the difference of number of parameters
    vuong -= log(n)*(dfdiff)/2
    vuong /= sqrt(n*v)
    vuong
end

function calvuongts!(kk::Integer,tseq::AbstractVector,logl::AbstractVector,
    singlelogl::AbstractMatrix,logbw::AbstractVector,fwdict::AbstractDict,
    bvpair::AbstractVector, priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, deriveddose::AbstractDict,
    doseerrls::AbstractVector, polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    caldistance::Bool=false)
    isnothing(snporder) && (snporder = 1:size(chrdose,1))
    # fwdict is modified inside of function
    dfdiff = -1
    isologl = singlelogl[tseq[kk],:]
    b=.!isnan.(isologl)
    if kk==1
        dataprobls = calsitedataprobls(tseq[2],doseerrls[tseq[2]],bvpair,
            deriveddose,chrdose,priorspace,polygeno)
        propfwdict = calinitforward(tseq[2],dataprobls,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
        proplogl0 = calindlogl(propfwdict,logbw,tseq[2])
        proplogl = proplogl0[b] .+ isologl[b]
    elseif kk==length(tseq)
        proplogl0 = [i[tseq[kk-1]] for i=fwdict["fwlogl"]]
        proplogl = proplogl0[b] .+ isologl[b]
    else
        dataprobls = calsitedataprobls(tseq[kk+1],doseerrls[tseq[kk+1]],bvpair,
            deriveddose,chrdose,priorspace,polygeno)
        if caldistance
            bvkeyls = getbvkeyls(bvpair,priorspace,polygeno)
            tnowdis= first(calinterdis(tseq[kk-1],tseq[kk+1],priorprocess,bvkeyls,
                dataprobls,fwdict,logbw))
        else
            pri1=first(values(priorprocess))
            tnowdis=sum(pri1.markerdeltd[tseq[kk-1:kk]])
        end
        calnextforward!(fwdict,tseq[kk-1],tseq[kk+1],dataprobls,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder,tnowdis=tnowdis)
        proplogl0 = calindlogl(fwdict,logbw,tseq[kk+1])
        proplogl = proplogl0[b] .+ isologl[b]
        dataprobls = calsitedataprobls(tseq[kk],doseerrls[tseq[kk]],bvpair,
            deriveddose,chrdose,priorspace,polygeno)
        calnextforward!(fwdict,tseq[kk-1],tseq[kk],dataprobls,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
    end
    logldiff = proplogl .- logl[b]
    isempty(logldiff) ? 0 : calvuongts(logldiff,dfdiff)
end

function polymarkerdel!(fhaploindex::AbstractVector,fhaploset::AbstractVector,
    bvpair::AbstractVector, priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, doseerrls::AbstractVector, polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    delsiglevel::Real=0.05,
    caldistance::Bool=false)
    fhaplo= getfhaplo(fhaploindex, fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    isnothing(snporder) && (snporder = 1:size(chrdose,1))
    singlelogl= calsinglelogl(chrdose, deriveddose, bvpair,priorspace,priorprocess,
        doseerrls, polygeno,snporder=snporder)
    logbw = callogbackward(chrdose, deriveddose, bvpair,priorspace,priorprocess,
        doseerrls, polygeno; snporder=snporder)
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    # TODO: cal dataprobls
    dataprobls = calsitedataprobls(tseq[1],doseerrls[tseq[1]],bvpair,deriveddose,
        chrdose,priorspace,polygeno)
    fwdict = calinitforward(tseq[1],dataprobls,bvpair,priorspace,priorprocess,
        polygeno,snporder=snporder)
    logl = calindlogl(fwdict,logbw, tseq[1])
    res=Vector{Union{Missing,Float64}}(missing,length(markerincl))
    res[tseq] = [calvuongts!(kk,tseq,logl,singlelogl,logbw,fwdict,
        bvpair,priorspace,priorprocess,chrdose, deriveddose, doseerrls, polygeno,
        snporder=snporder,caldistance=caldistance) for kk=1:length(tseq)]
    threshold = abs(quantile(Normal(),delsiglevel))
    delsnp = tseq[res[tseq] .> threshold]
    markerincl[delsnp] .= false
    for (key, val) in priorprocess
         setmarkerincl!(val,markerincl)
    end
    delsnp
end

function polymarkerdel!(bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, doseerrls::AbstractVector, polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    delsiglevel::Real=0.05,
    caldistance::Bool=false,
    chrindex::Integer=1)
    fphase = polygeno.parentgeno[chrindex]
    fhaploset=[[[i] for i=j] for j=eachcol(fphase)]
    fhaploindex=[ones(Int, length(j)) for j=fhaploset]
    polymarkerdel!(fhaploindex,fhaploset,bvpair,
        priorspace, priorprocess,chrdose, doseerrls, polygeno;
        snporder,delsiglevel,caldistance)
end
