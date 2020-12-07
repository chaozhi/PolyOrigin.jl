

function calsinglelogl(fhaplophase::AbstractMatrix,dataprobset::AbstractVector,
    bvpair::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing)
    snporder==nothing && (snporder = 1:size(dataprobset[1],1))
    popidls = polygeno.designinfo[!,:population]
    markerincl = first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    noff = size(polygeno.offspringinfo,1)
    singlelogl = Matrix{Union{Missing,Float64}}(missing,length(markerincl),noff)
    for t=tseq
        parentphase = fhaplophase[snporder[t],:]
        for popid = popidls
            offls = findall(polygeno.offspringinfo[!,:population] .== popid)
            ismono=all([unique(parentphase[i]) in [[1],[2]] for i=priorspace[popid]["parentindex"]])
            if ismono
                singlelogl[t,offls] .= NaN
            else
                # ploidy is the same for individuals in a population
                ploidy = polygeno.offspringinfo[first(offls),:ploidy]
                condstates = priorspace[popid]["condstate"]
                keyls = priorspace[popid]["valentkey"]
                singlelogl[t,offls]=[begin
                    bv=bvpair[off]
                    pri = priorprocess[keyls[bv]]
                    dataprob = dataprobset[off][snporder[t],condstates[bv]]
                    startprob = pri.startprob
                    log(dot(startprob,dataprob))
                end for off=offls]
            end
        end
    end
    singlelogl
end

function callogbackward(dataprobset::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing)
    popidls = polygeno.designinfo[!,:population]
    res=Vector{Matrix{Union{Float64,Missing}}}(undef,size(polygeno.offspringinfo,1))
    nsnp = size(dataprobset[1],1)
    snporder == nothing && (snporder = 1:nsnp)
    for popid = popidls
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        condstates = priorspace[popid]["condstate"]
        keyls = priorspace[popid]["valentkey"]
        for off = offls
            bv=bvpair[off]
            pri = priorprocess[keyls[bv]]
            dataprob = dataprobset[off][snporder,:]
            incl = pri.markerincl
            dataprobseq = [dataprob[i,condstates[bv]] for i=findall(incl)]
            tranprobseq=Vector{Matrix{Float64}}(pri.tranprobseq[incl][1:end-1])
            bw0 = logbackward(tranprobseq,dataprobseq)
            bw=hcat(bw0...)'
            resbw = Matrix{Union{Float64,Missing}}(missing,nsnp,size(bw,2))
            resbw[pri.markerincl,:] = bw
            res[off]=resbw
        end
    end
    res
end


function getdataprobls(snp::Union{Integer,AbstractVector},
    dataprobset::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,phasedgeno::PolyGeno)
    if typeof(snp) <: Integer
        dataprobls = Vector{Vector{Float64}}(undef,length(dataprobset))
    elseif typeof(snp) <: AbstractVector
        dataprobls = Vector{Matrix{Float64}}(undef,length(dataprobset))
    else
        @error string("wrong snp=",snp)
    end
    popidls = phasedgeno.designinfo[!,:population]
    for popid=popidls
        offls = findall(phasedgeno.offspringinfo[!,:population] .== popid)
        condstates = priorspace[popid]["condstate"]
        for off=offls
            ss=condstates[bvpair[off]]
            dataprobls[off]=dataprobset[off][snp,ss]
        end
    end
    dataprobls
end


function calinitforward(tinit::Integer,dataprobinfo::AbstractVector,
    bvpair::AbstractVector, priorspace::AbstractDict,
    priorprocess::AbstractDict, phasedgeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing)
    nsnp = length(first(values(priorprocess)).markerid)
    snporder == nothing && (snporder=1:nsnp)
    bvkeyls = getbvkeyls(bvpair,priorspace,phasedgeno)
    startprobls = [get(priorprocess,i,missing).startprob for i=bvkeyls]
    if eltype(dataprobinfo) <: Matrix
        dataprobls = getdataprobls(snporder[tinit],dataprobinfo,bvpair,priorspace,phasedgeno)
    else
        dataprobls = dataprobinfo
    end
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
    dataprobinfo::AbstractVector,bvpair::AbstractVector,
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
        tranprobdict = Dict([strkey=>getzygotetran(tnowdis,pri.nvalent)
            for (strkey, pri) in priorprocess])
    end
    if eltype(dataprobinfo) <: Matrix
        dataprobls = getdataprobls(snporder[tnext],dataprobinfo,bvpair,priorspace,phasedgeno)
    else
        dataprobls = dataprobinfo
    end
    popidls = phasedgeno.designinfo[!,:population]
    for popid = popidls
        offls = findall(phasedgeno.offspringinfo[!,:population] .== popid)
        keyls = priorspace[popid]["valentkey"]
        for off = offls
            fwoff = fwdict["fwprob"][off]
            dataprob = dataprobls[off]
            tranprob = tranprobdict[keyls[bvpair[off]]]
            pnow =fwoff[tnow,:]
            ls= (tranprob' * pnow) .* dataprob
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
    dataprobset::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    caldistance::Bool=false)
    snporder==nothing && (snporder = 1:size(dataprobset[1],1))
    # fwdict is modified inside of function
    dfdiff = -1
    isologl = singlelogl[tseq[kk],:]
    b=.!isnan.(isologl)
    if kk==1
        propfwdict = calinitforward(tseq[2],dataprobset,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
        proplogl0 = calindlogl(propfwdict,logbw,tseq[2])
        proplogl = proplogl0[b] .+ isologl[b]
    elseif kk==length(tseq)
        proplogl0 = [i[tseq[kk-1]] for i=fwdict["fwlogl"]]
        proplogl = proplogl0[b] .+ isologl[b]
    else
        if caldistance
            bvkeyls = getbvkeyls(bvpair,priorspace,polygeno)
            dataprobls = getdataprobls(snporder[tseq[kk+1]],dataprobset,bvpair,
                priorspace,polygeno)
            tnowdis= first(calinterdis(tseq[kk-1],tseq[kk+1],priorprocess,bvkeyls,
                dataprobls,fwdict,logbw))
        else
            pri1=first(values(priorprocess))
            tnowdis=sum(pri1.markerdeltd[tseq[kk-1:kk]])
        end
        calnextforward!(fwdict,tseq[kk-1],tseq[kk+1],dataprobset,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder,tnowdis=tnowdis)
        proplogl0 = calindlogl(fwdict,logbw,tseq[kk+1])
        proplogl = proplogl0[b] .+ isologl[b]
        calnextforward!(fwdict,tseq[kk-1],tseq[kk],dataprobset,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
    end
    logldiff = proplogl .- logl[b]
    isempty(logldiff) ? 0 : calvuongts(logldiff,dfdiff)
end

function polymarkerdel!(fhaploindex::AbstractVector,fhaploset::AbstractVector,
    dataprobset::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    delsiglevel::Real=0.05,
    caldistance::Bool=false)
    fhaplophase = hcat([map((x,y)->ismissing(y) ? x[1] .* missing : x[y],fhaploset[i],fhaploindex[i])
        for i=1:length(fhaploset)]...)
    snporder==nothing && (snporder = 1:size(dataprobset[1],1))
    singlelogl= calsinglelogl(fhaplophase,dataprobset,bvpair,priorspace,
        priorprocess,polygeno,snporder=snporder)
    logbw=callogbackward(dataprobset,bvpair,priorspace,priorprocess,polygeno,
        snporder=snporder)
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    fwdict = calinitforward(tseq[1],dataprobset,bvpair,priorspace,priorprocess,
        polygeno,snporder=snporder)
    logl = calindlogl(fwdict,logbw, tseq[1])
    res=Vector{Union{Missing,Float64}}(missing,length(markerincl))
    res[tseq] = [calvuongts!(kk,tseq,logl,singlelogl,logbw,fwdict,
        dataprobset,bvpair,priorspace,priorprocess,polygeno,
        snporder=snporder,caldistance=caldistance) for kk=1:length(tseq)]
    threshold = abs(quantile(Normal(),delsiglevel))
    delsnp = tseq[res[tseq] .> threshold]
    markerincl[delsnp] .= false
    for (key, val) in priorprocess
         setmarkerincl!(val,markerincl)
    end
    delsnp
end

function polymarkerdel!(dataprobset::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    delsiglevel::Real=0.05,
    caldistance::Bool=false,
    chrindex::Integer=1)
    fphase = polygeno.parentgeno[chrindex]
    fhaploset=[[[i] for i=j] for j=eachcol(fphase)]
    fhaploindex=[ones(Int, length(j)) for j=fhaploset]
    polymarkerdel!(fhaploindex,fhaploset,dataprobset,bvpair,
        priorspace, priorprocess,polygeno,
        snporder=snporder,delsiglevel=delsiglevel,caldistance=caldistance)
end
