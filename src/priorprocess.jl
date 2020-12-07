
function getcharpairslist(pairslist::AbstractVector,set::AbstractVector)
    vcat([getchrpair(pairs,set) for pairs=pairslist]...)
end

function getchrpair(pairs::AbstractVector,set::AbstractVector)
    set2=setdiff(set,pairs...)
    more=getchrpair(set2)
    vcat(pairs,more[1])
    [[pairs...,i] for i=more]
end

function getchrpair(set::AbstractVector)
    c1=min(set...)
    set2=filter(x->x!=c1,set)
    [[c1,i] for i=set2]
end


function getouterlist(v1::AbstractVector,v2::AbstractVector,x...)
    res=[vcat(i,j) for i=v1 for j=v2]
    foldl(getouterlist,x,init=res)
end

# a zygote consists of two gamtes for each of two parents
# eg. n=4, nvalents =[2,2] for parent1, nvalents = [4,4] for parent2
# then nvalentslist=[[2,2],[4,4]]
# it allows that a zygote consists of a single gamete, e.g., nvalentslist=[[2,2]]
function getzygotetran(d::AbstractFloat,nvalentslist::AbstractVector{Array{Int64,1}})
    ls=[getgametetran(d,nvalents) for nvalents=nvalentslist]
    foldl(kron,ls)
end

# a gamate consists of n/2 bivalent chromosomes or one multivalent pairing
# e.g. n=6, nvalents = [2,2,2]  or [6,6,6]
function getgametetran(d::AbstractFloat,nvalents::AbstractVector{Int64})
    ls=[getvalenttran(d,n) for n=nvalents]
    foldl(kron,ls)
end

# transition along a single chromome resulting from bi- or multi-valent pairing
# d is in unit of Morgan
function getvalenttran(d::AbstractFloat,nvalent::Integer)
    tran=Matrix{Float64}(I,nvalent,nvalent)
    ratio=(nvalent-1.0)/nvalent
    r= ratio*(1-exp(-d/ratio))
    tran *= (1-r*(1+1/(nvalent-1)))
    tran .+ r/(nvalent-1)
end

# calculate bivalent formation
function getgametebivalent(set::AbstractVector)
    len = length(set)
    foldl(getcharpairslist,[set for i=1:(len ÷ 2-1)],init=getchrpair([],set))
end

# calculate chromosome pairings, under constraint multivalent <= maxvalent
function getgametevalent(ploidy::Integer,maxvalent::Integer)
    set = collect(1:ploidy)
    res=getgametebivalent(set)
    valentset = collect(2:2:ploidy)
    if !(maxvalent in valentset)
        @error string("wrong maxvalent = ",maxvalent, ", maxvalent must be ",valentset)
    end
    if maxvalent>=4
        ls =collect(powerset(set, 4, 4))
        if ploidy==4
            # 4-4
            ls2 = [[i,i] for i=ls]
            res = vcat(res,ls2)
        elseif ploidy==6
            # 4-4-2
            ls2 = [[i,i,setdiff(set,i)] for i=ls]
            res = vcat(res,ls2)
        elseif ploidy==8
            # 4-4-2-2
            ls =collect(powerset(set, 4, 4))
            ls2 = [[vcat([i,i],j) for j=getgametebivalent(setdiff(set,i))] for i=ls]
            ls3 = vcat(ls2...)
            res = vcat(res,ls3)
            # 4-4-4-4
            ls =collect(powerset(set, 4, 4))
            temp2 = [[i,setdiff(set,i)] for i=ls]
            temp2=unique(sort.(temp2))
            temp3=[[i[1],i[1],i[2],i[2]] for i=temp2]
            res = vcat(res,temp3)
        end
    end
    if maxvalent>=6
        ls =collect(powerset(set, 6, 6))
        if ploidy==6
            # 6-6-6
            ls2 = [[i,i,i] for i=ls]
            res = vcat(res,ls2)
        elseif ploidy==8
            # 6-6-6-2
            ls2 = [[i,i,i,setdiff(set,i)] for i=ls]
            res = vcat(res,ls2)
        end
    end
    res
end

function getzygotevalent(ploidypp::AbstractVector, maxvalentpp::AbstractVector,isselfing::Bool,
    chrpairing::Integer)
    if length(ploidypp)==2
        n1,n2=ploidypp
        b1,b2=maxvalentpp
        ss1=getgametevalent(n1,b1)
        ss2=getgametevalent(n2,b2)
        if isselfing
            if n1==n2
                res=[[i,j] for i=ss1, j=ss2]
            else
                @error string("selfing parent must have unique ploidy.")
            end
        else
            res=[[i,map(x->x .+ n1,j)] for i=ss1, j=ss2]
        end
    elseif length(ploidypp)==1
        n1=first(ploidypp)
        b1=first(maxvalentpp)
        res = [[i] for i=getgametevalent(n1,b1),j=1:1]
    else
        error(string("length of ploidypp must be 1 or 2!, ploidypp =",ploidypp))
    end
    if chrpairing > 22
        res = [begin
            len = sort(map(x->max(length.(x)...),i),rev=true)
            len[1]*10+len[2] > chrpairing ? missing : i
        end for i=res]
    end
    res
end


function getgameteneighbor(gametevalents::AbstractVector)
    valents2=unique.(gametevalents)
    rule= Dict(valents2 .=> 1:length(valents2))
    [begin
        v=valents2[s]
        vswap=[[begin
            v2=deepcopy(v)
            v2[i][k], v2[j][l]=v2[j][l],v2[i][k]
            sort(sort.(v2))
        end for k=1:length(v[i]) for l=1:length(v[j])] for i=1:length(v) for j=i+1:length(v)]
        vswap2=unique(map(x->sort(sort.(x)),vcat(vswap...)))
        vswap3=map(x->x[sortperm(length.(x),rev=true)],vswap2)
        nb=[get(rule,i,missing) for i=vswap3]
        push!(nb,s)
        sort(nb)
    end for s=1:length(valents2)]
end

function getzygoteneighbor(ploidypp::AbstractVector, maxvalentpp::AbstractVector,isselfing::Bool,
    chrpairing::Integer)
    if length(ploidypp)==2
        n1,n2=ploidypp
        if isselfing && n1!=n2
            @error string("selfing parents must have unique ploidy.")
        end
        b1,b2=maxvalentpp
        ss1=getgametevalent(n1,b1)
        ss2=getgametevalent(n2,b2)
        nb1=getgameteneighbor(ss1)
        nb2=getgameteneighbor(ss2)
        linear=LinearIndices((1:length(ss1),1:length(ss2)))
        res=[[sort(linear[nb1[i],j]),sort(linear[i,nb2[j]])] for i=1:length(nb1), j=1:length(nb2)]
    elseif length(ploidypp)==1
        n1=first(ploidypp)
        b1=first(maxvalentpp)
        ss1=getgametevalent(n1,b1)
        res=getgameteneighbor(ss1)
    else
        error(string("length of ploidypp must be 1 or 2!, ploidypp =",ploidypp))
    end
    if chrpairing > 22 && length(ploidypp)==2
        valents=[[i,j] for i=ss1, j=ss2]
        bool=[begin
            len = sort(map(x->max(length.(x)...),i),rev=true)
            len[1]*10+len[2] > chrpairing
        end for i=valents]
        res2=Matrix{Union{Missing,eltype(res)}}(res)
        res2[bool] .= missing
        res2
    else
        res
    end
end

# function getgametestate(ploidy::Integer, maxvalent::Integer)
#     setls=[collect(1:ploidy) for i=1:ploidy ÷ 2]
#     res = getouterlist(setls...)
#     if maxvalent==2
#         bool=[length(unique(i))==length(i) for i=res]
#         res=res[bool]
#     elseif iseven(maxvalent)
#         a = [max([sum(i.==j ) for j=unique(i)]...) for i=res]
#         res=res[a .<= (maxvalent ÷ 2)]
#     else
#         @error string("wrong maxvalent = ",maxvalent, ", maxvalent must be 2, 4, 6, or 8.")
#     end
#     res
# end

# function getzygotestate(ploidypp::AbstractVector, maxvalentpp::AbstractVector,isselfing::Bool)
#     if length(ploidypp)==2
#         n1,n2=ploidypp
#         b1,b2=maxvalentpp
#         ss1=getgametestate(n1,b1)
#         ss2=getgametestate(n2,b2)
#         if isselfing
#             if n1==n2
#                 res=[vcat(i,j) for i=ss1 for j=ss2]
#             else
#                 @error string("selfing parent must have unique ploidy.")
#             end
#         else
#             res=[vcat(i,j .+ n1) for i=ss1 for j=ss2]
#         end
#     elseif length(ploidypp)==1
#         n1=first(ploidypp)
#         b1=first(maxvalentpp)
#         res = getgametestate(n1,b1)
#     else
#         error(string("length of ploidypp must be 1 or 2!, ploidypp =",ploidypp))
#     end
#     res
# end

# function getgametegroupstate(ploidy::Integer, maxvalent::Integer)
#     res = getgametestate(ploidy,maxvalent)
#     unique(sort.(res))
# end
#
# function getzygotegroupstate(ploidypp::AbstractVector, maxvalentpp::AbstractVector,isselfing::Bool)
#     if length(ploidypp)==2
#         n1,n2=ploidypp
#         b1,b2=maxvalentpp
#         ss1=getgametegroupstate(n1,b1)
#         ss2=getgametegroupstate(n2,b2)
#         if isselfing
#             if n1==n2
#                 res=[vcat(i,j) for i=ss1 for j=ss2]
#             else
#                 @error string("selfing parent must have unique ploidy.")
#             end
#         else
#             res=[vcat(i,j .+ n1) for i=ss1 for j=ss2]
#         end
#     elseif length(ploidypp)==1
#         n1=first(ploidypp)
#         b1=first(maxvalentpp)
#         res = getgametegroupstate(n1,b1)
#     else
#         error(string("length of ploidypp must be 1 or 2!, ploidypp =",ploidypp))
#     end
#     res
# end

function getzygotestate(zygotevalent::AbstractMatrix)
    states = [ismissing(valent) ? missing : getouterlist(vcat(valent...)...) for valent = zygotevalent]
    states2 = collect(skipmissing(reshape(states,:)))
    sort(unique(vcat(states2...)))
end

function getzygotegroupstate(zygotestate::AbstractVector)
    sort(unique(sort.(zygotestate)))
end

function getzygotecondstate(zygotevalent::AbstractMatrix,zygotestate::AbstractVector)
    rule = Dict(zygotestate .=> 1:size(zygotestate,1))
    [begin
        if ismissing(valent)
            missing
        else
            condstate=getouterlist(vcat(valent...)...)
            [get(rule,i,0) for i = condstate]
        end
    end for valent = zygotevalent]
end

function getpriorstatespace(polygeno::PolyGeno,chrpairing::Integer)
    getpriorstatespace(polygeno.designinfo,polygeno.parentinfo,chrpairing)
end

function getpriorstatespace(designinfo::DataFrame,parentinfo::DataFrame,
    chrpairing::Integer)
    maxvalent = max(digits(chrpairing)...)
    ploidyls = parentinfo[!,:ploidy]
    maxvalentls=[min(maxvalent,ploidyls[i]) for i=1:length(ploidyls)]
    Dict([begin
        gametecount=Vector(designinfo[pop,2:end])
        parentindex =  findall(gametecount .> 0)
        ploidypp=parentinfo[parentindex,:ploidy]
        parentid=parentinfo[parentindex,:individual]
        maxvalentpp=maxvalentls[gametecount .>= 1]
        isselfing = 2 in gametecount
        if isselfing
            ploidypp=repeat(ploidypp,2)
            maxvalentpp=repeat(maxvalentpp,2)
        end
        zygotevalent=getzygotevalent(ploidypp,maxvalentpp,isselfing,chrpairing)
        zygotestate =getzygotestate(zygotevalent)
        groupstate =getzygotegroupstate(zygotestate)
        condstate =getzygotecondstate(zygotevalent,zygotestate)
        valentkey =tovalentkey(zygotevalent)
        valentneighbor = getzygoteneighbor(ploidypp,maxvalentpp,isselfing,chrpairing)
        popstate=Dict("parent"=>parentid,"parentindex"=>parentindex,"valent"=>zygotevalent,
            "valentkey"=>valentkey,"valentneighbor"=>valentneighbor,
            "condstate"=>condstate,"state"=>zygotestate,"groupstate"=>groupstate)
        popid=designinfo[pop,:population]
        popid => popstate
    end for pop = 1: size(designinfo,1)])
end

# function getrfseq(polygeno::PolyGeno,chr::Integer)
#     deltd=diff(polygeno.markermap[chr][!,:position]) ./ 100.0
#     # Haldane map function, here rfreq defined one chromsome pairing
#     # so that the map function is the same as the diploid case
#     rfseq =@. 0.5 - 0.5 * exp(-2 * deltd)
#     rfseq
# end

function getpriorprocess(polygeno::PolyGeno,chr::Integer,chrpairing::Integer)
    d = sort(digits(chrpairing))
    chrpairing=sum([d[k]*10^(k-1) for k=1:length(d)])
    maxvalent = max(d...)
    ploidyls = polygeno.parentinfo[!,:ploidy]
    maxvalentls=[min(maxvalent,ploidyls[i]) for i=1:length(ploidyls)]
    markermap= polygeno.markermap
    locseq = markermap[chr][!,:position] ./ 100
    markerid = markermap[chr][!,:marker]
    priordict = Dict{String,MarkovPrior}()
    for pop = 1: size(polygeno.designinfo,1)
        gametecount=Vector(polygeno.designinfo[pop,2:end])
        ploidypp=polygeno.parentinfo[gametecount .> 0,:ploidy]
        maxvalentpp=maxvalentls[gametecount .>= 1]
        isselfing = 2 in gametecount
        if isselfing
            ploidypp=repeat(ploidypp,2)
            maxvalentpp=repeat(maxvalentpp,2)
        end
        zygotevalent=getzygotevalent(ploidypp,maxvalentpp,isselfing,chrpairing)
        keyls = tovalentkey(zygotevalent)
        indexls = [findfirst(x->(!ismissing(x)) && x==i,keyls) for i=skipmissing(unique(keyls))]
        for k = indexls
            strkey=keyls[k]
            if  !haskey(priordict,strkey)
                val = valentkey2priorproces(strkey,markerid,locseq)
                push!(priordict,strkey => val)
            end
        end
    end
    priordict
end

mutable struct MarkovPrior
    # for a given linkage group
    startprob::Vector{Float64}
    tranprobseq::Union{Missing,Vector{Union{Missing,Matrix{Float64}}}}
    nvalent::Vector{Vector{Int64}}
    markerdeltd::Vector{Union{Missing,Float64}}
    markerincl::BitVector
    markerid::Vector{String}
    function MarkovPrior(startprob,tranprobseq,nvalent,markerdeltd,markerincl,markerid)
        dims = length.([markerdeltd,markerincl,markerid])
        ismissing(tranprobseq) || push!(dims,size(tranprobseq,1))
        if length(union(dims))!=1
            @error("inconsistent number of markers: ",dims)
        end
        # if !ismissing(tranprobseq)
        #     s = length(startprob)
        #     dims = union(size.(skipmissing(tranprobseq)),[(s,s)])
        #     if length(union(dims))!=1
        #         @error("inconsistent number of states: ",dims)
        #     end
        # end
        new(startprob,tranprobseq,nvalent,markerdeltd,markerincl,markerid)
    end
end

function valentkey2priorproces(valentkey::AbstractString,
    markerid::AbstractVector,locseq::AbstractVector)
    # Haldane map function, here rfreq defined one chromsome pairing
    # so that the map function is the same as the diploid case
    deltd = abs.(diff(locseq))
    nvalent = map(x->parse.(Int,x),split.(split(valentkey,"|"),"-"))
    if max(length.(nvalent)...) >=3 && max(vcat(nvalent...)...)>=4
        tranprobseq = missing
    else
        # tranprobseq[i] =tranprob between marker i and i+1
        tranprobseq =vcat([getzygotetran(i,nvalent) for i=deltd],missing)
    end
    # nstate=size(tranprob[1],1)
    nstate=foldr(*,vcat(nvalent...))
    startprob = [1.0/nstate for i=1:nstate]
    markerincl = trues(length(locseq))
    markerdeltd = vcat(deltd,[0])
    MarkovPrior(startprob,tranprobseq,nvalent,markerdeltd,markerincl,deepcopy(markerid))
end


function reverseprior!(pri::MarkovPrior)
    bool=pri.markerincl
    tran = similar(pri.tranprobseq)
    tran .= missing
    tran[bool]=circshift(pri.tranprobseq[bool],1)
    reverse!(tran)
    pri.tranprobseq = tran
    pri.markerdeltd=circshift(reverse(pri.markerdeltd),-1)
    reverse!(pri.markerincl)
    reverse!(pri.markerid)
    pri
end


function setmarkerincl!(pri::MarkovPrior,markerincl::BitVector)
    if length(pri.markerincl) != length(markerincl)
        @error("markerincl has wrong #markers: ",length(markerincl))
    end
    deltd = pri.markerdeltd
    pre=0
    for i=1:length(markerincl)-1
        if markerincl[i]
            pre = i
        else
            pre>0 && (deltd[pre] += deltd[i])
            deltd[i] = 0
        end
    end
    markerincl[end] || (deltd[pre]=0)
    tran = similar(pri.tranprobseq)
    tran .= missing
    ii = findall(markerincl)[1:end-1]
    tran[ii]=[getzygotetran(i,pri.nvalent) for i=deltd[ii]]
    pri.tranprobseq = tran
    pri.markerincl = markerincl
    pri.markerdeltd = deltd
    pri
end

function setsegrev!(priorprocess::AbstractDict,tt::AbstractVector)
    tt2=tt[1:end-1]
    for (strkey, pri) in priorprocess
        pri.markerdeltd[tt2] = pri.markerdeltd[reverse(tt2)]
        pri.tranprobseq[tt2] = pri.tranprobseq[reverse(tt2)]
        pri.markerid[tt] = pri.markerid[reverse(tt)]
    end
    priorprocess
end

function setdistanceat!(priorprocess::AbstractDict,tnow::Integer,tnowdis::Real)
    for (strkey, pri) in priorprocess
        pri.markerdeltd[tnow] = tnowdis
        pri.tranprobseq[tnow] = getzygotetran(tnowdis,pri.nvalent)
    end
    priorprocess
end

function tovalentkey(valents::AbstractVecOrMat)
    keyls0 =[ismissing(v) ? missing : map(x->length.(x),v) for v=valents]
    keyls =[ismissing(v) ? missing : stringjoin(stringjoin.(v,"-"),"|") for v=keyls0]
end

function getvalentdict(priorspace::AbstractDict)
    Dict([begin
        valents = priorspace[popid]["valent"]
        keyls =[string(vcat(map(x->length.(x),v)...)...) for v=valents]
        keyset = unique(keyls)
        keybool = [keyls .== i for i = keyset]
        popid => Dict(keyset .=> keybool)
    end for popid= keys(priorspace)])
end
