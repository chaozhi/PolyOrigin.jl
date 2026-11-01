
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



# a gamate consists of n/2 bivalent chromosomes or one multivalent pairing
# e.g. n=6, nvalents = [2,2,2]  or [6,6,6]
function getgametetran(d::AbstractFloat,nvalents::AbstractVector,bindis::Real)
    # d ≈ 0.0 && return nothing    
    d < bindis/100 && return nothing    # d in unit of Morgan, bindis in unit of cM
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
    unique(reduce(vcat, [collect(permutations(i)) for i in res]))
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
    [ismissing(i) ? missing : Vector{Vector{Vector{Int16}}}(i) for i in res]    
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


function get_priorstatespace(polyped::PolyPed, chrpairing::Union{Nothing,Integer})
    (;ped, designinfo, parentinfo, offspringinfo) = polyped;
    maxvalent = isnothing(chrpairing) ? 4 : max(digits(chrpairing)...)
    ploidyls = parentinfo[!,:ploidy]
    maxvalentls=[min(maxvalent,ploidyls[i]) for i=eachindex(ploidyls)]
    OrderedDict([begin                        
        popid=designinfo[pop,:family]        
        offspringindex = findall(offspringinfo[!,:family] .== popid)
        offspringid = offspringinfo[offspringindex, :individual]
        midparentindex = findall(parentinfo[!,:family] .== popid)
        midparentid = parentinfo[midparentindex,:individual]
        if isempty(offspringid)
            ploidyls = unique(parentinfo[midparentindex, :ploidy])
            length(ploidyls) > 1 && @error string("differnet ploidies in family=", pop, ": ",ploidyls)
            offploidy = only(ploidyls)
        else
            ploidyls = unique(offspringinfo[offspringindex, :ploidy])
            length(ploidyls) > 1 && @error string("differnet ploidies in family=", pop, ": ",ploidyls)
            offploidy = only(ploidyls)
        end        
        parentidls = parentinfo[!,:individual]
        nprogenyls = parentinfo[!,:immediate_progeny]
        if isempty(offspringid)
            mem = findfirst(==(first(midparentid)), ped.member)
        else
            mem = findfirst(==(first(offspringid)), ped.member)
        end        
        moth = ped.mother[mem]
        fath = ped.father[mem]
        mothindex = findfirst(==(moth), parentidls) 
        fathindex = findfirst(==(fath), parentidls)
        parentid = [moth, fath]
        parentindex =  [mothindex, fathindex] 
        ploidypp=parentinfo[parentindex,:ploidy]
        maxvalentpp=maxvalentls[parentindex]
        isselfing = mothindex == fathindex
        if isnothing(chrpairing)
            minnprogeny = minimum(nprogenyls[parentindex]) 
            maxploidy = maximum(ploidypp)
            if maxploidy <= 2
                chrpairing2 = 22
            elseif maxploidy <= 4
                if minnprogeny > 2*maxploidy^2
                    chrpairing2 = 42
                else
                    chrpairing2 = 44
                end
            else
                chrpairing2 = 42
            end
        else
            chrpairing2 = chrpairing
        end
        zygotevalent=getzygotevalent(ploidypp,maxvalentpp,isselfing,chrpairing2)
        zygotestate =getzygotestate(zygotevalent)
        groupstate =getzygotegroupstate(zygotestate)
        condstate =getzygotecondstate(zygotevalent,zygotestate)
        valentkey =tovalentkey(zygotevalent)
        valentneighbor = getzygoteneighbor(ploidypp,maxvalentpp,isselfing,chrpairing2)        
        
        # offspring and midparent are full sibs with common parent
        popstate=Dict("parent"=>parentid,"parentindex"=>parentindex,
            "parentploidy" => ploidypp,
            "midparent"=>midparentid,"midparentindex"=>midparentindex,
            "offspring"=>offspringid,"offspringindex"=>offspringindex,
            "offspringploidy" => offploidy, 
            "chrpairing" => chrpairing2,
            "valent"=>zygotevalent,
            "valentkey"=>valentkey,"valentneighbor"=>valentneighbor,
            "condstate"=>condstate,"state"=>zygotestate,"groupstate"=>groupstate)
        popid => popstate
    end for pop = 1:size(designinfo,1)])
end


function get_popsum(priorspace::AbstractDict)
    df = DataFrame("familiy"=>String[], "parents"=>String[], "size"=>Int[],"chrpairing_phase"=>Int[])
    for (pop,val) in priorspace
        parents = join(val["parent"],"|")
        popsize = length(val["offspring"])
        push!(df, [pop,parents,popsize,val["chrpairing"]])
    end
    df
end

function get_priorprocess(polygeno::PolyGeno,chr::Integer,maxvalent::Integer;bindis::Real=1e-3)
    markerid = polygeno.markermap[chr][!,:marker]
    locseq = polygeno.markermap[chr][!,:poscm] ./ 100
    priordict = Dict{String,MarkovPrior}()
    for p = 1:size(polygeno.polyped.parentinfo,1)
        ploidypp=polygeno.polyped.parentinfo[p,:ploidy]
        maxvalentpp=min(maxvalent,ploidypp)
        gammetevalent = getgametevalent(ploidypp,maxvalentpp)
        keyls = unique([join(length.(i),"-") for i=gammetevalent])
        for strkey = keyls
            if  !haskey(priordict,strkey)
                val = valentkey2priorproces(strkey,markerid,locseq,bindis)
                push!(priordict,strkey => val)
            end
        end
    end
    priordict
end

function strkey2prior(priorprocess::AbstractDict, strkey::AbstractString)
    ls=[priorprocess[i] for i=split(strkey,"|")]
    startprob = [i.startprob for i=ls]
    tranprobseq = [i.tranprobseq[i.markerincl][1:end-1] for i=ls]    
    (markerincl=ls[1].markerincl,startprob=startprob, tranprobseq=tranprobseq)
end


function get_startprobls(priorprocess::AbstractDict,bvkeyls::AbstractVector)
    dict= Dict([begin
            priorls=[priorprocess[i] for i=split(strkey,"|")]
            startprob = [i.startprob for i=priorls]
            strkey => kron(startprob...)
        end for strkey = unique(bvkeyls)])
    [dict[i] for i=bvkeyls]
end

function get_startprobls(priorprocess::AbstractDict,bvdict::AbstractDict)
    bvkeyls = reduce(vcat, last.(values(bvdict)))
    get_startprobls(priorprocess, bvkeyls)
end

function get_tranprobls(priorprocess::AbstractDict,bvkeyls::AbstractVector, t::Integer; calkron::Bool=false)
    dict= Dict([begin
            priorls=[priorprocess[i] for i=split(strkey,"|")]
            tranprob = [i.tranprobseq[t] for i=priorls]
            if calkron
                tranprob2 = any(isnothing.(tranprob)) ? nothing : kron(tranprob...)
                strkey => tranprob2 
            else
                strkey => tranprob                 
            end            
        end for strkey = unique(bvkeyls)])
    [dict[i] for i=bvkeyls]
end

function get_tranprobls(priorprocess::AbstractDict,bvdict::AbstractDict,t; calkron::Bool=false)
    bvkeyls = reduce(vcat, last.(values(bvdict)))
    get_tranprobls(priorprocess, bvkeyls,t; calkron)
end

function cal_tranprobls(priorprocess::AbstractDict,bvdict::AbstractDict,tnowdis::Real,bindis::Real)
    bvkeyls = reduce(vcat, last.(values(bvdict)))
    cal_tranprobls(priorprocess, bvkeyls, tnowdis,bindis)
end


function cal_tranprobls(priorprocess::AbstractDict,bvkeyls::AbstractVector,tnowdis::Real,bindis::Real)
    gametedict = Dict([strkey=>getgametetran(tnowdis,pri.nvalent,bindis)
        for (strkey, pri) in priorprocess])        
    dict= Dict([begin        
        tranprob=[gametedict[i] for i=split(strkey,"|")]
        strkey => tranprob
    end for strkey = unique(bvkeyls)])    
    [dict[i] for i=bvkeyls]
end


mutable struct MarkovPrior
    # for a given linkage group
    startprob::Vector{Float64}
    tranprobseq::Vector{Union{Nothing,Matrix{Float64}}}
    nvalent::Vector{Int}
    markerdeltd::Vector{Union{Nothing,Float64}}
    markerincl::BitVector
    markerid::Vector{String}
    function MarkovPrior(startprob,tranprobseq,nvalent,markerdeltd,markerincl,markerid)
        dims = length.([markerdeltd,markerincl,markerid])
        if length(union(dims))!=1
            @error("inconsistent number of markers: ",dims)
        end
        new(startprob,tranprobseq,nvalent,markerdeltd,markerincl,markerid)
    end
end

function valentkey2priorproces(valentkey::AbstractString,
    markerid::AbstractVector,locseq::AbstractVector,bindis::Real)
    deltd = abs.(diff(locseq))
    nvalent = parse.(Int,split(valentkey,"-"))
    tranprobseq =vcat([getgametetran(i,nvalent,bindis) for i=deltd],nothing)
    # nstate=size(tranprob[1],1)
    nstate=foldr(*,vcat(nvalent...))
    startprob = [1.0/nstate for i=1:nstate]
    markerincl = trues(length(locseq))
    markerdeltd = vcat(deltd,[nothing])
    MarkovPrior(startprob,tranprobseq,nvalent,markerdeltd,markerincl,deepcopy(markerid))
end

function reverseprior!(pri::MarkovPrior)
    incl = pri.markerincl
    deltd = circshift(reverse(pri.markerdeltd[incl]),-1)
    tran = circshift(reverse(pri.tranprobseq[incl]),-1)
    reverse!(pri.markerincl) # incl is revised too
    pri.markerdeltd[incl] .= deltd
    pri.markerdeltd[.!incl] .= nothing
    pri.tranprobseq[incl] .= tran
    pri.tranprobseq[.!incl] .= nothing
    reverse!(pri.markerid)
    pri
end


function setpriorprocess!(priorprocess::AbstractDict, snporder::AbstractVector,
    delsnps::AbstractVector,bindis::Real)
    dict = Dict(snporder .=> 1:length(snporder))
    deltt = [get(dict, snp, nothing) for snp in delsnps]    
    setpriorprocess!(priorprocess,deltt,bindis)
end

function setpriorprocess!(priorprocess::AbstractDict, deltt::AbstractVector,bindis)
    for (_, pri) in priorprocess
         setmarkerincl!(pri,deltt,bindis)
    end
    priorprocess
end

function setmarkerincl!(pri::MarkovPrior,deltt::AbstractVector,bindis::Real)
    markerincl = pri.markerincl
    markerincl[deltt] .= false
    # calculate new deltd
    deltd = copy(pri.markerdeltd)
    i0 = findfirst(markerincl)
    if isnothing(i0)
        pri.markerdeltd[.!markerincl] .= nothing        
        pri.tranprobseq[.!markerincl] .= nothing
        @warn string("less than 2 markers are included! deltt=", deltt)
    else
        pre = i0
        for i=i0:length(markerincl)-1
            if markerincl[i]
                pre = i
            else
                if !isnothing(deltd[i])
                    deltd[pre] += deltd[i]
                    deltd[i] = nothing
                end
            end
        end
        markerincl[end] || (deltd[pre]=0.0)
        # println("pre=",pre,";i0=",i0,",deltt=",deltt)
        # update pri.markerdeltd
        pos = findall(deltd[1:end-1] .!= pri.markerdeltd[1:end-1])
        pri.markerdeltd[pos] .= deltd[pos]
        pos = pos[markerincl[pos]]
        if !isempty(pos)
            pri.tranprobseq[pos] .= [getgametetran(i,pri.nvalent,bindis) for i in deltd[pos]]
        end
    end
    pri
end

function setdistanceat!(priorprocess::AbstractDict,tnow::Integer,tnowdis::Real,bindis::Real)
    for (strkey, pri) in priorprocess
        pri.markerdeltd[tnow] = tnowdis
        pri.tranprobseq[tnow] = getgametetran(tnowdis,pri.nvalent,bindis)
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
