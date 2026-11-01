"""
    Pedigree{T}

immutable struct that stores information of a pedigree. See also [`readped`](@ref).

The Pedigree fields include `nfounder`, `member`, `mother`, `father`, `family`, `gender`, 
`generation`, and `ploidy`. The mothers and fathers of founders must be denoted by 0 or `notapplicable`.

    Pedigree(nfounder,member,mother,father,family,gender,generation,ploidy)

inner constructor. The gender of each member must be
"notapplicable", "female", or "male".

    Pedigree(df::DataFrame)

extenral constructor. The dataframe must have columns `member`,
`mother`, and `father`. The columns `family`, `gender`, and `ploidy` are optional. The constructor calculates
the generation for each member, which is the max path length from founders to
the member. The constructor orders pedigree by calling [`orderped`](@ref).

"""
struct Pedigree{T}
    nfounder::Int
    member::Vector{T}    
    mother::Vector{T}
    father::Vector{T}
    family::Vector{T}
    gender::Vector{T}
    generation::Vector{Int}
    ploidy::Vector{Int}
    function Pedigree(nfounder,member,mother,father,family, gender,generation,ploidy)
        if !allunique(member)
            error("pedigree members are not unique. ", member)
        end
        new{typeof(first(member))}(nfounder,member,mother,father,family, gender,generation,ploidy)
    end
end

function Pedigree(df::DataFrame)
    colls = propertynames(df)
    if issubset([:member, :mother,:father],colls)
        newcolls = intersect([:member, :mother,:father,:family, :gender, :ploidy], colls)
        # in(:gender, colls) && push!(newcolls, :gender)
        # in(:family, colls) && insert!(newcolls, 2, :family)
        df = df[:,newcolls]
    else
        @error string("pedigree dataframe must have columns member, mother, father")
    end
    n=size(df,1)
    colls = propertynames(df)
    if in(:gender, colls)
        gender=string.(df[!,:gender])
        unknown=setdiff(union(gender),["notapplicable","female","male"])
        if !isempty(unknown)
            @error("unknown genders = ", unknown, "; gender must be female, male, or notapplicable")
        end
    else
        gender = repeat(["notapplicable"],n)
    end
    if in(:family, colls)
        family = string.(df[!,:family])
        # TODO: check 
    else
        family= repeat(["missing"],n)
    end    
    if in(:ploidy, colls)
        ploidy = df[!,:ploidy]
        # TODO: check 
    else
        ploidy= repeat([2],n) # default is diploid
    end
    member,mother,father=[String.(strip.(string.(df[!,i]))) for i=[:member, :mother,:father]]
    allunique(member) || @error("pedigree members are not unique")
    set_founder_missing!(mother)
    set_founder_missing!(father)    
    isfounder= get_ismissing_parent(mother) .&& get_ismissing_parent(father)     
    isnonfounder = .!isfounder
    d=setdiff(mother[isnonfounder],member)
    isempty(d) || @error("nonfounders have unknown mothers  ", d)
    d=setdiff(father[isnonfounder],member)
    isempty(d) || @error("nonfounders have unknown father  ", d)
    nfounder=Int(sum(isfounder))
    generation=Int.(isnonfounder)
    ped=Pedigree(nfounder,member,mother,father,family, gender,generation, ploidy)    
    ped = orderped(ped) # generation is reset after ordering ped
    in("missing",family[isnonfounder]) && set_family!(ped)        
    ped    
end


function get_ismissing_parent(parent::AbstractVector{T} where T <: Union{AbstractString, Integer})
    missingset = eltype(parent) <: AbstractString ? ["notapplicable", "0","NA","missing"] : [0]
    [in(i, missingset) for i in parent]
end

function set_founder_missing!(parent::AbstractVector{T} where T <: AbstractString)
    for i in eachindex(parent)
        if in(parent[i], ["0","NA","missing"])
            parent[i] = "notapplicable"
        end
    end
    parent
end

function set_family!(ped::Pedigree; famprefix="Fam")
    family = ped.family    
    parentls = join.(tuple.(ped.mother, ped.father), "//")
    parentset = unique(parentls)
    for i in eachindex(parentset)
        pp = parentset[i]
        family[parentls .== pp] .= string(famprefix,i-1)
    end
    nf = ped.nfounder
    family[1:nf] .= string.("Founder",1:nf)     
    ped
end

function set_family!(peddf::AbstractDataFrame; famprefix="Fam")
    family = peddf[!,:family]   
    parentls = join.(tuple.(peddf[!,:mother], peddf[!,:father]), "//")
    parentset = unique(parentls)
    for i in eachindex(parentset)
        pp = parentset[i]
        family[parentls .== pp] .= string(famprefix,i-1)
    end
    # nf = sum(peddf[!,:mother] .== "notapplicable")     
    isfounder =  get_ismissing_parent(peddf[!,:mother]) .&& get_ismissing_parent(peddf[!,:father]) 
    nf = sum(isfounder)
    family[1:nf] .= string.("Founder",1:nf)     
    peddf
end


"""
    orderped(ped::Pedigree)

order pedigree members so that parents always come first.

"""
function orderped(ped::Pedigree)
    member, mother, father = ped.member, ped.mother, ped.father
    gen = [findall(ped.generation .== 0)]
    n = length(member)
    while true
        indexbef =vcat(gen...)
        indexaft=setdiff(1:n,indexbef)
        indexaft ==[] && break
        popbef=member[indexbef]
        b=[(mother[i] in popbef) && (father[i] in popbef)  for i= indexaft]
        push!(gen,indexaft[b])
    end
    ii=vcat(gen...)
    generation=vcat([(i-1)*ones(Int,length(gen[i])) for i=eachindex(gen)]...)    
    Pedigree(ped.nfounder,member[ii],mother[ii],father[ii],ped.family[ii], ped.gender[ii],generation,ped.ploidy[ii])
end

function extract_compactped(ped::Pedigree)
    nf = ped.nfounder
    founders = ped.member[1:nf]
    intermidparent = setdiff(union(ped.mother,ped.father), founders, ped.mother[[1]])
    peddf = ped2df(ped)
    offdf = peddf[nf+1:end,:]
    gdf = groupby(offdf,:family)
    rules = founders .=> founders
    for df in gdf
        midp =intersect(intermidparent, df[!,:member])
        fam = string(df[1,:family], "(", size(df, 1), ")")
        pushfirst!(midp,fam)
        # newnode = join(midp, "\n")
        newnode = fam
        append!(rules, df[!, :member] .=> newnode)
    end
    for col in [:member,:mother,:father]
        peddf[!, col] .= replace(peddf[!, col], rules...)
    end
    unique!(peddf)
    Pedigree(peddf)        
end


# function extract_compactped2(ped::Pedigree)
#     df = ped2df(ped)
#     set_family!(df)
#     nf = ped.nfounder
#     parentls = unique(Matrix(df[nf+1:end,[:mother, :father]]))
#     b = [in(i, parentls) for i in df[!,:member]]
#     ls = findall(.!b)
#     ls = ls[first.(splitindex(df[ls,:family]))]
#     df[ls,:member] .= df[ls,:family]
#     ls2 = findall(b)
#     union!(ls, ls2)
#     sort!(ls)
#     df2 = df[ls,:]
#     # rename intermediators 
#     rowls = findall(df2[!,:member] .!= df2[!,:family])
#     setdiff!(rowls, 1:nf)
#     newmemls = string.("M", rowls)
#     dict= Dict(df2[rowls,:member] .=> newmemls)
#     isempty(intersect(newmemls, df2[1:nf, :member])) || @error string("non-unique renaming")
#     for col in [:member, :mother, :father]
#         df2[nf+1:end, col] .= [haskey(dict, i) ? dict[i] : i for i in df2[nf+1:end,col]]
#     end
#     #
#     Pedigree(df2)        
# end


"""
    readped(pedfile::AbstractString, delim=',', commentstring="##")

read a pedigree from pedfile, ignoring the lines beginning with commentstring.

The pedigree file  must contain columns: member, mother, and father. The mothers
and fathers of founders must be denoted by 0. The value in the optional gender
column must  be "notapplicable", "female" or "male".

"""
function readped(pedfile::AbstractString; delim=',', commentstring="##")
    df=CSV.File(pedfile, delim=delim,comment=string(commentstring)) |> DataFrame
    rename_ped_cols!(df)
    Pedigree(df)
end


function rename_ped_cols!(peddf::AbstractDataFrame)
    colnamels = lowercase.(strip.(names(peddf)))    
    for i in eachindex(colnamels)
        s = colnamels[i]
        if occursin(r"^ind", s) || occursin(r"^id", s)
            colnamels[i] = "member"
        elseif occursin(r"^pop", s) || occursin(r"^fam", s)
            colnamels[i] = "family"
        elseif occursin(r"^moth", s)
            colnamels[i] = "mother"
        elseif occursin(r"^fath", s)
            colnamels[i] = "father"
        elseif occursin(r"^parent", s)
            if occursin(r"1$",s)
                colnamels[i] = "mother"
            elseif occursin(r"2$",s)
                colnamels[i] = "father"
            else
                @warn string("Could not parse pedigree column = ",colnamels[i])
            end
        else
            0
        end
    end    
    rename!(peddf, colnamels)
    peddf
end



"""
    saveped(filename::AbstractString,ped::Pedigree)

save pedigree into a csv file. The file contains columns: member,
mother, father, gender, generation.

"""
function saveped(filename::AbstractString,ped::Pedigree)
    df = ped2df(ped)    
    if in("ploidy",names(df))
        if all(df[!,:ploidy] .== 2)
            select!(df, Not([:ploidy]))
        end
    end        
    if in("gender",names(df))
        if all(df[!,:gender] .== "notapplicable")
            select!(df, Not([:gender]))
        end
    end        
    if in("family",names(df))
        family = string.(df[!, :family])
        b = [occursin(r"^Founder",i) || occursin(r"^Fam",i) for i in family] # internal naming rule for family
        if all(b)
            select!(df, Not([:family]))
        end
    end    
    CSV.write(filename,df)
end

function ped2df(ped::Pedigree)    
    # first field is nfounder, a scalar, not vector
    title=[fieldnames(typeof(ped))[2:end]...]
    df=DataFrame([getfield(ped,i) for i=title],title)
    df
end

"""
    setfounderid!(ped::Pedigree,founders::AbstractVector)

change founder IDs in a pedigree

"""
function setfounderid!(ped::Pedigree{String},founders::AbstractVector)
    # isfounder =  get_ismissing_parent(mother) .&& get_ismissing_parent(father) 
    isfounder = ped.generation .== 0
    if sum(isfounder) != ped.nfounder
        error("founders does not match nfounder")
    end
    if length(founders) != ped.nfounder
        error("the number of input founders must be  ",ped.nfounder)
    end
    if !allunique(founders)
        error("input founders are not unique: ",founders)
    end
    founders=string.(founders)
    ls = intersect(founders,ped.member[.!isfounder])
    if !isempty(ls)
        error("input founders overlap with non-founders: ",ls)
    end
    old= ped.member[isfounder]
    fdict = Dict(map((x,y)->x => convert(String,y),old, founders))
    ped.member .= [get(fdict,i,i) for i in ped.member]
    ped.mother .= [get(fdict,i,i) for i in ped.mother]
    ped.father .= [get(fdict,i,i) for i in ped.father]    
    ped
end

function getfounderid(ped::Pedigree)
    nf = ped.nfounder
    ped.member[1:nf]
end

"""
    toindexped(ped::Pedigree)

return a pedigree where members are coded as integers staring from 1.
The input pedigree must be ordered. See [`orderped`](@ref).

"""
function toindexped(ped::Pedigree)
    genderset = union(ped.gender)
    gendertype=promote_type(typeof.(genderset)...)
    if gendertype <: AbstractString
        unknown=setdiff(genderset,["notapplicable","female","male"])
        if isempty(unknown)
            sexrule=Dict(["notapplicable"=>0,"female"=>1,"male"=>2])
            genderno=[get(sexrule,i,-1) for i = ped.gender]
        else
            error("unknown genders = ", unknown, "; gender must be notapplicable, female, or male")
        end
    elseif gendertype <: Integer
        unknown=setdiff(genderset,[0,1,2])
        if isempty(unknown)
            genderno= ped.gender
        else
            error("unknown genders = ", unknown, "; gender must be  0=notapplicable, 1=female or 2=male")
        end
    else
        error("unexpected type of gender: ", gendertype)
    end
    member=ped.member
    family = ped.family
    mother=ped.mother
    father=ped.father
    dict=Dict(member .=> collect(1:size(member,1)))
    memberno= [get(dict,i,0) for i in member]
    motherno= [get(dict,i,0) for i in mother]
    fatherno= [get(dict,i,0) for i in father]
    famls = unique(family)
    dict = Dict(famls .=> eachindex(famls))
    familyno = [dict[i] for i in family]
    Pedigree(ped.nfounder,memberno,motherno,fatherno,familyno, genderno,ped.generation,ped.ploidy)
end


function getparents(peddict::Dict,memls::AbstractVector)
    res = copy(memls)
    while true
        memls = unique(vcat([get(peddict,i,nothing) for i=memls]...))
        memls = memls[.!isnothing.(memls)]
        isempty(memls) && break
        push!(res,memls...)
    end
    res
end

function getsubped(pedigree::Pedigree{T},member::T) where {T <: Union{String,Integer}}
    nlast = findlast(x->x==member,pedigree.member)
    isnothing(nlast) && error(string(member, " is not in the pedigree"))
    nf = pedigree.nfounder
    dict = Dict([pedigree.member[i] => [pedigree.mother[i],pedigree.father[i]] for i=nlast:-1:nf+1])
    parents = getparents(dict,[member])
    b = [in(i,parents) for i=pedigree.member]
    mother, father = pedigree.mother[b],pedigree.father[b]
    gender, generation = pedigree.gender[b], pedigree.generation[b]
    # cculcate new nf
    # isfounder =  get_ismissing_parent(mother) .&& get_ismissing_parent(father)         
    Pedigree(sum(generation .== 0),pedigree.member[b], mother,father,pedigree.family[b],gender,generation,pedigree.ploidy[b])
end

# ordering of memberls does not matter
function getsubped(pedigree::Pedigree{T},memberls::Vector{T}) where {T <: Union{String,Integer}}
    d = setdiff(memberls,pedigree.member)
    isempty(d) || @error string("members not in pedigree: ",d)
    nf = pedigree.nfounder
    nlast = length(pedigree.member)
    dict = Dict([pedigree.member[i] => [pedigree.mother[i],pedigree.father[i]] for i=nlast:-1:nf+1])
    parents = getparents(dict, memberls)
    b = [in(i,parents) for i in pedigree.member]
    mother, father = pedigree.mother[b],pedigree.father[b]
    gender, generation = pedigree.gender[b], pedigree.generation[b]
    # cculcate new nf
    # isfounder =  get_ismissing_parent(mother) .&& get_ismissing_parent(father)     
    Pedigree(sum(generation .== 0),pedigree.member[b],mother,father,pedigree.family[b],gender,generation,pedigree.ploidy[b])
end

function pedadjmtx(pedigree::Pedigree)
    iped=toindexped(pedigree)
    n = length(iped.member)
    adj = zeros(Int,n,n)
    nf=iped.nfounder
    for i=nf+1:n
        adj[iped.mother[i],i] = 1
        adj[iped.father[i],i] = 1
    end
    # for i=1:n
    #     adj[i,1:i-1] = adj[1:i-1,i]
    # end
    adj
end

# function splitvec(f::Function,A::AbstractVector)
#     res=Vector{typeof(A)}()
#     i0=1
#     for i=2:size(A,1)
#         if !f(A[i-1],A[i])
#             push!(res,A[i0:i-1])
#             i0=i
#         end
#     end
#     push!(res,A[i0:end])
#     res
# end

# function splitvec(A::AbstractVector)
#     f(x,y)= x == y
#     splitvec(f,A)
# end


"""
    plot(ped::Pedigree)

plot a pedigree. 

"""
function plotped(pedigree::Pedigree;
    curves = true,
    nodename = nothing,
    nodecolor = nothing,
    plotsize=nothing,
    markersize = nothing,
    fontsize = 10,
    nodesize = 0.1,
    xscale = 1, 
    yscale = 1, 
    edgecolor = :black,
    edgestyle = :solid,    
    plotargs... 
    )
    adj = pedadjmtx(pedigree)
    n=size(adj,1)
    if issubset(["female","male"], unique(pedigree.gender))
        nodeshape = [i=="male" ? :rect : :circle for i= pedigree.gender]
    else
        nodeshape = :circle
    end
    nf  = pedigree.nfounder
    if isnothing(nodecolor)    
        founder_colors = distinguishable_colors(nf, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
        nodecolor = vcat(founder_colors,[RGB(0.9,0.9,0.9) for i=nf+1:n])
    end
    nodename = isnothing(nodename) ? string.(pedigree.member) : (isempty(nodename) ? ["" for _ in eachindex(pedigree.member)] : string.(nodename))    
    xls, yls = get_ped_xy(pedigree)    
    xls .*= xscale
    yls .*= yscale
    xmin, xmax = extrema(xls)
    if isnothing(plotsize)
        plotsize = (max(120,120*(xmax-xmin)),max(120,100*first(yls)))
    end
    # @info "plotting pedigree with size = $(plotsize)..."

    if isnothing(markersize)
        lenls = length.(splitvec(pedigree.generation))
        maxlen = maximum(lenls)
        markersize = min(maxlen<=4 ? 2 : 5, max(length(lenls),maxlen)/3.0)
    end
    fig = graphplot(adj;
        x=xls, y=yls,
        nodeshape,
        nodecolor,        
        curves,
        markersize, 
        nodesize, 
        edgecolor, 
        edgestyle,
        size = plotsize, 
        plotargs...
    )
    labls = [(xls[i],yls[i],Plots.text(nodename[i],fontsize,:black)) for i=eachindex(xls, yls, nodename)]
    annotate!(fig,labls)
    fig
end

function get_ped_xy(ped::Pedigree; minxspace=nothing)
    gen = ped.generation
    gen2 = splitvec(gen)        
    gensize = reduce(vcat,[length(i)*ones(Int,length(i)) for i in gen2])      
    x_gen=reduce(vcat,[begin
        c = length(i)
        xx = (1:c) .- ((c+1)/2)        
    end for i=gen2])
    nf = ped.nfounder
    iped = toindexped(ped)    
    x = zeros(length(iped.member))    
    x[1:nf] .= (1:nf) .- ((nf+1)/2)
    for i in nf+1:length(iped.member)                
        x[i] = (x[iped.mother[i]] + x[iped.father[i]])/2
    end
    b = gensize .> 1
    @. x[b] = (x[b] + x_gen[b])/2    
    maxlen  = max(length.(gen2)...)
    xscale = 1.5*maxlen / length(gen2[1])
    x .*= xscale
    iils = splitindex(gensize)
    wls = [abs(-(extrema(x[ii])...)) for ii in iils]
    xscalels =maximum(wls) ./ wls
    for i in eachindex(iils,xscalels)
        length(iils[i]) <= 2 && continue
        m = mean(x[iils[i]])
        x[iils[i]] .*= xscalels[i]    
        x[iils[i]] .+= m - mean(x[iils[i]])
    end
    y = (last(gen)+1.0) .- gen    
    xmin,xmax = extrema(x)
    yscale = max(1.15, 1.5*(xmax-xmin)/first(y))
    y .*= yscale
    # enlarge founder populations
    b = gen .== 0
    x[b] .*= 1.1
    y[b] .+= 1.0

    minxspace = isnothing(minxspace) ? xscale : 1
    # set minxspace 
    iils = splitindex(gen)
    for ii in iils        
        subx = view(x, ii)
        oo = sortperm(subx)
        d = diff(subx[oo])
        d .= [i < minxspace ? minxspace : i for i in d]
        d .= accumulate(+, d)
        pushfirst!(d,0)
        start = subx[oo[1]] - (d[end] - (subx[oo[end]] - subx[oo[1]]))/2
        subx[oo] .= d .+ start        
    end
    x, y
end


function ped_nodename(pedigree::Pedigree)
    isstr_ped = eltype(pedigree.member) <: AbstractString
    if isstr_ped
        nf = pedigree.nfounder
        fls = pedigree.member[1:nf]
        offls = nf+1:length(pedigree.member)
        offls2 = [in(i,fls) ? string("O_",i) : i for i in offls]
        nodename = vcat(fls,offls2)
    else
        nodename = pedigree.member
    end
    nodename
end

function rename_nonfounder(ped::Pedigree;
    suffix::AbstractString=Random.randstring(),
    keepmembers::AbstractVector=[])
    b = ped.generation .> 0
    old=ped.member[b]
    if !isempty(keepmembers)   
        d =setdiff(keepmembers,ped.member)
        isempty(d) || @warn string(d," are not members of pedigree")
    end
    setdiff!(old,keepmembers)
    new = [string(i,"_",suffix) for i in old]
    fdict=Dict(old .=> new)
    member=[get(fdict,i,i) for i in ped.member]
    mother=[get(fdict,i,i) for i in ped.mother]
    father=[get(fdict,i,i) for i in ped.father]
    Pedigree(ped.nfounder,member,mother,father,ped.family, ped.gender,ped.generation,ped.ploidy)
end

function mergeped(pedigrees::Vector{Pedigree{String}}; 
    isfounderinbred::Bool=true,
    keeplast::Bool=true)
    pedls = [rename_nonfounder(ped; 
        suffix=string(i),
        keepmembers = (keeplast ? [ped.member[end]] : [])) for (i,ped) in enumerate(pedigrees)]
    peddf = unique(reduce(vcat,ped2df.(pedls)))
    sort!(peddf,:generation)
    if isfounderinbred
        f1df = peddf[peddf[!,:generation] .== 1,:]
        dict = Dict()
        for row in eachrow(f1df)            
            push!(dict,row[:member]=>join(row[[:mother,:father]],"/"))
        end
        peddf[!,[:member,:mother,:father]] .= [get(dict,i,i) for i in Matrix(peddf[!,[:member,:mother,:father]])]
    end
    unique!(peddf)
    # sort founders
    b = @. peddf[!,:mother] == "0" && peddf[!,:father]  == "0"
    peddf[b,:] .= sort(peddf[b,:], :member)
    Pedigree(peddf)
end

function renameped(ped)
    peddf = PolyOrigin.ped2df(ped)
    nf = ped.nfounder
    founderls = peddf[1:nf,:member]
    offls = peddf[(nf+1):end,:member]
    parentls = unique(Matrix(peddf[(nf+1):end,[:mother, :father]]))
    midpls = intersect(offls,parentls)
    leafls = setdiff(offls,midpls)
    dict = Dict(vcat(founderls .=> string.("P",eachindex(founderls)), midpls .=> string.("M",eachindex(midpls)), 
        leafls .=> string.("O",eachindex(leafls))))    
    peddf[!,:member] .= [dict[x] for x in peddf[:,:member]]
    for col in [:father,:mother]
        peddf[(nf+1):end,col] .= [dict[x] for x in peddf[(nf+1):end,col]]
    end
    PolyOrigin.Pedigree(peddf)
end

