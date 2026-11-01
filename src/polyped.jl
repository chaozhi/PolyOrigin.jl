"""
    PolyPed

struct that stores pedigree information

    PolyPed(ped,designinfo,parentinfo,offspringinfo) 
    PolyPed(ped)

constructor from a pedigree. 

# Fields

`parentinfo::DataFrame`: parent information with columns [:individual, :ploidy].

`offspringinfo::DataFrame`: offspring information with columns
[:individual,:family,:isoutlier].

`designinfo::DataFrame`: design information with columns
[:family, :parent1, :parent2, ...]. The dataframe element `n_ij` denotes
the nubmer of gametes contributed to each offspring in i-th population by j-th parent.
`n_ij`=2 means that i-th population is produced by self-fertilization of
j-th parent. Row sum must be 2.

"""
struct PolyPed
    ped::Pedigree    
    designinfo::DataFrame    
    parentinfo::DataFrame    
    offspringinfo::DataFrame
    function PolyPed(ped, designinfo,parentinfo,offspringinfo)
        # TODO: consistency test
        new(ped, designinfo,parentinfo,offspringinfo)
    end
end

function PolyPed(ped::Pedigree)
    parentinfo = get_parentinfo(ped)
    offspringinfo = get_offspringinfo(ped, parentinfo)
    designinfo = get_designinfo(ped, parentinfo)
    PolyPed(ped, designinfo,parentinfo,offspringinfo)
end

function readpolyped(pedfile::AbstractString; 
    delim::AbstractChar=',',
    commentstring::AbstractString="##", 
    workdir::AbstractString=pwd())
    pedfile2 = getabsfile(workdir, pedfile)
    isfile(pedfile2) || @error string(pedfile2, " does not exit")    
    ped = readped(pedfile2; delim, commentstring)
    PolyPed(ped)
end


function savepolyped(sink::Union{IO,AbstractString},polyped::PolyPed;
    delim::AbstractChar=',',
    workdir::AbstractString = pwd())
    if typeof(sink) <: AbstractString
        outputfile2 = getabsfile(workdir,sink)
        io = open(outputfile2, "w")
    elseif typeof(sink) <: IO
        io = sink
    else
        error(string("unknown sink: ",sink))
    end    
    ped=polyped.ped
    df = ped2df(ped)    
    CSV.write(io,df; header=true,append=true,delim)    
    flush(io)
    typeof(sink) <: AbstractString && close(io)
    sink
end

function get_designinfo(ped::Pedigree,parentinfo::DataFrame)
    nf = ped.nfounder
    famls = ped.family[nf+1:end] 
    famset = unique(famls)
    parentls = parentinfo[!,:individual]
    parentdict = Dict(parentls .=> eachindex(parentls))
    res = zeros(Int,length(famset),length(parentls))
    for f in eachindex(famset)
        fam = famset[f]
        ii = findfirst(==(fam), famls) + nf
        moth = parentdict[ped.mother[ii]]
        fath = parentdict[ped.father[ii]]
        res[f, moth] += 1
        res[f, fath] += 1
    end
    df=DataFrame(res,parentls)
    insertcols!(df, 1, :family=> convert.(String,famset))
    df
end

function get_offspringinfo(ped::Pedigree,parentinfo::DataFrame)    
    parentls = parentinfo[!,:individual]
    offls = setdiff(ped.member, parentls)
    dict = Dict(ped.member .=> eachindex(ped.member))
    ii = [dict[i] for i in offls]
    DataFrame("individual"=>ped.member[ii],
        "family" => ped.family[ii], "generation"=> ped.generation[ii], 
        "ploidy"=>ped.ploidy[ii], "isoutlier" => missing)
end

function get_parentinfo(ped::Pedigree)
    nf = ped.nfounder
    parentls = unique(vcat(ped.mother[nf+1:end],ped.father[nf+1:end]))
    union!(parentls, ped.member[1:nf])
    parentls = intersect(ped.member,parentls) # to maintain the eorder of parents  
    dict = Dict(ped.member .=> eachindex(ped.member))
    ii = [dict[i] for i in parentls]
    isfounder = ped.generation[ii] .== 0
    ls = tuple.(ped.mother,ped.father)[nf+1:end]    
    nprogenyls = [sum(in(p, j) for j in ls; init=0) for p in parentls]
    length(ls) * 2 - sum(allequal.(ls)) == sum(nprogenyls) || @error "unexpected nprogenyls"
    DataFrame("individual" => ped.member[ii], "isfounder" => isfounder,  
        "family" => ped.family[ii], "generation" => ped.generation[ii],
        "ploidy" => ped.ploidy[ii], "immediate_progeny" => nprogenyls
    )
end



function get_parent_2_immediateoff(polyped::PolyPed)
    (;ped, parentinfo) = polyped
    parentls  = parentinfo[!,:individual]
    nf = ped.nfounder
    ls = tuple.(ped.mother,ped.father)[nf+1:end]
    offlsls = [ped.member[findall([in(p, j) for j in ls] ) .+ nf] for p in parentls]
    Dict(parentls .=> offlsls)
end


function filter_polyped(polyped::PolyPed; minprogeny::Integer=12,keepparent::Bool=false)    
    dells = String[]    
    parent2off = get_parent_2_immediateoff(polyped)    
    parentls = collect(keys(parent2off))
    offlsls = collect(values(parent2off))        
    while true
        ls = findall(length.(offlsls) .< minprogeny)
        if keepparent
            b = [isempty(intersect(parentls,offlsls[i])) for i in ls]
            keepat!(ls, b)
        else
            if !isempty(ls)
                while true
                    midls = unique(reduce(vcat,[intersect(parentls,offlsls[i]) for i in ls]))
                    setdiff!(midls, parentls[ls])
                    isempty(midls) && break
                    append!(ls, [findfirst(==(i), parentls) for i in midls])
                    unique!(ls)
                end
                sort!(ls)
            end
        end        
        if isempty(ls)
            break
        else
            append!(dells,offlsls[ls]...)
            append!(dells,parentls[ls])
            unique!(dells)
            deleteat!(offlsls,ls)        
            deleteat!(parentls,ls)        
            for i in offlsls
                setdiff!(i,dells)
            end                        
        end                
    end
    keepls = setdiff(polyped.ped.member,dells)
    subped = PolyOrigin.getsubped(polyped.ped,keepls)
    PolyPed(subped)
end


# function readDesign(pedfile::AbstractString;
#     delimchar::AbstractChar=',',
#     commentstring::AbstractString="##")
#     # not considered: the case of offspring with only one parent, that is, gamete instead of zygote.
#     if !isfile(pedfile)
#         error("pedfile: ", pedfile, " does not exist")
#     end
#     design=CSV.read(pedfile,DataFrame; delim=delimchar,comment=commentstring)
#     if size(design,2)<5
#         error("At least 5 columns are required in ",pedfile)
#     else
#         design=design[:,1:5]
#         rename!(design,[:individual,:family,:mother,:father,:ploidy])
#         for i=1:4
#             design[!,i]=string.(strip.(string.(design[!,i])))
#         end
#     end
#     allunique(design[!,1]) || @error("individual IDs are not unique")
#     isparent = [design[i,3] == "0" && design[i,4] == "0" for i=1:size(design,1)]
#     isparent .+= [design[i,3] == "NA" && design[i,4] == "NA" for i=1:size(design,1)]
#     isoff  = .!isparent
#     d=setdiff(design[isoff,:mother],design[!,:individual])
#     isempty(d) || @error("offspring have unknown mothers  ", d)
#     d=setdiff(design[isoff,:father],design[!,:individual])
#     isempty(d) || @error("offspring have unknown father  ", d)
#     poplist=design[isoff,:family]
#     popindexlist=[findall(poplist .== i) for i=unique(poplist)]
#     parentinfo=design[isparent,[:individual,:ploidy]]
#     offspringinfo=design[isoff,[:individual,:family,:ploidy]]
#     outlier = Vector{Union{Missing,Bool}}(missing,size(offspringinfo,1))
#     offspringinfo[!,:isoutlier] = outlier
#     dict=Dict(parentinfo[!,:individual] .=> parentinfo[!,:ploidy])
#     for rg = popindexlist
#         subpop = unique(design[isoff, :][rg,2:5])
#         size(subpop,1) == 1 || @error string("parents and ploidy levels in a subpopulation are inconsistent: ", subpop)
#         parentploidy = [get(dict, i, missing) for i=Vector(subpop[1,2:3])]
#         offploidy = sum(parentploidy) ÷ length(parentploidy)
#         offploidy == subpop[1,end] || @error string("offspring ploidy = ",offploidy, " is consisstent with those of parents = ",parentploidy)
#     end
#     popparent=design[isoff,[:family,:mother,:father]][first.(popindexlist),:]
#     designinfo=getdesigninfor(popparent,parentinfo)
#     (parentinfo=parentinfo,offspringinfo=offspringinfo,designinfo=designinfo)
# end

# function getdesigninfor(popparent::DataFrame,parentinfo::DataFrame;
#     missingstring::String="NA")
#     npop=size(popparent,1)
#     nfounder = size(parentinfo,1)
#     rule=Dict(parentinfo[!,:individual] .=> 1:nfounder)
#     crossparent=setdiff(Matrix(popparent[!,2:3])[:],[missingstring])
#     d=setdiff(crossparent,parentinfo[!,:individual])
#     if !isempty(d)
#         error(string("parents that are used in crosses but not listed:", d, ",pp=",parentinfo[!,:individual]))
#     end
#     d=setdiff(parentinfo[!,:individual],crossparent)
#     if !isempty(d)
#         error(string("parents that are listed in ",pedfile,"  but not used in crosses:", d))
#     end
#     res = zeros(Int,npop,nfounder)
#     for i=1:npop,j =2:3
#         c=get(rule,popparent[i,j],0)
#         c>0 && (res[i,c]+=1)
#     end
#     df=DataFrame(res,Symbol.(parentinfo[!,:individual]))
#     insertcols!(df, 1, :family=> convert.(String,popparent[!,:family]))
#     df
# end
