function readparenthaplo(refhapfile::AbstractString, parentinfo::DataFrame; 
    workdir::AbstractString=pwd(),
    missingstring::AbstractString="NA",    
    commentstring::AbstractString="##")
    refhapfile2=getabsfile(workdir,refhapfile)
    isfile(refhapfile2) || @error(string(refhapfile," does not exist in workdir = ",workdir))    
    ext = last(split_allext(refhapfile))
    delim = occursin(r"^.vcf",ext) ? '\t' : ','
    refhap=CSV.read(refhapfile2,DataFrame; delim, comment=commentstring,missingstring=missingstring)
    parentls = parentinfo[!,:individual]    
    rename!(refhap, Dict([i => strip(i,[' ', '"']) for i in names(refhap)]))    
    pcolls = [findfirst(==(i),names(refhap)) for i in parentls]
    b = isnothing.(pcolls)
    if any(b)
        msg = string("parents in parentinfo but not in refhapfile: ", parentls[b])
        error(msg)
    end    
    if occursin(r"^.csv",ext)        
        _col_1stp = minimum(pcolls)
        if _col_1stp > 5
            mapdf = refhap[!,1:5]
        elseif _col_1stp == 4
            mapdf = refhap[!,1:3]            
            insertcols!(mapdf, 4, "physchrom"=> missing)
            insertcols!(mapdf, 5, "physposbp"=> missing)
        else 
            error(string("pcolls=",pcolls))
        end
        # parentdf
        parentdf = refhap[!, pcolls]
    elseif occursin(r"^.vcf",ext)
        mapdf = refhap[!,[3,1,2]]
        # [newinfo, linkagegroup, poscm, doseerror, baseerror, allelicbias]
        infols = parse_vcf_info.(refhap[!,:INFO])
        chrls = [i[2] == "NA" ? missing : i[2] for i in infols]
        posls = [i[3] == "NA" ? missing : parse(Float64,i[3]) for i in infols]
        insertcols!(mapdf, 2, "linkagegroup"=> chrls)
        insertcols!(mapdf, 3, "poscm"=> posls)
        # parentdf
        formatpriority = ["GT"]
        missingset = [".","./.", ".|."]    
        formatls = refhap[!, :FORMAT]
        parentdf = refhap[!, pcolls]
        parentgeno = Matrix{Any}(parentdf)
        for i in eachindex(formatls)
            format = split(formatls[i],":")
            formatcode = [in(i, formatpriority) ? 1 : -1 for i in format] # -1 denotes missing fromat        
            for j in 1:length(pcolls)
                parentgeno[i,j] = first(parse_vcfcell(parentgeno[i,j], formatcode, formatpriority, missingset, true) )
            end
        end
        # TODO: allow for missing in parentgeno
        parentdf .= parentgeno
    else
        error(string("refhaplfile with unknown file extension: ",ext))
    end
    res = hcat(mapdf, parentdf)
    if occursin(r"^.vcf",ext)        
        dict = Dict(names(res)[1:5] .=> ["marker", "linkagegroup", "poscm","physchrom","physposbp"])
        rename!(res, dict)
        res[!,2] .= missing
        res[!,3] .= missing
    end
    for col in [1,2,4]
        res[!,col] .= [isa(i, AbstractString) ? String(strip(i, [' ','"'])) : i for i in res[!,col]]
    end
    res    
end


function readparenthaplo(refhapfile::AbstractString, parentinfo::DataFrame,markermap::Vector{DataFrame};
    workdir::AbstractString=pwd(),
    missingstring::AbstractString="NA",    
    commentstring::AbstractString="##")    
    refhap = readparenthaplo(refhapfile, parentinfo; workdir, missingstring, commentstring)
    ploidyls = parentinfo[!,:ploidy]
    # added missing markers to refhap
    vcatmap = vcat(markermap...)
    dsnp = setdiff(vcatmap[!,1],refhap[!,1])    
    dsnpdf = copy(refhap[ones(Int,length(dsnp)),:])
    dsnpdf[!,2] .= "NA"
    dsnpdf[!,3] .= 0.0
    dsnpdf[!,5] .= 0
    markerls = vcatmap[!,1]
    for i in eachindex(dsnp)
        i2 = findfirst(==(dsnp[i]), markerls)
        dsnpdf[i,1:5] .= Vector(vcatmap[i2,1:5])
    end        
    # "NA" denotes missing used in parseinputgeno
    for j = 6:size(dsnpdf,2)       
        dsnpdf[!,j] .= String.(dsnpdf[!,j]) 
        dsnpdf[!,j] .= join([-1 for _ in 1:ploidyls[j-5]],"|")
    end
    refhap = vcat(dsnpdf, refhap)
    # # align markers  of refhap with input markermap    
    snpdict =Dict(refhap[!,1] .=> 1:size(refhap,1))
    snpindex = [get(snpdict,i, nothing) for i = vcatmap[!,1]]
    b = isnothing.(snpindex)    
    if any(b) 
        @error string("unexpected missing markers!")
        deleteat!(snpindex,b)
    end
    refhap = refhap[snpindex,:]
    # cal fhaplo2
    if all(.!ismissing.(refhap[!,2])) 
        groupcol = 2
    elseif all(.!ismissing.(refhap[!,5])) 
        groupcol = 4
    else
        error(string("Could not find groupign column without missing value"))
    end
    groupref = groupby(refhap,groupcol)
    fhaplo = [Matrix(groupref[ch][!,6:end]) for ch=1:length(groupref)]
    fhaplo2=parseinputgeno(fhaplo; missingstring="NA")
    # check if it is haplotypes
    hhcount = [unique(eachrow(length.(fhaplo2[i]))) for i=1:length(fhaplo2)]
    all(length.(hhcount) .==1) || @error("number of haplotypes inconsistent accross markers")
    # all([i[1] ==ploidy for i=hhcount]) || @error("number of haplotypes mismatch ploidy")
    parenttype1 = [eltype(fhaplo2[i]) for i=1:length(fhaplo2)]
    parenttype2 = unique(parenttype1)
    if length(parenttype2) != 1
        @error(string("inconsistent data types of parental genodata: ",parenttype1))
    end
    if !(eltype(parenttype2[1]) <: Union{Missing,Integer})
        @error string("unknown data types of parental genodata: ",parenttype2[1])
    end
    fhaplo2
end

function calnphaseerr(refhapch::AbstractVector,esthapch::AbstractVector)
    map((x,y)-> sum(skipmissing(sum((x .== y) .=== false,dims=2) .> 0)),refhapch,esthapch)
end

function calndoseerr(refhapch::AbstractVector,esthapch::AbstractVector)
    map((x,y)->sum(skipmissing(abs.(sign.((sum(x.-1,dims=2)-sum(y.-1,dims=2)))))),refhapch,esthapch)
end

function calpermhomolog(refhapch::AbstractVector,esthapch::AbstractVector)
    ploidy=size.(refhapch,2)
    nparent = length(ploidy)
    permesthapch = Vector{Matrix}(undef,nparent)
    permhomolog = Vector{Vector}(undef,nparent)
    for p=1:nparent
        perm=[argmin([sum(abs.(skipmissing(i .- col)))  for col =eachcol(esthapch[p])])
            for i=eachcol(refhapch[p])]
        if length(unique(perm)) < ploidy[p]
            perms = collect(permutations(1:ploidy[p]))
            dis = [sum(abs.(skipmissing(esthapch[p][:,perm] .- refhapch[p]))) for perm=perms]
            perm = perms[argmin(dis)]
        end
        permhomolog[p] = perm
        permesthapch[p] = esthapch[p][:,perm]
    end
    nphaseerr = calnphaseerr(refhapch,permesthapch)
    permhomolog, permesthapch, nphaseerr
end

function calabsolutehap(estparentgeno::AbstractVector,refparenthap::AbstractVector, polyped::PolyPed)
    refhap = [[Matrix(reduce(hcat,hh[:,i])') for i=1:size(hh,2)] for hh = refparenthap]
    hhls = getchrparentgeno.(estparentgeno)
    esthap = [[Matrix(reduce(hcat,hh[:,i])') for i=1:size(hh,2)] for hh = hhls]    
    ccparent = connected_parents(polyped)
    nch = length(esthap)
    nparent = size(polyped.parentinfo[!,:ploidy],1)
    rulehomologls = Vector(undef,nch)
    newesthap = Vector(undef,nch)
    for ch=1:nch
        # rulehomolog[ref_i,:] = [est_parentindex, est_permutation of homologs]
        # rulehomolog[2,:] = [1, [4,3,2,1]]: # refhaploch[2] matches esthaploch[1][4,3,2,1]
        rulehomolog = Matrix(undef,nparent,2)
        newesthapch = similar(esthap[ch])
        for pp = ccparent
            permhomolog, permesthapch,nphaseerr=calpermhomolog(refhap[ch][pp],esthap[ch][pp])
            rulehomolog[pp,:] = [pp permhomolog]
            newesthapch[pp] = permesthapch
            if length(pp)==2            
                # Given only unphased genotypic data of offspring,
                # ordering of the two parents in a simple F1 cross is indistinguishable
                permhomolog2, permesthapch2,nphaseerr2=calpermhomolog(refhap[ch][pp],esthap[ch][reverse(pp)])
                if sum(nphaseerr2) < sum(nphaseerr)
                    rulehomolog[pp,:] = [reverse(pp) permhomolog2]
                    newesthapch[pp] = permesthapch2
                end
            end
        end
        newesthap[ch] = newesthapch
        rulehomologls[ch] = rulehomolog
    end    
    # inverse map of parentindex
    # homologmapls[chr][est_i] = [ref_parentindex, est_permutation of homologs]
    # homologmapls[chr][2] = [1, [4,3,2,1]]: # esthaploch[2][4,3,2,1] matches refhaploch[1]
    homologmapls = [begin
        homologmap = Matrix(undef,nparent,2)
        for esti=1:nparent
            refi=findfirst(rulehomolog[:,1] .== esti)
            homologmap[esti,:]=[refi,rulehomolog[refi,2]]
        end
        homologmap
    end for rulehomolog = rulehomologls]
    # rulehomologls[chr]: matrix with size #parents x two columns(parentindex, permutation)
    # newesthap[chr][p]: phased genotypes for parent p, matrix of size #markers x ploidy 
    # newestparentgeno[chr]: phase genotyep matrix of size #markers x #parents
    newestparentgeno = [reduce(hcat,eachrow.(hh)) for hh=newesthap]    
    homologmapls,newestparentgeno
end

function gethomologdict(homologmap::AbstractMatrix,popstatespace::AbstractDict)
    parentindex = popstatespace["parentindex"]
    # groupstate = popstatespace["groupstate"]
    if length(parentindex)==1
        ls = homologmap[parentindex[1],2]
    elseif length(parentindex)==2
        if parentindex[1] < parentindex[2]
            hmin = homologmap[parentindex[1],2]
            ls =vcat(hmin,homologmap[parentindex[2],2] .+ length(hmin))
        else
            hmin = homologmap[parentindex[2],2]
            ls =vcat(hmin,homologmap[parentindex[1],2] .+ length(hmin))
        end
    else
        @error("wrong parentindex = ",parentindex)
    end
    Dict(ls .=> 1:length(ls))
end

function getstateorder(homologmap::AbstractMatrix,popstatespace::AbstractDict)
    homologdict = gethomologdict(homologmap,popstatespace)
    diplostate = popstatespace["state"]
    diplodstate2 = [[get(homologdict,i,missing) for i=j] for j=diplostate]
    rule = Dict(diplodstate2 .=> 1:length(diplodstate2))
    [get(rule,i,missing) for i=diplostate]
end

function ordercondprob!(polyancestry::PolyAncestry,homologmapls::AbstractVector)
    # homologmapls[chr]: matrix size nparent x 2 where 1st col is parent index, and 2nd col is the estimated permutation of homologs]
    # homologmapls[chr][p,:] = [p, perm], [i[perm] for i in polyest.parentgeno[chr][:,p]] == refparenthap[chr][:,p]
    nch = length(polyancestry.markermap)
    for pop = 1:size(polyancestry.polyped.designinfo,1)
        popid = polyancestry.polyped.designinfo[pop,:family]
        popstatespace = polyancestry.statespace[popid]        
        offls = findall(polyancestry.polyped.offspringinfo[!,:family] .== popid)
        for ch=1:nch
            neworder = getstateorder(homologmapls[ch],popstatespace)
            for off=offls
                polyancestry.diploprob[ch][off] = polyancestry.diploprob[ch][off][:,neworder]
            end
        end
    end
    setgenoprob!(polyancestry)
    sethaploprob!(polyancestry)
end

function ordercorrecion!(correctdf::DataFrame,homologmapls::AbstractVector,
    chridls::AbstractVector,ppidls::AbstractVector)
    chridrule = Dict(chridls .=> 1:length(chridls))
    ppidrule = Dict(ppidls .=> 1:length(ppidls))
    for i=1:size(correctdf,1)
        chr = get(chridrule,correctdf[i,:linkagegroup],missing)
        pp = get(ppidrule,correctdf[i,:parent],missing)        
        pp2, perm = homologmapls[chr][pp,:]
        correctdf[i,:parent] = ppidls[pp2]
        correctdf[i,:old_genotype] = join(split(correctdf[i,:old_genotype],"|")[perm],"|")
        correctdf[i,:new_genotype] = join(split(correctdf[i,:new_genotype],"|")[perm],"|")
    end
end

function getvalentorder(homologmap::AbstractMatrix,popstatespace::AbstractDict)
    homologdict = gethomologdict(homologmap,popstatespace)
    valentstate = popstatespace["valent"]
    states = [sort.(i) for i=valentstate[:]]
    staterule = Dict(states.=> 1:length(states))
    [begin
        vv2=[sort(map(x->sort([get(homologdict,i,missing) for i=x]),v)) for v=vv]
        get(staterule,vv2,missing)
    end for vv=states]
end

function ordervalentprob!(polyancestry::PolyAncestry,homologmapls::AbstractVector)
    nchr = length(polyancestry.markermap)
    for pop = 1:size(polyancestry.polyped.designinfo,1)
        popid = polyancestry.polyped.designinfo[pop,:family]
        popstatespace = polyancestry.statespace[popid]
        offls = findall(polyancestry.polyped.offspringinfo[!,:family] .== popid)
        for chr=1:nchr            
            neworder=getvalentorder(homologmapls[chr],popstatespace)
            chrvvprob = polyancestry.valentprob[chr]
            for off=offls
                chrvvprob[off][:,1] = neworder[chrvvprob[off][:,1]]
            end
        end
    end
end

"""
    setabsphase!(refhapfile,polyest::Union{PolyGeno, PolyAncestry}; workdir) 

set absolute parental phases, based on the reference haplotype file.

# Positional arguments

`refhapfile::AbstractString`: reference haplotype file has the same format
as the input genofile, except that parental genotypes are phased and offspring
genotypes are ignored if they exist.

`polyest::Union{PolyGeno, PolyAncestry}`: a struct with phased parent genotypes .

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for reading `refhapfile`.

"""
function setabsphase!(refhapfile::AbstractString, polyest::Union{PolyGeno, PolyAncestry};
    workdir::AbstractString=pwd(),    
    commentstring::AbstractString="##")
    refparenthap= readparenthaplo(refhapfile,polyest.polyped.parentinfo,
        polyest.markermap;  workdir, commentstring)        
    setabsphase0!(refparenthap,polyest)
end

function setabsphase0!(refparenthap::AbstractVector,polyest::Union{PolyGeno, PolyAncestry})    
    # homologmapls[chr]: matrix size nparent x 2 where 1st col is parent index, and 2nd col is the estimated permutation of homologs]
    # homologmapls[chr][p,:] = [p, perm], [i[perm] for i in polyest.parentgeno[chr][:,p]] == refparenthap[chr][:,p]
    # chridls =  [i[1,:linkagegroup] for i=polyancestry.markermap]
    homologmapls,newestparenthap = calabsolutehap(polyest.parentgeno,refparenthap,polyest.polyped)    
    # for (i,homologmap) in enumerate(homologmapls)
    #     @info string("chr=",i) homologmap
    # end
    polyest.parentgeno = newestparenthap 
    if isa(polyest, PolyAncestry)
        orderviterbipath!(polyest,homologmapls)
        ordercondprob!(polyest,homologmapls)
        isnothing(polyest.valentprob) || ordervalentprob!(polyest,homologmapls)
    end
    # if haskey(polyest.misc, "correction")        
    #     ppidls = polyest.polyped.parentinfo[!,:individual]
    #     ordercorrecion!(polyest.misc["correction"],homologmapls,chridls,ppidls)
    # end   
    refparenthap
end 

function orderviterbipath!(polyancestry::PolyAncestry,homologmapls::AbstractVector)
    # homologmapls[chr][est_i] = [ref_parentindex, est_permutation of homologs]
    # homologmapls[chr][2] = [1, [4,3,2,1]]: # esthaploch[2][4,3,2,1] matches refhaploch[1]    
    nch = length(polyancestry.markermap)
    for pop = 1:size(polyancestry.polyped.designinfo,1)
        popid = polyancestry.polyped.designinfo[pop,:family]
        popstatespace = polyancestry.statespace[popid]
        diplostate = popstatespace["state"]
        staterule = Dict(diplostate .=> 1:length(diplostate))
        offls = findall(polyancestry.polyped.offspringinfo[!,:family] .== popid)
        for ch=1:nch            
            homologdict = gethomologdict(homologmapls[ch],popstatespace)
            neworder = [staterule[[homologdict[i] for i=j]] for j=diplostate]
            polyancestry.viterbipath[ch][:,offls] .= neworder[polyancestry.viterbipath[ch][:,offls]] 
        end
    end
    nothing
end

function savephaseerr(phaseerrdf::DataFrame; 
    outid::AbstractString="", 
    missingstring::AbstractString="NA",
    io::Union{IO,Nothing}=nothing, verbose::Bool=true)
    ismissing(phaseerrdf[end,"allchr"]) && return nothing
    keyls = keys(phaseerrdf[end,"allchr"])
    startcol = findfirst(==("allchr"), names(phaseerrdf))    
    for isdose in [true,false]        
        if isdose
            df = deepcopy(phaseerrdf)            
            rg = 1:2
            for j in startcol:size(df,2)                
                df[!,j] .= [ismissing(i) ? i : join(values(i)[rg],"|") for i in df[!,j]]
            end            
            msg_phase = string("PolyOrigin,",outid*"unphasederror",",posterior_unphased_error=",join(keyls[rg],"|"))            
        else
            df = deepcopy(phaseerrdf)
            rg = 3:length(keyls)
            for j in startcol:size(df,2)                
                df[!,j] .= [ismissing(i) ? i : join(values(i)[rg],"|") for i in df[!,j]]
            end            
            msg_phase = string("PolyOrigin,",outid*"phasingerror",",posterior_phasing_error=",join(keyls[rg],"|"))            
        end
        if verbose
            @info msg_phase 
            show(df, show_row_number=true, eltypes=false, summary=false)
            println()         
        end
        if !isnothing(io) 
            printconsole(io, false, msg_phase)
            CSV.write(io, df; append=true, header=true,missingstring)
        end
    end
    nothing
end

function savephaseerr(refhapfile::AbstractString,polyest::Union{PolyGeno,PolyAncestry,JoinedPolyAncestry};     
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem")
    refparenthap = setabsphase!(refhapfile,polyest; workdir)
    phaseerrdf = calphaseerr(refparenthap,polyest)    
    outfile = getabsfile(workdir, outstem * "_phaseerr.csv")   
    open(outfile, "w") do io
        savephaseerr(phaseerrdf; io,verbose=false)
    end
    phaseerrdf, outfile
end    

function calancestryerror(trueorig::PolyGeno, joinedancestry::JoinedPolyAncestry; 
    ismidparent::Bool=true)
    polyest = joinedancestry;
    chridls =  [i[1,:linkagegroup] for i=polyest.markermap]
    nchr = length(chridls)
    if ismidparent
        # [!,[:individual,:isfounder,:family,:generation, :immediate_progeny]]
        infodf = polyest.polyped.parentinfo
        trueorigls = trueorig.parentgeno
        estorigls = polyest.parentviterbipath
    else
        # [!,[:individual, :family, :generation, :isoutlier]]
        infodf = polyest.polyped.offspringinfo
        trueorigls = trueorig.offspringgeno
        estorigls = polyest.viterbipath
    end
    nind = size(infodf,1)
    res = [begin 
        chrtruehap = trueorigls[chr]
        chresthap= estorigls[chr]
        [calphaseerr(chrtruehap[:,p],chresthap[:,p]) for p in 1:nind]
    end for chr in 1:nchr]
    phaseres2df(res, chridls, infodf; isancestry=true)    
end


function saveancestryerror(trueorigfile::AbstractString,joinedancestry::JoinedPolyAncestry;     
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem")
    trueorig = readtruegeno!(trueorigfile,joinedancestry; workdir,verbose=false)        
    ancestryerr_midp = calancestryerror(trueorig,joinedancestry; ismidparent = true)
    ancestryerr_off = calancestryerror(trueorig,joinedancestry; ismidparent = false)
    outfile = abspath(workdir, outstem*"_ancestryerr.csv")
    open(outfile, "w") do io    
        savephaseerr(ancestryerr_midp; outid = "midparent_", io, verbose=false)
        savephaseerr(ancestryerr_off; outid = "offspring_", io, verbose=false)
    end
    outfile
end


function phaseres2df(phaseres::AbstractVector, chridls::AbstractVector, infodf::DataFrame; isancestry::Bool=false)
    keyls = keys(phaseres[1][1])
    if isancestry
        isparent = in("isfounder",names(infodf))
        if isparent 
            fls = findall(infodf[!,:isfounder])
            res2 = [[in(i,fls) ? missing : chrres[i] for i in eachindex(chrres)] for chrres in phaseres]
        else
            res2 = phaseres
        end
    else
        res2 = phaseres
    end
    df = DataFrame(chridls .=> res2)    
    sumls = [begin 
        valls = skipmissing([ismissing(j) ? j : [values(j)...] for j in row])    
        isempty(valls) ? missing : (;zip(keyls,sum(valls))...)
    end for row in eachrow(df)]
    insertcols!(df, 1, "allchr"=>sumls) 
    df = hcat(infodf,df)
    # push a summary row    
    push!(df,last(df))    
    for j in setdiff(names(infodf),[:individual,:family])
        df[!,j] .= vcat(df[1:end-1,j],[missing])  
    end    
    df[end,:individual] = "sum_inds"
    df[end,:family] = "sum_families"
    for j in (size(infodf,2)+1):size(df,2)
        valls = skipmissing([ismissing(k) ? k : [values(k)...] for k in df[1:end-1,j]])            
        df[end,j] = isempty(valls) ? missing : (;zip(keyls,sum(valls))...)         
    end
    df
end

function calphaseerr(refparenthap::AbstractVector,
    polyest::Union{PolyGeno, PolyAncestry,JoinedPolyAncestry})
    chridls =  [i[1,:linkagegroup] for i=polyest.markermap]
    nchr = length(chridls)    
    parentinfo = polyest.polyped.parentinfo
    nparent = size(parentinfo,1)
    res = [begin 
        chrtruehap = refparenthap[chr]
        chresthap= getchrparentgeno(polyest.parentgeno[chr])
        [calphaseerr(chrtruehap[:,p],chresthap[:,p]) for p in 1:nparent]
    end for chr in 1:nchr]
    phaseres2df(res, chridls, parentinfo)    
end



function calphaseerr(true_phasedgeno::AbstractVector,est_phasedgeno::AbstractVector)    
    perm = get_absphase(true_phasedgeno,est_phasedgeno)
    est_phasedgeno_perm = [[i[perm]...] for i in est_phasedgeno]    
    b = map((x,y)-> sort(x) == sort(y), true_phasedgeno,est_phasedgeno_perm) 
    ngeno = length(b) - sum(ismissing.(b))
    nunphasederr = ngeno - sum(skipmissing(b))
    ishetero = map((x,y) -> ismissing(x) ? false : x && !allequal(y), b, true_phasedgeno)
    true_phasedgeno2 = true_phasedgeno[ishetero]
    est_phasedgeno2 = est_phasedgeno_perm[ishetero]
    nhetero = sum(ishetero)
    nheteroerr = nhetero == 0 ? 0 : sum(true_phasedgeno2 .!= est_phasedgeno2)
    npermerr = 0 # default,initial value
    if nheteroerr > 0 
        # calculate permutation error
        true_heterols = true_phasedgeno2
        est_heterols = copy.(est_phasedgeno2)        
        ploidy = length(first(true_heterols))
        permls = collect(permutations(1:ploidy))
        heterols = tuple.(true_heterols,est_heterols)
        segls = splitindex(heterols)    
        if length(segls) > 1
            nmarker = length(true_heterols)
            for i in 1:(length(segls)-1)
                now = first(segls[i])
                next = first(segls[i+1])
                trueg = true_heterols[now]
                estg = est_heterols[now]
                if trueg != estg
                    b = [estg[p] == trueg for p in permls]
                    permls2 = permls[b]
                    true_heterols2 = view(true_heterols,next:nmarker)
                    dls = [begin             
                        est_heterols2 = [est_heterols[j][perm] for j in next:nmarker]
                        sum(true_heterols2 .!= est_heterols2) 
                    end for perm in permls2]                    
                    perm = permls2[argmin(dls)]        
                    i > 1 && (npermerr += 1)                    
                    # println("i=",i, ", perm=", perm, ", trueg=", trueg, ", estg=", estg, ", dls=",dls, ",npermerr=",npermerr)
                    for m in now:nmarker 
                        est_heterols[m] .= est_heterols[m][perm]
                    end                    
                end
            end
            now = length(segls)
            trueg = true_heterols[now]
            estg = est_heterols[now]
            trueg == estg || (npermerr += 1)
        end
    end
    (nunphasederr=nunphasederr,ngeno=ngeno, nheteroerr=nheteroerr,npermerr=npermerr, npermoppty= max(0,nhetero-1))
end

function phasedgeno_vector2mtx(est_phasedgeno::AbstractVector)
    if isa(first(est_phasedgeno), Tuple)
        reduce(hcat, [[i...] for i in est_phasedgeno])'
    elseif isa(first(est_phasedgeno), AbstractVector)        
        permutedims(reduce(hcat,est_phasedgeno))
    else
        error(string("uknown est_phasedgeno: ",est_phasedgeno))
    end
end


function get_absphase(true_phasedgeno::AbstractVector,est_phasedgeno::AbstractVector)    
    truedata = phasedgeno_vector2mtx(true_phasedgeno) 
    estdata = phasedgeno_vector2mtx(est_phasedgeno)
    perm=[argmin([sum(skipmissing(i .!= col))  for col =eachcol(estdata)])  for i=eachcol(truedata)]
    if !allunique(perm)
        ploidy = length(perm)
        perms = collect(permutations(1:ploidy))
        dis = [sum(skipmissing(estdata[:,perm] .!= truedata)) for perm=perms]
        perm = perms[argmin(dis)]
    end
    perm
end

