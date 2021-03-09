
function readparenthaplo(refhapfile::AbstractString, parentinfo::DataFrame,
    markermap::Vector{DataFrame};
    workdir::AbstractString=pwd(),
    missingstring::AbstractString="NA",
    delim=',',
    comment::AbstractString="#")
    refhapfile2=getabsfile(workdir,refhapfile)
    isfile(refhapfile2) || @error(string(refhapfile," does not exist in workdir = ",workdir))
    ploidy = parentinfo[!,:ploidy]
    refhap=CSV.read(refhapfile2,DataFrame; delim=delim,comment=comment,missingstring=missingstring)
    refcols = strip.(names(refhap))
    refdict = Dict(refcols[4:end] .=> 4:length(refcols))
    refii = [get(refdict, i, nothing) for i = parentinfo[!,:individual]]
    if in(nothing, refii)
        string("parents in polyancestray but not refhapfile: ",parentinfo[isnothing.(refii),:individual])
    end
    refhap = refhap[!, vcat(1:3,refii)]
    for i=union(1:2,4:size(refhap,2))
        refhap[!,i] = string.(strip.(string.(refhap[!,i])))
    end
    vcatmap = vcat(markermap...)
    dsnp = setdiff(vcatmap[!,1],refhap[!,1])
    if length(dsnp)>0
        dsnpdf = copy(refhap[1:length(dsnp),:])
        # "0" denotes missing
        strcodemiss = replace.(replace.(Vector(refhap[1,4:end]),"1"=>"-1"),"2"=>"-1")
        for i=1:size(dsnpdf,1)
            for j=1:3
                dsnpdf[[i],j] = vcatmap[vcatmap[!,1] .== dsnp[i],j]
            end
            dsnpdf[i,4:end] .= strcodemiss
        end
        refhap=vcat(dsnpdf,refhap)
    end
    snpdict =Dict(refhap[!,1] .=> 1:size(refhap,1))
    snpindex = [get(snpdict,i,missing) for i = vcatmap[!,1]]
    refhap = refhap[snpindex,:]
    # check marker and chr
    parentid = strip.(string.(names(refhap)[4:end]),' ')
    if parentid != parentinfo[!,:individual]
        @error string("parent IDs are inconsistent")
    end
    refmap = refhap[!, 1:2]
    rename!(refmap,[:marker,:chromosome])
    for i=1:2
        refmap[!,i]=string.(strip.(string.(refmap[!,i])))
    end
    if refmap !=  strip.(vcatmap[!,1:2])
        error("marker and chromosome IDs do not match")
    end
    # cal fhaplo2
    groupref = groupby(refhap,2)
    fhaplo = [Matrix(groupref[ch][!,4:end]) for ch=1:length(groupref)]
    fhaplo2=parseinputgeno(fhaplo)
    # check if it is haplotypes
    hhcount = [unique(eachrow(length.(fhaplo2[i]))) for i=1:length(fhaplo2)]
    all(length.(hhcount) .==1) || @error("number of haplotypes inconsistent accross markers")
    all([i[1] ==ploidy for i=hhcount]) || @error("number of haplotypes mismatch ploidy")
    parenttype1 = [eltype(fhaplo2[i]) for i=1:length(fhaplo2)]
    parenttype2 = unique(parenttype1)
    if length(parenttype2) != 1
        @error(string("inconsistent data types of parental genodata: ",parenttype1))
    end
    if !(eltype(parenttype2[1]) <: Union{Missing,Integer})
        @error string("unknow data types of parental genodata: ",parenttype2[1])
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

function calabsolutehap(estparentgeno::AbstractVector,refparenthap::AbstractVector,
    parentinfo::DataFrame,designinfo::DataFrame;
    io::Union{Nothing,IOStream}=nothing,verbose::Bool=true)
    refhap = [[Matrix(hcat(hh[:,i]...)') for i=1:size(hh,2)] for hh = refparenthap]
    esthap = [[Matrix(hcat(hh[:,i]...)') for i=1:size(hh,2)] for hh = estparentgeno]
    ploidy = parentinfo[!,:ploidy]
    ccparent = connected_parents(designinfo)
    nch = length(esthap)
    nparent = size(parentinfo[!,:ploidy],1)
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
    if io !==nothing
        nphaseerr = [calnphaseerr(refhap[ch],newesthap[ch]) for ch=1:length(refhap)]
        nphaseerr2 = hcat(nphaseerr...)'
        msg = string("comparison: #mismatch phases between esthaplo and refhaplo = ",nphaseerr2)
        printconsole(io,verbose,msg)
        ndoseerr = [calndoseerr(refhap[ch],newesthap[ch]) for ch=1:length(refhap)]
        ndoseerr2 = hcat(ndoseerr...)'
        msg = string("comparison: #mismatch dosages between esthaplo and refhaplo = ", ndoseerr2)
        printconsole(io,verbose,msg)
    end
    # rulehomolog[ref_i,:] = [est_parentindex, est_permutation of homologs]
    # rulehomolog[2,:] = [1, [4,3,2,1]]: # refhaploch[2] matches esthaploch[1][4,3,2,1]
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
    homologmapls,newesthap
end

function gethomologdict(homologmap::AbstractMatrix,popstatespace::AbstractDict)
    parentindex = popstatespace["parentindex"]
    groupstate = popstatespace["groupstate"]
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
    groupstate = popstatespace["groupstate"]
    newgroupstate = [sort([get(homologdict,i,missing) for i=j]) for j=groupstate]
    rule = Dict(newgroupstate .=> 1:length(newgroupstate))
    [get(rule,i,missing) for i=groupstate]
end

function ordercondprob!(polyancestry::PolyAncestry,homologmapls::AbstractVector)
    # homologmapls[chr][est_i] = [ref_parentindex, est_permutation of homologs]
    # homologmapls[chr][2] = [1, [4,3,2,1]]: # esthaploch[2][4,3,2,1] matches refhaploch[1]
    nparent = size(polyancestry.parentinfo,1)
    nch = length(polyancestry.markermap)
    for pop = 1:size(polyancestry.designinfo,1)
        popid = polyancestry.designinfo[pop,:population]
        popstatespace = polyancestry.statespace[popid]
        offls = findall(polyancestry.offspringinfo[!,:population] .== popid)
        for ch=1:nch
            neworder = getstateorder(homologmapls[ch],popstatespace)
            for off=offls
                polyancestry.genoprob[ch][off] = polyancestry.genoprob[ch][off][:,neworder]
            end
        end
    end
    sethaploprob!(polyancestry)
end

function ordercorrecion!(correctdf::DataFrame,homologmapls::AbstractVector,
    chridls::AbstractVector,ppidls::AbstractVector)
    chridrule = Dict(chridls .=> 1:length(chridls))
    ppidrule = Dict(ppidls .=> 1:length(ppidls))
    for i=1:size(correctdf,1)
        chr = get(chridrule,correctdf[i,:chromosome],missing)
        pp = get(ppidrule,correctdf[i,:parent],missing)
        # homologmapls[chr][est_i] = [ref_parentindex, est_permutation of homologs]
        # homologmapls[chr][2] = [1, [4,3,2,1]]: # esthaploch[2][4,3,2,1] matches refhaploch[1]
        pp2, perm = homologmapls[chr][pp,:]
        correctdf[i,:parent] = ppidls[pp2]
        correctdf[i,:old_genotype] = join(split(correctdf[i,:old_genotype],"|")[perm],"|")
        correctdf[i,:new_genotype] = join(split(correctdf[i,:new_genotype],"|")[perm],"|")
    end
end

function getvalentorder(homologmap::AbstractMatrix,popstatespace::AbstractDict)
    homologdict = gethomologdict(homologmap,popstatespace)
    valentstate = popstatespace["valent"]
    valentrule = Dict(valentstate[:] .=> 1:length(valentstate))
    [begin
        vv2=[sort(map(x->sort([get(homologdict,i,missing) for i=x]),v)) for v=vv]
        get(valentrule,vv2,missing)
    end for vv=valentstate[:]]
end

function ordervalentprob!(polyancestry::PolyAncestry,homologmapls::AbstractVector)
    nchr = length(polyancestry.markermap)
    for pop = 1:size(polyancestry.designinfo,1)
        popid = polyancestry.designinfo[pop,:population]
        popstatespace = polyancestry.statespace[popid]
        offls = findall(polyancestry.offspringinfo[!,:population] .== popid)
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
    setAbsPhase!(refhapfile,phasedgeno::PolyGeno,io=nothing,verbose=true)

set absolute parental phases, based on the reference haplotype file.

# Positional arguments

`refhapfile::AbstractString`: reference haplotype file has the same format
as the input genofile, except that parental genotypes are phased and offspring
genotypes are ignored if they exist.

`phasedgeno::PolyGeno`: polygeno struct with phased parent genotypes .

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for reading `refhapfile`.

`io::Union{Nothing,IOStream}`: stream for writing log.

`verbose::Bool=true`: true if print messages on console.

"""
function setAbsPhase!(refhapfile::AbstractString, phasedgeno::PolyGeno;
    workdir::AbstractString=pwd(),
    delim=',',
    comment::AbstractString="#",
    io::Union{Nothing,IOStream}=nothing,
    verbose::Bool=true)
    refparenthap= readparenthaplo(refhapfile,
        phasedgeno.parentinfo,phasedgeno.markermap,
        workdir=workdir,delim=delim,comment=comment)
    homologmapls,newesthap = calabsolutehap(phasedgeno.parentgeno,refparenthap,
        phasedgeno.parentinfo,phasedgeno.designinfo,io=io,verbose=verbose)
    phasedgeno.parentgeno = [hcat([[h[i,:] for i=1:size(h,1)] for h=hh]...) for hh=newesthap]
    if size(phasedgeno.correction,1)>0
        chridls =  [i[1,:chromosome] for i=phasedgeno.markermap]
        ppidls = phasedgeno.parentinfo[!,:individual]
        ordercorrecion!(phasedgeno.correction,homologmapls,chridls,ppidls)
    end
    homologmapls
end

"""
    setAbsPhase!(refhapfile, polyancestry,io=nothing,verbose=true)

set absolute parental phases and consistent polyancestry.genoprob and
polyancestry.haploprob, based on the reference haplotype file.

# Positional arguments

`refhapfile::AbstractString`: reference haplotype file has the same format
as the input genofile, except that parental genotypes are phased and offspring
genotypes are ignored if they exist.

`polyancestry::PolyAncestry`: polyancestry struct with phased parent genotypes and inference genoprob.

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for reading `trueancestryfile`.

`io::Union{Nothing,IOStream}`: stream for writing log.

`verbose::Bool=true`: true if print messages on console.

"""
function setAbsPhase!(refhapfile::AbstractString,polyancestry::PolyAncestry;
    workdir::AbstractString=pwd(),
    delim=',',
    comment::AbstractString="#",
    io::Union{Nothing,IOStream}=nothing,
    verbose::Bool=true)
    refparenthap= readparenthaplo(refhapfile,
        polyancestry.parentinfo,polyancestry.markermap,
        workdir=workdir,delim=delim,comment=comment)
    setabsphase0!(refparenthap,polyancestry,io=io,verbose=verbose)
end

function setabsphase0!(refparenthaplo::AbstractVector,polyancestry::PolyAncestry;
    io::Union{Nothing,IOStream}=nothing,
    verbose::Bool=true)
    # homologmapls[est_i] = [ref_parentindex, est_permutation of homologs]
    # homologmapls[2] = [1, [4,3,2,1]]: # esthaploch[2][4,3,2,1] matches refhaploch[1]
    homologmapls,newesthap = calabsolutehap(polyancestry.parentgeno,refparenthaplo,
        polyancestry.parentinfo, polyancestry.designinfo,io=io,verbose=verbose)
    polyancestry.parentgeno = [hcat([[h[i,:] for i=1:size(h,1)] for h=hh]...) for hh=newesthap]
    ordercondprob!(polyancestry,homologmapls)
    if size(polyancestry.correction,1)>0
        chridls =  [i[1,:chromosome] for i=polyancestry.markermap]
        ppidls = polyancestry.parentinfo[!,:individual]
        ordercorrecion!(polyancestry.correction,homologmapls,chridls,ppidls)
    end
    polyancestry.valentprob ===nothing || ordervalentprob!(polyancestry,homologmapls)
    polyancestry
end
