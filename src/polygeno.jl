"""
    PolyGeno

mutable struct that stores genotypic data and pedigree information

    PolyGeno(markermap,parentgeno,offspringgeno,parentinfo,offspringinfo,designinfo,delmarker,correction)

constructor from genotypic data and pedigree information.

# Fields

`markermap::Vector{DataFrame}`: marker map for each chromosome. `markermap[c]`
gives the dataframe of chromosome c with columns [:marker, :chromosome, :position].

`parentgeno::Vector{Matrix}`: genotypic data in parents. `parentgeno[c][m,i]`
gives the genotype of parent i at marker m of chromosome c.

`offspringgeno::Vector{Matrix}`: genotypic data in offspring. `offspringgeno[c][m,i]`
 gives the genotype of offspring i at marker m of chomosome c.

`parentinfo::DataFrame`: parent information with columns [:individual, :ploidy].

`offspringinfo::DataFrame`: offspring information with columns
[:individual,:population,:isoutlier].

`designinfo::DataFrame`: design information with columns
[:population, :parent1, :parent2, ...]. The dataframe element `n_ij` denotes
the nubmer of gametes contributed to each offspring in i-th population by j-th parent.
`n_ij`=2 means that i-th population is produced by self-fertilization of
j-th parent. Row sum must be 2.

`delmarker::DataFrame`: dataframe for collecting markers that were removed from markermap,
in the stages of parental phasing or map refinement. It has columns
[:marker, :chromosome, :position].

`correction::DataFrame`: dataframe for collecting parental genotype corrections.
It has columns [:round,:marker,:chromosome,:parent,:old_genotype,:new_genotype,:old_nerr,:new_nerr].

!!! note
    A genotype in `offspringgeno` must be one of following formats 1-3, and
    a genotype in `parentgeno` must be one of following formats 1-4:

    1. `dosage`: ranges from 0, 1, ..., ploidy, and missing for missing dosage;
    2. `readcount`:  [c1, c2] where c1 and c2 are the number of reads for
        alleles 1 and 2, respectively. Missing genotypes are denoted by [0,0].
    3. `probability`: [p0, p1, ...], where p_i denotes the probability of observed
       data given dosage i = 0, ..., ploidy, and the probabilities are normalized so that
       their sum is 1.
    4. `phasedgeno`: [g1, g2,...],  where g_i=1 or 2 for i=1,..., ploidy.

"""
mutable struct PolyGeno
    # markermap  [chr][:marker,:chromosome,:position]
    markermap::Vector{DataFrame}
    # parentgeno: [chr][marker,parent]
    parentgeno::Vector{Matrix}
    # offspringgeno: [chr][marker,offspring]
    offspringgeno::Vector{Matrix}
    #parentinfo [:individual,:ploidy]
    parentinfo::DataFrame
    #offspringinfo [:individual, :population]
    offspringinfo::DataFrame
    #designinfo [:population, :parent1, :parent2, ...]
    designinfo::DataFrame
    # delmarker  [:marker,:chromosome,:position]
    delmarker::DataFrame
    # correction: [:marker,:chromosome,:parent,:old_genotype,:new_genotype, :old_nerr,:new_nerr]
    correction::DataFrame
    function PolyGeno(markermap,parentgeno,offspringgeno,parentinfo,
            offspringinfo,designinfo,delmarker,correction)
        checkgenodim(parentgeno) || @error("inconsistent parental genotypes")
        checkgenodim(offspringgeno) || @error("inconsistent offspring genotypes")
        dims = [length(markermap),size(parentgeno,1),size(offspringgeno,1)]
        if length(union(dims))!=1
            error("inconsistent number of linkage groups ",dims)
        end
        dims = [size.(markermap,1),size.(parentgeno,1),size.(offspringgeno,1)]
        if length(union(dims))!=1
            error("inconsistent number of markers within linkage groups ",dims)
        end
        dims = [union(size.(offspringgeno,2)),[size(offspringinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of offspring ",dims)
        end
        dims = [union(size.(parentgeno,2)),[size(parentinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of parents ",dims)
        end
        new(markermap,parentgeno,offspringgeno,parentinfo,
            offspringinfo,designinfo,delmarker,correction)
    end
end

"""
    readPolyGeno(genofile, pedfile, keyargs...)

reads input files and returns polygeno::PolyGeno.

# Positional arguments

`genofile::AbstractString`: filename for genotypic data. For example, a CSV formatted
genofile looks like

```
    marker, chromosome, pos, ind1,ind2, ind3, ...
    snp1, 1, 0.14, 0, 2, 1, 4, ...
    snp2, 1, 0.16, 4, 0, NA, 2, ...
    snp3, 1, 0.21, NA, 3, 0, 1, ...
```
!!! note
    * The first three columns specify the genetic map. Marker IDs in
      column 1 must be unique, chromosome IDs in column 2 must be
      consecutive, and positions (in unit of centi-Morgan or base pair)
      of markers in column 3 must be non-descreasing within a chromosome.

    * The rest of columns give the genotypes of sampled individuals. The indvidual
      IDs must be unique. The genotypes of all parents must be represented by
      one of the following formats 1-4, and the genotypes of all offspring must
      be represented by one of the following formats 1-3.

      1. `dosage`: ranges from 0, 1, ..., ploidy, and NA for missing dosage;
      2. `readcount`: c1|c2, where c1 and c2 are the number of reads for alleles
         1 and 2, respectively. Missing genotypesare given by 0|0
      3. `probability`: p(0)|p(1)|...|p(ploidy), where p(i) denotes the probability
          of observed data given dosage i = 0, ..., ploidy, and the probabilities
          are normalized so that their sum is 1.
      4. `phasedgeno`: g(1)|g(2)|...|g(ploidy),  where g(i)=1 or 2 for i=1,..., ploidy.

    * All individuals must be in the `pedfile`.

`pedfile::AbstractString`: filename for pedigree information. For example, a CSV formatted
pedfile looks like
```
    individual, population, motherid, fatherid, ploidy
    P1, 0, 0, 0, 4
    P3, 0, 0, 0, 4
    P3, 0, 0, 0, 4
    offspring1, pop1, P1, P2, 4
    offspring2, pop1, P1, P2, 4
    offspring3, pop2, P1, P3, 4
    offspring4, pop2, P1, P3, 4
    offspring5, pop3, P2, P3, 4
    offspring6, pop4, P3, P3, 4
```
!!! note
    * The pedigree contains three founders (parents), two offspring from the cross
      beween parents 1 and 2, two offspring from the cross between parents 1 and 3,
      one offspring from the cross between parents 2 and 3, and one offspring from the selfing of
      parent 3.

    * All individual IDs in column 1 must be unique,
      column 2 denotes ID for the founder population and IDs for each F1 cross or selfing,
      columns 3 and 4 denotes the parents of each sub-population (motherID and fatherID
      of founders are set to 0), and column 5 denotes the ploidy level.

    * All parents and all offspring must be in the `genofile`.

# Keyword arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`delimchar::AbstractChar=','`:  text delimiter.

`commentstring::AbstractString="#"`: rows that begins with commentstring will be ignored.

`isphysmap::Bool=false`: true if input markermap is physical map, where
marker locations are in unit of base pair(bp).

`recomrate::Real=1` recombination rate in unit of 1 cM/Mbp (centiMorgan per million base pair).

`workdir::AbstractString = pwd()`: directory for reading genofile and pedfile

"""
function readPolyGeno(inpugenofile::AbstractString, inpupedfile::AbstractString;
    isphysmap::Bool=false,
    recomrate::Real=1, # 1 cM per Mbp
    missingstring::AbstractString="NA",
    delimchar::AbstractChar=',',
    commentstring::AbstractString="#",
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    genofile = getabsfile(workdir,inpugenofile)
    isfile(genofile) || @error(string(inpugenofile," does not exist in workdir = ",workdir))
    pedfile = getabsfile(workdir,inpupedfile)
    isfile(pedfile) || @error(string(inpupedfile," does not exist in workdir = ",workdir))
    # imort pedfile
    parentinfo,offspringinfo,designinfo=readDesign(pedfile,
            delimchar=delimchar,
            commentstring=commentstring)
    parentid=parentinfo[!,:individual]
    offid=offspringinfo[!,:individual]
    nparent = length(parentid)
    noff = length(offid)
    # imort genofile
    geno=CSV.read(genofile, DataFrame; delim=delimchar,comment=commentstring,missingstring=missingstring)
    if size(geno,2) != 3+nparent+noff
        error(string(genofile, " must contain #column = 3+#parent+#offspring = ",3+nparent+noff))
    end
    markermap = geno[!, 1:3]
    parsemarkermap!(markermap)
    typedind = string.(strip.(string.(names(geno)[4:end])))
    allunique(typedind) || @error ("individual IDs are not unique")
    d=setdiff(parentid,typedind)
    isempty(d) || @error string("parents in ", pedfile,  " but not genotyped: ", d)
    d=setdiff(offid,typedind)
    isempty(d) || @error string("offspring in ", pedfile,  " but not genotyped: ", d)
    d=setdiff(typedind,parentid,offid)
    isempty(d) || @error string("offspring genoytped but not in ", pedfile,  ": ", d)
    snpindexlist=splitindex(markermap[!,:chromosome])
    rule= Dict(typedind .=>4:(3+nparent+noff))
    parentcols=[get(rule,i,missing) for i=parentid]
    parentgeno = [Matrix(geno[i,parentcols]) for i=snpindexlist]
    parentgeno2=parseinputgeno(parentgeno,missingstring=missingstring)
    offcols=[get(rule,i,missing) for i=offid]
    offspringgeno = [Matrix(geno[i,offcols]) for i=snpindexlist]
    offspringgeno2=parseinputgeno(offspringgeno,missingstring=missingstring)
    if isphysmap
        maxloc = max(markermap[!,:position]...)
        if maxloc<10^6.0
            error(string("marker positions in physical map must be in base pair"))
        end
        markermap[!,:position] .*= recomrate*10^(-6.0)
    else
        if eltype(markermap[!,:position]) <: Real
            maxloc = max(markermap[!,:position]...)
            if maxloc>10^6.0
                error(string("marker positions in genetic map must be in centiMorgan"))
            end
        end        
    end
    markermap2=[DataFrame(i) for i=groupby(markermap,:chromosome)]
    verbose && @info string("data: #pop=", size(designinfo,1),
        ", #parent=",nparent, ", #offspring=",noff,
        ", #chr=",length(markermap2), ", #marker=",size(markermap,1))
    delmarker = getemptydf([String,String,Float64],[:marker,:chromosome,:position])
    colid = [:round,:marker,:chromosome,:parent,:old_genotype,:new_genotype, :old_nerr,:new_nerr]
    coltype = [Integer,String,String,String,String,String,Integer,Integer]
    correction = getemptydf(coltype,colid)
    PolyGeno(markermap2,parentgeno2,offspringgeno2,parentinfo,offspringinfo,
        designinfo,delmarker,correction)
end

function parsemarkermap!(markermap::DataFrame)
    rename!(markermap,[:marker,:chromosome,:position])
    for i=1:2
        markermap[!,i]=string.(strip.(string.(markermap[!,i])))
    end
    # check unqiue of marker id
    allunique(markermap[:,1]) || @error ("marker IDs are not unique")
    # check unqiue of chrid
    snpindexlist=splitindex(markermap[:,2])
    chridls = [markermap[first(i),2] for i=snpindexlist]
    if !allunique(chridls)
        @error ("some markers in the same chromosome are not consecutive")
    end
    # set chrid with same string length
    # idlen = max(length.(markermap[!,2])...)
    # markermap[!,2]=[lpad(i,idlen) for i=markermap[!,2]]
    # check position
    if eltype(markermap[!,3]) <: Real
        markermap[!,3]=Float64.(markermap[!,3])
        bool=[begin
            A=markermap[i,3]
            all(A[1:end-1] .<= A[2:end])
        end for i=snpindexlist]
        if !all(bool)
            wrongchrs = chridls[.!bool]
            @error string("markers positions in chrs=",wrongchrs, " are not in non-descreasing order")
        end
    end
    markermap
end

# TODO: more testing for genodata.
function parseinputgeno(geno::AbstractVector;missingstring::AbstractString="NA")
    kind = [eltype(geno[i]) for i=1:length(geno)]
    if length(unique(kind)) != 1
        @error(string("inconsistent genotypes among chromosomes: ",kind))
    end
    kind2 = kind[1]
    if kind2 <: AbstractString
        geno2=[strip.(i) for i=geno]
        isfloat= any([occursin(".",string(i...)) for i = geno2])
        if isfloat
            # missing is not allowed in probability
            type = Float64
            geno4 = [[map(x->parse.(type,x),i) for i=split.(j,"|")] for j=geno2]
        else
            # missing is allowed in read count
            type = Int
            geno3 = [[map(x->parse.(type,replace(x,missingstring=>"-1")),i) for i=split.(j,"|")] for j=geno2]
            geno4 = [replace.(i,-1=>missing) for i=geno3]
        end
        return geno4
    elseif kind2 <: Union{Missing,Integer}
        return geno
    else
        @error(string("unknown genotype data type: ",kind2))
    end
end

function checkgenodim(geno::AbstractVector)
    kind = eltype(geno[1])
    if kind <: Union{Missing,Integer}
        true
    elseif kind <: AbstractVector
        # missing is not allowed
        len = [unique(length.(i)) for i=geno]
        length.(unique(len)) == [1]
    else
        @error(string("unknown genotype data type: ",kind))
    end
end

function kindofgeno(geno::AbstractVector)
    type = eltype(geno[1])
    if type <: Union{Missing,Integer}
        return "dosage"
    elseif type <: AbstractVector
        subtype = eltype(type)
        if subtype <:  Union{Missing,Integer}
            bool = unique(length.(geno[1][1,:])) == [2]
            bool2 = issubset(skipmissing(unique(vcat(geno[1]...))),[1,2])
            return  (bool && !bool2) ? "readcount" : "phasedgeno"
        elseif subtype <: AbstractFloat
            return "probability"
        else
            @error string("unknow geno data type: ", type)
        end
    else
        @error string("unknow geno data type: ", type)
    end
end


function readDesign(pedfile::AbstractString;
    delimchar::AbstractChar=',',
    commentstring::AbstractString="#")
    # not considered: the case of offspring with only one parent, that is, gamete instead of zygote.
    if !isfile(pedfile)
        error("pedfile: ", pedfile, " does not exist")
    end
    design=CSV.read(pedfile,DataFrame; delim=delimchar,comment=commentstring)
    if size(design,2)<5
        error("At least 5 columns are required in ",pedfile)
    else
        design=design[:,1:5]
        rename!(design,[:individual,:population,:mother,:father,:ploidy])
        for i=1:4
            design[!,i]=string.(strip.(string.(design[!,i])))
        end
    end
    allunique(design[!,1]) || @error("individual IDs are not unique")
    isparent = [design[i,3] == "0" && design[i,4] == "0" for i=1:size(design,1)]
    isparent .+= [design[i,3] == "NA" && design[i,4] == "NA" for i=1:size(design,1)]
    isoff  = .!isparent
    d=setdiff(design[isoff,:mother],design[!,:individual])
    isempty(d) || @error("offspring have unknown mothers  ", d)
    d=setdiff(design[isoff,:father],design[!,:individual])
    isempty(d) || @error("offspring have unknown father  ", d)
    poplist=design[isoff,:population]
    popindexlist=[findall(poplist .== i) for i=unique(poplist)]
    parentinfo=design[isparent,[:individual,:ploidy]]
    offspringinfo=design[isoff,[:individual,:population,:ploidy]]
    outlier = Vector{Union{Missing,Bool}}(missing,size(offspringinfo,1))
    offspringinfo[!,:isoutlier] = outlier
    dict=Dict(parentinfo[!,:individual] .=> parentinfo[!,:ploidy])
    for rg = popindexlist
        subpop = unique(design[isoff, :][rg,2:5])
        size(subpop,1) == 1 || @error string("parents and ploidy levels in a subpopulation are inconsistent: ", subpop)
        parentploidy = [get(dict, i, missing) for i=Vector(subpop[1,2:3])]
        offploidy = sum(parentploidy) ÷ length(parentploidy)
        offploidy == subpop[1,end] || @error string("offspring ploidy = ",offploidy, " is consisstent with those of parents = ",parentploidy)
    end
    popparent=design[isoff,[:population,:mother,:father]][first.(popindexlist),:]
    designinfo=getdesigninfor(popparent,parentinfo)
    (parentinfo=parentinfo,offspringinfo=offspringinfo,designinfo=designinfo)
end

function getdesigninfor(popparent::DataFrame,parentinfo::DataFrame;
    missingstring::String="NA")
    npop=size(popparent,1)
    nfounder = size(parentinfo,1)
    rule=Dict(parentinfo[!,:individual] .=> 1:nfounder)
    crossparent=setdiff(Matrix(popparent[!,2:3])[:],[missingstring])
    d=setdiff(crossparent,parentinfo[!,:individual])
    if !isempty(d)
        error(string("parents that are used in crosses but not listed:", d, ",pp=",parentinfo[!,:individual]))
    end
    d=setdiff(parentinfo[!,:individual],crossparent)
    if !isempty(d)
        error(string("parents that are listed in ",pedfile,"  but not used in crosses:", d))
    end
    res = zeros(Int,npop,nfounder)
    for i=1:npop,j =2:3
        c=get(rule,popparent[i,j],0)
        c>0 && (res[i,c]+=1)
    end
    df=DataFrame(res,Symbol.(parentinfo[!,:individual]))
    insertcols!(df, 1, :population=> convert.(String,popparent[!,:population]))
    df
end

function getploidy(polygeno::PolyGeno, pop::Integer)
    ff=polygeno.designinfo[pop,2:end]
    popfounder = convert.(Bool,ff)
    nn=polygeno.parentinfo[:ploidy][popfounder]
    ploidy=Int(mean(nn))
    sum(ff) == 1 ? (ploidy ÷ 2) : ploidy
end

getploidy(polygeno::PolyGeno) =
    [getploidy(polygeno,pop) for pop=1:size(polygeno.designinfo,1)]

function getGenodf(polygeno::PolyGeno;missingstring::AbstractString="NA")
    markermap = Matrix(vcat(polygeno.markermap...))
    nch = length(polygeno.parentgeno)
    kind = kindofgeno(polygeno.parentgeno)
    if  kind in  ["phasedgeno", "readcount", "probability"]
        # fgeno =vcat([stringjoin.(polygeno.parentgeno[ch],"|") for ch=1:nch]...)
        fgeno = vcat([begin
            a=join.(polygeno.parentgeno[ch],"|")
            ii=findall(occursin.("missing",a))
            for i=ii
                a[i] = replace(a[i],"missing"=>missingstring)
            end
            a
        end for ch=1:nch]...)
    elseif kind in  ["dosage"]
        fgeno =vcat(polygeno.parentgeno...)
    else
        @error string("unkown data type of parental genodata: ",kind)
    end
    kind = kindofgeno(polygeno.offspringgeno)
    if  kind in  ["phasedgeno", "readcount", "probability"]
        offgeno =vcat([stringjoin.(polygeno.offspringgeno[ch],"|") for ch=1:nch]...)
    elseif kind in  ["dosage"]
        offgeno =vcat(polygeno.offspringgeno...)
    else
        @error string("unkown data type of offspring genodata: ",kind)
    end
    geno=hcat(markermap,fgeno,offgeno)
    parentid = Symbol.(polygeno.parentinfo[!,:individual])
    offid = Symbol.(polygeno.offspringinfo[!,:individual])
    rowid = vcat(propertynames(polygeno.markermap[1]),parentid,offid)
    DataFrame(geno,rowid)
end

"""
    savegenodata(outfile,polygeno,missingstring="NA",workdir=pwd())

saves genotypic data of the struct polygeno into outfile

# Positional arguments

`outfile::AbstractString`: filename for saving results.

`polygeno::PolyGeno`: a struct returned by [`readPolyGeno`](@ref).

# Keyward arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for writing outfile.

"""
function savegenodata(outfile::AbstractString,polygeno::PolyGeno;
    missingstring::AbstractString="NA",
    workdir::AbstractString = pwd())
    genodf =getGenodf(polygeno,missingstring=missingstring)
    outputfile2 = getabsfile(workdir,outfile)
    CSV.write(outputfile2, genodf; missingstring=missingstring,
        header=true,append=false)
end

function rawDoseCall!(polygeno::PolyGeno; seqerr::Real=0.001)
    errdigit0 = split(last(split(@sprintf("%.f",seqerr),".")),"")
    errdigit2=findfirst(x->x!="0",errdigit0)
    errdigit = isnothing(errdigit2) ? 6 : max(4,errdigit2+1)
    if kindofgeno(polygeno.parentgeno) == "readcount"
        ploidyls = polygeno.parentinfo[!,:ploidy]
        polygeno.parentgeno = [calreadpostprob(i,ploidyls,seqerr,digits=errdigit) for i=polygeno.parentgeno]
    end
    if kindofgeno(polygeno.offspringgeno) == "readcount"
        ploidyls = polygeno.offspringinfo[!,:ploidy]
        polygeno.offspringgeno = [calreadpostprob(i,ploidyls,seqerr,digits=errdigit+2) for i=polygeno.offspringgeno]
    end
    polygeno
end

function getchrid(polygeno::PolyGeno)
    [i[1,:chromosome] for i=polygeno.markermap]
end

function getsubPolyGeno!(polygeno::PolyGeno;
        chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
        snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing)
    markermap = polygeno.markermap
    nchr=length(markermap)
    nsnp = max(size.(polygeno.markermap,1)...)
    chrsubset2 = isnothing(chrsubset) ? (1:nchr) : chrsubset[chrsubset.<= nchr]
    snpsubset2 = isnothing(snpsubset) ? (1:nsnp) : snpsubset
    snpsubsetlist = [snpsubset2[snpsubset2 .<= size(markermap[ch],1)] for ch=chrsubset2]
    polygeno.markermap= map((x,y)->markermap[x][y,:],chrsubset2,snpsubsetlist)
    polygeno.parentgeno=map((x,y)->polygeno.parentgeno[x][y,:],chrsubset2,snpsubsetlist)
    polygeno.offspringgeno=map((x,y)->polygeno.offspringgeno[x][y,:],chrsubset2,snpsubsetlist)
    polygeno
end

function getsubPolyGeno(polygeno::PolyGeno;
        chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
        snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing)
    markermap = polygeno.markermap
    nchr=length(markermap)
    nsnp = max(size.(polygeno.markermap,1)...)
    chrsubset2 = isnothing(chrsubset) ? (1:nchr) : chrsubset[chrsubset.<= nchr]
    snpsubset2 = isnothing(snpsubset) ? (1:nsnp) : snpsubset
    snpsubsetlist = [snpsubset2[snpsubset2 .<= size(markermap[ch],1)] for ch=chrsubset2]
    markermap= map((x,y)->markermap[x][y,:],chrsubset2,snpsubsetlist)
    parentgeno=map((x,y)->polygeno.parentgeno[x][y,:],chrsubset2,snpsubsetlist)
    offgeno=map((x,y)->polygeno.offspringgeno[x][y,:],chrsubset2,snpsubsetlist)
    parentinfo=deepcopy(polygeno.parentinfo)
    offinfo = deepcopy(polygeno.offspringinfo)
    designinfo = deepcopy(polygeno.designinfo)
    PolyGeno(markermap,parentgeno,offgeno,parentinfo,offinfo,designinfo,
        deepcopy(polygeno.delmarker),deepcopy(polygeno.correction))
end
