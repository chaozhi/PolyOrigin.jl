

"""
    ErrorSetting

keyword-based struct for the setting of an error parameter. 

# Fields 

`name::String`: name of an error parameter. It is required and must be "doseerror", "allelicbias", "overdispersion", or "baseerror". 

`isinfer::Bool=true`: infer the error parameter only if isrefineerror && isinfer is true. 

`prior::Distribution`: prior distribution of an error rate. 

`init::Real`: initial error rate value. 

`softthresh::Real`: soft threshold for marker deletion. 
For "doseerror", "baseerror", or "overdispersion", a marker is deleted if 
error_rate > softthresh and it is an outlier among markers. 
For "allelicbias", a marker is deleted if f(error_rate) > f(softthresh) 
and it is an outlier. Here f(x) = max(x, 1/x). 

`hardthresh::Real`: hard threshold for marker deletion. 
Similar to softthresh except that a marker is deleted in spite of being outlier. 

"""
Base.@kwdef struct ErrorSetting
    name::String # "doseerror", "allelicbias", "overdispersion", "baseerror"
    isinfer::Bool =  seterrorfield(name, (true,true,true,false))   
    prior::Distribution = seterrorfield(name, (Beta(1,9), LogNormal(0,1), Exponential(0.5), Beta(1,49)))
    init::Union{Nothing,Real} = seterrorfield(name,(0.01,1,0.0,0.001)) # if doseerror.init == nothing, it is set to 0.001 if format is AD and otherwise 0.01 
    softthresh::Real = seterrorfield(name,(0.05,2,0.1,0.01))
    hardthresh::Real = seterrorfield(name,(0.2,10,1.0,0.05))
end

function seterrorfield(name::AbstractString, valuels::Tuple) 
    # valuels is in the order of "doseerror", "allelicbias", "overdispersion", "baseerror"
    length(valuels) == 4 || error("valuels must be of length 4")
    if name == "doseerror"
        valuels[1]
    elseif name == "allelicbias"
        valuels[2]
    elseif name == "overdispersion"
        valuels[3]
    elseif name == "baseerror"
        valuels[4]
    else
        error("unknown error name: ",name)
    end
end

# """
#     BiasSetting

# keyword-based struct for the setting of allelicbias. 
# The allelic bias is defined as the relative enrichment (>0) of alternative alleles in PCR sequencing; 
# the enrichment of reference allele is set to 1. 
# # Fields 

# `prior::Distribution = LogNormal(0,2)`: prior distribution of allelic bais that is used  
# in the inference in the 2nd step of refinement. 

# `init::Real = 1.0`: initial value in the inference. The default 1.0 denotes no allelic bias. 

# `softthresh::Real = 2`: a marker is deleted if max(allelicbias,1/allelicbias) > softthresh  AND 
# it is an outlier among markers. 

# """
# Base.@kwdef struct BiasSetting
#     prior::Distribution = LogNormal(0,2)    
#     init::Real = 1.0
#     softthresh::Real = 2
# end

"""
    PolyGeno

mutable struct that stores genotypic data and pedigree information

    PolyGeno(markermap,parentgeno,offspringgeno,parentinfo,offspringinfo,designinfo,isdelvuong,correction)

constructor from genotypic data and pedigree information.

# Fields

`polyped::PolyPed`: breeding pedigree information. See also [`PolyPed`](@ref).

`markermap::Vector{DataFrame}`: marker map for each chromosome. `markermap[c]`
gives the dataframe of chromosome c. 

`parentgeno::Vector{Matrix}`: genotypic data in parents. `parentgeno[c][m,i]`
gives the genotype of parent i at marker m of chromosome c.

`offspringgeno::Vector{Matrix}`: genotypic data in offspring. `offspringgeno[c][m,i]`
 gives the genotype of offspring i at marker m of chomosome c.

`misc::Dict{String, DataFrame}`: contains information such as (1) "deletion" dataframe
for markers that were removed from markermap, (2) "correction" dataframe for
parental error correction. 

"""
mutable struct PolyGeno
    polyped::PolyPed
    # markermap  [chr][:marker,:linkagegroup,:poscm]
    markermap::Vector{DataFrame}
    # parentgeno: [chr][marker,parent]
    parentgeno::Union{Nothing,Vector{Matrix}}
    # offspringgeno: [chr][marker,offspring]
    offspringgeno::Union{Nothing,Vector{Matrix}}
    #parentinfo [:individual,:ploidy]    
    misc::Dict{String, DataFrame}        

    function PolyGeno(polyped, markermap,parentgeno,offspringgeno, misc)
        # check genodata        
        dims = [length(markermap)]
        isnothing(parentgeno) || push!(dims, length(parentgeno))
        isnothing(offspringgeno) || push!(dims, length(offspringgeno))
        if length(union(dims))!=1
            error("inconsistent number of linkage groups ",dims)
        end
        dims = [size.(markermap,1)]
        isnothing(parentgeno) || push!(dims, size.(parentgeno,1))
        isnothing(offspringgeno) || push!(dims, size.(offspringgeno,1))
        if length(union(dims))!=1
            error("inconsistent number of markers within linkage groups ",dims)
        end
        # check dimension consistency between polyped and genodata
        if !isnothing(parentgeno)
            dims = size.(parentgeno,2)
            push!(dims,size(polyped.parentinfo,1))
            if length(unique(dims))!=1
                error("inconsistent number of parents ",dims)
            end
        end
        if !isnothing(offspringgeno)
            dims = size.(offspringgeno,2)
            push!(dims,size(polyped.offspringinfo,1))
            if length(unique(dims))!=1
                error("inconsistent number of offspring ",dims)
            end
        end
        new(polyped, markermap,parentgeno,offspringgeno,misc)
    end
end


"""
    formpolygeno(genofile, pedfile, keyargs...)

reads input files and returns polygeno::PolyGeno.

# Positional arguments

`genofile::AbstractString`: filename for genotypic data in a vcf format. 

`pedfile::AbstractString`: filename for pedigree information in a CSV format. 

# Keyword arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`delimchar::AbstractChar=','`:  text delimiter.

`commentstring::AbstractString="##"`: rows that begins with commentstring will be ignored.

`isphysmap::Bool=false`: true if input markermap is physical map, where
marker locations are in unit of base pair(bp).

`recomrate::Real=1` recombination rate in unit of 1 cM/Mbp (centiMorgan per million base pair).

`workdir::AbstractString = pwd()`: directory for reading genofile and pedfile

"""
function formpolygeno(genofile::AbstractString, pedfile::AbstractString;
    isphysmap::Bool=false,
    recomrate::Real=1, # 1 cM per Mbp
    formatpriority::AbstractVector=["AD","GP", "GT"],
    delmultiallelic::Bool=true,    
    missingstring=["NA","missing"],
    delimchar::AbstractChar=',',
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    # import pedfile
    polyped = readpolyped(pedfile; delim=delimchar, commentstring, workdir)    
    # import genofile
    genodf,commentlines = readgenodf(genofile,pedfile; formatpriority,
        delmultiallelic,commentstring, missingstring, workdir)    
    # formpoygeo
    polygeno = formpolygeno!(genodf, polyped; isphysmap, recomrate)
    if !isempty(commentlines)
        filecomment = DataFrame(sourcefile=getabsfile(workdir,genofile),
            commentlines=commentlines)
        pushmisc!(polygeno,"filecomment"=>filecomment)
    end
    polygeno
end

function pushmisc!(polygeno::PolyGeno,miscpair::Pair{String,DataFrame})
    strkey, df = miscpair
    if haskey(polygeno.misc,strkey)
        polygeno.misc[strkey] = vcat(polygeno.misc[strkey],df)
    else
        push!(polygeno.misc, miscpair)
    end
end

#  colunns in parse dataframe: ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat",
#     "multiallele", "alleles", "quality", "filter", "info", "doseerror","baseerror","allelicbias", "overdispersion"]

const _col_1sterr = 13

const _col_1stsample = 17

function formpolygeno!(genodf::DataFrame,polyped::PolyPed;    
    isphysmap::Bool=false,
    recomrate::Real=1.0)
    isfounder = polyped.parentinfo[!, :isfounder]
    founderid = polyped.parentinfo[isfounder,:individual]
    parentid = polyped.parentinfo[.!isfounder,:individual]
    offspringid = polyped.offspringinfo[!,:individual]    
    indls = names(genodf)[_col_1stsample:end]
    d= setdiff(indls,founderid, parentid,offspringid)
    if !isempty(d)
        msg = string(length(d), " individuals in genofile by not in pedfile: ",d)
        @warn msg
    end   
    # fgenodf has colunns consistent with newparentid; offgenodif has columns consistent with newoffspringid    
    mapdf, fgenodf, newfounderid, pgenodf, newparentid, offgenodf,newoffspringid = splitgenodf(genodf,
        founderid, parentid,offspringid; isphysmap,recomrate)       
    oldparents = vcat(founderid, parentid)
    newparents = vcat(newfounderid, newparentid)
    if newoffspringid != offspringid  || oldparents != newparents       
        peddf = ped2df(polyped.ped)
        filter!(row -> in(row[:member], vcat(newparents,newoffspringid)), peddf)
        polyped = PolyPed(Pedigree(peddf))
    end 
    set_format_dosage!(mapdf, fgenodf,pgenodf, offgenodf,polyped)        
    markermap, parentgeno, offspringgeno = parsegenodf!(mapdf, fgenodf,pgenodf, offgenodf)        
    misc = Dict{String, DataFrame}()    
    polygeno = PolyGeno(polyped,markermap,parentgeno,offspringgeno,misc)
    isphysmap && sortmarker!(polygeno; byphysmap=true)
    polygeno
end

function sortmarker!(polygeno::PolyGeno; byphysmap::Bool=false)
    if byphysmap 
        chrcol,poscol = :physchrom, :physposbp
    else
        chrcol,poscol = :linkagegroup, :poscm
    end
    b = [any(ismissing.(i[!,chrcol]) .|| ismissing.(i[!,poscol])) for i in polygeno.markermap]
    if any(b)
        @warn string("physical map is missing for some markers; sortmarker! is not perform")
    else
        merge_chromosome!(polygeno)
        chridls = sort(unique(polygeno.markermap[1][!,chrcol]))
        ls = [begin 
            posls = findall(polygeno.markermap[1][!,chrcol] .== chrid)
            posls[sortperm(polygeno.markermap[1][posls,poscol])]
            # all(markermap[1][posls,:physchrom] .== chrid)
        end for chrid in chridls]
        polygeno.markermap = [polygeno.markermap[1][i,:] for i in ls]
        if !isnothing(polygeno.parentgeno)
            polygeno.parentgeno =  [polygeno.parentgeno[1][i,:] for i in ls]
        end
        if !isnothing(polygeno.offspringgeno)
            polygeno.offspringgeno = [polygeno.offspringgeno[1][i,:] for i in ls]        
        end
    end
    polygeno
end

function merge_chromosome!(polygeno::PolyGeno)
    polygeno.markermap = [reduce(vcat,polygeno.markermap)]
    if !isnothing(polygeno.parentgeno)
        polygeno.parentgeno=[reduce(vcat,polygeno.parentgeno)]
    end
    if !isnothing(polygeno.offspringgeno)
        polygeno.offspringgeno=[reduce(vcat,polygeno.offspringgeno)]
    end
    polygeno
end

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
    elseif kind in  ["GT_dosage"]
        fgeno =vcat(polygeno.parentgeno...)
    else
        @error string("unkown data type of parental genodata: ",kind)
    end
    kind = kindofgeno(polygeno.offspringgeno)
    if  kind in  ["phasedgeno", "readcount", "probability"]
        offgeno =vcat([stringjoin.(polygeno.offspringgeno[ch],"|") for ch=1:nch]...)
    elseif kind in  ["GT_dosage"]
        offgeno =vcat(polygeno.offspringgeno...)
    else
        @error string("unkown data type of offspring genodata: ",kind)
    end
    geno=hcat(markermap,fgeno,offgeno)
    parentid = Symbol.(polygeno.polyped.parentinfo[!,:individual])
    offid = Symbol.(polygeno.polyped.offspringinfo[!,:individual])
    rowid = vcat(propertynames(polygeno.markermap[1]),parentid,offid)
    DataFrame(geno,rowid)
end


"""
    savegenodata(outfile,polygeno; missingstring="NA",workdir=pwd(),delim=',')

save genotypic data of polygeno::PolyGeno into outfile

# Positional arguments

`outfile::AbstractString`: output genofile for saving genotypic data.

`polygeno::PolyGeno`: a struct returned by [`formpolygeno`](@ref).

# Keyward arguments

`commentstring::AbstractString="##"`: lines starting with commentstrig are comment lines. 

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for writing outfile.

`delim::AbstractChar=','`: delimitor character. For vcf outfile, set delim = '\t'. 

"""
function savegenodata(sink::Union{IO,AbstractString},polygeno::PolyGeno;
    target::AbstractString="all",
    keepcomment::Bool=true,
    workdir::AbstractString = pwd(),
    delim::AbstractChar = ',',
    commentstring::AbstractString="##",
    missingstring::AbstractString="NA")    
    in(target, ["all","founder","parent","offspring"]) || @error string("unknown target=",
        target, ", target must be in [all, founder, parent, offspring]")
    delim2 = delim
    if typeof(sink) <: AbstractString
        _, ext = split_allext(sink)
        myopen = ext in [".vcf.gz",".csv.gz"] ? GZip.open : open
        if ext in [".vcf",".vcf.gz"]
            delim2 = '\t'
        elseif ext in [".csv",".csv.gz"]
            delim2 = ','
        else
            @error string("outfile extension = ", ext, ", not .vcf, .vcf.gz, .csv, or .csv.gz")
        end
        outputfile2 = getabsfile(workdir,sink)
        io = myopen(outputfile2, "w")
    elseif typeof(sink) <: IO
        io = sink
    else
        error(string("unknown sink: ",sink))
    end
    isvcf = delim2 == '\t'
    commentlines = haskey(polygeno.misc,"filecomment") ? split(polygeno.misc["filecomment"][1,:commentlines], "\n") : String[]    
    # @info "bef" commentlines
    if isempty(commentlines)
        commentlines = [commentstring*"fileformat=VCFv4.3", commentstring*"source=PolyOrigin"]
    else
        if !startswith(first(commentlines), commentstring*"fileformat")            
            pushfirst!(commentlines, commentstring*"fileformat=VCFv4.3")
        end
        if !any([startswith(i, commentstring*"source") for i in commentlines])
            insert!(commentlines, 2, commentstring*"source=PolyOrigin")
        end
    end       
    if !any([startswith(i, commentstring*"INFO=<ID=LINKAGEGROUP") for i in commentlines])
        push!(commentlines, commentstring*"INFO=<ID=LINKAGEGROUP,Number=1,Type=String,Description=\"Linkage group in genetic map\">")
    end
    if !any([startswith(i, commentstring*"INFO=<ID=POSCM") for i in commentlines])
        push!(commentlines, commentstring*"INFO=<ID=POSCM,Number=1,Type=Float,Description=\"Genetic marker position in centiMorgan\">")
    end    
    errnamels = ["DOSEERROR","BASEERROR","ALLELICBIAS","OVERDISPERSION"]
    errnamedesc = ["Dosage error rate","Sequencing base error rate","Sequencing allelic bias","Sequencing overdispersion"]
    markerdf = isempty(polygeno.markermap) ? nothing : polygeno.markermap[1]
    for (errorname, errordesc) in zip(errnamels, errnamedesc)
        colname = lowercase(errorname)
        ismiss = isnothing(markerdf) ? true : all(ismissing.(markerdf[!,colname]))
        if !ismiss && !any([startswith(i, commentstring*"INFO=<ID=$(errorname)") for i in commentlines])
            push!(commentlines, commentstring*"INFO=<ID=$(errorname),Number=1,Type=Float,Description=\"$(errordesc)\">")
        end
    end
    if isempty(polygeno.markermap) 
        formatls = String[]
    else
        formatls = reduce(vcat, [vcat(unique(df[!,:parentformat]), unique(df[!,:offspringformat])) for df in polygeno.markermap])
        formatls = unique([i[1:2]  for i in formatls])
    end
    for format in formatls
        if !any([startswith(i, commentstring*"FORMAT=<ID=$(format)") for i in commentlines])
            if format == "GT"
               des = "Genotype"
               num = 1
               ty = "String"
            elseif format == "AD"
                des = "Allelic depths for the reference and alternate alleles in the order listed"
                num = "R"
                ty = "Integer"
            elseif format == "GP"
                des = "Genotype probabilities"
                num = "G"
                ty = "Float"
            else 
                des = "Genotype format $(format)"
                num = "."
                ty = "String"
                @warn string("unknown genotype format: ",format,", FORMAT=<ID=$(format),Number=$(num),Type=$(ty),Description=\"$(des)\"> is added to commentlines")
            end               
            push!(commentlines, commentstring*"FORMAT=<ID=$(format),Number=$(num),Type=$(ty),Description=\"$(des)\">")
        end
    end
    # @info "aft" commentlines
    keepcomment && write(io,string(join(commentlines,"\n"), "\n"))
    isunphased_2dosage = !isvcf && all([unique(i[!,:multiallele]) == [2] for i in polygeno.markermap])
    if isvcf 
        nmapcol_keep = nothing    
    else        
        nmapcol_keep = get_nmapcol_keep(polygeno)    
    end
    colnames = geno_colnames(polygeno,isvcf,nmapcol_keep, target)    
    write(io, join(colnames,delim2),"\n")    
    nchr = length(polygeno.markermap)
    for chr in 1:nchr
        mtx = togenomtx(polygeno, chr; isvcf, isunphased_2dosage, nmapcol_keep, missingstring, target)
        writedlm(io,mtx,delim2)
        flush(io)
    end
    typeof(sink) <: AbstractString && close(io)
    sink
end


function get_nmapcol_keep(polygeno::PolyGeno)
    colnames = names(polygeno.markermap[1])
    isdelcols = [if col in ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat"]
        false
    elseif col in ["alleles","doublereduction"]
        true 
    elseif col in ["multiallele"]
        unique(reduce(vcat, [unique(i[!,col]) for i in polygeno.markermap])) == [2]
    else
        all([all(ismissing.(i[!,col])) for i in polygeno.markermap])
    end for col in colnames]
    findlast(.!isdelcols)
end

function geno_colnames(polygeno::PolyGeno,isvcf::Bool,nmapcol_keep, target::AbstractString)
    if isvcf
        colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
    else
        colnames = ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat", 
            "multiallele", "alleles", "quality", "filter", "info",
            "doseerror","baseerror","allelicbias", "overdispersion", "doublereduction"]     
        if !isnothing(nmapcol_keep)  
            colnames = colnames[1:nmapcol_keep]
        end
    end    
    if !isnothing(polygeno.parentgeno) && in(target,["all","parent","founder"])
        founderid = copy(polygeno.polyped.parentinfo[!,:individual])
        if target == "founder"
            b = polygeno.polyped.parentinfo[!,:isfounder]
            keepat!(founderid,b)
        end
        append!(colnames, founderid)
    end
    if !isnothing(polygeno.offspringgeno) && in(target,["all","offspring"])
        offid = polygeno.polyped.offspringinfo[!,:individual]
        append!(colnames, offid)
    end
    colnames
end

function togenomtx(polygeno::PolyGeno,chr::Integer;
    isvcf::Bool,
    isunphased_2dosage::Bool, 
    nmapcol_keep, 
    missingstring::AbstractString="NA",
    target::AbstractString="all")
    in(target, ["all","founder", "parent","offspring"]) || @error string("unknown target=",
        target, ", target must be in [all, founder,parent, offspring]")    
    parentmtx = nothing
    if in(target, ["all","parent","founder"])
        if !isnothing(polygeno.parentgeno)
            parentmtx = Matrix{Any}(polygeno.parentgeno[chr])
            if target == "founder"
                b = polygeno.polyped.parentinfo[!,:isfounder]
                parentmtx = parentmtx[:,b]
            end
            formatls = polygeno.markermap[chr][!,:parentformat]                        
            inverse_parsegenomtx!(parentmtx,formatls,isvcf,isunphased_2dosage, missingstring)
        end        
    end
    offmtx = nothing
    calledoffmtx = nothing        
    if in(target, ["all","offspring"])
        if !isnothing(polygeno.offspringgeno)
            offmtx = Matrix{Any}( polygeno.offspringgeno[chr])
            formatvec = polygeno.markermap[chr][!,:offspringformat]
            if unique(formatvec) == ["GP"]
                for i in offmtx
                    i .= round.(i,digits=4)
                end
                calledoffmtx = biallele_gp2gt.(offmtx)
            end
            inverse_parsegenomtx!(offmtx,formatvec,isvcf,isunphased_2dosage, missingstring)            
        end
    end
    #  ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat", 
    #   "multiallele", "alleles", "quality", "filter", "info",
    #   "doseerror","baseerror","allelicbias","overdispersion"]     
    mapmtx = Matrix{Any}(polygeno.markermap[chr]);   
    submapmtx = view(mapmtx,:, vcat([3],_col_1sterr:_col_1stsample-1))    
    submapmtx .= [ismissing(i) ? i : round(i,digits=8) for i in submapmtx]
    if isvcf        
        mapmtx =  tovcfgenomtx!(mapmtx, parentmtx, offmtx,calledoffmtx)
    else
        submapmtx = view(mapmtx,:, vcat(1:_col_1stsample-1))
        submapmtx[ismissing.(submapmtx)] .= missingstring
        if target == "parent"
            mapmtx[:, 7] .= missingstring  # offspringformat            
            mapmtx[:, _col_1sterr:_col_1stsample-1] .= missingstring  # doseerror...overdispersion
        elseif target == "offspring"            
            mapmtx[:, 6] .= missingstring  # parentformat            
        end
    end    
    if !isvcf
        # recover GT_unphased format for input polygeno
        if isunphased_2dosage
            for col in [:parentformat, :offspringformat]
                polygeno.markermap[chr][!,col] .= replace.(polygeno.markermap[chr][!,col], "GT_dosage" =>"GT_unphased")
            end
        end
        if !isnothing(nmapcol_keep)
            mapmtx = mapmtx[:, 1:nmapcol_keep]
        end
    end
    b = [isnothing(parentmtx),isnothing(offmtx)]
    if b == [true,true]
        mapmtx
    elseif b == [true,false]
        hcat(mapmtx,offmtx)
    elseif b == [false,true]
        hcat(mapmtx,parentmtx)
    else
        hcat(mapmtx, parentmtx, offmtx)
    end
end


function biallele_gp2gt(gp::AbstractVector)
    i = argmax(gp) -1
    gt = ones(Int, length(gp)-1) 
    gt[i+1:end].= 0
    join(reverse(gt),"/")
end


# modify parentmtx and offspingmtx, return new vcfmapmtx
function tovcfgenomtx!(mapmtx::AbstractMatrix,
    parentmtx::Union{Nothing,AbstractMatrix},
    offspringmtx::Union{Nothing,AbstractMatrix},
    calledoffspringmtx::Union{Nothing,AbstractMatrix},
    )    
    # input mapmtx columns: 
    #  ["marker", "linkagegroup","poscm","physchrom", "physposbp","parentformat","offspringformat", 
    #   "multiallele", "alleles", "quality", "filter", "info",
    #   "doseerror","baseerror","allelicbias","overdispersion"]      
    # vcfmap colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
    # mid: col4-col9
    nsnp = size(mapmtx,1)
    vcfmapmtx = Matrix(undef, nsnp, 9)
    vcfmapmtx[:,1:3] .= mapmtx[:,[4,5,1]]
    vcfmapmtx[:,6] .= [ismissing(i) ? "." : i for i in mapmtx[:,10]] # QUAL
    vcfmapmtx[:,7] .= [ismissing(i) ? "." : i for i in mapmtx[:,11]]  # FILTER       
    submapmtx = view(vcfmapmtx,:, 1:3)
    submapmtx[ismissing.(submapmtx)] .= "."
    isnofounder = isnothing(parentmtx)
    isnooff = isnothing(offspringmtx)    
    for snp in 1: nsnp        
        # set vcf 4-8 columns: "REF","ALT","QUAL","FILTER","INFO"        
        if ismissing(mapmtx[snp, 9]) # 9 of mapmtx = "alleles"
            vcfmapmtx[snp,4] = "."
            vcfmapmtx[snp,5] = "."
        else            
            alleles = split(mapmtx[snp, 9],"|") # alleles
            if !ismissing(mapmtx[snp,8])
                multiallele = mapmtx[snp,8]
                if isa(multiallele, AbstractString) 
                    multiallele = parse(Int, multiallele)
                end
                if length(alleles) != multiallele 
                    @warn string("inconsisent between multiallele=", multiallele, " and alleles=",alleles) maxlog=10
                end
            end
            ref = popfirst!(alleles)
            alt = join(alleles,",")
            vcfmapmtx[snp,4] = ref 
            vcfmapmtx[snp,5] = alt == "NA" ? "." : alt
        end
        info = ismissing(mapmtx[snp,12]) ? "." : mapmtx[snp,12]
        linkagegroup = mapmtx[snp,2]
        poscm = mapmtx[snp,3]
        errls = mapmtx[snp, _col_1sterr:_col_1stsample-1]        
        vcfmapmtx[snp,8] = combine_2_vcfinfo(info, linkagegroup,poscm, errls)        
        #  set vcf 9th column: "FORMAT"
        if isnofounder && isnooff
            vcfmapmtx[snp,9] = "."
        else
            # col 6: parentformat, col 7: offspringformat
            offspringformat = sort(unique([i[1:2] for i in split(mapmtx[snp,7],"&")]))
            length(offspringformat) > 1 && @warn string("multiple offspringformat=",offspringformat) maxlog=10    
            if isnofounder
                # col 7: offspringformat    
                vcfmapmtx[snp,9] = join(offspringformat,":")
            else
                parentformat = sort(unique([i[1:2] for i in split(mapmtx[snp,6],"&")]))
                genoformat = isnooff ? parentformat : sort(unique(vcat(parentformat,offspringformat)))        
                length(genoformat) == 2 && genoformat[2] == "GT" && reverse!(genoformat) 
                vcfmapmtx[snp,9] = join(genoformat,":")
                # reset genodata
                nformat = length(genoformat)
                if nformat > 1 
                    if !isnooff                        
                        formatindex = findfirst(==(only(offspringformat)),genoformat)
                        if isnothing(calledoffspringmtx)
                            offgenols = view(offspringmtx, snp,:)
                            for j in eachindex(offgenols)
                                offgenols[j] = join([i==formatindex ? offgenols[j] : "." for i in 1:nformat],":") 
                            end
                        else
                            offgenols = view(offspringmtx, snp,:)
                            calledoffgenols = view(calledoffspringmtx, snp,:)
                            for j in eachindex(offgenols,calledoffgenols)
                                offgenols[j] = join([if i==formatindex 
                                    offgenols[j] 
                                elseif genoformat[i] == "GT" 
                                    calledoffgenols[j] 
                                else
                                    "." 
                                end for i in 1:nformat],":") 
                            end
                        end                        
                    end
                    parentgenols = view(parentmtx, snp,:)
                    if length(genoformat) == 2 
                        if genoformat[1] == "GT" 
                            parentgenols .= [(occursin("|",i) || occursin("/",i)) ? string(i,":.") : string(".:",i) for i in parentgenols]
                        elseif genoformat[2] == "GT" 
                            parentgenols .= [(occursin("|",i) || occursin("/",i)) ? string(".:",i) : string(i,":.") for i in parentgenols]
                        elseif genoformat[1] == "AD" && in(genoformat[2],["GP","GL"])
                            parentgenols .= [(occursin(",",i) && occursin(".",i)) ? string(".:",i) : string(i,":.") for i in parentgenols]
                        else
                            @error string("TODO for genoformat=",genoformat)
                        end
                    else
                        @error string("TODO for genoformat=",genoformat)
                    end

                end
            end
        end
    end
    # VCFv4.3 requires that REF must be one of A, T, G, C, N
    vcfmapmtx[vcfmapmtx[:,4] .== ".",4] .= "N" 
    vcfmapmtx
end


function combine_2_vcfinfo(info::AbstractString, linkagegroup, poscm, errls::AbstractVector)    
    oldinfo = ""
    if !ismissing(linkagegroup)
        oldinfo *= string(";LINKAGEGROUP=",linkagegroup)
    end
    if !ismissing(poscm)
        oldinfo *= string(";POSCM=",poscm)
    end
    errnamels = ["DOSEERROR", "BASEERROR", "ALLELICBIAS","OVERDISPERSION"]
    for j in eachindex(errls,errnamels)
        if !ismissing(errls[j])
            oldinfo *= string(";", errnamels[j], "=", errls[j])
        end
    end
    oldinfo =lstrip(oldinfo, ';')
    if isempty(oldinfo)
        info
    else
        info == "." ? oldinfo : string(oldinfo,";", info)
    end    
end

function rawDoseCall!(polygeno::PolyGeno; baseerror::Real=0.001)
    errdigit0 = split(last(split(@sprintf("%.f",baseerror),".")),"")
    errdigit2=findfirst(x->x!="0",errdigit0)
    errdigit = isnothing(errdigit2) ? 6 : max(4,errdigit2+1)
    if kindofgeno(polygeno.parentgeno) == "readcount"
        ploidyls = polygeno.polyped.parentinfo[!,:ploidy]
        polygeno.parentgeno = [calpostprob_AD(i,ploidyls,baseerror,digits=errdigit) for i=polygeno.parentgeno]
    end
    if kindofgeno(polygeno.offspringgeno) == "readcount"
        ploidyls = polygeno.polyped.offspringinfo[!,:ploidy]
        polygeno.offspringgeno = [calpostprob_AD(i,ploidyls,baseerror,digits=errdigit+2) for i=polygeno.offspringgeno]
    end
    polygeno
end

function getsubpolygeno!(polygeno::PolyGeno;
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

function getsubpolygeno(polygeno::PolyGeno;
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
    PolyGeno(deepcopy(polygeno.polyped), markermap,parentgeno,offgeno,deepcopy(polygeno.misc))
end

function filterdepth_AD!(polygeno::PolyGeno; 
    mindepth_parent::Integer=10,
    mindepth_offspring::Integer=10,
    maxmiss_parent::Real=0.5,
    maxmiss_offspring::Real=0.5)
    nchr = length(polygeno.markermap)
    for chr in 1:nchr
        # parents
        chrdepths = [ismissing(i) ? 0 : sum(i) for i in polygeno.parentgeno[chr]];
        sumdepthls = sum(chrdepths,dims=2)[:,1] 
        nind = size(chrdepths,2)        
        bls_parent = (sumdepthls ./ nind) .>= mindepth_parent
        bls_parent_miss = mean(sumdepthls .== 0) .<= maxmiss_parent
        # offspring
        chrdepths = [ismissing(i) ? 0 : sum(i) for i in polygeno.offspringgeno[chr]]
        sumdepthls = sum(chrdepths,dims=2)[:,1] 
        nind = size(chrdepths,2)        
        bls_off = (sumdepthls ./ nind) .>= mindepth_offspring
        bls_off_miss = mean(sumdepthls .== 0) .<= maxmiss_offspring
        bls = bls_parent .&& bls_off .&& bls_parent_miss .&& bls_off_miss
        polygeno.markermap[chr] = polygeno.markermap[chr][bls,:]
        polygeno.parentgeno[chr] = polygeno.parentgeno[chr][bls,:]
        polygeno.offspringgeno[chr] = polygeno.offspringgeno[chr][bls,:]
    end
    polygeno
end

function filtermissing_GT!(polygeno::PolyGeno;     
    maxmiss_parent::Real=0.5,
    maxmiss_offspring::Real=0.5)
    nchr = length(polygeno.markermap)
    for chr in 1:nchr
        # parents
        chrmissing = [ismissing(i) ? 1 : occursin(".",i) ? 1 : 0 for i in polygeno.parentgeno[chr]];
        nind = size(chrmissing,2)        
        missls = sum(chrmissing,dims=2)[:,1] ./ nind
        bls_parent_miss = missls .<= maxmiss_parent
        # offspring
        chrmissing = [ismissing(i) ? 1 : occursin(".",i) ? 1 : 0 for i in polygeno.offspringgeno[chr]];
        nind = size(chrmissing,2)        
        missls = sum(chrmissing,dims=2)[:,1] ./ nind
        bls_off_miss = missls .<= maxmiss_offspring
        # 
        bls = bls_parent_miss .&& bls_off_miss
        polygeno.markermap[chr] = polygeno.markermap[chr][bls,:]
        polygeno.parentgeno[chr] = polygeno.parentgeno[chr][bls,:]
        polygeno.offspringgeno[chr] = polygeno.offspringgeno[chr][bls,:]
    end
    polygeno
end

function delmonomarker!(polygeno::PolyGeno)
    isparentphased = all([unique(i[!,:parentformat]) == ["GT_phased"] for i in polygeno.markermap])
    if !isparentphased 
        @warn "some parental genotypes are not phased; filtermono! works only when all parent are phased"
        return polygeno
    end
    for chr in eachindex(polygeno.markermap)
        # parents
        ispolyls = [length(unique(reduce(vcat,skipmissing(i)))) > 1 for i  in eachrow(polygeno.parentgeno[chr])]
        if !any(ispolyls)
            @warn "no polymorphic markers!"
            @info "" chr polygeno.parentgeno[chr]
        end
        polygeno.markermap[chr] = polygeno.markermap[chr][ispolyls,:]
        polygeno.parentgeno[chr] = polygeno.parentgeno[chr][ispolyls,:]        
        if !isnothing(polygeno.offspringgeno)
            polygeno.offspringgeno[chr] = polygeno.offspringgeno[chr][ispolyls,:]
        end
    end
    polygeno
end

function readpolygeno(polygenofile::AbstractString;
    delim::AbstractChar=',',    
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    polygenofile2=getabsfile(workdir,polygenofile)    
    # issparse = true, CSV.read is used to parse string, when cells with too long string will be truncated
    res = readmultitable(polygenofile2; delim,commentstring,isparse=true)    
    haskey(res, "genodata") || error(string("key genodata is not in ", keys(res)))
    genodf = res["genodata"]            
    haskey(res, "pedigree") || error(string("key genodata is not in ", keys(res)))  
    ped = Pedigree(res["pedigree"])        
    polyped = PolyPed(ped)
    polygeno = formpolygeno!(genodf,polyped; isphysmap=false)
    misc = Dict{String, DataFrame}()
    for (key,val) in res
        if key[1:5]=="misc|"                         
            push!(misc, key[6:end]=>val)            
        end
    end
    polygeno.misc = misc
    polygeno
end


"""
    savepolygeno(outfile,polygeno,workdir=pwd(),delim=',')

save genotypic data of polygeno::PolyGeno into outfile

# Positional arguments

`outfile::AbstractString`: filename for saving results.

`polygeno::PolyGeno`: a struct returned by [`formpolygeno`](@ref).

# Keyward arguments

`missingstring::AbstractString = "NA"`: string for missing values
`workdir::AbstractString = pwd()`: directory for writing outfile.

`delim::AbstractChar=','`: delimitor character

"""
function savepolygeno(sink::Union{IO,AbstractString},polygeno::PolyGeno;
    missingstring::AbstractString = "NA",    
    workdir::AbstractString = pwd(),
    delim::AbstractChar = ',')
    if typeof(sink) <: AbstractString
        ext = last(splitext(sink))
        myopen = ext == ".gz" ? GZip.open : open
        outputfile2 = getabsfile(workdir,sink)
        io = myopen(outputfile2, "w")
    elseif typeof(sink) <: IO
        io = sink
    else
        error(string("unknown sink: ",sink))
    end
    initial = "PolyOrigin"
    write(io,initial,delim,"pedigree","\n")
    savepolyped(io,polygeno.polyped; workdir,delim)
    write(io,initial,delim,"genodata","\n")
    savegenodata(io,polygeno; keepcomment=false,missingstring, workdir,delim)
    misc = polygeno.misc
    if !isempty(misc)
        for (key,val) in misc
            key == "filecomment" && continue         
            write(io,initial,delim,string("misc|",key),"\n")
            CSV.write(io,val; header=true,append=true,delim)    
        end
    end
    flush(io)
    typeof(sink) <: AbstractString && close(io)
    sink
end



function fill_parent_missing!(polygeno::Union{PolyGeno, PolyAncestry})
    for chr in eachindex(polygeno.markermap)
        fill_parent_missing!(polygeno,chr)
    end
    nothing
end

function fill_parent_missing!(polygeno::Union{PolyGeno, PolyAncestry},chr::Integer)    
    ploidyls = polygeno.polyped.parentinfo[!,:ploidy]
    chrpgeno = polygeno.parentgeno[chr]
    # chrmultiallele = polygeno.markermap[chr][!,:multiallele]
    # ["GT_dosage","GT","GT_unphased","GT_phased","AD","GL","GP","AD&GT","GP&GT","GL&GT"]
    chrparentformat = polygeno.markermap[chr][!,:parentformat]
    formatls = ["GT", "GT_unphased","GT_phased","AD","AD&GT","GP&GT","GL&GT"]    
    intersect!(formatls, unique(chrparentformat))
    for p in eachindex(ploidyls)
        ploidy = ploidyls[p]
        for format in formatls
            b = chrparentformat .== format .&& ismissing.(chrpgeno[:, p])
            any(b) || continue
            if in(format, ["GT_unphased","GT","AD&GT","GP&GT","GL&GT"])
                misscode = join(repeat(["N"],ploidy),"/")
                chrpgeno[b, p] .= misscode
            elseif format == "GT_phased"
                misscode = join(repeat(["N"],ploidy),"|")                
                chrpgeno[b, p] .= [misscode for _ in 1:sum(b)]            
            else
                0
            end
        end        
    end
    nothing
end



function info_polygeno(polygeno::PolyGeno;
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    genolen = sum([i[end,:poscm] - i[1,:poscm] for i in polygeno.markermap])
    physlen = sum([i[end,:physposbp] - i[1,:physposbp] for i in polygeno.markermap])
    isfounder = polygeno.polyped.parentinfo[!,:isfounder]
    msg = string("#pop=", size(polygeno.polyped.designinfo,1))
    if all(isfounder)
        msg *= string(", #parent=",sum(isfounder))
    else
        msg *= string(", #founder=",sum(isfounder),", #midparent=",sum(.!isfounder))
    end
    msg *= string(", #offspring=",size(polygeno.polyped.offspringinfo,1),                
            ", #chr=",length(polygeno.markermap),        
            ", #marker=",sum(size.(polygeno.markermap,1)))
    ismissing(genolen) || (msg *= string(", genolen=", round(genolen, digits=1), "cM"))
    ismissing(physlen) || (msg *= string(", physlen=", round(1e-6 * physlen, digits=1), "Mbp"))
    printconsole(io,verbose,msg)
end

function get_genoformat(polygeno::PolyGeno; isparentformat::Bool)
    col = isparentformat ? :parentformat : :offspringformat
    unique(reduce(vcat, [unique(i[!, col]) for i in polygeno.markermap]))
end

function rename_samples(genofile::AbstractString,pedfile::AbstractString; 
    issavefhaplo::Bool=false,   
    workdir::AbstractString=pwd(), outstem::AbstractString="outstem")
    polygeno = formpolygeno(genofile,pedfile)
    ped2 = renameped(polygeno.polyped.ped)
    polygeno.polyped = PolyPed(ped2)    
    if issavefhaplo
        genofile2 = abspath(workdir, outstem*"_rename_fhaplo.vcf.gz")
    else
        genofile2 = abspath(workdir, outstem*"_rename_geno.vcf.gz")
    end    
    pedfile2 = abspath(workdir, outstem*"_rename_ped.csv")    
    saveped(pedfile2, polygeno.polyped.ped)
    savegenodata(genofile2, polygeno; target = issavefhaplo ? "founder" : "all")
    genofile2,pedfile2
end




function disturbmap!(polygeno::PolyGeno; disturbjump::Integer=0, disturbstd::Real=0)
    disturbjump == 0 && disturbstd ≈ 0.0 && return polygeno
    for chr = eachindex(polygeno.markermap)
        df = polygeno.markermap[chr]
        nmarker = size(df, 1)
        bef = collect(1:nmarker)
        if disturbstd ≈ 0.0
            aft = bef
        else
            noise = rand(Normal(0,disturbstd),nmarker)
            aft = sortperm(bef .+ noise)
        end
        if disturbjump > 0             
            nhalf = div(nmarker,2)
            nhalf <= disturbjump && @warn string("disturbjump=",disturbjump, ", greater than nhalf=",nhalf, " in chr=",chr)
            disturbjump2 = nhalf <= disturbjump ? nhalf : disturbjump
            leftls = sample(1:nhalf,disturbjump2; replace=false)
            rightls = sample((nhalf+1):nmarker,disturbjump2; replace=false)
            leftls2 = copy(aft[leftls])
            rightls2 = copy(aft[rightls])
            aft[leftls] .= rightls2
            aft[rightls] .= leftls2
        end
        # change markermap
        poscmls = copy(df[!,:poscm])
        df2 = df[aft,:]
        df2[!,:poscm] .= poscmls
        polygeno.markermap[chr] = df2
        # change geno
        polygeno.parentgeno[chr] .= polygeno.parentgeno[chr][aft, :]
        polygeno.offspringgeno[chr] .= polygeno.offspringgeno[chr][aft, :]
    end
    polygeno
end