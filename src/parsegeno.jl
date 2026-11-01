
function readgenodf(genofile::AbstractString,pedfile::AbstractString;    
    formatpriority::AbstractVector=["AD","GT"],
    delmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    logfile::Union{Nothing, AbstractString,IO} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=false)
    ext = last(split_allext(genofile))
    ext in [".csv",".csv.gz", ".vcf",".vcf.gz"] || @error string("genofile ext must be .csv, .csv.gz, .vcf or .vcf.gz")
    genofile2=getabsfile(workdir,genofile)
    isfile(genofile2) || @error string(genofile2, " does not exit")
    isvcf = ext in [".vcf",".vcf.gz"]    
    lastcomment = findlastcomment(genofile2; commentstring)
    if isvcf
        lastcomment == 0 && @warn string("no comment lines for vcfile=",genofile,", commentstring=",commentstring)
        lastcomment -= 1 # vcf last comment line: #CHROM, ..., is the headrow of dataframe
    end
    commentlines = read_headlines(genofile2, lastcomment)    
    outstem = tempname(workdir,cleanup=false)
    outext = ".csv.gz"
    tempgenofile = outstem*outext
    try 
        if isvcf            
            parsevcf(genofile2, pedfile; formatpriority,delmultiallelic,
                commentstring, outstem,outext, workdir, logfile, verbose
            ) # "NA" as misingstring in output of parsevcf            
        end        
        if isvcf || ext == ".csv.gz"
            genofile3 = isvcf ? tempgenofile : genofile2            
            genodf = readdataframe(genofile3;delim=',',missingstring,commentstring)
        else
            # ext == ".csv"
            genodf = CSV.read(genofile2,DataFrame;delim=',',comment = commentstring,missingstring)
        end
        rename!(genodf, string.(strip.(names(genodf))))
        genodf,commentlines
    finally
        isvcf && rm(tempgenofile; force=true)        
    end
end

function splitgenodf(genodf::DataFrame,
    founderid::AbstractVector,
    parentid::AbstractVector,
    offspringid::AbstractVector;
    isphysmap::Bool=false,
    recomrate::Real=1.0)    
    colnames = names(genodf)
    coldict = Dict(colnames .=> eachindex(colnames))    
    mincol_sample = minimum([haskey(coldict,i) ? coldict[i] : length(colnames) + 1 for i in vcat(founderid, parentid,offspringid)])
    # @info string("1st column for samples: ", mincol_sample)
    if mincol_sample < _col_1stsample
        # make it compatitable with the input csv genofile for the ealier version of polyorigin
        mapcols = ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat",
            "multiallele", "alleles", "quality", "filter", "info","doseerror","baseerror","allelicbias", "overdispersion","doublereduction"]
        mapcolnames = replace(lowercase.(colnames[1:mincol_sample-1]), 
            "markers" => "marker", 
            "id" => "marker", 
            "chromosome" => "linkagegroup", "chrom" => "linkagegroup", "chr" => "linkagegroup", 
            "position"=>"poscm","pos"=>"poscm", "cm" => "poscm"
        ) 
        setdiff!(mapcolnames,mapcols)
        if !isempty(mapcolnames)
            @error string("unknown column names = ", mapcolnames, " before the 1st sample")
        end
        inferredcol_1stsample = mincol_sample
    else
        inferredcol_1stsample = _col_1stsample
    end
    sampleid = colnames[inferredcol_1stsample:end]
    # println("names(genodf)=",names(genodf), ";parentid=",parentid,",offspringid=",offspringid)    
    d = setdiff(sampleid,founderid, parentid,offspringid)
    if !isempty(d)
        msg = string("genotyped individuals in genofile but not in pedfile: ",d)
        @warn msg
    end
    # mapdf
    mapdf = genodf[!,1:inferredcol_1stsample-1]
    if size(mapdf,2) == 3
        if isphysmap
            insertcols!(mapdf, 2, "linkagegroup"=>missing)        
            insertcols!(mapdf, 3, "poscm"=>missing)        
        else
            insertcols!(mapdf, size(mapdf,2)+1, "physchrom"=>missing)        
            insertcols!(mapdf, size(mapdf,2)+1, "physposbp"=>missing)                
        end
    end        
    parsemarkermap!(mapdf) # mapdf is expanded to have #cols == _col_1stsample    
    iscmmiss= ismissing.(mapdf[!,:poscm])
    if isphysmap
        b = ismissing.(mapdf[!,:physposbp])
        if all(b)
            msg = string("physical map is missing")
            @warn msg
        else
            if any(b)
                msg = string(sum(b), " out of ", length(iscmmiss), " markers' physical positions are missing")
                @warn msg
            end
            maxloc = maximum(skipmissing(mapdf[!,:physposbp]))
            if maxloc<10^6.0
                @warn string("marker positions in physical map must be in base pair")
            end
            all(iscmmiss) || @warn string("genetic positions are overwritten by transforming physical position using recomrate =",recomrate," cM/Mbp")
            mapdf[!,:poscm] .= [ismissing(i) ? i : i*recomrate*10^(-6.0) for i in mapdf[!,:physposbp]]
            mapdf[!,:linkagegroup] .= mapdf[!,:physchrom]            
        end
    else
        if all(iscmmiss)            
            msg = string("genetic map is missing")
            @info msg
        else
            if any(iscmmiss)
                msg = string(sum(iscmmiss), " out of ", length(iscmmiss), " markers' genetic positions are missing")
                @warn msg
            end
            maxloc = maximum(skipmissing(mapdf[!,:poscm]))
            if maxloc>10^6.0
                msg = "marker positions in genetic map must be in centiMorgan. "
                if eltype(mapdf[!,:poscm]) <: Union{Missing, Integer}
                    msg *= "Try to set isphysmap = true. "
                end
                @warn msg 
            end
        end        
    end
    coldict = Dict(names(genodf)[inferredcol_1stsample:end] .=> inferredcol_1stsample:size(genodf,2))
    # fgenodf    
    newfounderid = copy(founderid)
    if isempty(founderid)
        fgenodf = nothing
    else
        cols = [get(coldict,i,nothing) for i=founderid]
        missind = founderid[isnothing.(cols)]
        if !isempty(missind)
            msg = string(length(missind), " founders in pedfile but not in genofile: ",missind)
            @warn msg
            b = .!isnothing.(cols)
            cols = cols[b]
            newfounderid = founderid[b]
        end
        fgenodf = genodf[!, cols]
    end

    # pgenodf    
    newparentid = copy(parentid)
    if isempty(parentid)
        pgenodf = nothing
    else
        cols = [get(coldict,i,nothing) for i=parentid]
        missind = parentid[isnothing.(cols)]
        if !isempty(missind)
            # msg = string(length(missind), " parents in pedfile but not in genofile: ",missind)
            # @warn msg
            b = .!isnothing.(cols)
            cols = cols[b]
            newparentid = parentid[b]
        end
        pgenodf = genodf[!, cols]
    end

    # offgenodf
    newoffspringid = copy(offspringid)
    if isempty(offspringid)
        offgenodf = nothing
    else
        cols = [get(coldict,i,nothing) for i=offspringid]
        missind = offspringid[isnothing.(cols)]
        if !isempty(missind)
            # msg = string(length(missind), " offspring in pedfile but not in genofile: ",missind)
            # @warn msg
            b = .!isnothing.(cols)
            cols = cols[b]
            newoffspringid = offspringid[b]
        end
        offgenodf = genodf[!, cols]
    end    
    mapdf, fgenodf, newfounderid, pgenodf, newparentid, offgenodf,newoffspringid
end

function parsemarkermap!(markermap::DataFrame)    
    mapcols = ["marker", "linkagegroup","poscm","physchrom", "physposbp", "parentformat","offspringformat",
        "multiallele", "alleles", "quality", "filter", "info", "doseerror","baseerror","allelicbias","overdispersion"]
    nmapcol = length(mapcols)
    ncol = size(markermap,2)
    if ncol < 3
        @warn string("number of columns from markermap = ", ncol, ",< 3")
    elseif ncol < nmapcol
        for i in ncol+1:nmapcol
            insertcols!(markermap, i, mapcols[i]=>missing)        
        end
        markermap[!,:multiallele] .= 2
    elseif (ncol == nmapcol+1) && in("doublereduction", names(markermap))
        push!(mapcols, "doublereduction")        
    elseif ncol > nmapcol
        @warn string("number of markermap columns = ", ncol,  " > 12")
        select!(markermap, names(markermap)[1:nmapcol])
    end
    rename!(markermap,mapcols)
    markerls = Vector{Any}(markermap[!,:marker])
    b = ismissing.(markerls)
    if any(b)
        markerls[b].= string.("snp",findall(b))
        markermap[!,:marker] .= string.(markerls)
    end            
    allunique(markermap[!,:marker]) || @error ("marker IDs are not unique")
    for j = [:marker, :linkagegroup]        
        markermap[!,j] .= [ismissing(i) ? i : string(strip(string(i))) for i in markermap[!,j]]
    end
    b = .!ismissing.(markermap[!,:linkagegroup])
    if any(b)
        snpindexlist=splitindex(markermap[b,:linkagegroup])
        chridls = [markermap[first(i),:linkagegroup] for i=snpindexlist]
        if !allunique(chridls)
            @warn ("markers in the same chromosome are not consecutive")
        end
    else
        snpindexlist = [1:size(markermap,1)]
        chridls = markermap[1:1,:linkagegroup]
    end
    if !any(ismissing.(markermap[!,:poscm]))
        poscmls = [isa(i,AbstractString) ? tryparse(Float64,i) : i for i in markermap[!,:poscm]]
        b = isnothing.(poscmls)
        if any(b)
            msg = string("delete ", sum(b), " markers with poscm: ", unique(markermap[b,:poscm]))
            @warn msg                        
            deleteat!(markermap, b)
            markermap[!,:poscm] = Vector{Float64}(poscmls[.!b])
            snpindexlist=splitindex(markermap[!,:linkagegroup])
        end
        for chr in eachindex(snpindexlist)
            ii = snpindexlist[chr]
            A=markermap[ii,:poscm]
            b = A[1:end-1] .<= A[2:end]
            if !all(b)
                posls = findall(.!b)
                union!(posls,posls .+ 1)
                sort!(posls)
                ncol = min(5, size(markermap,2))
                msg = string("Non-descreasing genetic positions in linkagegroup=",chridls[chr], "\n", markermap[ii[posls],1:ncol])
                @error msg
            end
        end        
    end    
    markermap
end

function get_kindofgeno(genodatadf::AbstractDataFrame, ploidyset::AbstractVector)
    coltypes = unique([eltype(genodatadf[!,i]) for i in names(genodatadf)])
    missingcodes = ["NA", "missing",""]
    if all(coltypes .<: Union{Missing, Integer})
        "GT_dosage"
    elseif all(coltypes .<: Union{Missing, AbstractString})        
        ls = strip.(unique(Matrix(genodatadf)))
        b = [ismissing(i) || in(i, missingcodes) for i in ls]
        deleteat!(ls, b)
        isempty(ls) && return "GT_unphased" # default when all data  are missing         
        if all([occursin("/",i) for i = ls])
            "GT_unphased"
        elseif all([occursin("&",i) for i = ls])
            # AD, GP, GL
            ls1 = split.(ls,"&")
            lenls = unique(length.(ls1)) 
            maximum(lenls) >= 2 || @error string("unexpected dimensions for genodata = ",ls)  # number of alleles >= 2
            if maximum(ploidyset) > 2 && maximum(lenls) == 2                 
                "AD"   # biallelic markers for ploidy > 2
            else
                ls2 = unique(reduce(vcat,ls1))
                filter!(x->!in(x,missingcodes),ls2)
                ls3 = tryparse.(Int, ls2)
                if in(nothing, ls3) 
                    # GP or GL
                    if length(ploidyset) == 1 && (only(lenls) != only(ploidyset) + 1)
                        @error string("unexpected dimensions for genodata = ",ls[1:min(100,length(ls))])  
                    end                    
                    parse_GP_GL(ls2)                                        
                else
                    maximum(ls3) < 0 && error(string("Cound not parse genodata with negative elements: ", ls3))                        
                    "AD"
                end
            end
        elseif all([occursin("|",i) for i = ls])
            # GT_phased, AD, GP, GL
            ls1 = split.(ls,"|")
            lenls = unique(length.(ls1))            
            maximum(lenls) >= 2 || @error string("unexpected dimensions for genodata = ",ls)  # number of alleles >= 2            
            ls2 = unique(reduce(vcat,ls1))
            filter!(x->!in(x,missingcodes),ls2)
            ls3 = tryparse.(Int, ls2)
            if in(nothing, ls3) 
                # GP or GL                  
                if length(ploidyset) == 1 && (only(lenls) != only(ploidyset) + 1)
                    @error string("unexpected dimensions for genodata = ",ls[1:min(100,length(ls))])  
                end                                        
                parse_GP_GL(ls2)  
            else
                if maximum(ploidyset) > 2 && maximum(lenls) == 2 
                    "AD"   # biallelic markers for ploidy > 2
                elseif ploidyset == [2] && maximum(lenls) > 2 
                    "AD" # multallelic markers for diploid
                else                    
                    # GT_phased or AD; ls3 is a list of different alleles for GT_phased; ls3 is a list of read counts for AD                                        
                    maxls = maximum(ls3) 
                    maximum(ls3) < 0 && error(string("Cound not parse genodata with negative elements: ", ls3))                        
                    if maxls <= 2
                        if maxls == 2
                            @warn string("check if alleles = ",ls3, ". In case of allelic depths, transform into a VCF-format genofile")
                        end
                        "GT_phased"                    
                    else 
                        if length(ploidyset) == 1 && length(lenls) > 1 
                            "AD" # multiallelic markers
                        elseif issetequal(ls3, 0:length(ls3)-1)        
                            if maximum(ls3) > 50
                                @warn string("check if read counts = ",ls3, ". In case of phased genotypes, transform into a VCF-format genofile")
                                "AD" # multiallelic markers
                            else
                                @warn string("check if multialleles = ",ls3, ". In case of allelic depths, transform into a VCF-format genofile")
                                "GT_phased"
                            end   
                        else
                            @warn string("check if read counts = ",ls3, ". In case of phased genotypes, transform into a VCF-format genofile")
                            "AD"
                        end                        
                    end                    
                end
            end            
        else
            error(string("Cound not parse genodata = ", ls[1:min(100,length(ls))]))
        end    
    else
        error(string("Cound not parse genodata with column types = ", coltypes))
    end
end

# ls is a list of sring probablities (GP) or loglike(GL)
function parse_GP_GL(ls::AbstractVector)
    # GP or GL
    ls4 = tryparse.(Float32, ls)
    in(nothing, ls4) && error(string("Cound not parse genodata = ", ls4))                    
    if all(ls4 .<= 0)
        "GL"
    elseif all(0 .<= ls4 .<= 1)
        "GP"
    else
        error(string("Cound not parse genodata = ", ls4))                        
    end    
end
   

function get_dose2unphased(dose::Union{Missing, Integer}, ploidy)
    ismissing(dose) && return missing
    join(vcat(["0" for _ in 1:ploidy-dose], ["1" for _ in 1:dose]),"/")
end


function set_format_dosage!(mapdf::DataFrame, 
    fgenodf::Union{Nothing,DataFrame}, 
    pgenodf::Union{Nothing,DataFrame}, 
    offgenodf::Union{Nothing,DataFrame},polyped::PolyPed)    
    errmsg = string("dosage format does not apply to multiallelic markers")
    # founders and intermeiate parents
    isfounder = polyped.parentinfo[!,:isfounder]
    col = :parentformat
    mapdf[!,col] .= [if ismissing(i) 
            Vector{Any}([missing, missing]) 
        else
            ls = String.(split(i,"|")) 
            length(ls) == 1 ? [ls[1], ls[1]] : ls
        end for i in mapdf[!,col]]
    # founders
    ploidyls = polyped.parentinfo[isfounder,:ploidy]
    nind = length(ploidyls)
    nind == size(fgenodf, 2) || error(string("inconsistent #founders=",nind))
    formatls = Vector{Any}(first.(mapdf[!,col]))
    formatset = unique(formatls)
    for format in formatset
        # format and formatls will be updated if format is missing 
        if ismissing(format)
            mls = findall(ismissing.(formatls))
            format = get_kindofgeno(fgenodf[mls,:], unique(ploidyls)) 
            for m in mls
                mapdf[m, col][1] = format
                formatls[m] = format
            end                
        end
        if format == "GT_dosage"                        
            mls = findall(formatls .== format)  
            multialleleset = unique(mapdf[mls,:multiallele]) 
            multialleleset == [2] || error(string(errmsg,". multialleleset = ", multialleleset))           
            for j in eachindex(ploidyls)
                ls = Vector{Any}(fgenodf[!,j])
                ls[mls] .= get_dose2unphased.(fgenodf[mls,j],ploidyls[j])
                fgenodf[!,j] .= ls
            end
            for m in mls
                mapdf[m, col][1] = "GT_unphased"
            end
        end
    end
    # change 1 based alleli coding into 0 based coding: input/output files from polyorigin v1
    formatls = first.(mapdf[!,col])    
    markls = findall(formatls .== "GT_phased")    
    if !isempty(markls)
        fgenomtx = Matrix(fgenodf[markls,:])
        alleleset = sort(unique(reduce(vcat,split.(reduce(vcat,unique(fgenomtx)),"|"))))
        if alleleset == ["1","2"]
            fgenodf[markls,:] .= replace.(fgenomtx,"1"=>"0","2"=>"1")
        end
    end

    # intermeiate parents
    if isnothing(pgenodf)
        mapdf[!,col] .= first.(mapdf[!,col])
    else
        ploidyls = polyped.parentinfo[.!isfounder,:ploidy]
        nind = length(ploidyls)
        nind == size(pgenodf, 2) || error(string("inconsistent #intermidparent parents=",nind))
        formatls = Vector{Any}(last.(mapdf[!,col]))
        formatset = unique(formatls)
        for format in formatset
            if ismissing(format)
                mls = findall(ismissing.(formatls))
                format = get_kindofgeno(pgenodf[mls,:], unique(ploidyls)) 
                for m in mls
                    mapdf[m, col][end] = format
                    formatls[m] = format
                end
            end
            if format == "GT_dosage"        
                mls = findall(formatls .=== format)  
                unique(mapdf[mls,:multiallele]) == [2] || @error errmsg            
                for j in eachindex(ploidyls)
                    ls = Vector{Any}(pgenodf[!,j])
                    ls[mls] .= get_dose2unphased.(pgenodf[mls,j],ploidyls[j])
                    pgenodf[!,j] .= ls
                end
                for m in mls
                    mapdf[m, col][end] = "GT_unphased"
                end
            end
        end      
        mapdf[!,col] .= [allequal(i) ? i[1] : join(i,"|") for i in mapdf[!,col]]  
    end

    # offspring 
    ploidyls = polyped.offspringinfo[!,:ploidy]
    col = :offspringformat
    formatset = unique(mapdf[!,col])
    for format in formatset
        newformat = format
        if ismissing(format)
            b = ismissing.(mapdf[!,col])
            newformat = get_kindofgeno(offgenodf[b,:],unique(ploidyls)) 
            ls = Vector{Any}(mapdf[!,col])
            ls[b] .= newformat
            mapdf[!,col] .= ls
        end
        if newformat == "GT_dosage"      
            b = mapdf[!,col] .=== newformat  # missing === "GT_dosage" results in false            
            unique(mapdf[b,:multiallele]) == [2] || @error errmsg            
            for j in eachindex(ploidyls)
                ls = Vector{Any}(offgenodf[!,j])
                ls[b] .= get_dose2unphased.(offgenodf[b,j],ploidyls[j])
                offgenodf[!,j] .= ls
            end             
            mapdf[b,col] .= "GT_unphased"         
        end
    end
    nothing
end

function parsegenodf!(mapdf::DataFrame, 
    fgenodf::Union{Nothing,DataFrame}, 
    pgenodf::Union{Nothing,DataFrame}, 
    offgenodf::Union{Nothing,DataFrame};
    missingstring=["NA","."]) 
    # split into linagegroups    
    markerls =splitindex(mapdf[!,:linkagegroup])
    markermap = [mapdf[snps,:] for snps in markerls]    
    for df in markermap
        df[!,:parentformat] .= [occursin("|",i) ? String.(split(i, "|")) : [i,i] for i in df[!,:parentformat]]
    end    
    # founders
    if isnothing(fgenodf)
        fgeno = nothing
    else        
        fgeno = [Matrix{Any}(fgenodf[snps,:]) for snps in markerls];    
        empty!(fgenodf) 
        for chr in eachindex(markerls)
            chrgeno = fgeno[chr]
            b = [isa(i, AbstractString) && in(i,missingstring) for i in chrgeno]
            chrgeno[b] .= missing
            chrformatls = first.(markermap[chr][!,:parentformat])
            d = setdiff(unique(chrformatls),["GT_dosage","GT","GT_unphased","GT_phased","AD","GL","GP","AD&GT","GP&GT","GL&GT"])
            if !isempty(d)
                @info "test" unique(markermap[chr][!,:parentformat])
                @error string("unknown format ",d," in chr=",chr,". formatset=",unique(chrformatls))
            end
            isnonmiss = .!ismissing.(chrgeno);
            chrgeno[isnonmiss] .= string.(chrgeno[isnonmiss])                                    
            parsegenomtx!(chrgeno,chrformatls)            
            for i in eachindex(chrformatls)
                markermap[chr][i,:parentformat][1] = chrformatls[i]
            end
        end                 
    end    
    # intermeiate parents
    if isnothing(pgenodf) || isempty(pgenodf)
        pgeno = nothing
    else
        pgeno = [Matrix{Any}(pgenodf[snps,:]) for snps in markerls];    
        empty!(pgenodf) 
        for chr in eachindex(markerls)
            chrgeno = pgeno[chr]
            b = [isa(i, AbstractString) && in(i,missingstring) for i in chrgeno]
            chrgeno[b] .= missing            
            chrformatls = last.(markermap[chr][!,:parentformat])
            d = setdiff(unique(chrformatls),["GT_dosage","GT","GT_unphased","GT_phased","AD","GL","GP","AD&GT","GP&GT","GL&GT"])
            if !isempty(d)
                @info "test" unique(markermap[chr][!,:parentformat])
                @error string("unknown format ",d," in chr=",chr,". formatset=",unique(chrformatls))
            end
            
            isnonmiss = .!ismissing.(chrgeno);
            chrgeno[isnonmiss] .= string.(chrgeno[isnonmiss])                                    
            parsegenomtx!(chrgeno,chrformatls)                        
            for i in eachindex(chrformatls)
                markermap[chr][i,:parentformat][end] = chrformatls[i]
            end
        end                 
    end    
    for df in markermap
        df[!,:parentformat] .= [allequal(i) ? first(i) : join(i, "|") for i in df[!,:parentformat]]
    end    
    if isnothing(fgeno)
        parentgeno = pgeno
    elseif isnothing(pgeno)
        parentgeno = fgeno
    else
        parentgeno = [[fgeno[chr] pgeno[chr]] for chr in eachindex(fgeno,pgeno)]
    end
    # offspring
    if isnothing(offgenodf) || isempty(offgenodf)
        offspringgeno = nothing
    else         
        formatls = mapdf[:,:offspringformat]
        offspringgeno = [Matrix{Any}(offgenodf[snps,:]) for snps in markerls];     
        empty!(offgenodf)   
        GC.gc()
        for chr in eachindex(markerls)
            chroffgeno = offspringgeno[chr]
            b = [isa(i, AbstractString) && in(i,missingstring) for i in chroffgeno]
            chroffgeno[b] .= missing
            chrformatls = formatls[markerls[chr]]
            isnonmiss = .!ismissing.(chroffgeno);
            chroffgeno[isnonmiss] .= string.(chroffgeno[isnonmiss])
            parsegenomtx!(chroffgeno,chrformatls) # takes lots of time for GP                        
        end        
    end        
    empty!(mapdf)
    markermap, parentgeno, offspringgeno
end


function parsegenomtx!(genomtx::AbstractMatrix, formatvec::AbstractVector; 
    missingstring=[".", "NA"], geno_delimiter="&") 
    res = genomtx    
    formatset = unique(formatvec)
    if any(ismissing.(formatvec))
        println("there exist unknown genoformat: ",formatset)
    end
    for format = formatset
        ii = findall(formatvec .== format)
        subres =view(res, ii,:)
        nonmissls = findall(.!ismissing.(subres))            
         if format in ["AD&GT","GP&GT","GL&GT"]       
            for i in nonmissls
                g = subres[i]
                if occursin("|",g) 
                    delim= "|"                                         
                elseif occursin("/", g)
                    delim = "/"                    
                elseif occursin("&", g)
                    delim = "&"                    
                else
                    # delim = ""
                    error(string("format=",format, ",g=",g))
                end            
                ls = string.(strip.(split(g,delim)))
                b = [in(i,missingstring) for i in ls]
                if all(b)
                    subres[i] = missing
                else
                    if delim == "&"
                        if any(b) 
                            subres[i] = missing
                        else
                            if format == "AD&GT"
                                subres[i] = parse.(Int, ls)            
                            else
                                subres[i] = parse.(Float32, ls)
                            end
                        end
                    else
                        ls[b] .= "N"
                        subres[i] = delim == "|" ? ls : join(ls, delim)
                    end
                end
            end 
        elseif format in ["AD", "GL","GP"]
            type = format == "AD" ? Int : Float32                        
            for i in nonmissls
                if occursin(geno_delimiter,subres[i]) 
                    ls = string.(strip.(split(subres[i],geno_delimiter)))
                elseif occursin("|", subres[i])
                    ls = string.(strip.(split(subres[i],"|")))
                else
                    # delim = ""
                    error(string("format=",format, ",g=",subres[i]))
                end                
                subres[i] = isempty(intersect(ls, missingstring)) ? parse.(type, ls) : missing
            end            
        elseif format in ["GT","GT_unphased","GT_phased"]
            for i in nonmissls
                g = subres[i]
                if occursin("|",g) 
                    delim= "|"                                         
                elseif occursin("/", g)
                    delim = "/"                    
                else
                    # delim = ""
                    error(string("format=",format, ",g=",g))
                end            
                ls = string.(strip.(split(g,delim)))
                b = [in(i,missingstring) for i in ls]
                if all(b)
                    subres[i] = missing
                else
                    ls[b] .= "N"
                    subres[i] = delim == "|" ? ls : join(ls, delim)
                end    
            end        
        elseif format == "discretefgl"
            # no missing for the truevalue resulting from simulation
            # discretefgl must be stored in csv file
            subres .= [tryparse.(Int,i) for i=split.(subres,"|")]        
        else
            error(string("unknown geno format: ",format))
        end        
    end
    res
end


function inverse_parsegenomtx!(genomtx::AbstractMatrix, formatvec::AbstractVector,
    isvcf::Bool,isunphased_2dosage::Bool, 
    missingstring::AbstractString)
    formatset = unique(formatvec)
    res = genomtx
    if isvcf
        for format = formatset
            ii = formatvec .== format
            subres = view(res, ii,:)     
            if format in ["AD&GT","GP&GT","GL&GT"]       
                for i in eachindex(subres)
                    if ismissing(subres[i]) 
                        subres[i] =  "." 
                    elseif isa(subres[i], AbstractVector)
                        if eltype(subres[i]) <: AbstractString
                            subres[i] = join(replace(subres[i], "N"=>"."), "|")  # GT_phased                      
                        else
                            subres[i] = join(subres[i], ",")       # AD, GL, or GP                
                        end
                    elseif isa(subres[i], AbstractString)
                        subres[i] = replace(subres[i], "N"=>".") # GT_unphased
                    else 
                        @error string("unknonw geno=",subres[i]," with format=",format) maxlog=10
                    end
                end    
            elseif format in ["GT_unphased","GT_phased","GT"]                                      
                for i in eachindex(subres)
                    if ismissing(subres[i]) 
                        subres[i] =  "." 
                    elseif isa(subres[i], AbstractVector)
                        subres[i] = join(replace(subres[i], "N"=>"."), "|")                        
                    else
                        subres[i] = replace(subres[i], "N"=>".")
                    end
                end            
            elseif format in ["AD", "GL", "GP"]
                subdelimiter = ","
                b = ismissing.(subres)
                if any(b)                    
                    subres[b] .= "."
                    subres[.!b] .= join.(subres[.!b],subdelimiter)
                else
                    subres .= join.(subres,subdelimiter)
                end
            elseif format in ["discretefgl"]
                error(string(format," must be stored in csv file"))
            else
                error(string("unknown genoformat = ",format))
            end
        end
    else
        for format = formatset
            if format in ["GT", "GT_unphased", "GT_phased","AD", "GL", "GP","AD&GT","GP&GT","GL&GT", "discretefgl"]
                ii = formatvec .== format
                subres = view(res, ii,:)
                if format in ["GT", "AD&GT","GP&GT","GL&GT"]       
                    for i in eachindex(subres)
                        if ismissing(subres[i]) 
                            subres[i] =  missingstring
                        elseif isa(subres[i], AbstractVector)
                            if eltype(subres[i]) <: AbstractString
                                subres[i] = join(subres[i], "|")  # GT_phased                      
                            else
                                subres[i] = join(subres[i], "&")       # AD, GL, or GP                
                            end
                        elseif isa(subres[i], AbstractString) # GT_unphased
                            0
                            # do nothing 
                            # if isunphased_2dosage                        
                            #     subres[i] = occursin("N",subres[i]) ? missingstring : count('1',subres[i])                                
                            # end                            
                        else 
                            @error string("unknonw geno=",subres[i]," with format=",format) maxlog=10
                        end
                    end    
                elseif format in ["AD", "GL", "GP"]   
                    geno_delimiter = "&"                  
                    b = ismissing.(subres)
                    if any(b)                    
                        subres[b] .= missingstring                        
                        subres[.!b] .= join.(subres[.!b],geno_delimiter)
                    else
                        subres .= join.(subres,geno_delimiter)
                    end
                elseif format in ["GT_phased","discretefgl"]                    
                    subres .= join.(subres,"|") # no missing allowed
                elseif format in ["GT_unphased"]                    
                    if isunphased_2dosage                        
                        b = [ismissing(i) || occursin("N",i) for i in subres]                        
                        subres[b] .= missingstring
                        subres[.!b] .= count.('1',subres[.!b])
                        formatvec[ii] .= "GT_dosage"
                    else
                        b = ismissing.(subres)
                        subres[b] .= missingstring
                    end                                                 
                end
            else
                error(string("unknown genoformat = ",format))
            end
        end
    end
    res
end

# TODO: to update;  more testing for genodata.
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
            type = Float32
            geno4 = [[map(x->parse.(type,x),i) for i=split.(j,"|")] for j=geno2]
        else
            # missing is allowed in read count or phased geno
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


# function kindofgeno(geno::AbstractVector)
#     type = eltype(geno[1])
#     if type <: Union{Missing,Integer}
#         return "GT_dosage"
#     elseif type <: AbstractVector
#         subtype = eltype(type)
#         if subtype <:  Union{Missing,Integer}
#             bool = unique(length.(geno[1][1,:])) == [2]
#             bool2 = issubset(skipmissing(unique(vcat(geno[1]...))),[1,2])
#             return  (bool && !bool2) ? "readcount" : "phasedgeno"
#         elseif subtype <: AbstractFloat
#             return "probability"
#         else
#             error(string("unknown geno data type: ", type))
#         end
#     else
#         # @error string("unknown geno data type: ", type)
#         error(string("unknown geno data type: ", type))
#     end
# end