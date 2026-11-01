
function tocsvgeno(genofile::AbstractString,pedfile::AbstractString;
    outstem::AbstractString = first(split_allext(genofile)),
    outext::AbstractString = ".csv.gz",    
    formatpriority::AbstractVector=["AD","GT"],
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    last(split_allext(genofile)) in [".vcf",".vcf.gz"] || @error string("inpu genofile ext must .vcf or .vcf.gz")    
    parsevcf(genofile, pedfile;
            outstem, outext, formatpriority, 
            commentstring, workdir, logfile = nothing, verbose=false)
    outstem*outext
end

function tovcfgeno(genofile::AbstractString,pedfile::AbstractString;
    outstem::AbstractString = first(split_allext(genofile)),
    outext::AbstractString = ".vcf.gz",    
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    last(split_allext(genofile)) in [".csv",".csv.gz"] || @error string("input genofile ext must .csv or .csv.gz")
    magicgeno = formpolygeno(genofile, pedfile;
        isphysmap=false,
        commentstring, missingstring=["NA","missing"], workdir)
    outfile = outstem*outext    
    savegenodata(outfile,magicgeno;
        delim='\t',commentstring,workdir)
    outfile
end

function parsevcf(genofile::AbstractString,pedfile::AbstractString;    
    formatpriority::AbstractVector=["AD","GP", "GT"],
    delmultiallelic::Bool=true,    
    commentstring::AbstractString="##",
    outstem::Union{IO,AbstractString} = "outstem",
    outext::AbstractString = ".csv.gz",
    logfile::Union{Nothing, AbstractString,IO} = (isa(outstem, AbstractString) ? outstem : first(split_allext(genofile)))*"_parsevcf.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = set_logfile_begin(logfile, workdir, "parsevcf"; verbose,delim="-")
    inext = last(split_allext(genofile))
    inext in [".vcf",".vcf.gz"] || @error string("genofile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    msg = string("list of args: \n",
        "genofile = ", genofile, "\n",
        "pedfile = ", pedfile, "\n",        
        "formatpriority = ", formatpriority, "\n",
        "delmultiallelic = ", delmultiallelic, "\n",
        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",
        "outext = ", outext, "\n",
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    genofile2=getabsfile(workdir,genofile)
    lastcomment = findlastcomment(genofile2; commentstring)
    nheader = lastcomment + 1
    outext in [".csv",".csv.gz",".vcf",".vcf.gz"] || @error string("outfile ext must be in [.csv,.csv.gz,.vcf,.vcf.gz]")
    keepvcf = outext in [".vcf",".vcf.gz"]
    in_open(genofile2,"r") do inio
        if isa(outstem, AbstractString)
            out_open = outext in [".vcf.gz",".csv.gz"] ? GZip.open : open
            outfile = getabsfile(workdir,outstem*outext)
            out_open(outfile, "w") do outio
                parsevcf_io(inio, outio, logio, pedfile,
                    formatpriority,delmultiallelic, commentstring, nheader, keepvcf, workdir, verbose)
            end
            msg = string("save parsed genofile: ", outfile)
            printconsole(logio,verbose,msg)
        else
            # isa(outstem, IO)
            parsevcf_io(inio, outstem, logio, pedfile, 
                formatpriority,delmultiallelic,commentstring, nheader, keepvcf, workdir, verbose)
        end
    end
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, logio, tused,"parvcf"; verbose,delim="-")
    return
end

function parsevcf_io(inio::IO,outio::IO,logio::Union{Nothing, IO},
    pedfile::AbstractString, 
    formatpriority::AbstractVector,
    delmultiallelic::Bool,    
    commentstring::AbstractString,
    nheader::Integer, 
    keepvcf::Bool,    
    workdir::AbstractString, verbose::Bool)
    filedate = string(commentstring*"fileDate=",replace(string(Date(now())),"-"=>""))
    write(outio, filedate,"\n")
    write(outio, commentstring*"source= modified by polyorigin\n")    
    for _ in 1:nheader-1
        line = readline(inio,keep=true)                
        write(outio, line)
    end
    # parse title row
    titlerow = split(readline(inio,keep=false),"\t")
    ncolumn = length(titlerow)
    strip(titlerow[1]) == "#CHROM" || @error string("unexpected head row =",titlerow)
    _, fcols, pcols, offcols,newtitlerow = parse_titlerow(titlerow, pedfile;
        keepvcf, commentstring,workdir,logio)
    write(outio,newtitlerow,"\n")
    # parse rest rows
    msg = string("begin, #founders=",length(fcols), ", #midparents=", length(pcols), ", #offspring=",length(offcols))
    printconsole(logio,verbose,msg)
    fgeno = Vector(undef, length(fcols))
    pgeno = Vector(undef, length(pcols)) #  intermidparent parents
    offgeno = Vector(undef, length(offcols))
    keys = Symbol.(formatpriority)
    values = 1:length(formatpriority)
    prioritytuple = (; zip(keys, values)...)
    missingset = [".","./.", ".|.","./././.",".|.|.|."]    
    nmarker = 0
    nmultia = 0
    outdelim = keepvcf ? '\t' : ','
    startt = time()
    missingstring = "NA"
    while !eof(inio)
        # rem(nmarker, 1000) == 0 && (startt = time())
        nmarker += 1
        rowgeno = split(readline(inio,keep=false),"\t")        
        if length(rowgeno) <= 9
            # error(string("#columns < 9 for rowgeno=",rowgeno, " at marker index=", nmarker))
            continue
        end
        if length(rowgeno) != ncolumn
            error(string("marker index=", nmarker, ": #columns = ",  length(rowgeno),  " unequal to expected #columns = ", ncolumn))
        end
        # filter for multiallelic
        alt = rowgeno[5] # col5=alternative                
        ismultia =length(split(alt,",")) > 1                                
        if ismultia 
            nmultia += 1
            delmultiallelic && continue
        end          # 
        format = Symbol.(split(rowgeno[9],":"))
        formatcode = [get(prioritytuple,i,-1) for i in format] # -1 denotes missing fromat        
        founderformat = parse_rowgeno!(fgeno,rowgeno[fcols],
            formatcode, formatpriority, missingset,keepvcf)        
        if isempty(pcols)
            parentformat = founderformat
        else
            parentformat = parse_rowgeno!(pgeno,rowgeno[pcols],
                formatcode, formatpriority, missingset,keepvcf)        
        end
        offspringformat = parse_rowgeno!(offgeno, rowgeno[offcols],
            formatcode, formatpriority, missingset,keepvcf)
        # if nmarker == 1
        #     @info "" founderformat parentformat offspringformat rowgeno[fcols] rowgeno[pcols]
        # end
        if keepvcf
            if founderformat == parentformat 
                if parentformat == offspringformat
                    newformat = offspringformat
                else
                    newformat = string(parentformat,":", offspringformat)
                    fgeno .= string.(fgeno, ":.")
                    pgeno .= string.(pgeno, ":.")
                    offgeno .= string.(".:",offgeno)
                end
            else                               
                if founderformat == offspringformat
                    newformat = string.(founderformat,":", parentformat)
                    fgeno .= string.(fgeno, ":.")
                    pgeno .= string.(".:", pgeno)
                    offgeno .= string.(offgeno,":.")
                elseif parentformat == offspringformat                    
                    newformat = string(founderformat,":", parentformat)
                    fgeno .= string.(fgeno, ":.")
                    pgeno .= string.(".:", pgeno)
                    offgeno .= string.(".:", offgeno)
                else
                    newformat = string(founderformat, ":", parentformat,":", offspringformat)
                    fgeno .= string.(fgeno, ":.:.")
                    pgeno .= string.(".:", pgeno, ":.")
                    offgeno .= string.(".:.:",offgeno)
                end
            end
            # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
            res = vcat(rowgeno[1:8], [newformat], fgeno, pgeno, offgeno)
        else
            # keep information of some vcf columns            
            # parse vcf info column
            info = rowgeno[8]        
            newinfo, linkagegroup, poscm, doseerror, baseerror, allelicbias,overdispersion = parse_vcf_info(info)            
            newinfo = strip(newinfo,';')    
            isempty(newinfo) && (newinfo = missingstring)
            physchrom,physposbp,snpid,ref, alt, quality, filter = [i== "." ? missingstring : i for i in rowgeno[1:7]] 
            allelels = split(alt, ",")
            pushfirst!(allelels, ref)
            multiallele = length(allelels)            
            if founderformat != parentformat  
                if occursin(r"^GT",founderformat) 
                    if occursin(r"^GT",parentformat) 
                        parentformat = "GT"
                    else
                        parentformat = string(parentformat, "&GT")
                    end
                else
                    if occursin(r"^GT",parentformat) 
                        parentformat = string(founderformat, "&GT")
                    else
                        parentformat = string(founderformat, "&",parentformat)
                    end
                end
            end
            res = vcat([snpid,linkagegroup,poscm,physchrom,physposbp,parentformat, offspringformat], 
                [multiallele, join(allelels,"|"), quality, filter, newinfo],
                [doseerror,baseerror,allelicbias,overdispersion], fgeno,pgeno, offgeno)
            0
        end
        write(outio,join(res,outdelim),"\n")
        if rem(nmarker, 1000) == 0
            msg = string("marker index = ", nmarker, ", tused = ", 
                round(time()-startt,digits=1),"s")
            printconsole(logio,verbose,msg)
        end
    end    
    msg = string("end, #founders=",length(fcols), ", #midparents=", length(pcols),
        ", #offspring=",length(offcols), ", #markers=",nmarker, 
        ",#multiallelic=",nmultia, delmultiallelic ? "(del)" : "(keep)")
    printconsole(logio,verbose,msg)
end

function parse_vcf_info(info::AbstractString)    
    newinfo = ""
    linkagegroup = "NA"          
    poscm = "NA"    
    doseerror = "NA"
    baseerror = "NA"
    allelicbias = "NA"            
    overdispersion = "NA"    
    if info != "."
        infols = split.(split(info,";"),"=")
        otherinfo = infols[length.(infols) .== 1]
        if !isempty(otherinfo)
            newinfo *=";"*join(reduce(vcat,otherinfo),";")
        end
        infols = infols[length.(infols) .== 2]
        infols = [strip.(i) for i in infols]                
        infodict = OrderedDict([Pair(i...) for i in infols])
        if haskey(infodict,"LINKAGEGROUP")
            linkagegroup = infodict["LINKAGEGROUP"]
            delete!(infodict,"LINKAGEGROUP")
        end
        if haskey(infodict,"POSCM")
            poscm = infodict["POSCM"]
            delete!(infodict,"POSCM")
        end        
        if haskey(infodict,"DOSEERROR")
            doseerror = infodict["DOSEERROR"]
            delete!(infodict,"DOSEERROR")
        end
        if haskey(infodict,"BASEERROR")
            baseerror = infodict["BASEERROR"]
            delete!(infodict,"BASEERROR")
        end
        if haskey(infodict,"ALLELICBIAS")
            allelicbias = infodict["ALLELICBIAS"]
            delete!(infodict,"ALLELICBIAS")
        end        
        if haskey(infodict,"OVERDISPERSION")
            overdispersion = infodict["OVERDISPERSION"]
            delete!(infodict,"OVERDISPERSION")
        end        
        if !isempty(infodict)
            newinfo *=";"*join(map(i->string(i[1],"=",replace(i[2],","=>"|")),collect(infodict)),";")
        end
    end
    newinfo = strip(newinfo,';')    
    [newinfo, linkagegroup, poscm, doseerror, baseerror, allelicbias, overdispersion]
end

function parse_titlerow(titlerow::AbstractVector,pedfile::AbstractString;
    keepvcf::Bool=true,
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd(),    
    logio::Union{Nothing, IO})
    polyped = readpolyped(pedfile; commentstring,workdir)
    isfounder = polyped.parentinfo[!,:isfounder]
    founderid = polyped.parentinfo[isfounder,:individual]
    midparentid = polyped.parentinfo[.!isfounder,:individual]
    offspringid = polyped.offspringinfo[!,:individual]        
    dict = Dict(titlerow .=> 1:length(titlerow))    

    # founders and midparents
    fcols = [get(dict,i,nothing) for i in founderid]
    pcols = [get(dict,i,nothing) for i in midparentid]
    b = isnothing.(fcols)
    if any(b)
        msg = string(sum(b)," parents in pedfile but not in genofie: ",founderid[b])
        @warn msg
        printconsole(logio, false, "Warning: "*msg)
        fcols = fcols[.!b]
        founderid = founderid[.!b]
    end
    b = isnothing.(pcols)
    if any(b)
        msg = string(sum(b)," parents in pedfile but not in genofie: ",midparentid[b])
        @warn msg
        printconsole(logio, false, "Warning: "*msg)
        pcols = pcols[.!b]
        midparentid = midparentid[.!b]
    end
    # offspring
    offcols = [get(dict,i,nothing) for i in offspringid]
    b = isnothing.(offcols)
    if any(b)
        msg = string(sum(b)," offspring in pedfile but not in genofie: ",offspringid[b])
        @warn msg
        printconsole(logio, false, "Warning: "*msg)
        offcols = offcols[.!b]
    end
    if keepvcf
        res = titlerow[vcat(1:9, fcols, pcols, offcols)]
    else
        res = ["marker", "linkagegroup","poscm","physchrom", "physposbp", 
            "parentformat","offspringformat", 
            "multiallele", "alleles", "quality", "filter", "info",
            "doseerror","baseerror","allelicbias","overdispersion"]
        res = vcat(res,titlerow[vcat(fcols, pcols, offcols)])
    end
    delim = keepvcf ? '\t' : ','    
    newtitlerow = join(res,delim)
    polyped, fcols, pcols, offcols,newtitlerow
end

function parse_rowgeno!(resgeno::AbstractVector,
    rowgeno::AbstractVector,formatcode::AbstractVector,
    formatpriority::AbstractVector, missingset::AbstractVector,keepvcf::Bool)
    ncell = length(rowgeno)
    # ncell == length(resgeno) || error(string("unexpected length of rowgeno=",rowgeno))
    formatls = Vector(undef,ncell)
    if all(formatcode .< 0) #-1 denotes missing fromat
        resgeno .= missing
        formatls .= missing
    else
        # 
        ThreadsX.foreach(eachindex(resgeno,formatls,rowgeno)) do i
            @inbounds resgeno[i], formatls[i] = parse_vcfcell(rowgeno[i], formatcode,
                formatpriority, missingset,keepvcf)
        end
    end
    commonformat, missval = parse_format(formatls,formatpriority; isvcfformat = keepvcf)       
    resgeno[ismissing.(resgeno)] .= missval
    if ismissing(commonformat)
        if keepvcf
            commonformat = "GT"
            resgeno .= "."
        else
            commonformat = "GT_unphased"
            resgeno .= missing
        end    
    end
    commonformat
end

function parse_format(formatls::AbstractVector,formatpriority::AbstractVector; isvcfformat::Bool)
    ls = skipmissing(formatls)
    if isempty(ls)
        (missing,missing)
    else
        # commonformat = commonest(ls)               
        # assuming that only frist 2 characters of formats  are retained in formatpriority
        formatls2 = unique(ls)
        if length(formatls2) == 1
            commonformat = formatls2[1]            
        else
            commonformat = join(sort(unique([i[1:2] for i in formatls2])),"&")
        end
        # posls = [findfirst(formatpriority .== i[1:2]) for i in formatls2]
        # commonformat =  formatls2[argmin(posls)]
        if isvcfformat
            missval = "."
        else
            missval = "NA"
        end
        commonformat,missval
    end
end

function parse_vcfcell(cellgeno::AbstractString,formatcode::AbstractVector,
    formatpriority::AbstractVector, missingset::AbstractVector,keepvcf::Bool)
    n = length(formatcode)
    geno = split(cellgeno,":")
    ng = length(geno)
    if ng < n
        geno = vcat(geno, repeat(["."],n-ng))
    elseif ng > n
        geno = geno[1:n]
    end
    b = [(geno[i] in missingset) || formatcode[i] < 0 for i in 1:n]  # -1 of formatcode denotes missing fromat
    if all(b)
        genoset = geno[formatcode .>= 1]
        if in("./.", genoset)
            g = "./."
            f = "GT" 
        elseif in(".|.", genoset)
            g = ".|."
            f = "GT" 
        else
            g = missing
            f = missing
        end
    else
        # the less the formatcode element, the higher the priority
        upbound = length(formatpriority)+10
        i = argmin(@. formatcode + upbound * b)             
        g = geno[i]
        f = formatpriority[formatcode[i]]
    end
    if keepvcf || ismissing(f)
        g, f
    else
        if f == "GT"
            if occursin("|", g)
                f2 = "GT_phased"
            elseif occursin("/", g)
                f2 = "GT_unphased"
            else
                error(string("unknown GT geno: ",g))
            end
            g2 = replace(g, "."=>"N")
        else
            g2 = replace(g, "," => "&") # seperator & for AD, GL, and GP
            f2 = f
        end
        g2, f2
    end
end
