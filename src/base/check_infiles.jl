function check_infiles(genofile::AbstractString,pedfile::AbstractString;     
    commentstring = "##",
    workdir::AbstractString=pwd(),
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    # check genofile extension and exist
    fileext = last(split_allext(genofile))
    if in(fileext,[".vcf",".vcf.gz"])         
        msg = string("check genofile extension = ",fileext, ": OK")            
        printconsole(io,verbose,msg)        
    elseif in(fileext,[".csv",".csv.gz"]) 
        msg = string("check genofile extension = ",fileext, ": VCF-format genofile is recommended")
        printconsole(io,verbose,msg)        
    else
        msg = string("check genofile extension = ",fileext, ": FAILED. Must be .vcf or .vcf.gz")
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    genofile2 = getabsfile(workdir,genofile)
    if isfile(genofile2) 
        msg =  string("check genofile existence: YES")
        printconsole(io,verbose,msg)        
    else
        msg = string(genofile2, " does not exist")
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    # check pedfile and exist    
    fileext = last(split_allext(pedfile))    
    # pedigree file
    if in(fileext,[".csv",".csv.gz"]) 
        msg = string("check pedfile extension = ",fileext, ": OK")            
        printconsole(io,verbose,msg)        
    else
        msg = string("check pedfile extension = ",fileext, ": FAILED. Must be .csv or .csv.gz")
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    pedfile2 = getabsfile(workdir,pedfile)
    if isfile(pedfile2) 
        msg = string("check pedfile existence: YES")
        printconsole(io,verbose,msg)        
    else
        printconsole(io,false,"Error: "*msg)
        error(string(pedfile2, " does not exist"))
    end            
    # check consistency betwen sampleid of genofile and pedfile
    if in(last(split_allext(genofile)),[".csv",".csv.gz"]) 
        samples = csv_get_samples(genofile; commentstring, workdir)    
    else
        samples = vcf_get_samples(genofile; commentstring, workdir)    
    end
    check_allunique(samples,"samples";io, verbose)    
    polyped = readpolyped(pedfile; commentstring, workdir)
    parents = polyped.parentinfo[!,:individual]
    check_allunique(parents,"parents";io, verbose)    
    d = setdiff(parents,samples)            
    if isempty(d)                     
        msg = string("check parents: OK. All ", length(parents), " parents = ", parents, " in pedfile are also in samples in genofile")
        printconsole(io,verbose,msg)        
    else
        if length(d) == length(parents)                    
            msg = string("check parents: FAILED. All ", length(d), " parents = ", d, " in pedfile but not in samples in genofile")
            printconsole(io,false,"Error: "*msg)
            error(msg)
        else
            msg = string("check parents: PASS. Among ", length(parents), " parents, ", length(d), " parents = ", d, " in pedfile but not in samples in genofile. Ignore these parentrs and their offspring.")                    
            printconsole(io,false,"Warning: "*msg)
            @warn msg
        end
    end
    offspring = polyped.offspringinfo[!,:individual]
    check_allunique(offspring,"offspring";io, verbose)    
    d = setdiff(offspring,samples)                        
    if isempty(d)                     
        msg = string("check offspring: OK. All ", length(offspring), " offspring in pedfile are also in samples in genofile")
        printconsole(io,verbose,msg)        
    else
        msg = string("check offspring: PASS. Among ", length(offspring), " offspring, ignore ", length(d), " offspring = ", d, " in pedfile but not in samples in genofile")
        printconsole(io,false,"Warning: "*msg)
        @warn msg
    end
    d = setdiff(samples, parents, offspring)
    if isempty(d)                     
        msg = string("check samples: OK. All ", length(samples), " samples in genofile are also in pedfile")
        printconsole(io,verbose,msg)        
    else
        msg = string("check samples: PASS. Among ", length(samples), " samples, ignore ", length(d), " samples in genofile but not in pedfile")        
        printconsole(io,false,"Warning: "*msg)
        @warn msg
    end            
    msg = string("Pass all checks for input genofile and pedfile")
    printconsole(io,verbose,msg)        
end

function check_allunique(samples::AbstractVector, varname::AbstractString;
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    if allunique(samples)        
        msg = string("check ", varname, " allunique: YES")            
        printconsole(io,verbose,msg)
    else
        ls = unique(samples)
        isdupe = [count(==(i)) > 1 for i in ls]
        msg = string("check ", varname, " allunique: FAILED. Duplicated ", varname, " = ",ls[isdupe])            
        printconsole(io,verbose,msg)
    end
end


function csv_get_samples(csvgenofile::AbstractString; commentstring = "##",workdir::AbstractString=pwd())
    inext = last(split_allext(csvgenofile))
    inext in [".csv",".csv.gz"] || @error string("csvgenofile ext must be .csv or .csv.gz")
    in_open = inext in [".csv.gz"] ? GZip.open : open
    csvgenofile2 = getabsfile(workdir,csvgenofile)
    lastcomment = findlastcomment(csvgenofile2; commentstring)
    in_open(csvgenofile2,"r") do io
        for _ in 1:lastcomment
            readline(io)
        end
        # parse title row
        ls = split(readline(io,keep=false),",")[4:end] 
        res = [strip(string(i),[' ','"']) for i in ls]
        setdiff(res, ["physchrom", "physposbp", "parentformat", "offspringformat", "multiallele", "alleles", "quality", "filter", "info", "doseerror", "baseerror", "allelicbias", "overdispersion"])
    end
end


function vcf_get_samples(vcffile::AbstractString; commentstring = "##",workdir::AbstractString=pwd())
    inext = last(split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    vcffile2 = getabsfile(workdir,vcffile)
    isfile(vcffile2) || error(string(vcffile2, " does not exist"))
    lastcomment = findlastcomment(vcffile2; commentstring)    
    in_open(vcffile2,"r") do io
        for _ in 1:lastcomment
            readline(io)
        end
        # parse title row
        ls = split(readline(io,keep=false),"\t")[10:end]        
        [strip(string(i),[' ','"']) for i in ls]
    end
end

