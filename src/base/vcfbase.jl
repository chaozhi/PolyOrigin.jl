

"""
    vcffilter(vcffile; keyargs...)

filter markers line by line. 

# Positional arguments

`vcffile::AbstractString`: vcf genofile. 

# Keyword arguments

`setmarkerid::Union{Nothing,Bool}=nothing`: if true, set markerid. If it is nothing, setmarkerid = true only if markerid is missing

`delsamples::Union{Nothing,AbstractVector}=nothing`: list of sample IDs to be deleted. If it is nothing, no deletion of samples 

`deldupe::Bool=false`: if true, delete sucessive markers that have exactly duplicated genotypes in format of GT

`isdelmultiallelic::Bool=true`: if true, delete markers with >2 alleles. 

`isdelmonomorphic::Bool=true`: if true, delete markers with single allele. 

`seqstretch::Integer=0`: delete non-initial markers in a sequence stretch of length <= seqstretch (in bp), assuming marker are ordered by physical positions. If it is not positive, no filtering for short streches.

`maxmiss::Real = 0.99`: delete markers with missing fraction > maxmiss

`minmaf::Real = 0.01`: delete markers with minor allele frequency < minmaf

`commentstring::AbstractString="##"`: the lines beginning with  are ignored

`outstem::AbstractString=popid`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

`logfile::AbstractString = outstem*"_vcffilter.log"`: log filename. 

`verbose::Bool=true`: if true, print details on the stdout.

"""
function vcffilter(vcffile::AbstractString;
    setmarkerid::Union{Nothing,Bool}=nothing,  # if it is nothign, set it in format of CHROM_POS only if markerid is missing
    delsamples::Union{Nothing,AbstractVector}=nothing,     
    deldupe::Bool=false, 
    isdelmultiallelic::Bool=true,    
    isdelmonomorphic::Bool=true,    
    keeponlypass::Bool=true, 
    seqstretch::Integer=0, # 64 bp for length of sequence tags
    maxmiss::Real = 0.99,
    minmaf::Real = 0.01,
    commentstring::AbstractString="##",
    outstem::AbstractString = first(split_allext(basename(vcffile))),      
    workdir::AbstractString=pwd(),
    logfile::AbstractString = outstem*"_vcffilter.log",
    verbose::Bool=true)    
    starttime = time()
    logio = set_logfile_begin(logfile, workdir, "vcffilter"; verbose,delim="-")
    inext = last(split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    msg = string("list of args: \n",
        "vcffile = ", vcffile, "\n",        
        "setmarkerid = ", setmarkerid, "\n",        
        "delsamples = ", delsamples, "\n",        
        "deldupe = ", deldupe, " \n",        
        "isdelmonomorphic = ", isdelmonomorphic, "\n",        
        "isdelmultiallelic = ", isdelmultiallelic, "\n",                
        "keeponlypass = ", keeponlypass, "\n",                
        "seqstretch = ", seqstretch, " \n",        
        "maxmiss = ", maxmiss, "\n",        
        "minmaf = ", minmaf, "\n",        
        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)    
    msg = deldupe ? "delete sucessive duplicate markers\n" : ""
    seqstretch > 0 && (msg *= string("delete noninitial markers in a sequence stretch of length < ", seqstretch, "\n"))
    msg *= string("filter sequentially for biallelic, non-monomorphic, missing <= ", maxmiss, ", and maf >= ",minmaf)
    printconsole(logio,verbose,msg)        
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    vcffile2=getabsfile(workdir,vcffile)
    lastcomment = findlastcomment(vcffile2; commentstring)    
    outext = inext
    in_open(vcffile2,"r") do inio        
        out_open = outext in [".vcf.gz"] ? GZip.open : open
        outfile = outstem*"_vcffilter"*outext
        outfile2 = getabsfile(workdir,outfile)
        out_open(outfile2, "w") do outio
            # commentstring = "##" for vcf file
            filedate = string(commentstring*"fileDate=",replace(string(Date(now())),"-"=>""))
            write(outio, filedate,"\n")
            write(outio, commentstring*"source=RABBIT\n")
            for _ in 1:lastcomment
                line = readline(inio,keep=true)                
                write(outio, line)
            end            
            # parse title row    
            titlerow = split(readline(inio,keep=false),"\t")
            samples = string.(strip.(titlerow[10:end]))
            if isnothing(delsamples)
                delsamples2 = []
            elseif eltype(delsamples) <: Integer 
                delsamples2 = samples[delsamples]
            elseif eltype(delsamples) <: AbstractString
                delsamples2 = string.(strip.(delsamples))
                d = setdiff(delsamples2, samples)
                if !isempty(d)
                    @warn string("delsamples that are not in vcffile: ",d)
                end
                intersect!(delsamples2,samples)
            else
                @error string("unknown type=",eltype(delsamples), " for input delsamples=",delsamples)
            end               
            dict = Dict(samples .=> 1:length(samples))
            cols = [get(dict,old,nothing) for old in delsamples2]
            any(isnothing.(cols)) && error(string("unexpected cols=",cols))
            cols .+= 9
            deleteat!(titlerow, cols)
            write(outio, join(titlerow,"\t"),"\n")
            nind = length(titlerow)-9
            msg = string("#samples=", nind, ",sample (up to 10): ", join(titlerow[10:min(19,length(titlerow))],","))
            printconsole(logio,verbose,msg)                      

            # read marker line
            startt = time()            
            nmarker_incl = 0
            nmarker = 0
            ndupe = 0
            nnotpass = 0
            nstretch = 0            
            nmultia = 0            
            nmono = 0
            nmiss = 0
            nmaf = 0
            binchr = ""
            binpos = -10*seqstretch
            bingeno = ["" for _ in 1:nind]
            while !eof(inio)
                # rem(nmarker, 1000) == 0 && (startt = time())
                nmarker += 1
                rowgeno = split(readline(inio,keep=false),"\t")
                deleteat!(rowgeno, cols)           
                # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
                formatls = split(rowgeno[9],":")  # col9 = FORAMT
                gt = findfirst(x->x=="GT",formatls)
                if isnothing(gt) 
                    @error string("GT format does not exist! marker_row=",rowgeno) maxlog = 10
                end
                genols = [i[gt] for i in split.(rowgeno[10:end],":")]                
                nowchr = rowgeno[1]
                nowpos = tryparse(Int, rowgeno[2])                                
                nowgeno = genols
                ismultia =length(split(rowgeno[5],",")) > 1  # col5=alternative      
                if keeponlypass && (rowgeno[7] != "." && !occursin("PASS", uppercase(rowgeno[7]))) #col7 = FILTER
                    isdel = true
                    nnotpass += 1                    
                elseif nowchr == binchr && all(nowgeno .== bingeno) && deldupe
                    isdel = true 
                    ndupe += 1
                else                    
                    if nowchr == binchr && !isnothing(binpos) && !isnothing(nowpos) && abs(nowpos - binpos)+1 <= seqstretch
                        isdel = true 
                        nstretch += 1
                    else
                        binchr = nowchr
                        binpos = nowpos
                        bingeno = nowgeno
                        if ismultia 
                            nmultia += 1
                            isdel = isdelmultiallelic                    
                        else
                            isdel = false
                        end            
                    end                
                end                
                if !isdel
                    # filter for non-monomorphic        
                    # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]                                               
                    alleles = replace(string(unique(genols)...),"/"=>"","|"=>"","."=>"")
                    alleleset = unique(split(alleles,""))
                    ismono = isempty(alleles) || length(alleleset) == 1    # isempty(alleles) == 100% missing
                    if ismono
                        nmono += 1
                        isdel = isdelmonomorphic
                    end                
                    if length(alleleset) >=3  && !ismultia && !deldupe
                        msg = string("inconsistent multi-allelic: marker index = ", nmarker, ", allele_set=",alleleset, ", ALT=",alt)
                        @warn msg maxlog = 20
                        printconsole(logio, false, msg)                            
                    end                        
                    if !isdel                        
                        # filter missing                        
                        if ismono
                            if minmaf > 0 && maxmiss < 1.0
                                nmaf += 1
                                isdel = true
                            end
                        else
                            freqmiss = mean([in(i,[".","./.",".|."]) for i in genols])
                            if freqmiss > maxmiss 
                                nmiss += 1
                                isdel = true
                            else
                                # filter for maf
                                if minmaf > 0 
                                    alleles = replace(string(genols...),"/"=>"","|"=>"","."=>"")
                                    freqls = [count(i,alleles)/length(alleles) for i in alleleset]
                                    maf = 1.0- max(freqls...)
                                    if maf < minmaf || maf > 1-minmaf
                                        nmaf += 1
                                        isdel = true
                                    end
                                end
                            end
                        end
                    end
                end
                if !isdel 
                    # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]                                               
                    isimputeid = setmarkerid
                    if isnothing(setmarkerid) && strip(rowgeno[3]) == "."
                        isimputeid = true
                    end
                    if  !isnothing(isimputeid) && isimputeid
                        physchr, physpos = strip.(rowgeno[1:2])
                        if physchr == "." || physpos == "."
                            rowgeno[3] = string("snp",nmarker)
                        else
                            rowgeno[3] = string(physchr,"_", physpos)
                        end
                    end
                    write(outio, join(rowgeno,"\t"),"\n")                    
                    nmarker_incl += 1
                end
                if rem(nmarker, 1000) == 0
                    msg = string("#marker=", nmarker, 
                        ", #marker(keep)=", nmarker_incl, 
                        keeponlypass ? string(", #del_notpass=", nnotpass) : "",                         
                        deldupe ? string(", #del_dupe=", ndupe) : "",                         
                        seqstretch > 0 ? string(", #del_stretch=", nstretch) : "",                         
                        string(", #multi", isdelmultiallelic ? "=" : "(keep)=",nmultia), 
                        string(", #mono", isdelmonomorphic ? "=" : "(keep)=",nmono), 
                        ", #>maxmiss=", nmiss, ", #<minmaf=", nmaf,                         
                        ", tused=", round(time()-startt,digits=1),"s")
                    printconsole(logio,verbose,msg)
                    flush(outio)
                end
            end            
            msg = string("#individuals=", nind, ", #markers=",nmarker, 
                ", #marker(keep)=", nmarker_incl, 
                keeponlypass ? string(", #del_notpass=", nnotpass) : "",                         
                deldupe ? string(", #del_dupe=", ndupe) : "",                         
                seqstretch > 0 ? string(", #del_stretch=", nstretch) : "",                         
                string(", #multi", isdelmultiallelic ? "=" : "(keep)=",nmultia), 
                string(", #mono", isdelmonomorphic ? "=" : "(keep)=",nmono), 
                ", #>maxmiss=", nmiss, ", #<minmaf=", nmaf,                                         
            )
            printconsole(logio,verbose,msg)
        end
        msg = string("output vcffile: ", outfile)
        printconsole(logio,verbose, msg)
    end
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, logio, tused,"vcffilter"; verbose,delim="-")
    return 
end

