function exportall(mod)
    for n in names(mod, all = true)
        if Base.isidentifier(n) && n ∉ (Symbol(mod), :eval)
            @eval mod export $n
        end
    end
end

function rmfile(file::AbstractString;
    workdir::AbstractString=pwd())
    olddir = pwd()
    a = splitpath(file)
    if length(a)==1
        cd(workdir)
        rm(file,force=true)
    else
        workdir2,file2=splitdir(file)
        cd(workdir2)
        rm(file2,force=true)
    end
    cd(olddir)
end

function riffle(A::AbstractVector, x)
    n=length(A)
    res = permutedims(hcat(A, repeat([x],n)))
    reshape(res,length(res))[1:end-1]
end

function stringjoin(A::AbstractVector, x)
    string(riffle(A,x)...)
end


function sortsplit(A::AbstractMatrix,col::Integer)
   mtx=A[sortperm(A[:, col]), :]
   v = mtx[:,col]
   [mtx[v .== i,:] for i=unique(v)]
end

function splitindex(f::Function,A::AbstractVector)
    size(A,1)==1 && return [1:1]
    res=Vector{typeof(1:1)}()
    i0=1
    for i=2:size(A,1)
        if !f(A[i-1],A[i])
            push!(res,i0:i-1)
            i0=i
        end
    end
    push!(res,i0:size(A,1))
    res
end

function splitindex(A::AbstractVector)
    f(x,y)= x == y
    splitindex(f,A)
end

function splitvec(f::Function,A::AbstractVector)
    res=Vector{typeof(A)}()
    i0=1
    for i=2:size(A,1)
        if !f(A[i-1],A[i])
            push!(res,A[i0:i-1])
            i0=i
        end
    end
    push!(res,A[i0:end])
    res
end

function splitvec(A::AbstractVector)
    f(x,y)= x == y
    splitvec(f,A)
end


# using DelimitedFiles, DataStructures
function readdlm2dict(filename::AbstractString;
    delim::AbstractChar=',', comment_char::AbstractChar='#',
    isordereddict::Bool=false)
    if !isfile(filename)
        error(filename," does not exist")
    end
    ls = readdlm(filename,delim,comment_char=comment_char,comments=true)
    regionkey = ls[1,1]
    seg = findall(x->x==regionkey,ls[:,1])
    push!(seg,size(ls,1)+1)
    datadict = isordereddict ? OrderedDict : Dict
    datadict([begin
        # println("i=",i)
        i1,i2=seg[i:i+1]
        key=ls[i1,2]
        cols = [findfirst(x->x=="",ls[j,:]) for j=i1+1:i2-1]
        # cols = cols[cols .!= nothing]
        cols = cols[.!isnothing.(cols)]
        ncol = cols == [] ?  length(ls[i1+1,:]) : max(cols...) -1
        segmtx=ls[i1+1:i2-1,1:ncol]
        df=DataFrame(segmtx[2:end,:],Symbol.(segmtx[1,:]))
        for j=1:size(df,2)
            ty = typejoin(unique(typeof.(df[!,j]))...)
            df[!,j]=convert.(ty, df[!,j])
        end
        key => df
    end for i = 1:length(seg)-1])
end


function savedict2dlm(outfile::Union{AbstractString,IOStream},
    dict::AbstractDict;
    regionkey::AbstractString = "PolyOrigin",
    delim::AbstractChar=',')
    io = typeof(outfile) <: AbstractString ? open(outfile, "w+") : outfile
    for (key, value) in dict
        write(io,string(regionkey,string(delim), key,"\n"))
        if isa(value,AbstractDataFrame)
            if isempty(value)
                write(io, join(names(value),delim),"\n")
            else
                CSV.write(io, value; delim=delim,
                    header=true,append=true)
            end
            flush(io)
        else
            # error in writing nothing
            writedlm(io,value,delim)
        end
    end
    typeof(outfile) <: AbstractString ? close(io) : flush(io)
end

function readmultitable(filename::AbstractString;
    delim::AbstractChar=',',
    commentstring::AbstractString="##",
    missingstring::AbstractString="NA", # not adopted when isparse is false
    isparse::Bool= true, 
    isordereddict::Bool = true)
    ext = last(splitext(filename))
    myopen = ext == ".gz" ? GZip.open : open
    myopen(filename, "r") do io
        lines = readlines(io,keep=isparse)        
        b = startswith.(lines,commentstring)
        deleteat!(lines,b)
        length(lines) >= 1 || @error "empty file"
        regionkey = first(split(lines[1],delim))
        seg = findall(startswith.(lines,regionkey))
        push!(seg,length(lines)+1)
        mydict = isordereddict ? OrderedDict : Dict    
        mydict([begin            
            i1=seg[i]
            i2 =seg[i+1]
            tabtitle = split(lines[i1],delim)            
            length(tabtitle) >= 2 || error("tabe name not exist")
            tabkey = strip(chomp(tabtitle[2]),'\"')
            tabval = view(lines, i1+1:i2-1)                        
            if isparse
                # bug in CSV.read, a cell with value of  a very very long string will be truncated!!
                # for example, sparse probability matrix
                val = CSV.read(IOBuffer(join(tabval)),DataFrame;
                    missingstring,dateformat = nothing)          # prevent parsing ancestralgenotype/state into date
            else
                strval = reduce(hcat, split.(chomp.(tabval),delim))          
                mtx =  permutedims(strip.(strval[:,2:end],'\"'))
                cols = strip.(strval[:,1],'\"')
                val = DataFrame(mtx,cols)
            end
            tabval .= ""                        
            tabkey => val
        end for i=1:length(seg)-1])
    end
end


function metroplos_norm(f::Function,xnow;
    maxiter=maxiter,temperature::Real=0,
    initlogl::Union{Nothing,Real}=nothing,propstd::Real=1)
    logl = isnothing(initlogl) ? f(xnow) : initlogl
    acceptls=zeros(maxiter)
    for it=1:maxiter
        xprop = rand(Normal(xnow,propstd))
        proplogl = f(xprop)
        if temperature <=0
            isaccept = proplogl >= logl
        else
            isaccept = rand() < min(1,exp((proplogl - logl)/temperature))
        end
        if isaccept
            xnow, logl = xprop, proplogl
        end
        acceptls[it] = isaccept
    end
    xnow,logl,mean(acceptls)
end


function gridpartition(xls::AbstractVector,gridsize::Real)
    xls .-= xls[1]
    grid=collect(0:gridsize:xls[end])
    grid[end] += gridsize+0.1
    res=Vector{Vector}()
    b=i=1
    bin=[]
    while true
        if grid[b]<=xls[i]<grid[b+1]
            push!(bin,i)
            i+=1
        else
            push!(res,bin)
            b+=1
            bin=[]
        end
        i>length(xls) && break
    end
    push!(res,bin)
    res
end

function getemptydf(column_eltypes::AbstractVector{T}, cnames::Vector{Symbol})::DataFrame where T<:Type
    columns = [Vector{elty}(undef,0) for elty = column_eltypes]
    return DataFrame(columns, cnames, copycols=false)
end

function set_logfile_begin(logfile::Union{Nothing,AbstractString,IO},
    workdir::AbstractString,
    funname::AbstractString;
    verbose::Bool=true,
    delim::Union{Nothing,AbstractString}="-")
    if isnothing(logfile)
        io=nothing
    else
        if isa(logfile, AbstractString)
            io=open(getabsfile(workdir,logfile), "w")
            msg = string(funname, ", logfile=", logfile,
                ", ", round(Dates.now(),Dates.Second))
            printconsole(io,verbose,msg)
        else
            # isa(io, IO)
            io=logfile
        end
    end
    if !isnothing(delim)
        nleft = (76 - length(funname)) ÷ 2
        nright = 76 - nleft -length(funname)
        printconsole(io,verbose,string(repeat(delim,nleft),
            funname,repeat(delim,nright)))
    end
    io
end

function set_logfile_end(logfile::Union{Nothing,AbstractString,IO},
    io::Union{Nothing,AbstractString,IO},
    tused::Real,
    funname::AbstractString;
    verbose::Bool=true,
    delim::Union{Nothing,AbstractString}="-")
    msg= string("End, ", round(Dates.now(),Dates.Second),", tused = ",
        tused, " seconds by ",funname)
    printconsole(io,verbose,msg)
    isnothing(delim) || printconsole(io,verbose,repeat(delim,76))
    isa(logfile, AbstractString) && close(io)
end

function split_allext(filename::AbstractString)
    list_rabbitext = [".vcf",".csv",".tar",".gz",".png",".log",".jl",".txt",".tsv",".xlsx",".pdf",".zip",".7z"]
    filebase,ext = splitext(filename)
    isempty(ext) && return (filebase,ext)
    while true    
        filebase,ext2 = splitext(filebase)
        if isempty(ext2)  
            break
        else
            if ext2 in list_rabbitext
                ext = ext2 * ext
            else
                filebase *=ext2
                break
            end
        end        
    end
    filebase,ext
end

function create_tar(dir::AbstractString,outfile_tar::AbstractString)
    if !occursin(r".tar$",outfile_tar)
        @warn string(outfile_tar, " is not .tar")
    end
    isdir(dir) || @error string(dir, " is not a directory")
    Tar.create(dir, outfile_tar)
end

function create_targz(dir::AbstractString,outfile_targz::AbstractString)
    if !occursin(r".tar.gz$",outfile_targz)
        @warn string(outfile_targz, " is not .tar.gz")
    end
    isdir(dir) || @error string(dir, " is not a directory")
    open(outfile_targz, "w") do tar_gz    
        tar = GzipCompressorStream(tar_gz)
        try
            Tar.create(dir, tar)
        finally
            close(tar)
        end    
    end
    outfile_targz
end

function extract_tar(tarfile::AbstractString, outdir::AbstractString)
    if !occursin(r".tar$",tarfile)
        @warn string(tarfile, " is not .tar")
    end
    isdir(outdir) || @error string(outdir, " is not a directory")
    Tar.extract(tarfile,outdir)
end

function extract_targz(targzfile::AbstractString, outdir::AbstractString)
    if !occursin(r".tar.gz$",targzfile)
        @warn string(targzfile, " is not .tar.gz")
    end
    isdir(outdir) || @error string(outdir, " is not a directory")
    tar_gz = open(targzfile)
    try
        tar = GzipDecompressorStream(tar_gz)
        try
            Tar.extract(tar,outdir)
        finally
            close(tar)
        end
    finally
        close(tar_gz)
    end
end

function read_headlines(file::AbstractString, nline::Integer)
    nline <= 0 && return ""
    isgz = last(splitext(file)) == ".gz"
    myopen = isgz ? GZip.open : open
    myopen(file,"r") do io
        lines = ""
        for i in 1:nline
            eof(io) && break
            line = readline(io,keep=true)
            lines *= line
        end
        rstrip(lines,['\n'])
    end
end

function findlastcomment(file::AbstractString; commentstring::AbstractString="##")
    isgz = last(splitext(file)) == ".gz"
    myopen = isgz ? GZip.open : open
    myopen(file,"r") do io
       commentline = 0
       while !eof(io)
           word = readline(io)
           if isempty(word) 
                break
           else
                startswith(word,commentstring) ? (commentline+=1) : break
           end
       end
       commentline
    end
end

function readdataframe(filename::AbstractString;
    delim::AbstractChar= ',',
    commentstring::AbstractString="##",    
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    ext = last(splitext(filename))
    myopen = ext == ".gz" ? GZip.open : open
    myopen(getabsfile(workdir,filename), "r") do io
        lines = readlines(io,keep=true)
        b = startswith.(lines,commentstring)
        deleteat!(lines,b)        
        length(lines) >= 1 || @error "empty file"        
        df = CSV.read(IOBuffer(join(lines)),DataFrame;delim, missingstring)                
        df
    end
end

function printpkgst(io::Union{Nothing,IO},verbose::Bool,pkgname::AbstractString)
    verbose && Pkg.status(pkgname)
    if !isnothing(io)
        redirect_stdout(io) do
            Pkg.status(pkgname)
        end
    end
end

function getabsfile(workdir::AbstractString,file::AbstractString)  # equivalent to abspath(path1,path2)
    if isabspath(file)
        file
    else
        abspath(workdir,file)
    end
end

function printconsole(io::Union{Nothing,IO},verbose::Bool,msg::AbstractString)
    verbose && @info(msg)
    if !isnothing(io)
        write(io,string(msg,"\n"))
        flush(io)
    end
end


function logsumexp(a::AbstractVector)
    amax, ind = findmax(a)
    isinf(amax) && amax<0 && (return -Inf)
    b = exp.(a .- amax)
    b[ind] -= 1.0    
    # log1p(x): log(1+x)
    log1p(sum(b))+amax
end

function logsumexp!(a::AbstractVector)
    scale = logsumexp(a)
    a .-= scale
    scale
end


function normprobvec!(probvec::AbstractVector)
    a = sum(probvec)    
    probvec ./= a
    a    
end


function metroplis1d(loglfun::Function;
    xstart::Real=0.0, temperature::Real=1.0,
    constraint::Function = x->true,
    stepsize::Real=1.0, nstep::Integer=1)
    xnow = xstart
    logl = loglfun(xnow)
    chain = [[xnow, logl,NaN]]
    for it in 1:nstep
        xprop = rand(Normal(xnow, stepsize))
        if constraint(xprop)
            proplogl = loglfun(xprop)
            if temperature < 10^(-5.0)
                cond = proplogl > logl
            else
                cond = rand() < min(1.0, exp((proplogl-logl)/temperature))
            end
            if cond
                xnow, logl = xprop, proplogl
            end
        else
            cond = false
        end
        push!(chain, [xnow, logl, cond])
    end
    chain
end