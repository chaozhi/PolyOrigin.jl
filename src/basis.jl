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

# """
#     getabsfile(workdir,file)
#
# get an absolute path by adding the work directory to file if necessary.
#
# """
function getabsfile(workdir::AbstractString,file::AbstractString)
    if isabspath(file)
        file
    else
        abspath(workdir,file)
    end
end

function printpkgst(io::Union{Nothing,IO},verbose::Bool,pkgname::AbstractString)
    verbose && Pkg.status(pkgname)
    if io !== nothing
        redirect_stdout(io) do
            Pkg.status(pkgname)
        end
    end
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
        cols = cols[cols .!= nothing]
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
    regionkey::AbstractString = "juliaDict",
    delim::AbstractChar=',')
    io = typeof(outfile) <: AbstractString ? open(outfile, "w+") : outfile
    for (key, value) in dict
        write(io,string(regionkey,string(delim), key,"\n"))
        if isa(value,AbstractDataFrame)
            CSV.write(io, value; delim=delim,
                header=true,append=true)
            flush(io)
        else
            # error in writing nothing
            writedlm(io,value,delim)
        end
    end
    typeof(outfile) <: AbstractString ? close(io) : flush(io)
end


function metroplos_norm(f::Function,xnow;
    maxiter=maxiter,temperature::Real=0,
    initlogl::Union{Nothing,Real}=nothing,propstd::Real=1)
    logl = initlogl==nothing ? f(xnow) : initlogl
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
    # ls=[Vector{Float64}(xls[i]) for i=res]
    # [grid[2:end] ls]
    res
end

function logsumexp(a::AbstractVector)
    off, ind = findmax(a)
    isinf(off) && off<0 && (return -Inf)
    b = exp.(a .- off)
    b[ind] -= 1.0
    # using KahanSummation
    # log1p(sum_kbn(b))+off
    log1p(sum(b))+off
end

function getemptydf(column_eltypes::AbstractVector{T}, cnames::Vector{Symbol})::DataFrame where T<:Type
    columns = [Vector{elty}(undef,0) for elty = column_eltypes]
    return DataFrame(columns, cnames, copycols=false)
end
