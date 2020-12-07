"""
    PolyAncestry

mutable struct that stores phased polygeno and ancestral inference.

    PolyAncestry(markermap,parentgeno,parentinfo,offspringinfo,designinfo,
        statespace,valentprob,genoprob)

Constructor from phased polygeno and ancestral inference.

# Fields

`markermap::Vector{DataFrame}`: marker map for each chromosome. `markermap[c]`
gives the dataframe of chromosome c with columns [:marker, :chromosome, :position].

`parentgeno::Vector{Matrix}`: genotypic data in parents. `parentgeno[c][m,i]`
gives the genotype of parent i at marker m of chromosome c.

`parentinfo::DataFrame`: parent information with columns [:individual, :ploidy].

`offspringinfo::DataFrame`: offspring information with columns
[:individual,:population].

`designinfo::DataFrame`: design information with columns
[:population, :parent1, :parent2, ...], where element D\\_{ij} denotes the nubmer of gametes contributed to each offspring in i-th population by j-th parent.
D\\_{ij}=2 means that i-th population is produced by self-fertilization of j-th parent.

`delmarker::DataFrame`: dataframe for collecting markers that were removed from markermap,
in the stages of parental phasing or map refinement. It has columns
[:marker, :chromosome, :position].

`correction::DataFrame`: dataframe for collecting parental genotype correction.
It has columns [:marker,:chromosome,:parent,:old_genotype,:new_genotype,:old_nerr,:new_nerr].

`statespace::Dict`: specify valents and ancestral genotypes for each population.
Each key is a population id, and its value is in turn a dict with keys:
"parent", "parentindex", "valent", "valentneighbor", and "groupstate".

`valentprob::Union{Nothing,Vector}`, posterior valent probability for each offspring in each chromosome.
valentprob[c][o] for offspring o in chromsome c is a matrix with three columns being
valentindex, loglike, and posterior probability. Here posterior is calcuated from loglike,
assuming a discrete uniform prior distribution of valents.

`genoprob::Vector`: posterior ancestral genotype probability for each offspring
in each chromosome. genoprob[c][o] for offspring o in chromsome c is a sparse matrix,
with element (m,s) being the probability of genotype s at marker m.

`haploprob::Union{Nothing,Vector}`: posterior ancestral haplotype probability for each offspring
in each chromosome. haploprob[c][o] for offspring o in chromsome c is a sparse matrix,
with element (m,s) being the probability of haplotype s at marker m.

"""
mutable struct PolyAncestry
    # markermap [:marker,:chromosome,:position]
    markermap::Vector{DataFrame}
    # parentgeno: [chr][marker,parent]
    parentgeno::Vector{Matrix}
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
    # keys is a list of sub-population IDs.
    # each value is Dict again with keys = valent, parent, ancestralstate
    statespace::Dict{String,Dict{String,Array}}
    # valentprob[chr]: a list of matrix;
    # matrix has 3 columns: valent_index, valent_logl, valent_posteriorprob
    valentprob::Union{Nothing,Vector{Vector{Matrix}}}
    # genoprob[chr][offspring]: a sparse matrix with dimension nsnp x nstate
    genoprob::Vector{Vector{SparseMatrixCSC}}
    # haploprob[chr][offspring]: a sparse matrix with dimension nsnp x nstate
    haploprob::Union{Nothing,Vector{Vector{SparseMatrixCSC}}}
    function PolyAncestry(markermap,parentgeno,parentinfo,offspringinfo,designinfo,
        delmarker,correction,statespace,valentprob,genoprob,haploprob)
        dims = [length(markermap),size(parentgeno,1),size(genoprob,1)]
        if length(union(dims))!=1
            error("inconsistent number of linkage groups ",dims)
        end
        dims = [size.(markermap,1),size.(parentgeno,1),size.(first.(genoprob),1)]
        if length(union(dims))!=1
            error("inconsistent number of markers within linkage groups ",dims)
        end
        dims = [union(length.(genoprob)),[size(offspringinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of offspring ",dims)
        end
        dims = [union(size.(parentgeno,2)),[size(parentinfo,1)]]
        if length(union(dims))!=1
            error("inconsistent number of parents ",dims)
        end
        new(markermap,parentgeno,parentinfo,offspringinfo,designinfo,
            delmarker,correction,statespace,valentprob,genoprob,haploprob)
    end
end

"""
    savePolyAncestry(outfile,polyancestry,missingstring="NA",workdir=pwd())

saves polyancestry into outfile in CSV format.

# Positional arguments

`outfile::AbstractString`: file for saving polyancestry in the directory workdir.
The file consists of several tables: the first row of each table has two cells:
[PolyOrigin-PolyAncestry, nameoftable], and the rest rows denote a dataframe with column names.
The following is a list of table names and their descriptions.

1. `designinfo`: design information with columns being [population, parent1, parent2, ...].
   Matrix element D\\_{ij} denotes the nubmer of gametes contributed to each
   offspring in i-th population by j-th parent. D\\_{ij}=2 means that i-th population
   is produced by self-fertilization of j-th parent.

2. `parentinfo`: parent information with columns being [individual, ploidy].

3. `offspringinfo`: offspring information with columns being [individual, population, ploidy].

4. `valentlist`: list of possible bi- or multi-valent formations (chromosome pairings)
   for each population. The dataframe has columns [population, parent, valentindex, valent].

5. `valentprob`: posterior valent probability for each offspring in each chromosome.
   The dataframe has columns
   [chromosome, individual, population, valentindex, valent, loglike, valentprob].

6. `parentgeno`: phased parental genotypes. The dataframe has columns
   [marker, chromosome, position, parent1, parent2, ...].

7. `ancestralgenotype`: list of ancestral genotypes for each population. The dataframe
   has columns [population, parent, stateindex, state]

8. `genoprob`: marginal posterior probabilities for ancestral genotypes. For each offspring
   at each marker, the posterior probability vector is represented by a sparse vector
   in the form of I=>V: state index vector I and non-zero probability vector V such
   that I[k]=V[k] for k=1, ..., K. The elements of vector I or V are delimited by "|".

`polyancestry::PolyAncestry`: results of ancestral inference returned from [`polyOrigin`](@ref).

# Keyward arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for writing outfile.

"""
function savePolyAncestry(outfile::AbstractString,polyancestry::PolyAncestry;
    missingstring::AbstractString="NA",
    workdir::AbstractString = pwd())
    keyls = ["designinfo", "parentinfo", "offspringinfo", "delmarker","correction"]
    res = OrderedDict(keyls .=>
        [polyancestry.designinfo,polyancestry.parentinfo,polyancestry.offspringinfo,
         polyancestry.delmarker,polyancestry.correction])
    dfvalents,dfstates = statespace2df(polyancestry,missingstring)
    push!(res,"valentlist"=>dfvalents)
    push!(res,"valentprob"=>valentprob2df(polyancestry))
    push!(res,"parentgeno"=>parentgeno2df(polyancestry))
    push!(res,"ancestralgenotype"=>dfstates)
    dfprob = condprob2df(polyancestry,ishaploprob=false)
    push!(res,"genoprob"=>dfprob)
    # dfprob = condprob2df(polyancestry,ishaploprob=true)
    # push!(res,"haploprob"=>dfprob)
    outputfile2 = getabsfile(workdir,outfile)
    savedict2dlm(outputfile2,res,regionkey="PolyOrigin-PolyAncestry")
end

function parentgeno2df(polyancestry::PolyAncestry)
    nch = length(polyancestry.parentgeno)
    fhaplo=vcat([stringjoin.(polyancestry.parentgeno[ch],"|") for ch=1:nch]...)
    vcatmap = Matrix(vcat(polyancestry.markermap...))
    fhaplo2=hcat(vcatmap,fhaplo)
    parentid = polyancestry.parentinfo[!,:individual]
    parentinfo = polyancestry.parentinfo
    DataFrame(fhaplo2,[propertynames(polyancestry.markermap[1]); Symbol.(parentinfo[!,1])])
end

function statespace2df(polyancestry::PolyAncestry,missingstring::AbstractString)
    states=polyancestry.statespace
    popls = polyancestry.designinfo[!,:population]
    resvalents = []
    resstates = []
    for i=popls
        # valents0=states[i]["valent"]
        valents0 = [[unique(j) for j=i] for i = states[i]["valent"]]
        a=[string(i,"|",j) for i=1:size(valents0,1), j=1:size(valents0,2)]
        valents = vec(valents0)
        v = [ismissing(k) ? missingstring : stringjoin([stringjoin(stringjoin.(j,":"),"-") for j=k],"|") for k=valents]
        parent = stringjoin(states[i]["parent"],"|")
        parentindex = stringjoin(states[i]["parentindex"],"|")
        resv = hcat(repeat([i parentindex parent],length(v)),reshape(a,:),v)
        push!(resvalents,resv)
        s = stringjoin.(states[i]["groupstate"],"-")
        ns = length(s)
        ress = hcat(repeat([i parentindex parent],ns),1:ns,s)
        push!(resstates,ress)
    end
    dfvalents = DataFrame(vcat(resvalents...),Symbol.(["population","parentindex", "parent", "valentindex", "valent"]))
    dfstates = DataFrame(vcat(resstates...),Symbol.(["population","parentindex", "parent", "stateindex", "state"]))
    dfvalents,dfstates
end

function condprob2df_old(polyancestry::PolyAncestry;ishaploprob::Bool=true)
    condprob = ishaploprob ? polyancestry.haploprob : polyancestry.genoprob
    offls = polyancestry.offspringinfo[!,:individual]
    markermap= polyancestry.markermap
    chrls = [markermap[i][1,2] for i=1:length(markermap)]
    rowid=["chromosome", "individual", "nmarker", "ngenotype", "markerindex", "stateindex", "genoprob"]
    if ishaploprob
        rowid[4]="nhaplotype"
        rowid[7]="haploprob"
    end
    noff = length(offls)
    res = []
    for ch=1:length(chrls)
        for ind=1:noff
            prob = condprob[ch][ind]
            m,n = size(prob)
            I, J, V = findnz(prob)
            resrow = vcat([chrls[ch],offls[ind]],m,n,stringjoin.([I,J], "|"),stringjoin(V, "|"))
            push!(res,resrow)
        end
    end
    res2= permutedims(hcat(res...))
    DataFrame(res2,Symbol.(rowid))
end


function condprob2df(polyancestry::PolyAncestry;ishaploprob::Bool=true)
    condprob = ishaploprob ? polyancestry.haploprob : polyancestry.genoprob
    offls = polyancestry.offspringinfo[!,:individual]
    noff = length(offls)
    vcatmap = vcat(polyancestry.markermap...)
    nsnp = size(vcatmap,1)
    res= Matrix(undef,nsnp,noff+3)
    res[:,1:3]=Matrix(vcatmap)
    nchr=length(condprob)
    for ind=1:noff
        indpr=[begin
            prob = condprob[ch][ind]
            [begin
                I, V = findnz(prob[m,:])
                join([join(I,"|"),join(V,"|")],"=>")
            end for m = 1:size(prob,1)]
        end for ch=1:nchr]
        res[:,ind+3] =  vcat(indpr...)
    end
    rowid = [propertynames(vcatmap); Symbol.(offls)]
    DataFrame(res,rowid)
end

function viterbipath2df(polyancestry::PolyAncestry)
    if polyancestry.viterbipath == nothing
        return reshape(["nothing"],1,:)
    end
    res = [[HMM.toStringpath(HMM.toJumppath(i)) for i=eachcol(j)] for j=polyancestry.viterbipath]
    res2=[join(i,"|") for i=eachrow(hcat(res...))]
    offls = polyancestry.offspringinfo[!,:individual]
    DataFrame([offls res2],[:individual,:viterbipath])
end

function valentprob2df(polyancestry::PolyAncestry)
    statespace = polyancestry.statespace
    valentprob = polyancestry.valentprob
    offls = polyancestry.offspringinfo[!,:individual]
    offpopls = polyancestry.offspringinfo[!,:population]
    markermap= polyancestry.markermap
    chrls = [markermap[i][1,2] for i=1:length(markermap)]
    rowid=["chromosome", "individual", "population", "valentindex", "valent", "loglike", "valentprob"]
    noff = length(offls)
    res = []
    for chr=1:length(chrls)
        for off=1:noff
            vv = valentprob[chr][off]
            valentindex, logl, postprob = [stringjoin(vv[:,j],"|") for j=1:3]
            valent00 = statespace[offpopls[off]]["valent"]
            valent0=  [[unique(j) for j=i] for i =valent00]
            valent = valent0[vv[:,1]]
            valent2 = [stringjoin([stringjoin([stringjoin(j,":") for j=k],"-") for k=i],"&") for i=valent]
            valent3=stringjoin(valent2,"|")
            resrow=[chrls[chr],offls[off],offpopls[off],valentindex,valent3,logl,postprob]
            push!(res,resrow)
        end
    end
    res2= permutedims(hcat(res...))
    DataFrame(res2,Symbol.(rowid))
end

"""
    readPolyAncestry(genoprobfile, pedfile, missingstring="NA",workdir=pwd())

return a struct varaible with type PolyAncestry from  genoprobfile and pedfile
    in the directory workdir.

# Positional argument

`genoprobfile`: file storing conditional genotype probability.

`pedfile::AbstractString`: filename for pedigree information.

# Keyword arguments

`chrpairing::Integer=44`: chromosome pairing that has been used in producing
genoprobfile. chrpairing = 22 indicates only bivalent formations and 44 for
bivalent and quadrivalent formations.

`missingstring::AbstractString="NA"`: string code for missing value.

`delimchar::AbstractChar=','`:  text delimiter.

`commentstring::AbstractString="#"`: rows that begins with commentstring will be ignored.

`workdir::AbstractString = pwd()`: directory for reading genoprobfile and pedfile

"""
function readPolyAncestry(genoprobfile::AbstractString,
    pedfile::AbstractString;
    chrpairing::Integer=44,
    missingstring::AbstractString="NA",
    delimchar::AbstractChar=',',
    commentstring::AbstractString="#",
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    genoprobfile = getabsfile(workdir,genoprobfile)
    isfile(genoprobfile) || @error(string(genoprobfile," does not exist in workdir = ",workdir))
    pedfile = getabsfile(workdir,pedfile)
    isfile(pedfile) || @error(string(pedfile," does not exist in workdir = ",workdir))
    # imort pedfile
    parentinfo,offspringinfo,designinfo=readDesign(pedfile,
            delimchar=delimchar,
            commentstring=commentstring)
    parentid=parentinfo[!,:individual]
    offid=offspringinfo[!,:individual]
    nparent = length(parentid)
    noff = length(offid)
    # imort genoprobfile
    genoprob=CSV.read(genoprobfile,DataFrame; delim=delimchar,comment=commentstring,
        missingstring=missingstring)
    if size(genoprob,2) != 3+nparent+noff
        error(string(genoprobfile, " must contain #column = 3+#parent+#offspring = ",3+nparent+noff))
    end
    markermap = genoprob[!, 1:3]
    parsemarkermap!(markermap)
    typedind = string.(strip.(string.(names(genoprob)[4:end])))
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
    parentgeno = [Matrix(genoprob[i,parentcols]) for i=snpindexlist]
    parentgeno2=parseinputgeno(parentgeno,missingstring=missingstring)
    offcols=[get(rule,i,missing) for i=offid]
    offgenoprob = [Matrix(genoprob[i,offcols]) for i=snpindexlist]
    priorspace = getpriorstatespace(designinfo,parentinfo,chrpairing)
    statespace = getstatespace(priorspace)
    rule=Dict([key=>length(statespace[key]["groupstate"]) for key=keys(statespace)])
    nstatels=[get(rule,i,missing) for i=offspringinfo[!,:population]]
    offgenoprob2=[begin
        mtx=offgenoprob[ch]
        [sparse(hcat(parseprobcell.(mtx[:,off],nstatels[off])...)') for off=1:length(nstatels)]
    end for ch=1:length(offgenoprob)]
    markermap2=[DataFrame(i) for i=groupby(markermap,:chromosome)]
    verbose && @info string("data: #pop=", size(designinfo,1),
        ", #parent=",nparent, ", #offspring=",noff,
        ", #chr=",length(markermap2), ", #marker=",size(markermap,1))
    delmarker = getemptydf([String,String,Float64],[:marker,:chromosome,:position])
    colid = [:round,:marker,:chromosome,:parent,:old_genotype,:new_genotype, :old_nerr,:new_nerr]
    coltype = [Integer,String,String,String,String,String,Integer,Integer]
    correction = DataFrame(coltype,colid)
    polyancestry=PolyAncestry(markermap2,parentgeno2, parentinfo,offspringinfo,
        designinfo,delmarker,correction,statespace,nothing,offgenoprob2,nothing)
end

"""
    readPolyAncestry(ancestryfile, missingstring="NA",workdir=pwd())

return a struct varaible with type PolyAncestry from  ancestryfile in the
    directory workdir.

# Positional argument

`ancestryfile`: file storing polyancestry that is generated by [`savePolyAncestry`](@ref)
It is one of output files by [`polyOrigin`](@ref).

# Keyward arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for reading ancestryfile,

"""
function readPolyAncestry(ancestryfile::AbstractString;
    missingstring::AbstractString="NA",
    workdir::AbstractString=pwd())
    ancestryfile2 = getabsfile(workdir,ancestryfile)
    isfile(ancestryfile2) || @error(string(ancestryfile," does not exist in workdir = ",workdir))
    res = readdlm2dict(ancestryfile2)
    markermap, parentgeno = parseparentgeno(res["parentgeno"])
    parentinfo = res["parentinfo"]
    parentinfo[!,:individual]=string.(parentinfo[!,:individual])
    offspringinfo = res["offspringinfo"]
    offspringinfo[!,:individual]=string.(offspringinfo[!,:individual])
    offspringinfo[!,:population]=string.(offspringinfo[!,:population])
    designinfo = res["designinfo"]
    designinfo[!,:population]=string.(designinfo[!,:population])
    delmarker = res["delmarker"]
    for i=1:2
        delmarker[!,i]=string.(delmarker[!,i])
    end
    correction = res["correction"]
    for i=2:6
        correction[!,i]=string.(correction[!,i])
    end
    statespace=parsestatespace(res["ancestralgenotype"],res["valentlist"],missingstring)
    # priorspace = getpriorstatespace(designinfo,parentinfo,chrpairing)
    # statespace = getstatespace(priorspace)
    valentprob = parsevalentprob(res["valentprob"])
    # ancestralgenotype is named internally as groupstate
    rule=Dict([key=>length(statespace[key]["groupstate"]) for key=keys(statespace)])
    nstatels=[get(rule,i,missing) for i=offspringinfo[!,:population]]
    genoprob=parsecondprob(res["genoprob"],nstatels)
    polyancestry=PolyAncestry(markermap,parentgeno, parentinfo,offspringinfo,
        designinfo,delmarker,correction,statespace,valentprob,genoprob,nothing)
    sethaploprob!(polyancestry)
    polyancestry
end


function parseparentgeno(parentgeno)
    markermap = parentgeno[:,1:3]
    for i=1:2
        markermap[!,i]=string.(markermap[!,i])
    end
    groupdf = groupby(parentgeno, :chromosome)
    parentgeno = [begin
        haplo = split.(Matrix(i[!,4:end]),"|")
        map(x->parse.(Int,x),haplo)
    end for i = groupdf]
    markermap2=[DataFrame(i) for i=groupby(markermap,:chromosome)]
    markermap2, parentgeno
end

function parsecondprob_old(condprob)
    groupdf = groupby(condprob,1)
    [begin
        mtx=Matrix(groupdf[ch][!,3:7])
        mtx[:,3:end]=split.(mtx[:,3:end],"|")
        [sparse(parse.(Int,i[3]),parse.(Int,i[4]),parse.(Float64,i[5]),i[1],i[2]) for i=eachrow(mtx)]
    end for ch=1:length(groupdf)]
end

function parseprobcell(x::AbstractString,nstate::Integer)
    ls=split.(split(x,"=>"),"|")
    sparsevec(parse.(Int,ls[1]),parse.(Float64,ls[2]),nstate)
end


function parsecondprob(condprobdf::DataFrame,nstatels::AbstractVector)
    groupdf = groupby(condprobdf,2)
    noff = length(nstatels)
    [begin
        mtx=Matrix(groupdf[ch][:,4:end])
        [sparse(hcat(parseprobcell.(mtx[:,off],nstatels[off])...)') for off=1:noff]
    end for ch=1:length(groupdf)]
end

function parsevalentprob(valentprob)
    groupdf = groupby(valentprob,1)
    [begin
        # col 4,6,7= valuentindex, loglike,valentprob
        mtx=Matrix(groupdf[ch][!,[4,6,7]])
        mtx2=Matrix{Any}(split.(string.(mtx),"|"))
        for j=1:3
            mtx2[:,j]=map(x->parse.(Float64,x),mtx2[:,j])
        end
        [begin
            vp=Matrix{Real}(hcat(i...))
            vp[:,1]=Int.(vp[:,1])
            vp
        end for i=eachrow(mtx2)]
    end for ch=1:length(groupdf)]
end




function repeatmultivalent(zygotevalent::AbstractVector)
    [vcat([repeat([j],div(length(j), 2)) for j=i]...) for i=zygotevalent]
end

function parsestatespace(ancestralstate,valent,missingstring::AbstractString)
    # dfvalent:["population","parentindex", "parent", "valentindex", "valent"]))
    # dfstate: ["population","parentindex", "parent", "stateindex", "state"]))
    dfstate=groupby(ancestralstate,1)
    dfvalent=groupby(valent,1)
    Dict([begin
        popid = string(dfstate[pop][1,1])
        parentindex= map(x->parse.(Int,x),split.(string(dfstate[pop][1,2]),"|"))
        parent= split(dfstate[pop][1,3],"|")
        state = map(x->parse.(Int,x),split.(dfstate[pop][!,5],"-"))
        vv0=dfvalent[pop][:,5]
        bool = vv0 .== missingstring
        vv0[bool] .= "0"
        vv1 = map(x->split.(x,"-"),split.(vv0,"|"))
        vv2=[map(x->map(y->parse.(Int,y),split.(x,":")),v) for v=vv1]
        if 1 in bool
            vv2 = Vector{Union{Missing,eltype(vv2)}}(vv2)
            vv2[bool] .= missing
        end
        vv3=[ismissing(i) ? i : repeatmultivalent(i) for i=vv2]
        index=hcat(map(x->parse.(Int,x),split.(dfvalent[pop][!,4],"|"))...)
        vv4=reshape(vv3,max(index[1,:]...),:)
        popid=>Dict(["valent"=>vv4,"parent"=>parent,"parentindex"=>parentindex,"groupstate"=>state])
    end for pop = 1:length(dfstate)])
end

function sethaploprob!(polyancestry::PolyAncestry)
    dict=Dict{String,SparseMatrixCSC}()
    polyancestry.haploprob = [begin
        [begin
            offpop = polyancestry.offspringinfo[off,:population]
            if !haskey(dict,offpop)
                genotypes = polyancestry.statespace[offpop]["groupstate"]
                g2h=geno2haplo(genotypes)
                push!(dict,offpop =>sparse(g2h))
            end
            mtx = dict[offpop]
            chgeno[off] * mtx
        end for off=1:length(chgeno)]
    end for chgeno = polyancestry.genoprob]
end

function geno2haplo(genotypes::AbstractVector)
    nrow =size(genotypes,1)
    ncol = max(vcat(genotypes...)...)
    mtx = zeros(nrow,ncol)
    ploidy = length(genotypes[1])
    for i=1:nrow for j=1:ploidy
            mtx[i,genotypes[i][j]] += 1
        end
    end
    mtx ./= ploidy
    sparse(mtx)
end


"""
    savegenoprob(outfile,polyancestry,missingstring="NA",workdir=pwd())

saves genoprob of the struct polyancestry into outfile, phased parental
    genotypes being also included.

# Positional arguments

`outfile::AbstractString`: filename for saving results.

`polyancestry::PolyAncestry`: a struct returned by [`polyOrigin`](@ref).

# Keyward arguments

`missingstring::AbstractString="NA"`: string code for missing value.

`workdir::AbstractString = pwd()`: directory for writing outfile.

"""
function savegenoprob(outfile::AbstractString,polyancestry::PolyAncestry;
    missingstring::AbstractString="NA",
    workdir::AbstractString = pwd())
    parentgenodf = parentgeno2df(polyancestry)
    genoprobdf = condprob2df(polyancestry,ishaploprob=false)
    res=hcat(parentgenodf,genoprobdf[!,4:end])
    outputfile2 = getabsfile(workdir,outfile)
    CSV.write(outputfile2,res; missingstring=missingstring,
        header=true,append=false)
end
