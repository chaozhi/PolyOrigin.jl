"""
    plotCondprob(polyancestry,offspring=nothing,ishaploprob=true,
        colorgradient = ColorGradient([:white,:blue,:red]),
        boundaryline = (1.5,:dot,:black),
        truemarker=(:star, 5, 0.5,:gray,stroke(:gray)),
        truegeno=nothing)

plot heatmap for conditional probability.

# Positional arguments

`polyancestry::PolyAncestry`: polyancestry returned from [`polyOrigin`](@ref).

# Keyword arguments

`offspring::Union{Nothing,Integer}=nothing`: offsprign index and random offspring
by default.

`colorgradient::ColorGradient=ColorGradient([:white,:blue,:red])`: color gradient
for heatmap

`left_margin = :match`: Specifies the extra padding to the left of the subplot

`bottom_margin = :match`: Specifies the extra padding to the bottom of the subplot

`boundaryline=(1.5,:dot,:black)`: vertical lines for chromosome boundaries.

`truemarker=(:star, 5, 0.5,:gray,stroke(:gray))`: scatter markers for
true ancestral states.

`truegeno::Union{Nothing,NamedTuple}`: provides true parental origins. It is generated
by [`readTruegeno!`](@ref).

"""
function plotCondprob(polyancestry::PolyAncestry;
    offspring::Union{Nothing,Integer}=nothing,
    ishaploprob::Bool=true,
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),
    left_margin = :match,
    bottom_margin = :match,
    boundaryline = (1.5,:dot,:black),
    truemarker=(:star, 5, 0.5,:gray,stroke(:gray)),
    truegeno::Union{Nothing,NamedTuple}=nothing)
    chridls= [i[1,:chromosome] for i=polyancestry.markermap]
    condprob = ishaploprob ? polyancestry.haploprob : polyancestry.genoprob
    if offspring==nothing
        noff = size(condprob[1],1)
        off = rand(1:noff)
    else
        off=offspring
    end
    offid = polyancestry.offspringinfo[off,:individual]
    offid2 = string(offid,polyancestry.offspringinfo[off,:isoutlier] === true ? "(outlier)" : "")
    prob = [i[off] for i=condprob]
    heat=Matrix(vcat(prob...)')
    ibdmap=heatmap(heat; c=colorgradient,left_margin, bottom_margin)
    len = size.(prob,1)
    nstate= size(heat,1)
    x0=accumulate(+,len[1:end])
    pushfirst!(x0,0)
    xy=reduce(hcat,[[[i,0] [i,nstate+0.75] [NaN,NaN]] for i=x0])
    chrlabx=xy[1,1:3:end]
    chrlabx=(chrlabx[1:end-1] .+ chrlabx[2:end]) ./ 2
    annotate!(ibdmap,[(chrlabx[i],ishaploprob ? nstate + 0.75 : nstate+1.5,Plots.text(chridls[i],12,:black)) for i=1:length(chrlabx)])
    ylabel = ishaploprob ? "Haplotype" : "Genotype"
    titlestem = ishaploprob ? "Haplotype probablility" : "Genotype probablility"
    plot!(ibdmap,xy[1,:],xy[2,:],
        line = boundaryline,
        legend=false,
        xlabel="SNP index",
        ylabel=ylabel,
        title=string(titlestem, " for ",off,"th offspring =",offid2),
        titlefont = font("Times", 12)
    )
    if !isnothing(truegeno)
        offtrue=vcat([i[:,off] for i=truegeno.offspringgeno]...)
        if ishaploprob
            xy=vcat([[[i,j] for j=offtrue[i]] for i=1:length(offtrue)]...)
            xy2=reduce(hcat,xy)
            scatter!(ibdmap,xy2[1,:],xy2[2,:],
                marker=truemarker,
                legend=false)
        else
            @error string("TODO")
        end
    end
    if ishaploprob
        yticks = ([1:nstate;], [string("h",i) for i=1:nstate])
        plot!(ibdmap,yticks = yticks)
    end
    ibdmap
end


function saveProbPlot(
    polyancestry::PolyAncestry;
    ishaploprob::Bool = true,
    nplot_subpop::Integer = 10, 
    colorgradient::ColorGradient = cgrad([:white, :blue, :red]),
    left_margin = :match,
    bottom_margin = :match,
    boundaryline = (1.5, :dot, :black),
    truemarker = (:star, 5, 0.5, :gray, stroke(:gray)),
    truegeno::Union{Nothing,NamedTuple} = nothing,
    outstem::AbstractString = "outstem",
    workdir::AbstractString = pwd(),
)
    figdir = joinpath(workdir, outstem * "_plots")
    isdir(figdir) || mkdir(figdir)
    offinfo = polyancestry.offspringinfo
    popls = offinfo[!,:population]
    offls = sort(reduce(vcat, [begin 
        suboffls = findall(popls .== i )
        if nplot_subpop >= length(suboffls)
            suboffls
        else
            suboffls[1:nplot_subpop]
        end
    end for i in unique(popls)]))
    iiout = findall(offinfo[!, :isoutlier])
    iinonout = setdiff(offls, iiout)
    for k = 1:2
        ii = k == 1 ? iiout : iinonout
        outid = k == 1 ? "probplot_outlier_" : "probplot_"
        for i in ii
            offid = polyancestry.offspringinfo[i, :individual]
            offid2 = replace(replace(offid,"/"=>"_"),"\\"=>"_")
            fn = joinpath(
                figdir,
                string(outid, i, "th_offspring_", offid2, ".png"),
            )
            plotCondprob(polyancestry; offspring = i,
                ishaploprob,colorgradient,left_margin,bottom_margin,
                boundaryline,truemarker,truegeno)
            savefig(fn)
        end
    end
    figdir
end

"""
    animCondprob(polyancestry,fps=1,outfile=nothing,kewargs...)

animation for plots of conditional probability.

# Positional arguments

`polyancestry::PolyAncestry`: polyancestry returned from [`polyOrigin`](@ref).

# Keyword arguments

`fps::Real=1`: number of frames per seconds.

`outfile::AbstractString="condprob.gif"`: output file for saving animation.

see [`plotCondprob`](@ref) for keyargs.

"""
function animCondprob(polyancestry::PolyAncestry;
    fps::Real=1,    
    ishaploprob::Bool=true,
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),
    left_margin = :match,
    bottom_margin = :match,
    boundaryline = (1.5,:dot,:black),
    truemarker=(:xcross, 3, 0.5,:gray,stroke(:gray)),
    truegeno::Union{Nothing,NamedTuple}=nothing,
    outfile::Union{Nothing,AbstractString}=nothing)
    noff=size(polyancestry.offspringinfo,1)
    anim = Animation()
    n1 = ceil(Int,fps)
    n2 = fps == 0 ? 1 : ceil(Int,1/fps)
    @time for off=1:noff
        plotCondprob(polyancestry; offspring=off,
            ishaploprob, colorgradient,
            left_margin, bottom_margin,
            boundaryline, truemarker,truegeno
        )
        for i=1:n2
            frame(anim)
        end
    end
    if outfile == nothing
        gif(anim,fps=n1)
    else
        gif(anim,outfile,fps=n1)
    end
end

#######################################################################################

"""
    plotMapComp(mapx, mapy;
        boundaryline = (1.5,:dot,:black),
        plotmarker=(:star, 1, 0.5,:blue,stroke(:blue)),
        isannotate= true,
        maplabels=["mapx position", "mapy position"])

plot postions of mapx vs those of mapx.

# Positional arguments

`mapx::Vector{DataFrame}`: marker map for all chromosomes.

`mapy::Vector{DataFrame}`: comparing map for all chromosomes.

# Keyword arguments

`boundaryline=(1.5,:dot,:black)`: vertical lines for chromosome boundaries.

`isannotate::Bool=true`: if ture,  annotate chromosome ID and kendall correlation.

`mmaplabels::Union{Nothing,AbstractString}=nothing`: labels of comparing marker maps.

"""
function plotMapComp(mapx::Vector{DataFrame},mapy::Vector{DataFrame};
    boundaryline = (1.5,:dot,:black),    
    markersize = 3, 
    markershape = :circle, 
    fontsize = 12,     
    isannotate::Bool=true,
    cordigits::Integer=3,
    maplabels::AbstractVector=["map1 position", "map2 position"],
    plotargs...)
    # reordereing linkage groups of mapX
    truelg = findtruelg(mapx,mapy)
    mapx = deepcopy(mapx)
    mapx = [reduce(vcat,mapx[i]) for i in truelg]
    nchr=length(mapx)
    snps = [intersect(mapx[chr][!,:marker],mapy[chr][!,:marker]) for chr=1:nchr]
    #
    res=[begin
        df=filter(x->in(x[1],snps[chr]),mapx[chr])
        rule = Dict(mapy[chr][!,:marker] .=> 1:size(mapy[chr],1))
        pos = [get(rule,i,missing) for i=df[!,:marker]]
        isnonmiss= .!ismissing.(pos)
        posnon = pos[isnonmiss]
        xls = df[isnonmiss,:position]
        yls =mapy[chr][posnon,:position]
        tau=round(corkendall(posnon,collect(1:length(posnon))),digits=cordigits)
        [xls,yls,tau]
    end for chr=1:nchr]
    len = [i[1][end] for i=res]
    len = accumulate(+,len)
    for i=2:nchr
        res[i][1] .+= len[i-1]
    end
    ymax = max([max(i[2]...) for i=res]...)
    chridls=[i[1,:chromosome] for i =mapx]
    marker_colors = distinguishable_colors(nchr, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    fig = plot([0,0],[0,ymax]; line = boundaryline)
    for chr=1:nchr
        c = marker_colors[chr]
        plot!(fig,res[chr][1],res[chr][2];            
            marker=(markershape, markersize, 0.5, c,stroke(c)),            
            legend=false,
        )
        if isannotate
            labx=mean(res[chr][1][[1,end]])
            lab0 = chr==1 ? "\nr=" : "\n"
            lab = string(chridls[chr],lab0,res[chr][3])
            annotate!(fig,[(labx,ymax,text(lab,:left,font(fontsize)))])
        end
        boundaryx=res[chr][1][end]
        plot!(fig,[boundaryx,boundaryx],[0,ymax];
            line = boundaryline)
    end
    plot!(fig; 
        xlabel=maplabels[1],ylabel=maplabels[2],
        titlefontsize = fontsize,          
        tickfontsize = fontsize, 
        guidefontsize = fontsize, 
        legendfontsize= fontsize,
        left_margin = 10Plots.mm, 
        bottom_margin = 10Plots.mm, 
        size = (1200, 400), 
        plotargs...
    )
end


function findtruelg(truemarkermap::Vector{DataFrame}, estmarkermap::Vector{DataFrame})
    snpinter = [length(intersect(i[!,:marker],j[!,:marker]))
        for i=truemarkermap, j in estmarkermap]
    rowmax = maximum(snpinter,dims=2)[:,1]
    truelg = [begin
        pos = findall(snpinter[:,chr] .> 1)
        pos2 = pos[snpinter[pos,chr] .>= rowmax[pos]]
        pos2 = union(pos2,argmax(snpinter[:,chr]))
    end for chr in 1:size(snpinter,2)]
    pos = findall(length.(truelg) .> 1)
    if !isempty(pos)
        ids = [i[1,:chromosome] for i in truemarkermap[pos]]
        @warn string("1st markermap chromosomes merged: ", truelg[pos])
    end
    ls = vcat(truelg...)
    for chr in 1:length(truemarkermap)
        rep = sum(ls .== chr)
        if rep>1
            id = truemarkermap[chr][1,:chromosome]
            msg =string("1st markermap ", id)
            msg *= string(" is splited into ", rep, " chrs in 2nd markermap")
            @warn msg
        end
    end
    pos = findall(isempty.(truelg))
    if !isempty(pos)
        ids = [i[1,:chromosome] for i in estmarkermap[pos]]
        @error string("2nd markermap chromosomes not founder: ", ids)
    end
    truelg
end

function plotMapComp(mapxfile::AbstractString,mapyfile::AbstractString,
    moremapfiles::Vararg{AbstractString};
    boundaryline = (1.5,:dot,:black),
    plotmarker=(:star, 2, 0.5,:blue,stroke(:blue)),
    isannotate::Bool=true,
    cordigits::Integer=3,
    maplabels::Union{Nothing,AbstractVector}=nothing,
    delimchar::AbstractChar=',',
    missingstring::AbstractString="NA",
    commentstring::AbstractString="#",
    workdir::AbstractString = pwd())
    filels = [mapxfile,mapyfile, moremapfiles...]
    filels = [getabsfile(workdir,mapfile) for mapfile = filels]
    b = isfile.(filels)
    all(b) || error(string("files do not exist: ",filels[.!b]))
    mapls = [begin
        mapdf=CSV.read(mapfile, DataFrame; delim=delimchar,
            comment=commentstring,missingstring)
        if size(mapdf,2) < 3
            error(string(mapfile, " must contain #column >= 3"))
        end
        b = .!ismissing.(mapdf[!,2])
        mapdf = mapdf[b,1:3]
        parsemarkermap!(mapdf)
        [DataFrame(i) for i=groupby(mapdf,:chromosome)]
    end for mapfile = filels]
    labells = [string("map",i, " position") for i=1:length(mapls)]
    if !isnothing(maplabels)
        n = min(length(maplabels),length(mapls))
        labells[1:n] .= maplabels[1:n]
    end
    n = length(mapls)
    res = [plotMapComp(mapls[i],mapls[j];
                boundaryline, plotmarker, isannotate,cordigits,
                maplabels=labells[[i,j]])
                for i=1:n for j=(i+1):n]
    length(res)==1 ? res[1] : plot(res...,layout=(n,1))
end

function plotdesign(polygeno::PolyGeno;
    curves = false,
    nodecolor = 1,
    plotsize = nothing,
    fontsize = 7,
    method = :circular,
    nodesize = nothing,
    self_edge_size = 0.1,
    edgecolor = :black,
    edgestyle = :solid
    )
    offinfo = polygeno.offspringinfo
    designinfo = polygeno.designinfo
    npop, nf= size(designinfo)
    nf -= 1
    adj = zeros(Int,nf,nf)
    # edgelabel_dict = Dict()
    edgelabel_mtx = ["" for i=1:nf, j=1:nf]
    for i=1:size(designinfo,1)
        jj = findall(Vector(designinfo[i,2:end]) .> 0)
        pop = designinfo[i,1]
        popsize = sum(offinfo[!,:population] .== pop)
        lab = string("N(",pop, ")=",popsize)
        if length(jj)==1
            s = jj[1]
            adj[s,s] =1
            # edgelabel_dict[(s,s)] = lab
            edgelabel_mtx[s,s]=lab
        elseif length(jj)==2
            s, d = jj
            adj[s,d] =1
            adj[d,s] =1
            # edgelabel_dict[(s,d)] = lab
            # edgelabel_dict[(d,s)] = lab
            edgelabel_mtx[s,d]=lab
            # edgelabel_mtx[d,s]=lab
        else
            error(string("wrong parentage: ",designinfo[i,:]))
        end
    end
    if isnothing(nodecolor)
        nodecolor = vcat([3 for i=1:nf],[1 for i=nf+1:n])
    end
    nodename = names(designinfo)[2:end]
    isnothing(plotsize) && (plotsize=(200,200) .* sqrt(nf))
    if isnothing(nodesize)
        nodename = names(designinfo)[2:end]
        labsize=sum(length.(nodename))/nf
        nodesize = 0.2/sqrt(labsize)
    end
    graphplot(adj;
        names=nodename,
        nodeshape=:circle,
        edgelabel=edgelabel_mtx,
        curves,self_edge_size,
        fontsize, method,
        nodesize, nodecolor,
        edgecolor, edgestyle
        )
    plot!(size=plotsize)
end

#######################################################################################

function cal_valentcount(polyancestry::PolyAncestry; valentthresh::Real=0.2)
    parentls = polyancestry.parentinfo[!,:individual]
    chrls = [i[1,:chromosome] for i=polyancestry.markermap]
    popls = keys(polyancestry.statespace)
    inddict = Dict([pop =>findall(polyancestry.offspringinfo[!,:population] .== pop) for pop=popls])
    # TODO: parent specific ploidy/valentname
    ploidy =  polyancestry.offspringinfo[1,:ploidy]
    valents = polyancestry.statespace[first(popls)]["valent"]
    # valentname = [i[1][1] == 1:ploidy ? join(i[1][1],":") : join(join.(i[1],":"),"-") for i=valents[:,1]]
    valentname = [i[1][1] == 1:ploidy ? join(i[1][1],":") : join(join.(unique(i[1]),":"),"-") for i=valents[:,1]]
    valentindices = collect(eachindex(IndexCartesian(),valents))
    resls =[begin
        chrvalentprob = polyancestry.valentprob[chr]
        res = zeros(length(parentls),length(valentname))
        for pop in popls
            offls = inddict[pop]
            for off in offls
                offprob = chrvalentprob[off]
                ls = offprob[:,3] # 3rd column = poster prob
                b = ls .>= valentthresh*maximum(ls)
                vvls = Tuple.(valentindices[offprob[b, 1]])
                wwls = normalize(offprob[b,3],1)
                pp = polyancestry.statespace[pop]["parentindex"]
                if length(pp)==1
                    p1=p2=only(pp)
                elseif length(pp)==2
                    p1,p2 = pp
                else
                    error(string("wrong parents ",pp))
                end
                for i in eachindex(vvls, wwls)
                    res[p1,vvls[i][1]] +=wwls[i]
                    res[p2,vvls[i][2]] +=wwls[i]    
                end
            end
        end
        # reduce(hcat,[i ./ sum(i) for i in eachrow(res)])'
        res
    end for chr in eachindex(polyancestry.valentprob)]
    vvfrac = round.(vcat(resls...),digits=6)
    resdf = DataFrame(vvfrac,string.("valent", valentname))
    df0 = DataFrame(["chromosome"=>repeat(chrls,inner=length(parentls)),
        "parent"=>repeat(parentls,outer=length(chrls))])
    resdf2 = hcat(df0,resdf)
    resdf2
end


# valentthresh has little effect when chromosome length >= 200cM
# valentthresh = 0.2 ~ 0.4  has similar effects. 
# when chrlength = 100cM, valentthresh = 0.1 increases quadrivalent eestimtion by around 15%, 
# when chrlength = 100cM, valentthresh = 0.5 deccrease quadrivalent eestimtion by around 15%, 
# when chrlength = 100cM (M=100, N=200), quadrivalent=0.0: estimated quadrivalent = 0.085, 0.022 when valentthresh = 0.1, 0.2
# when chrlength = 200cM (M=100, N=200), quadrivalent=0.0: estimated quadrivalent = 0.023, 0.014 when valentthresh = 0.1, 0.2

function cal_valentsum(polyancestry::PolyAncestry; 
    valentthresh::Real=0.2, prefstep::Real=0.001)
    valentcount = PolyOrigin.cal_valentcount(polyancestry; valentthresh)
    col_1st_valent = 3
    valentnames = replace.(names(valentcount)[col_1st_valent:end],"valent"=>"")
    valentls = [[parse.(Int, j) for j in split.(split(i,"-"),":")] for i in valentnames]
    # group
    dfgroup = valentcount[:,1:col_1st_valent-1]
    groupls = [string("valentgroup", join(length.(i),"-")) for i in valentls]
    groupset = unique(groupls)
    for g = groupset
        cols = findall(groupls .== g) .+ (col_1st_valent - 1)
        insertcols!(dfgroup, size(dfgroup, 2)+1, g => round.(sum.(eachrow(valentcount[:,cols])),digits=6))
    end    
    # pair
    dfpair = valentcount[:,1:col_1st_valent-1]
    maxvalent = [maximum(length.(i)) for i in valentls]
    pairls = unique(reduce(vcat, valentls[maxvalent .== 2]))
    for p = pairls
        cols = findall([in(p, i) && all(length.(i) .== 2) for i in valentls]) .+ (col_1st_valent - 1)
        pname = string("pair", join(p, ":"))
        insertcols!(dfpair, size(dfpair, 2)+1, pname => round.(sum.(eachrow(valentcount[:,cols])),digits=6))
    end    
    prefpair = infer_prefpair_bivalentonly(valentcount; pstep=prefstep)
    OrderedDict("valentcount"=>valentcount, "prefpairing" => prefpair, 
        "valentfreq"=>normalize_valent(valentcount), 
      "valentgroupcount"=>dfgroup, "valentgroupfreq"=>normalize_valent(dfgroup), 
      "paircount"=>dfpair,"pairfreq"=>normalize_valent(dfpair))
end

function normalize_valent(dffreq::AbstractDataFrame)
    df = copy(dffreq)
    for row in eachrow(df)
        # first two columns: chromosome and parent
        row[3:end] .= round.(Vector(row[3:end]) ./ sum(row[3:end]),digits=6)
    end
    df
end


function get_marker_shapels(nshape::Integer; seed=1234)
    mshapels = Plots.supported_markers()
    setdiff!(mshapels, [:auto, :none,:x,:pixel])
    rng = MersenneTwister(seed)
    rand(rng, mshapels,nshape)
end


function plot_valentfreq(valentfreq::AbstractDataFrame; 
    fontsize = 14, 
    markersize = 8,
    markerseed = 1234,  
    plotkeyargs...)
    df = copy(valentfreq)
    col_1st_valent = 3
    valentnames = replace.(names(df)[col_1st_valent:end],"valent"=>"")
    valentls = [[parse.(Int, j) for j in split.(split(i,"-"),":")] for i in valentnames]
    groupls = [string("valentgroup", join(length.(i),"-")) for i in valentls]
    groupset = unique(groupls)
    keepat!(groupset, occursin.("-",groupset))
    for vgname in groupset
        cols = findall(groupls .== vgname) .+ (col_1st_valent - 1)
        cols = vcat(1:col_1st_valent-1,cols)
        subdf= df[!,cols]
        for row in eachrow(subdf)        
            row[col_1st_valent:end] .= Vector(row[col_1st_valent:end]) ./ sum(row[col_1st_valent:end])
        end    
    end
    ymax = maximum(Matrix(df[:, col_1st_valent:end]))    
    nchr = length(unique(df[!,1]))
    mshapels = get_marker_shapels(nchr; seed=markerseed)
    res = [begin     
        cols = findall(groupls .== vgname) .+ (col_1st_valent - 1)
        cols = vcat(1:col_1st_valent-1,cols)        
        subdf = df[!,cols]
        ctgls = replace.(names(subdf)[col_1st_valent:end],"valent"=>"")
        tickls = (collect(eachindex(ctgls)),ctgls)    
        gdf = groupby(subdf,:parent)
        [begin 
            legendls=gdf[i][!,1]            
            yls = collect.(eachrow(Matrix(gdf[i][!,col_1st_valent:end])))        
            title = string(PolyOrigin.tansform_valentgroupname(vgname),", parent = ",  gdf[i][1,:parent])
            colorls = distinguishable_colors(length(yls), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)                        
            fig = plot(;yrange=(0,1.05*ymax))
            for chr in eachindex(yls)
                c = colorls[chr]    
                plot!(fig, yls[chr];
                    labels=legendls[chr],
                    legend = :outerright, 
                    linecolor = c, 
                    title,                     
                    xticks = tickls, 
                    xrotation = 90, 
                    marker=(mshapels[chr], markersize, 0.5, c,stroke(c)),           
                    titlefontsize = fontsize,          
                    tickfontsize = fontsize, 
                    guidefontsize = fontsize, 
                    legendfontsize= fontsize,  
                    xlabel = "", 
                    ylabel = "Frequency", 
                    plotkeyargs...           
                )    
            end        
            plot!(x->1/length(ctgls),line=(:gray, :dash),labels="")
            fig
        end for i in eachindex(gdf)]
    end for vgname in groupset]
    res
end

function tansform_valentgroupname(vgname::AbstractString)    
    if vgname == "valentgroup2-2" 
        vgname = "double-bivalent"
    elseif vgname == "valentgroup4" 
        vgname = "quadrivalent"
    elseif vgname == "valentgroup2-2-2" 
        vgname = "triple-bivalent"    
    elseif vgname == "valentgroup4-2" 
        vgname = "quadrivalent+bivalent"
    elseif vgname == "valentgroup6"
        vgname = "hexavalent"
    else
        @error string("unexpected valenggroup=",vgname)
    end
end

function plot_valentgroupfreq(valentgroupfreq::AbstractDataFrame; 
    fontsize = 14, 
    plotkeyargs...)
    df = valentgroupfreq    
    gdf = groupby(df,:parent)
    parentls = [i[1,:parent] for i in gdf]
    chrls = gdf[1][!,:chromosome]
    ctg = repeat(parentls, inner = length(chrls))
    name = repeat(chrls, outer = length(parentls))
    [begin
        y = reduce(hcat,[df[!,col] for df in gdf])
        vgname = names(valentgroupfreq)[col]
        vgname = tansform_valentgroupname(vgname)
        ylabel = string("Freq of ", vgname)
        ymax = col == 3 ? maximum(Vector(df[!, 3])) : maximum(Matrix(df[!, 4:end]))
        barargs = (yrange = (0,1.1*ymax),
            title = string("<",vgname,">=",round(mean(y),digits=4)),
            lw = 0,                                   
            framestyle = :box, 
            titlefontsize = fontsize,          
            tickfontsize = fontsize, 
            guidefontsize = fontsize, 
            legendfontsize= fontsize)      
        if length(unique(name)) == 1
            xls = ctg
            yls = vec(y)
            bar(xls, yls;         
                labels=first(name),
                xlabel = "Parent", 
                ylabel,    
                barargs...,
                plotkeyargs...,            
            )
        else
            groupedbar(name, y; 
                group = ctg, 
                xlabel = "Linkage group", 
                ylabel,         
                barargs...,      
                plotkeyargs...,            
            )
        end
    end for col in 3:size(valentgroupfreq,2)]
end

function plot_pairfreq(pairfreq::AbstractDataFrame; 
    fontsize = 14, 
    markersize = 8, 
    markerseed=1234,
    plotkeyargs...)
    df = pairfreq
    col_1st_valent = 3        
    ymin = minimum(Matrix(df[:, col_1st_valent:end]))
    ymax = maximum(Matrix(df[:, col_1st_valent:end]))
    ctgls = names(df)[col_1st_valent:end]
    tickls = (collect(eachindex(ctgls)),ctgls)
    npair = length(ctgls)
    groupname = npair == 6 ? "double-bivalent" : (npair == 15 ? "triple-bivalent" : "")
    gdf = groupby(df,:parent)
    nchr = length(unique(df[!,1]))
    mshapels = get_marker_shapels(nchr; seed=markerseed)
    [begin 
        legendls=gdf[i][!,1]        
        yls = collect.(eachrow(Matrix(gdf[i][!,col_1st_valent:end])))        
        title = string(groupname, ", parent = ",  gdf[i][1,:parent])
        colorls = distinguishable_colors(length(yls), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)                        
        fig = plot(; 
            xticks = tickls,
            xrotation = 90,
            yrange = (0.5*ymin, 1.1*ymax),
            title,
            titlefontsize = fontsize+2,          
            tickfontsize = fontsize, 
            guidefontsize = fontsize, 
            legendfontsize= fontsize,  
            xlabel = "", 
            ylabel = "Frequency", 
            legend = :outerright, 
            left_margin=10Plots.mm,
            bottom_margin=10Plots.mm,
            plotkeyargs...
        )
        for chr in eachindex(yls)
            c = colorls[chr]    
            plot!(fig, yls[chr], 
                labels=legendls[chr],
                linecolor = c, 
                marker=(mshapels[chr], markersize, 0.5, c,stroke(c)),                                           
            )    
        end 
        plot!(x->1/length(ctgls),line=(:gray, :dash),labels="")
    end for i in eachindex(gdf)]
end


function plot_prefpair_posterior(prefpairdf::AbstractDataFrame; 
    fontsize = 14, 
    linewidth = 2, 
    plotkeyargs...)
    namels = names(prefpairdf)
    col_1st_post = findfirst(x->occursin(r"^posteriorPDF",x), namels)
    xls = parse.(Float64,replace.(namels[col_1st_post:end],"posteriorPDF"=>""))
    xmaxls = []
    gdf = groupby(prefpairdf,:chromosome)
    gls = [begin 
        fig = plot()
        colorls = distinguishable_colors(size(df,1), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)                        
        xmax = 0.0
        for i in 1:size(df, 1)
            yls = Vector(df[i, col_1st_post:end])            
            plot!(fig, xls, yls; linecolor=colorls[i], linewidth,
                label=df[i,:parent]
            )
            reverse!(yls)
            k = findfirst(>=(1e-2),yls)
            k2 = 1+length(yls) - k 
            xls[k2] > xmax && (xmax = xls[k2])
        end
        push!(xmaxls,xmax)
        plot!(fig, 
            xlabel = string("Pref pairing (chr=", df[1,:chromosome],")"), 
            ylabel = "PDF",
            titlefontsize = fontsize,          
            tickfontsize = fontsize, 
            guidefontsize = fontsize, 
            legendfontsize= fontsize,  
            plotkeyargs...
        )
    end for df in gdf]
    xmax = max(0.5,maximum(xmaxls))
    for g in gls
        plot!(g, xrange=(0, xmax))
    end
    gls
end

function plot_prefpair(prefpairdf::AbstractDataFrame; 
    fontsize = 14, 
    linewidth = 2, 
    plotkeyargs...)
    namels = names(prefpairdf)
    colerr1 = findfirst(==("quantile0.025"),namels)
    colerr2 = findfirst(==("quantile0.975"),namels)
    gdf = groupby(prefpairdf,:parent)
    parentls = [i[1,:parent] for i in gdf]
    chrls = gdf[1][!,:chromosome]
    ctg = repeat(parentls, inner = length(chrls))
    name = repeat(chrls, outer = length(parentls))
    [begin
        y = reduce(hcat,[df[!,col] for df in gdf])
        if col == 3
            y .*= -1
            ylabel = " - log10 Pr(pref<=0.05)" 
            ymax = 1.1*maximum(y) 
            yerror = nothing
        else
            ylabel =  "Estimate of pref pairing"
            yerr1 = reduce(hcat,[df[!,colerr1] for df in gdf])
            yerr2 = reduce(hcat,[df[!,colerr2] for df in gdf])
            ymax = 1.1* maximum(yerr2)
            yerror = (vec(y .- yerr1), vec(yerr2 .- y))
        end    
        barargs = (markerstrokecolor = :black, 
            lw = linewidth,                                   
            framestyle = :box, 
            titlefontsize = fontsize,          
            tickfontsize = fontsize, 
            guidefontsize = fontsize, 
            legendfontsize= fontsize)                 
        if length(unique(name)) == 1
            xls = ctg
            yls = vec(y)
            bar(xls, yls;         
                labels=first(name),
                xlabel = "Parent", 
                ylabel,    
                yrange = (0,ymax),  
                yerror, 
                barargs...,
                plotkeyargs...,            
            )
        else
            groupedbar(name, y; 
                group = ctg, 
                bar_position = :dodge,
                xlabel = "Linkage group", 
                ylabel,         
                yrange = (0,ymax),          
                yerror,    
                barargs..., 
                plotkeyargs...,            
            )
        end
    end for col in [4,3]]
end

function plot_valentsum(valentsum::AbstractDict;     
    workdir::AbstractString=pwd(), 
    outstem::AbstractString="outstem",
    io::Union{IO, Nothing}=nothing, 
    verbose=true, 
    fontsize = 14, 
    linewidth = 3, 
    markersize = 8, 
    markerseed = 1234, 
    plotkeyargs...)
    # valentformat: multivalent frequency and  preferential pairing
    gls = plot_valentgroupfreq(valentsum["valentgroupfreq"]; fontsize, markersize, plotkeyargs...)
    gls2 = plot_prefpair(valentsum["prefpairing"]; fontsize,linewidth=linewidth-1, plotkeyargs...)    
    append!(gls,gls2)
    fig1 = plot(gls...; 
        layout = (length(gls), 1),
        size=(1200,400*length(gls)),                
        left_margin = 20Plots.mm,
        bottom_margin = 10Plots.mm
    )
    outfile = string(outstem,"_valent_formation.png")
    savefig(fig1, getabsfile(workdir,outfile))
    msg = string("valentformation plot: ", outfile)
    printconsole(io,verbose,msg)    

    # valent frequency
    gls = plot_valentfreq(valentsum["valentfreq"]; 
        fontsize, markersize, markerseed, plotkeyargs...)
    gls_pair = PolyOrigin.plot_pairfreq(valentsum["pairfreq"]; 
        fontsize, markersize, markerseed, plotkeyargs...)
    pushfirst!(gls, gls_pair)
    for i in eachindex(gls)
        if i < length(gls)
            for k in gls[i]
                plot!(k; legend=false)
            end
        end
    end
    kcol = length(gls)
    krow = length(gls_pair)
    gls2 = vec(permutedims(reduce(hcat, gls)))
    fig2 = plot(gls2...; 
        layout = (krow,kcol),
        size=(700*kcol,450*krow),        
        left_margin = 20Plots.mm,        
        bottom_margin = 30Plots.mm
    )
    outfile = string(outstem,"_valent_freq.png")
    savefig(fig2, getabsfile(workdir,outfile))
    msg = string("valentfreq plot: ", outfile)
    printconsole(io,verbose,msg)



    gls = plot_prefpair_posterior(valentsum["prefpairing"]; fontsize,linewidth, plotkeyargs...)
    nplot = length(gls)
    krow = round(Int, sqrt(nplot))
    kcol = ceil(Int,nplot/krow)
    fig3 = plot(gls...;
        layout = (krow, kcol),
        size = (kcol*500,krow*400), 
        left_margin = 10Plots.mm,    
        bottom_margin = 10Plots.mm,     
    )
    outfile = string(outstem,"_valent_prefpairing.png")
    savefig(fig3, getabsfile(workdir,outfile))
    msg = string("prefpairing plot: ", outfile)
    printconsole(io,verbose,msg)
    # 
    fig1,fig2,fig3
end

#######################################################################################

function cal_doublereduction!(polyancestry::PolyAncestry; minprob::Real= 0.5)
    offinfo = polyancestry.offspringinfo
    isoutlier = offinfo[!,:isoutlier]
    if any(isoutlier)
        @info string("exclude ", sum(isoutlier)," outlier offspring")
    end
    offspringset = findall(.!isoutlier)    
    statespace = polyancestry.statespace 
    for (popid, val) in statespace
        offls = findall(offinfo[!,:population] .== popid)
        if length(val["parentindex"] )==1            
            setdiff!(offspringset, offls)
            @info string("calculation of double reduction excluded ", length(offls)," self-fertilized offspring in subpop = ",popid)
        else
            if offinfo[offls[1], :ploidy] == 2
                setdiff!(offspringset, offls)
                @info string("calculation of double reduction excluded ", length(offls)," diploid offspring in subpop = ",popid)
            end
        end
    end
    if isempty(offspringset)
        @warn string("all offspring are produced from selfing. Calculation of dobule reduction applies only to non-selfing offspring")
        return nothing
    end    
    parentinfo = polyancestry.parentinfo
    for chr in eachindex(polyancestry.markermap)
        markermap = polyancestry.markermap[chr]
        chrprob = polyancestry.genoprob[chr]
        nmarker = size(markermap,1)
        noff = length(chrprob)
        res = Matrix{Union{Missing, Float32}}(missing, nmarker, noff)
        for popid in keys(statespace)            
            groupstate = statespace[popid]["groupstate"]
            pls = statespace[popid]["parentindex"]
            length(pls) == 1 && continue
            d1 = div(parentinfo[pls[1], :ploidy],2)
            d2 = div(parentinfo[pls[2], :ploidy],2)
            unique(length.(groupstate)) == [d1+d2] || @error string("inconsistent beetween length of each state and parental ploidies")                        
            drstate = zeros(length(groupstate))
            if d1 > 1
                drstate .+= [!allunique(i[1:d1]) for i in groupstate] # check IBD for each gamete
            end
            if d2 > 1
                drstate .+= [!allunique(i[d1+1:d1+d2]) for i in groupstate] # check IBD for each gamete
            end      
            factor = Int((d1 > 1) + (d2> 1))
            drstate ./= factor
            # factor 2 denoting the fraction of dr per gamete, rather than per zygote
            offls = findall(offinfo[!,:population] .== popid)
            for off in offls
                res[:,off] .= [begin 
                    maxp, groupindex = findmax(prob)
                    maxp >= minprob ? drstate[groupindex] : missing
                end for prob in eachrow(chrprob[off])]
            end
        end
        drfreq = [all(ismissing.(i)) ? NaN : mean(skipmissing(i)) for i in eachrow(res[:, offspringset])]
        colname = "doublereduction" 
        if in(colname, names(markermap))
            markermap[!,colname] .= drfreq
        else
            insertcols!(markermap, size(markermap,2)+1, colname=>drfreq)
        end
    end
    nothing
end


# assuming parents have the same ploidy level
function cal_true_doublereduction(chrtrueancestry::AbstractMatrix)
    ploidy = length(first(chrtrueancestry))
    iseven(ploidy) || @error string("unexpected odd ploidy=",ploidy)
    halfploidy = div(ploidy, 2)
    # /2 accounts for DR per gammete instead of per zygote
    [mean([!allunique(i[1:halfploidy]) + !allunique(i[halfploidy+1:end])  for i in row])/2 for row in eachrow(chrtrueancestry)]
end

function plot_doublereduction(polyancestry::PolyAncestry; 
    boundaryline = (1.0,:dot,:gray),
    markersize::Real= 3,         
    fontsize = 14, 
    isannotate::Real=true, 
    plotkeyargs...)
    hasdr = all([in("doublereduction",names(i)) for i in polyancestry.markermap])
    if !hasdr
        msg = strubg("No column doublereduction in markermap; use cal_doublereduction! to calculate it")
        @error msg
        return nothing
    end
    yposls = [i[:, :doublereduction] for i in polyancestry.markermap]
    xposls = [i[:, :position] for i in polyancestry.markermap]
    xendls = last.(xposls)
    xaddls = accumulate(+,xendls)
    pushfirst!(xaddls, 0.0)
    pop!(xaddls)
    for i in eachindex(xposls, xaddls)
        xposls[i] .+= xaddls[i]
    end
    nchr = length(xposls)
    xmax = xposls[end][end]
    ymax = maximum([begin 
        b = .!isnan.(i)
        i2 =i[b]
        isempty(i2) ? 0.0 : maximum(i2)
    end for i in yposls])
    fig = plot(;         
        size = (max(1200,nchr*65),300),
        dpi = 600,     
        xrange = (-xmax*0.02,xmax*1.02),
        yrange = (-ymax*0.12,ymax*1.08),    
        xlabel = "Genetic position(cM)", 
        ylabel = "Double reduction", 
        legend = false,
        grid = false,                
        left_margin=15Plots.mm,         
        bottom_margin=15Plots.mm,
        top_margin = 5Plots.mm, 
        titlefontsize = fontsize,          
        tickfontsize = fontsize, 
        guidefontsize = fontsize, 
        legendfontsize= fontsize,             
        plotkeyargs...,
    )
    xminls = minimum.(xposls)
    push!(xminls, xmax)
    xy=hcat([[[i,0] [i,ymax] [NaN,NaN]] for i=xminls]...)'
    plot!(fig,xy[:,1],xy[:,2];
        line = boundaryline,
        legend=false,
    )
    if isannotate
        chrposls = (xminls[1:end-1] .+ xminls[2:end]) ./ 2    
        chridls = [i[1,:chromosome] for i in polyancestry.markermap]
        annotate!(fig,[(chrposls[i], ymax*1.1, text(chridls[i], fontsize)) for i in eachindex(chrposls)])
    end
    # scatter        
    marker_colors = distinguishable_colors(nchr, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    for chr=1:nchr
        xls = xposls[chr]
        yls = yposls[chr]
        c = marker_colors[chr]
        scatter!(fig,xls, yls, 
            legend=false,
            marker=(:circle, markersize, 0.5, c,stroke(c)),            
        )        
    end
    fig
end


#######################################################################################

function cal_posterior_recom(polyancestry::PolyAncestry;minprob::Real= 0.5)    
    bestancestry=ancestrycall(polyancestry.genoprob,minprob=minprob)
    nrecomls = calnumrecom.(bestancestry)
    colid = [string("posterior_recomnum_", i[1,:chromosome]) for i in polyancestry.markermap]
    df = DataFrame(colid .=> nrecomls)
    insertcols!(df,1, "individual" => polyancestry.offspringinfo[!,:individual])
    chrlenls = [i[end,:position]-i[1,:position] for i in polyancestry.markermap]
    chrlen = join(round.(chrlenls,digits=4),"|")
    insertcols!(df,size(df,2)+1,"chrlen_cM"=>chrlen)
    sum_nrecom = sum(nrecomls)
    insertcols!(df,size(df,2)+1,"posterior_recomnum"=>sum_nrecom)    
    genolen = 0.01 * sum(chrlenls)
    recomden = sum_nrecom./ genolen
    insertcols!(df,size(df,2)+1,"posterior_recomden"=>recomden)
    df
end

function cal_prior_recom(polyancestry::PolyAncestry)
    pinfo = polyancestry.parentinfo
    pls = pinfo[!,:individual]
    pdict = Dict(pinfo[!,:individual] .=> pinfo[!,:ploidy])
    res = Dict()
    for row in eachrow(polyancestry.designinfo)
        popid = row[:population]        
        ploidyls = [pdict[i] for i in pls[Vector(row[pls]) .> 0]]
        if length(ploidyls) == 1
            # selfing
            priden = only(ploidyls)
        else
            priden = div(sum(ploidyls),2) # assuming bivalent-only formation when producing S1 or F1 gammetes
        end
        push!(res, popid => priden)    
    end
    recomdenls = [res[i] for i in polyancestry.offspringinfo[!,:population]]
    DataFrame("individual" => copy(polyancestry.offspringinfo[!,:individual]), 
        "prior_recomden" => recomdenls)
end

function cal_posterior_prior_recom!(polyancestry::PolyAncestry;minprob::Real= 0.5)            
    recomdf = cal_posterior_prior_recom(polyancestry; minprob)
    polyancestry.offspringinfo = recomdf
    recomdf
end

function cal_posterior_prior_recom(polyancestry::PolyAncestry;minprob::Real= 0.5)        
    postnrecom = cal_posterior_recom(polyancestry; minprob)
    prinrecom = cal_prior_recom(polyancestry)
    ncol = min(4, size(polyancestry.offspringinfo,2))
    recomdf = innerjoin(polyancestry.offspringinfo[!,1:ncol],postnrecom, prinrecom, on= :individual)
    recomdf
end

function plot_posterior_recom(recomdf::AbstractDataFrame; 
    priorline = (2,:dot,:blue),
    fontsize = 12, 
    markersize = 5, 
    plotkeyargs...)        

    chrlen = 0.01 .* parse.(Float64, split(string(recomdf[1,:chrlen_cM]),"|"))
    genolen = sum(chrlen)
    isoutlier = recomdf[!,:isoutlier]
    # non-outlier
    nonoutlierdf = recomdf[.!isoutlier, :]
    colnames = names(nonoutlierdf)
    colls = occursin.("posterior_recomnum_", colnames)
    nrecom = Matrix(nonoutlierdf[:,colls])                
    nrecom_den = nrecom ./ permutedims(chrlen)
    chridls = replace.(colnames[colls], "posterior_recomnum_"=>"")
    xls = repeat(chridls,inner = size(nonoutlierdf,1))
    yls = vec(nrecom_den)
    grecom = violin(xls,yls, 
        side=:left, 
        linewidth=0,         
        label="posterior recomden",
        xdiscrete_values=chridls
    )
    # dotplot!(xls,yls, side=:right, linewidth=0, label="")
    prior_den = mean(nonoutlierdf[!,:prior_recomden])
    plot!(grecom, x->prior_den, 
        line = priorline, 
        label="prior density",                
    )    
    yls = vec(nrecom)
    violin!(grecom,xls,yls, 
        side=:right, 
        linewidth=0,         
        color = :red, 
        label="posterior recomnum"
    )    
    post_den = mean(nonoutlierdf[!,:posterior_recomnum] ./ genolen)
    plot!(grecom; 
        xlabel = "Linkage group", 
        ylabel = "Recombination breakpoints", 
        title = string("prior density=", round(prior_den,digits=3),", posterior density = ", round(post_den,digits=3), " per Morgan per offspring"),         
        size = (1000,500),        
        left_margin = 10Plots.mm,         
        bottom_margin = 10Plots.mm,               
        titlefontsize = fontsize,          
        tickfontsize = fontsize, 
        guidefontsize = fontsize, 
        legendfontsize= fontsize,             
        plotkeyargs...,
    )  
    # outlier
    if any(isoutlier)    
        outlierdf = recomdf[isoutlier, :]
        colnames = names(outlierdf)
        # colls and chridls are the same as that of nonoutlierdf
        nrecom = Matrix(outlierdf[:,colls])                
        xls = repeat(chridls,inner = size(outlierdf,1))
        yls = vec(nrecom)
        outlabel = string(size(outlierdf, 1), " outliers of ",size(recomdf,1)," progeny")
        dotplot!(grecom, xls, yls; 
            bar_width = 0.2, 
            marker=(:red,:xcross,markersize),
            label = outlabel
        )
    end  
    grecom    
end

function plot_valent_DR_recom!(polyancestry::PolyAncestry;
    isplot::Bool=true, 
    chrpairing::Integer= 44, 
    valentthresh = 0.2, 
    prefstep = 0.001,
    minprob::Real= 0.5, 
    priorline = (2,:dot,:blue),
    fontsize = 14,     
    workdir::AbstractString=pwd(), 
    outstem::AbstractString="outstem", 
    io::Union{Nothing,IO}=nothing, 
    verbose::Bool=true
    )
    if isplot
        figdir = joinpath(workdir, outstem * "_plots")
        isdir(figdir) || mkdir(figdir)
    end
    if chrpairing > 22
        try             
            valentsum = cal_valentsum(polyancestry; valentthresh, prefstep)        
            outfile = string(outstem,"_valentsummary.csv")
            msg = string("valentsummary file: ", outfile)
            printconsole(io,verbose,msg)
            savedict2dlm(getabsfile(workdir,outfile), valentsum)
            if isplot
                try
                    plot_valentsum(valentsum; workdir=figdir, outstem, io, verbose, 
                        fontsize, markersize=6)            
                catch err
                    @warn string("Failed to visualize valent summary!")
                    @error err
                end
            end
        catch err
            @warn string("Failed to calculate valent summary!")
            @error err
        end
        try 
            # double reduction 
            cal_doublereduction!(polyancestry; minprob)                
            if isplot 
                fig = plot_doublereduction(polyancestry; fontsize, markersize=3);           
                outfile = string(outstem,"_valent_DR.png")
                savefig(fig, getabsfile(figdir,outfile))
                msg = string("double reduction plot: ", outfile)
                printconsole(io,verbose,msg)
            end
        catch err
            @warn string("Failed to visualize double reduction!")
            @error err
        end
    end

    try 
        recomdf = cal_posterior_prior_recom!(polyancestry; minprob)
        if isplot 
            fig = plot_posterior_recom(recomdf; priorline, fontsize)
            outfile = string(outstem,"_posterior_nrecom.png")
            savefig(fig, getabsfile(figdir,outfile))
            msg = string("posterior_nrecom plot: ", outfile)
            printconsole(io,verbose,msg)
        end
    catch err
        @warn string("Failed to visualize posterior_recom")
        @error err
    end
    nothing
end

#######################################################################################