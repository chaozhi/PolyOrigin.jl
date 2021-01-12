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
    xy=hcat([[[i,0] [i,nstate+0.75] [NaN,NaN]] for i=x0]...)
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
    if truegeno != nothing
        offtrue=vcat([i[:,off] for i=truegeno.offspringgeno]...)
        if ishaploprob
            xy=vcat([[[i,j] for j=offtrue[i]] for i=1:length(offtrue)]...)
            xy2=hcat(xy...)
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
    iiout = findall(polyancestry.offspringinfo[!, :isoutlier])
    iinonout = setdiff(1:size(polyancestry.offspringinfo, 1), iiout)
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
    offspring::Union{Nothing,Integer}=nothing,
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


#
# function calnumrecom(offancestry::AbstractVector)
#     nrecom=[[length(MagicDataPrepare.splitindex(i))-1 for i=eachcol(a)]
#             for a=offancestry]
#     nrecom2 = hcat(nrecom...)
#     sum(nrecom2,dims=2)[:,1]
# end
#
# """
#     plotRecombreak(magicancestry,chrindex=1,
#         colorgradient = ColorGradient([:yellow,:blue,:red]),
#         truegeno=nothing)
#
# plot recombination changepoints.
#
# # Positional arguments
#
# `magicancestry::MagicAncestry`: magicancestry returned from [`magicReconstruct`](@ref).
#
# # Keyword arguments
#
# `chrindex::Integer=1`: chromosome index.
#
# `colorgradient::ColorGradient=ColorGradient([:yellow,:blue,:red])`: color gradient
# for heatmap
#
# `truegeno::Union{Nothing,NamedTuple}`: provides true ancestral origins. It is generated
# by [`readTruegeno!`](@ref).
#
# """
# function plotRecombreak(magicancestry::MagicAncestry;
#     chrindex::Integer=1,
#     colorgradient::ColorGradient=ColorGradient([:yellow,:blue,:red]),
#     truegeno::Union{Nothing,NamedTuple}=nothing)
#     if magicancestry.viterbipath == nothing
#         @error string("viterbipath is required")
#         return -1
#     end
#     a=magicancestry.viterbipath[chrindex]'
#     gmarkermap = groupby(magicancestry.markermap,:chromosome)
#     chrid = gmarkermap[chrindex][1,:chromosome]
#     g1=heatmap(a,c=colorgradient,
#         xlabel=string("SNP index in chr=", chrid),
#         ylabel="Offspring index")
#     estrecom = calnumrecom(magicancestry.viterbipath)
#     g2=histogram(estrecom,xlabel="Estimated #recom in all chrs",
#         ylabel="Frequency",legend=false)
#     if truegeno==nothing
#         plot(g1,g2,
#             title=["(a)" "(b)"],
#             titleloc = :right, titlefont = font(12))
#     else
#         truerecom = calnumrecom(truegeno.offspringgeno)
#         g3=scatter(truerecom,estrecom, smooth = true,
#             xlabel="True #recom in all chrs",
#             ylabel="Estimated #recom",legend=false)
#         plot!(g3,x->x)
#         l = @layout [a  [b
#                          c]]
#         plot(g1,g2,g3,layout=l,
#             title=["(a)" "(b)" "(c)"],
#             titleloc = :right, titlefont = font(12))
#     end
# end


"""
    plotMapComp(mapx, mapy,
        boundaryline = (1.5,:dot,:black),
        plotmarker=(:star, 5, 0.5,:blue,stroke(:blue)),
        isannotate= true,
        xlabel="mapx position",
        ylabel="mapy position")

plot postions of mapx vs those of mapx.

# Positional arguments

`mapx::Vector{DataFrame}`: marker map for all chromosomes.

`mapy::Vector{DataFrame}`: comparing map for all chromosomes.

# Keyword arguments

`boundaryline=(1.5,:dot,:black)`: vertical lines for chromosome boundaries.

`plotmarker=(:star, 5, 0.5,:blue,stroke(:blue))`: scatter markers.

`isannotate::Bool=true`: if ture,  annotate chromosome ID and kendall correlation.

`xlabel::AbstractString="mapx position"`: axis x label

`ylabel::AbstractString="mapy position"`: axis y label

"""
function plotMapComp(mapx::Vector{DataFrame},mapy::Vector{DataFrame};
    boundaryline = (1.5,:dot,:black),
    plotmarker=(:star, 5, 0.5,:blue,stroke(:blue)),
    isannotate::Bool=true,
    xlabel::AbstractString="mapx position",
    ylabel::AbstractString="mapy position"
    )
    nchr=length(mapx)
    nchr == length(mapy) || @error string("inconsistent number of chromosomes")
    snps = [intersect(mapx[chr][!,:marker],mapy[chr][!,:marker]) for chr=1:nchr]
    res=[begin
        df=filter(x->in(x[1],snps[chr]),mapx[chr])
        rule = Dict(mapy[chr][!,:marker] .=> 1:size(mapy[chr],1))
        pos = [get(rule,i,missing) for i=df[!,:marker]]
        isnonmiss= .!ismissing.(pos)
        posnon = pos[isnonmiss]
        xls = df[isnonmiss,:position]
        yls =mapy[chr][posnon,:position]
        tau=round(corkendall(posnon,collect(1:length(posnon))),digits=4)
        [xls,yls,tau]
    end for chr=1:nchr]
    len = [i[1][end] for i=res]
    len = accumulate(+,len)
    for i=2:nchr
        res[i][1] .+= len[i-1]
    end
    ymax = max([max(i[2]...) for i=res]...)
    chridls=[i[1,:chromosome] for i =mapx]
    fig = plot([0,0],[0,ymax],line = boundaryline)
    for chr=1:nchr
        scatter!(fig,res[chr][1],res[chr][2],
            marker=plotmarker,
            legend=false,
        )
        if isannotate
            labx=mean(res[chr][1][[1,end]])
            lab = string(chridls[chr],"\nr=",res[chr][3])
            annotate!(fig,[(labx,ymax,text(lab,:left,font(12)))])
        end
        boundaryx=res[chr][1][end]
        plot!(fig,[boundaryx,boundaryx],[0,ymax],line = boundaryline)
    end
    plot(fig,xlabel=xlabel,ylabel=ylabel)
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


function calvalentfreq(polyancestry::PolyAncestry)
    parentls = polyancestry.parentinfo[!,:individual]
    chrls = [i[1,:chromosome] for i=polyancestry.markermap]
    popls = keys(polyancestry.statespace)
    inddict = Dict([pop =>findall(polyancestry.offspringinfo[!,:population] .== pop) for pop=popls])
    # TODO: parent specific ploidy/valentname
    ploidy =  polyancestry.offspringinfo[1,:ploidy]
    valents = polyancestry.statespace[first(popls)]["valent"]
    valentname = [i[1][1] == 1:ploidy ? join(i[1][1],":") : join(join.(i[1],":"),"-") for i=valents[:,1]]
    resls =[begin
        res = zeros(length(parentls),length(valentname))
        maxvv = [divrem(i[1,1]-1,length(valentname)) .+ 1  for i=polyancestry.valentprob[chr]]
        p1v = [i[1] for i=maxvv]
        p2v = [i[2] for i=maxvv]
        for pop=popls
            pp = polyancestry.statespace[pop]["parentindex"]
            if length(pp)==1
                p1=p2=only(pp)
            elseif length(pp)==2
                p1,p2 = pp
            else
                error(string("wrong parents ",pp))
            end
            for off = inddict[pop]
                res[p1,p1v[off]] +=1
                res[p2,p2v[off]] +=1
            end
        end
        for i=1:size(res,1)
            res[i,:] ./= sum(res[i,:])
        end
        res
    end for chr=1:length(chrls)]
    vvfrac = round.(vcat(resls...),digits=4)
    resdf = DataFrame(vvfrac,Symbol.(valentname))
    df0 = DataFrame(["chromosome"=>repeat(chrls,inner=length(parentls)),
        "parent"=>repeat(parentls,outer=length(chrls))])
    resdf2 = hcat(df0,resdf)
    resdf2
end

function plotvalentfreq(valentfreq::DataFrame)
    valentname = names(valentfreq)[3:end]
    freqparent = combine(groupby(valentfreq,:parent),3:6 .=> sum)
    freqchr = combine(groupby(valentfreq,:chromosome),3:6 .=> sum)
    freqparent[:,2:end] ./= size(freqchr,1)
    freqchr[:,2:end] ./= size(freqparent,1)
    ctg = repeat(valentname, inner = size(freqparent,1))
    x = repeat(freqparent[!,:parent],outer=length(valentname))
    y = Matrix(freqparent[:,2:end])
    barparent = groupedbar(x, y, group = ctg,
        bar_position = :stack,
        xlabel = "Parents", ylabel = "Frequencies",
        lw = 0, framestyle = :box
    )
    ctg = repeat(valentname, inner = size(freqchr,1))
    x = repeat(freqchr[!,:chromosome],outer=length(valentname))
    y = Matrix(freqchr[:,2:end])
    barchr = groupedbar(x, y, group = ctg,
        bar_position = :stack,
        xlabel = "Chromosomes", ylabel = "Frequencies",
        lw = 0, framestyle = :box
    )
    plot(barparent,barchr)
end
