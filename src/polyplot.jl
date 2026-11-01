
"""
    plotcondprob(joinedancestry,offspring=nothing,
        colorgradient = ColorGradient([:white,:blue,:red]),
        boundaryline = (1.5,:dot,:black),
        truemarker=(:star, 5, 0.5,:gray,stroke(:gray)),
        truegeno=nothing,
        plotkeyargs...)

plot heatmap for conditional probability.

# Positional arguments

`joinedancestry::JoinedPolyAncestry`: joinedancestry returned from [`polyorigin`](@ref).

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

`trueorig::Union{Nothing,PolyGeno}`: provides true founder origins. 

"""
function plotcondprob(joinedancestry::JoinedPolyAncestry;    
    isviterbi::Bool=false,   
    ismidparent::Bool=false, 
    individual::Union{Nothing,Integer}=nothing,        
    fontsize = 12, 
    isindexchrlab::Bool=false, 
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),    
    boundaryline = (1.5,:dot,:black),
    truemarker=(:star, 5, 0.5,:gray,stroke(:gray)),
    trueorig::Union{Nothing,PolyGeno}=nothing,    
    plotkeyargs...)
    chridls= [i[1,:linkagegroup] for i=joinedancestry.markermap]    
    if isnothing(individual)
        if ismidparent
            midls = findall(.!joinedancestry.polyped.parentinf[!,:isfounder])
            if isempty(midls) 
                @error "no midparents"
                return nothing
            else
                ind = rand(midls)
            end
        else
            nind = size(joinedancestry.polyped.offspringinfo,1)
            ind = rand(1:nind)
        end    
    else
        ind = individual
    end        
    fglspace = reduce(vcat, [parse.(Int,split(i,"|")) for i in joinedancestry.joinedspace[:,:homologs]])
    nfgl = length(fglspace)
    fglspace == 1:nfgl || @error "fglspace must be 1:nfgl"
    if ismidparent
        indploidy = joinedancestry.polyped.parentinfo[ind,:ploidy]
        indid = joinedancestry.polyped.parentinfo[ind,:individual]
        indgen =  joinedancestry.polyped.parentinfo[ind,:generation]
        isfounder = joinedancestry.polyped.parentinfo[ind,:isfounder]
        if isfounder
            msg = string("no condprob for founder = ", indid)
            @warn msg
            return nothing
        end
    else
        indploidy = joinedancestry.polyped.offspringinfo[ind,:ploidy]
        indid = joinedancestry.polyped.offspringinfo[ind,:individual]        
        indid *= string(joinedancestry.polyped.offspringinfo[ind,:isoutlier] === true ? "(outlier)" : "")
        indgen =  joinedancestry.polyped.offspringinfo[ind,:generation]
    end
    if isviterbi 
        respath = ismidparent ? joinedancestry.parentviterbipath : joinedancestry.viterbipath
        indprob = [begin 
            nmarker = size(chrpath,1)
            prob = spzeros(nmarker,nfgl) 
            for m in 1:nmarker
                ss = chrpath[m,ind]
                v = 1/length(ss)
                for i in ss
                    prob[m,i] += v
                end                
            end
            prob
        end for chrpath in respath]
        # maxprob = maximum([maximum(unique([maximum(i[2]) for i in chrprob])) for chrprob in joinedancestry.haploprob])
        maxprob = maximum([maximum(i) for i in indprob])
        maxprob = round(round(Int,maxprob*indploidy)/indploidy,digits=2)      
    else
        condprob = ismidparent ? joinedancestry.parenthaploprob : joinedancestry.haploprob
        indprob = [begin 
            nmarker = size(chrprob,1)
            prob = spzeros(nmarker,nfgl) 
            for m in 1:nmarker
                i,v = chrprob[m,ind]
                prob[m,[i...]] .= v
            end
            prob
        end for chrprob in condprob]
        # use condprob instead of indprob: keep consisstent maxprob among individuals
        # maxprob = maximum([maximum(unique([maximum(i[2]) for i in chrprob])) for chrprob in condprob])
        maxprob = maximum([maximum(i) for i in indprob])        
        maxprob = round(round(Int,maxprob*indploidy)/indploidy,digits=2)      
    end
    
    heat= reduce(vcat,indprob)'
    if !isnothing(trueorig)
        if ismidparent
            indtrue=vcat([i[:,ind] for i=trueorig.parentgeno]...)    
        else
            indtrue=vcat([i[:,ind] for i=trueorig.offspringgeno]...)    
        end
    end
    nzfglls = findall(sum(heat,dims=2)[:,1] .> 0)
    if !isnothing(trueorig)
        union!(nzfglls,unique(reduce(vcat,indtrue)))
        sort!(nzfglls)
    end
    heat = heat[nzfglls,:]
    if !isnothing(trueorig)        
        dict = Dict(nzfglls .=> eachindex(nzfglls))
        indtrue = [[dict[i] for i in j] for j in indtrue]
    end        

    ibdmap=heatmap(heat;         
        c=colorgradient,
        clim=(0,maxprob),        
        size = (1500, max(400,length(nzfglls)*50)), 
        # dpi = 300,        
        xlabel="Marker index",
        ylabel = "Founder homolog" 
    )
    len = size.(indprob,1)
    nnzstate= size(heat,1)
    x0=accumulate(+,len[1:end])
    pushfirst!(x0,0)
    xy=reduce(hcat,[[[i,0] [i,nnzstate+0.75] [NaN,NaN]] for i=x0])
    chrlabx=xy[1,1:3:end]
    chrlabx=(chrlabx[1:end-1] .+ chrlabx[2:end]) ./ 2        
    if isindexchrlab
        chridls2 = string.(eachindex(chrlabx))
    else
        chridls2 = chridls
    end
    annotate!(ibdmap,[(chrlabx[i],nnzstate + 0.8,Plots.text(chridls2[i],fontsize,:black)) for i=1:length(chrlabx)])        
    
    titlestem = isviterbi ? "Viterbi path" : "Haplotype probablility" 
    yticks = ([eachindex(nzfglls);], [string("h",i) for i=nzfglls])
    plot!(ibdmap,xy[1,:],xy[2,:];         
        line = boundaryline,
        legend=false,         
        yrange = (0,nnzstate+1.1),       
        yticks = yticks,
        title=string(titlestem, " for ", ismidparent ? " midparent" : " offspring", ind," ",indid, " in generation G",indgen),
        # titlefont = font("Times", fontsize)
    )
    if !isnothing(trueorig)        
        xy=vcat([[[i,j] for j=indtrue[i]] for i=1:length(indtrue)]...)
        xy2=reduce(hcat,xy)
        scatter!(ibdmap,xy2[1,:],xy2[2,:],
            marker=truemarker,
            legend=false
        )                
    end
    plot!(ibdmap;         
        titlefontsize = fontsize,          
        tickfontsize = fontsize, 
        guidefontsize = fontsize, 
        legendfontsize= fontsize,
        left_margin = 12Plots.mm,
        bottom_margin = 12Plots.mm,               
        plotkeyargs...,
    )
end

function saveprobplot(
    joinedancestry::JoinedPolyAncestry;    
    isviterbi::Bool=false,     
    nplot_subpop::Integer=5, 
    colorgradient::ColorGradient = cgrad([:white, :blue, :red]),    
    boundaryline = (1.5, :dot, :black),
    truemarker = (:star, 5, 0.5, :gray, stroke(:gray)),
    trueorig::Union{Nothing,PolyGeno}=nothing,
    outstem::AbstractString = "outstem",
    workdir::AbstractString = pwd(),
    plotkeyargs...)
    figdir = getabsfile(workdir, outstem * "_plots")
    isdir(figdir) || mkdir(figdir)
    outid = isviterbi ? "viterbipath_" : "haploprob_" 
    # plot offspring
    offinfo = joinedancestry.polyped.offspringinfo
    popls = offinfo[!,:family]
    offls = sort(reduce(vcat, [begin 
        suboffls = findall(popls .== i )
        if nplot_subpop >= length(suboffls)
            suboffls
        else
            suboffls[1:nplot_subpop]
        end
    end for i in unique(popls)]))
    isoutls = offinfo[!, :isoutlier]
    iiout = any(ismissing.(isoutls)) ? [] : findall(isoutls)
    iinonout = setdiff(offls, iiout)
    for k = 1:2
        ii = k == 1 ? iiout : iinonout        
        outid2 = k == 1 ? outid*"outlier_" : outid
        for i in ii
            indid = joinedancestry.polyped.offspringinfo[i, :individual]
            indid2 = replace(replace(indid,"/"=>"_"),"\\"=>"_")
            indgen =  joinedancestry.polyped.offspringinfo[i,:generation]
            fn = joinpath(
                figdir,
                string(outid2, "offspring",i, "_", indid2, "_G", indgen,".png"),
            )
            plotcondprob(joinedancestry; 
                isviterbi, 
                ismidparent = false, 
                individual = i,
                colorgradient,
                boundaryline,
                truemarker,trueorig,
                plotkeyargs...)
            savefig(fn)
        end
    end
    # plot midparent
    midparentls = findall(.!joinedancestry.polyped.parentinfo[!,:isfounder])
    for midparent in midparentls
        indid = joinedancestry.polyped.parentinfo[midparent, :individual]
        indid2 = replace(replace(indid,"/"=>"_"),"\\"=>"_")
        indgen =  joinedancestry.polyped.parentinfo[midparent,:generation]
        fn = joinpath(
            figdir,
            string(outid, "midparent", midparent,"_", indid2, "_G", indgen,".png"),
        )
        plotcondprob(joinedancestry; 
            isviterbi, 
            ismidparent = true, 
            individual = midparent,
            colorgradient,
            boundaryline,
            truemarker,trueorig,
            plotkeyargs...)
        savefig(fn)
    end
    figdir
end


"""
    animcondprob(joinedancestry,fps=1,outfile=nothing,kewargs...)

animation for plots of conditional probability.

# Positional arguments

`joinedancestry::JoinedPolyAncestry`: joinedancestry returned from [`polyorigin`](@ref).

# Keyword arguments

`fps::Real=1`: number of frames per seconds.

`outfile::AbstractString="condprob.gif"`: output file for saving animation.

see [`plotcondprob`](@ref) for keyargs.

"""
function animcondprob(joinedancestry::JoinedPolyAncestry;
    isviterbi::Bool=false,     
    fps::Real=1,        
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),
    left_margin = :match,
    bottom_margin = :match,
    boundaryline = (1.5,:dot,:black),
    truemarker=(:xcross, 3, 0.5,:gray,stroke(:gray)),
    trueorig::Union{Nothing,PolyGeno}=nothing,
    outfile::Union{Nothing,AbstractString}=nothing)
    noff=size(joinedancestry.polyped.offspringinfo,1)
    anim = Animation()
    n1 = ceil(Int,fps)
    n2 = fps == 0 ? 1 : ceil(Int,1/fps)
    @time for ind=1:noff
        plotcondprob(joinedancestry; 
            isviterbi, 
            ismidparent = false, 
            individual=ind,
            colorgradient,
            left_margin, 
            bottom_margin,
            boundaryline, 
            truemarker,
            trueorig
        )
        for i=1:n2
            frame(anim)
        end
    end
    if isnothing(outfile)
        gif(anim,fps=n1)
    else
        gif(anim,outfile,fps=n1)
    end
end

#######################################################################################

"""
    plotmapcomp(mapx, mapy;
        isphysmap::AbstactVector=[false,false],
        boundaryline = (1.5,:dot,:black),        
        isannotate= true,
        maplabels=["mapx position", "mapy position"])

plot postions of mapx vs those of mapx.

# Positional arguments

`mapx::Vector{DataFrame}`: marker map for all chromosomes.

`mapy::Vector{DataFrame}`: comparing map for all chromosomes.

# Keyword arguments

`isphysmap::AbstactVector=[false,false]`: specicify if [mapx, mapy] is a physical map, respectively. If true, use physical position (bp) instead of genetic position (cM).

`boundaryline=(1.5,:dot,:black)`: vertical lines for chromosome boundaries.

`isannotate::Bool=true`: if ture,  annotate chromosome ID and kendall correlation.

`mmaplabels::Union{Nothing,AbstractString}=nothing`: labels of comparing marker maps.

"""
function plotmapcomp(mapx::Vector{DataFrame},mapy::Vector{DataFrame};
    isphysmap::AbstractVector=[false,false],
    boundaryline = (1.5,:dot,:black),    
    markersize = 3, 
    markershape = :circle, 
    fontsize = 12,     
    isannotate::Bool=true,
    isannotatetau::Bool=true,
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
    isphys_mapx,isphys_mapy = isphysmap
    res=[begin
        df=filter(x->in(x[1],snps[chr]),mapx[chr])
        rule = Dict(mapy[chr][!,:marker] .=> 1:size(mapy[chr],1))
        pos = [get(rule,i,missing) for i=df[!,:marker]]
        isnonmiss= .!ismissing.(pos)
        posnon = pos[isnonmiss]
        xls = isphys_mapx ? df[isnonmiss,:physposbp] .* 1e-6 : df[isnonmiss,:poscm] # Mbp or cM
        yls =isphys_mapy ? mapy[chr][posnon,:physposbp] .* 1e-6 : mapy[chr][posnon,:poscm]
        tau=round(corkendall(posnon,collect(1:length(posnon))),digits=cordigits)
        [xls,yls,tau]
    end for chr=1:nchr]
    len = [i[1][end] for i=res]
    len = accumulate(+,len)
    for i=2:nchr
        res[i][1] .+= len[i-1]
    end
    ymax = max([max(i[2]...) for i=res]...)
    chridls=[i[1,:linkagegroup] for i =mapx]
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
            if isannotatetau
                lab0 = chr==1 ? "\nr=" : "\n"
                lab = string(chridls[chr],lab0,res[chr][3])
            else
                lab = string(chridls[chr])
            end
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
        ids = [i[1,:linkagegroup] for i in truemarkermap[pos]]
        @warn string("1st markermap chromosomes merged: ", truelg[pos])
    end
    ls = vcat(truelg...)
    for chr in 1:length(truemarkermap)
        rep = sum(ls .== chr)
        if rep>1
            id = truemarkermap[chr][1,:linkagegroup]
            msg =string("1st markermap ", id)
            msg *= string(" is splited into ", rep, " chrs in 2nd markermap")
            @warn msg
        end
    end
    pos = findall(isempty.(truelg))
    if !isempty(pos)
        ids = [i[1,:linkagegroup] for i in estmarkermap[pos]]
        @error string("2nd markermap chromosomes not founder: ", ids)
    end
    truelg
end

function plotmapcomp(mapxfile::AbstractString,mapyfile::AbstractString,
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
        [DataFrame(i) for i=groupby(mapdf,:linkagegroup)]
    end for mapfile = filels]
    labells = [string("map",i, " position") for i=1:length(mapls)]
    if !isnothing(maplabels)
        n = min(length(maplabels),length(mapls))
        labells[1:n] .= maplabels[1:n]
    end
    n = length(mapls)
    res = [plotmapcomp(mapls[i],mapls[j];
                boundaryline, plotmarker, isannotate,cordigits,
                maplabels=labells[[i,j]])
                for i=1:n for j=(i+1):n]
    length(res)==1 ? res[1] : plot(res...,layout=(n,1))
end

function plotdesign(polyped::PolyPed;
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
    offinfo = polyped.offspringinfo
    designinfo = polyped.designinfo
    npop, nf= size(designinfo)
    nf -= 1
    adj = zeros(Int,nf,nf)
    # edgelabel_dict = Dict()
    edgelabel_mtx = ["" for i=1:nf, j=1:nf]
    for i=1:size(designinfo,1)
        jj = findall(Vector(designinfo[i,2:end]) .> 0)
        pop = designinfo[i,1]
        popsize = sum(offinfo[!,:family] .== pop)
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

# function plotmarkererror(mapfile::AbstractString;
#     tukeyfence::Union{Nothing,Real}=nothing.0,
#     missingstring=["NA","missing"],
#     commentstring::AbstractString="##",
#     workdir::AbstractString = pwd())
#     markermap = readmarkermap(mapfile; del_ungrouped = false,        
#         commentstring,missingstring, workdir)  
#     plotmarkererror(markermap;tukeyfence)
# end

function plotmarkererror(markermap::AbstractDataFrame;
    tukeyfence::Real=3.0, 
    fontsize=12)
    colls = [:doseerror,:baseerror,:allelicbias,:overdispersion]
    fencels = tukeyfence*ones(length(colls))
    b = [begin         
        isanymiss = any(ismissing.(markermap[!,i]))         
        if isanymiss 
            true
        else
            ls = unique(round.(unique(markermap[!,i]),digits=4))
            # length(ls) <=1 && println("col=",i, ",unqiue values=",ls)
            length(ls) <=1 
        end
    end for i in colls]
    deleteat!(colls,b)
    deleteat!(fencels,b)
    isempty(colls) && return nothing    
    gls = [begin     
        col = colls[c]        
        whisker_range = fencels[c]           
        yls = Vector(markermap[!, col])        
        yls[yls .== 0.0] .= 1e-6
        if in(col, [:allelicbias,:overdispersion])             
            @. yls = log10(yls)            
            ylabel = string("log10 ", col)
        else
            @. yls = log10(yls/(1-yls))            
            ylabel = string("logit10 ", col)
        end
        xls = string.(markermap[!,:linkagegroup])                
        g =violin(xls, yls;
            fillalpha=0.75, linewidth=0.5,                
            legend=false,            
            ylabel
        )              
        boxplot!(xls, yls;
            fillalpha=0.75, linewidth=2,
            whisker_range, 
            legend=false, 
            title = string(col, ": tukeyfence=",whisker_range)
        )   
        if col == :allelicbias
            plot!(x->0.0, line=(2,:dash,:gray), legend=false) # no allelic bias
        end
        g
    end for c in eachindex(colls)]
    nchr = length(unique(markermap[!,:linkagegroup]))
    plot(gls...;
        layout=(length(gls),1),
        size = (max(1000,nchr*60),300*length(gls)),
        titlefontsize = fontsize+1,          
        tickfontsize = fontsize, 
        guidefontsize = fontsize, 
        legendfontsize= fontsize,
        left_margin=15Plots.mm,         
        bottom_margin=10Plots.mm,
    )
end

function plot_phaseerr(accfile::AbstractString; 
    marker_midparent=(6,:diamond,:cyan),     
    marker_founder=(8,:xcross,:blue),
    label_midparent = "Intermediator",
    label_founder = "Founder",
    workdir = pwd(),
    )    
    accfile2 = abspath(workdir, accfile)    
    phaseerrdict = PolyOrigin.readmultitable(accfile2)
    [plot_phaseerr(phaseerrdict; 
        errname=i, marker_midparent, marker_founder, label_midparent, label_founder
    ) for i in ["doseerr","permerr","heteroerr"]] 
end



function plot_phaseerr(phaseerrdict::AbstractDict; 
    errname::AbstractString,
    marker_midparent=(6,:diamond,:brown),     
    marker_founder=(8,:xcross,:blue),    
    label_midparent = "Intermediator",
    label_founder = "Founder",
    )    
    if errname == "doseerr"        
        dosedf = phaseerrdict["unphasederror"]        
        isfounderls = [isa(i,AbstractString) ? parse(Bool,i) : i for i in dosedf[1:end-1,:isfounder]]
        nprogenyls = dosedf[1:end-1,:immediate_progeny]
        if typeof(nprogenyls[1]) <: AbstractString
            nprogenyls = parse.(Int, nprogenyls)
        end
        data = [parse.(Float64,split(i,"|")) for i in dosedf[1:end-1,:allchr]]
        errls = [i[1]/i[2] for i in data]
    else
        phasedf = phaseerrdict["phasingerror"]
        data = [parse.(Float64,split(i,"|")) for i in phasedf[1:end-1,:allchr]] # each element is [nheteroerror, npermerror, noppty]
        nprogenyls = phasedf[1:end-1,:immediate_progeny]        
        if typeof(nprogenyls[1]) <: AbstractString
            nprogenyls = parse.(Int, nprogenyls)
        end
        isfounderls = [isa(i,AbstractString) ? parse(Bool,i) : i for i in phasedf[1:end-1,:isfounder]]
        if errname == "permerr"
            errls = [i[2]/i[3] for i in data]
        elseif errname == "heteroerr"
            namels = names(phasedf)
            if in("nchr", namels)
                nchr = phasedf[1, "nchr"] # same for all individuals
            else
                col = findfirst(==("allchr"), namels)
                nchr = size(phasedf,2) - col
            end
            errls =  [i[1]/(nchr+i[3]) for i in data] 
        else
            error("errname must be doseerr, permerror or heteroerr")
        end
    end    
    maxpower = ceil(Int,log(2, maximum(nprogenyls)))
    xtickls = [2^i for i in 0:maxpower]
    # xtickls = [xtickls, string.(xtickls)]
    xls = nprogenyls[isfounderls]    
    yls = errls[isfounderls]
    gerr = scatter(xls,yls;
        marker=marker_founder,
        # markerstrokecolor = strokecolor_founder, 
        xscale = :log10,
        label = label_founder,
    )
    xls = nprogenyls[.!isfounderls]
    yls = errls[.!isfounderls]
    scatter!(gerr, xls,yls;
        marker=marker_midparent, 
        markerstrokewidth = 0.5,
        # markerstrokecolor = strokecolor_midparent,
        xscale = :log10, 
        label = label_midparent,
        xticks = (xtickls,string.(xtickls)),
        xrange = (0.8,1.1 * 2^maxpower),
        xlabel = "nprogeny", 
        ylabel = errname
    )    
    gerr
end


function write_avgphaseerr(accfilels::AbstractVector;     
    outstem::AbstractString="outstem",
    workdir::AbstractString=pwd())
    resls = [begin 
        accfile2 = abspath(workdir, accfile)
        phaseerrdict = PolyOrigin.readdlm2dict(accfile2)
        dosedf, phasedf = phaseerrdict["unphasederror"],phaseerrdict["phasingerror"]
        dosedf[!,:allchr] = [parse.(Int,i) for i in split.(dosedf[!,:allchr],"|")]
        phasedf[!,:allchr] = [parse.(Int,i) for i in split.(phasedf[!,:allchr],"|")]
        (dosedf, phasedf)
    end for accfile in accfilels]
    namels = names(resls[1][1])
    col = findfirst(==("allchr"), namels)
    dosedf = copy(resls[1][1][!,1:col])
    dosedf[!,col] .= join.(mean(i[1][!,col] for i in resls),"|")

    df = resls[1][2]
    namels = names(df)
    col = findfirst(==("allchr"), namels)
    nchr = size(df,2) - col
    phasedf = copy(df[!,1:col])
    phasedf[!,col] .= join.(mean(i[2][!,col] for i in resls),"|")
    insertcols!(phasedf, size(phasedf,2), :nchr => nchr)

    phaseerrdict = Dict("unphasederror"=>dosedf, "phasingerror"=> phasedf)
    outfile = abspath(workdir,string(outstem,"_avgphaseerr.csv"))
    savedict2dlm(outfile,phaseerrdict)
    outfile
end





function write_avgancestryerr(accfilels::AbstractVector;     
    outstem::AbstractString="outstem",
    workdir::AbstractString=pwd())
    resls = [begin 
        accfile2 = abspath(workdir, accfile)
        dict = PolyOrigin.readdlm2dict(accfile2)
        for (key,df) in dict
            namels = names(df)
            col = findfirst(==("allchr"), namels)
            nchr = size(df,2) - col
            df[!,col] .= [i=="NA" ? i : parse.(Int,split(i,"|")) for i in df[!,col]]
            df2 = copy(df[!,1:col])
            insertcols!(df2, size(df2,2), :nchr => nchr)
            dict[key] = df2
        end
        dict
    end for accfile in accfilels]
    keyls = keys(resls[1])
    resdict = Dict([begin     
        df = deepcopy(resls[1][key])
        namels = names(df)
        col = findfirst(==("allchr"), namels)
        nonmiss = df[!,col] .!= "NA"
        if all(nonmiss) 
            df[!,col] .= join.(mean(i[key][!,col] for i in resls),"|")
        else
            df[nonmiss,col] .= join.(mean(i[key][nonmiss,col] for i in resls),"|")
        end    
        key => df
    end for key in keyls])
    outfile = abspath(workdir,string(outstem,"_avgancestryerr.csv"))
    savedict2dlm(outfile,resdict)
    outfile
end

function plot_ancestryerr(accfile::AbstractString;     
    isinclaprent::Bool=true,
    workdir = pwd(),
    )    
    accfile2 = abspath(workdir, accfile)    
    phaseerrdict = PolyOrigin.readmultitable(accfile2)
    [plot_ancestryerr(phaseerrdict; errname=i, isinclaprent) for i in ["doseerr","permerr","heteroerr"]] 
end

function plot_ancestryerr(ancestryerr_dict::AbstractDict; errname::AbstractString, isinclaprent::Bool=true)
    if errname == "doseerr"                
        groupid = "offspring_"
        dosedf = ancestryerr_dict[groupid*"unphasederror"]   
        genls = dosedf[1:end-1,:generation]
        data = [parse.(Float64,split(i,"|")) for i in dosedf[1:end-1,:allchr]]
        errls = [i[1]/i[2] for i in data]         
        groupid = "parent_"
        if isinclaprent && haskey(ancestryerr_dict, groupid*"unphasederror")            
            dosedf = ancestryerr_dict[groupid*"unphasederror"]   
            nonmiss = findall(.!ismissing.(dosedf[1:end-1,:allchr]))
            append!(genls, dosedf[nonmiss,:generation])
            data = [parse.(Float64,split(i,"|")) for i in dosedf[nonmiss,:allchr]]
            append!(errls, [i[1]/i[2] for i in data])            
        end        
    else
        groupid = "offspring_"
        phasedf = ancestryerr_dict[groupid*"phasingerror"]        
        namels = names(phasedf)
        nchr = length(namels) - findfirst(==("allchr"), namels)
        genls = phasedf[1:end-1,:generation]
        data = [parse.(Float64,split(i,"|")) for i in phasedf[1:end-1,:allchr]]
        if errname == "permerr"
            errls = [i[2]/i[3] for i in data]
        elseif errname == "heteroerr"            
            errls =  [i[1]/(nchr+i[3]) for i in data] 
        else
            error("errname must be doseerr, permerror or heteroerr")
        end
        groupid = "parent_"
        if isinclaprent && haskey(ancestryerr_dict, groupid*"unphasederror")            
            phasedf = ancestryerr_dict[groupid*"phasingerror"]        
            nonmiss = findall(.!ismissing.(phasedf[1:end-1,:allchr]))
            append!(genls, phasedf[nonmiss,:generation])
            data = [parse.(Float64,split(i,"|")) for i in phasedf[nonmiss,:allchr]]
            if errname == "permerr"
                append!(errls, [i[2]/i[3] for i in data])
            elseif errname == "heteroerr"                
                append!(errls, [i[1]/(nchr+i[3]) for i in data])
            else
                error("errname must be doseerr, permerror or heteroerr")
            end            
        end
    end    
    xls = sort(unique(genls))
    avgerrls = [mean(errls[genls .== i]) for i in  xls]
    gerr = violin(genls, errls;    
        # linewidth=0,         
        # color = :cyan, 
        label="",
        xticks = (xls, string.("G", xls)),
        line = (0, :cyan),
    )    
    plot!(xls, avgerrls, label="", line=(2,:red))
    plot!(xlabel="Generation",ylabel=errname)
    gerr
end
