##########################################################################################


##########################################################################################

function get_nphased(ploidy::Integer, nallele::Integer)
    nallele^ploidy
end

function get_nunphased(ploidy::Integer, nallele::Integer)
    1+ sum(binomial(nallele+i-2,i) for i in 1:ploidy)
end


# assume nallele < typemax(Int8)
function get_phased_codeset(ploidy::Integer, nallele::Integer)
    res = [[i] for i in Int8.(0:nallele-1)]
    for _ in 2:ploidy
        res = [vcat([j],i) for i in res for j in Int8.(0:nallele-1)]
    end
    res
end

function get_unphased_codeset(ploidy::Integer, nallele::Integer)
    res = [[i] for i in Int8.(0:nallele-1)]
    for _ in 2:ploidy
        res = [vcat([j],i) for i in res for j in Int8.(0:first(i))]
    end    
    res
end

function unphased_genocode_2_index(genostring::AbstractString)        
    ls = split(genostring,"/")    
    genocode = in("N", ls) ? missing : sort(parse.(Int8, ls))
    unphased_genocode_2_index(genocode)    
end

function unphased_genocode_2_index(genocode::Union{Missing,AbstractVector})    
    ismissing(genocode) && return missing    
    index = sum(binomial(genocode[i]+i-1,i) for i in eachindex(genocode))  # starting from 0
    index += 1  # starting from 1    
    index 
end

function phased_genocode_2_index(genocode::Union{Missing,AbstractVector},nallele::Integer)    
    ismissing(genocode) && return missing    
    index = sum(nallele^(i-1)*genocode[i] for i in eachindex(genocode))     # starting from 0
    index += 1  # starting from 1    
    index     
end


##########################################################################################

function get_sitephase(snp::Integer, pindex::AbstractVector,
    phaploindex::AbstractVector,phaploset::AbstractVector,
    phaploweight::AbstractVector; minweight = 1e-8)    
    weightls = [phaploweight[i][snp] for i in pindex]
    phalols = [phaploset[i][snp] for i in pindex]
    if length(pindex) == 1
        sitephaseweight = weightls[1]
        sitephaseset = [[i] for i in eachindex(phalols[1])]        
        posls = eachindex(sitephaseweight,sitephaseset)
    else
        sitephaseweight = [*(i...) for i in getouterlist(weightls...)]    
        sitephaseset = getouterlist([collect(eachindex(i)) for i in phalols]...)        
        posls = findall(sitephaseweight .> minweight)
    end        
    hh = [phaploset[p][snp][phaploindex[p][snp]] for p in  eachindex(phaploindex,phaploset)]
    hhls = [begin 
        hh2 = copy.(hh)
        ls = sitephaseset[pos]
        for j in eachindex(pindex, phalols, ls)
            hh2[pindex[j]] .= phalols[j][ls[j]]
        end
        hh2
    end for pos = posls]
    sitephaplo = permutedims(reduce(hcat,hhls))
    sitephaplo, sitephaseset[posls], sitephaseweight[posls] 
end


function get_sitephaplo(snp::Integer, pindex::Integer,phaploindex::AbstractVector,phaploset::AbstractVector)
    hh = [phaploset[p][snp][phaploindex[p][snp]] for p in  eachindex(phaploindex,phaploset)]
    ls = phaploset[pindex][snp]
    hhls = [begin 
        hh2 = copy.(hh)
        hh2[pindex] = ls[i]
        hh2
    end for i in eachindex(ls)]
    permutedims(reduce(hcat,hhls))
end

function get_phaplo(phaploindex::AbstractVector,phaploset::AbstractVector)
    # return a matrix with size #amrkers x #parents, each element being a phased parental genotype
    hh0=[map((x,y)->ismissing(y) ? x[1] .* missing : x[y], phaploset[p],phaploindex[p]) for p = eachindex(phaploset,phaploindex)]
    reduce(hcat,hh0)
end

function cal_logl_parentdata(phaploindex::AbstractVector,phaploweight::AbstractVector)    
    [mapreduce((x,y)->log(max(1e-8,x[y])), +, phaploweight[p], phaploindex[p]) for p = eachindex(phaploweight,phaploindex)]
end

function cal_logl_parentdata(phaploindex::AbstractVector,phaploweight::AbstractVector, 
    priorprocess::AbstractDict; snporder=nothing)    
    markerincl = first(values(priorprocess)).markerincl
    snps = isnothing(snporder) ? findall(markerincl) : collect(eachindex(markerincl))
    [mapreduce((x,y)->log(max(1e-8, x[y])), +, phaploweight[p][snps], phaploindex[p][snps]) for p = eachindex(phaploweight,phaploindex)]
end


function init_phaploindex(phaploweight::AbstractVector;isrand::Bool=true)
    if isrand
        [Vector{Union{Missing,Int}}([rand(Categorical(i)) for i=j]) for j=phaploweight]
    else
        [Vector{Union{Missing,Int}}([randargmax(i) for i=j]) for j=phaploweight]
    end
end

function init_phaploprior(chrparentgeno::AbstractMatrix, 
    chrparentformat::AbstractVector, 
    chrmultiallele::AbstractVector, 
    polyped::PolyPed; 
    phaseerror::Real = 0.0,
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},    
    overdispersion::Union{Real,AbstractVector},
    onlyfounder::Bool=false)
    nseq = size(chrparentgeno,1)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nseq) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nseq) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nseq) : allelicbias    
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nseq) : overdispersion  
    isfounder = polyped.parentinfo[!,:isfounder]    
    phaploset = Matrix{Any}(nothing, size(chrparentgeno)...)
    phaploweight = Matrix{Any}(nothing, size(chrparentgeno)...)    
    chrparentformat2 = split.(chrparentformat,"|")
    issetfounderls = onlyfounder ? [true] : [true, false]
    for issetfounder in issetfounderls
        if issetfounder
            ppls = findall(isfounder)
            ploidyls = polyped.parentinfo[ppls,:ploidy]    
            formatls = first.(chrparentformat2)
        else
            ppls = findall(.!isfounder)
            ploidyls = polyped.parentinfo[ppls,:ploidy]    
            formatls = last.(chrparentformat2)            
        end        
        formatset = unique(formatls)
        d = setdiff(formatset, ["GT_unphased","GT", "GT_phased","AD","AD&GT"])
        isempty(d) || error(string("TODO for formats=",d))
        for m in eachindex(formatls, chrmultiallele, doseerrorls,baseerrorls,allelicbiasls,overdispersionls)
            format = formatls[m]
            nallele = chrmultiallele[m]     
            doseerror = doseerrorls[m]                    
            baseerror = baseerrorls[m]
            allelicbias = allelicbiasls[m]
            overdispersion = overdispersionls[m]
            for j in eachindex(ploidyls,ppls)
                ploidy = ploidyls[j]  
                pp = ppls[j]                            
                g = chrparentgeno[m, pp]  
                if isa(g, AbstractVector) && in("N",g) 
                    g = missing
                elseif isa(g, AbstractString) && occursin("N", g)
                    g = missing
                end
                if isa(g, AbstractVector) && eltype(g) <: AbstractString && occursin("GT", format) 
                    # phased parent genotypes                                                            
                    if phaseerror ≈ 0 && !in("N",g)                        
                        g2 = parse.(Int8, g) 
                        phaploset[m,pp] = [g2]
                        phaploweight[m,pp] = [1.0]
                    else          
                        genotypes = get_phased_codeset(ploidy, nallele)    
                        nphased = length(genotypes)                        
                        probs = ones(nphased) ./ nphased                        
                        phaploset[m,pp] = genotypes
                        phaploweight[m,pp] = probs
                    end
                    continue
                end
                genotypes = get_unphased_codeset(ploidy, nallele)    
                nunphased = length(genotypes)
                probs_prior = ones(nunphased) ./ nunphased              
                if ismissing(g)
                    probs = probs_prior
                elseif occursin("AD", format)  && isa(g, AbstractVector) && eltype(g) <: Integer 
                    probs = postprob_AD(g, ploidy; doseerror, baseerror,allelicbias, overdispersion, digits=6)                                            
                    if ismissing(probs)  
                        probs = probs_prior                                          
                    end
                elseif in(format, ["GT_unphased","GT","AD&GT","GP&GT"]) && isa(g, AbstractString)                    
                    if occursin("N", g)
                        probs = probs_prior
                    else    
                        gindex = unphased_genocode_2_index(g)
                        probs = (doseerror/(nunphased-1)) .* ones(nunphased)                            
                        try                             
                            probs[gindex] = (1-doseerror)                                
                        catch err
                            @info "" g gindex probs 
                            rethrow(err)
                        end
                    end                                
                else
                    # @info "" format g
                    @error string("TODO for format=",format, ", g=",g) maxlog=10
                end
                hls = []
                wls = Float64[]
                for k in eachindex(probs, genotypes)
                    if probs[k] > 0
                        gphased = get_unphased2phased(join(genotypes[k],"/"),nallele) 
                        ng = length(gphased)
                        append!(hls, gphased)
                        append!(wls, (probs[k]/ng) .* ones(ng))
                    end
                end
                wls .= normalize(round.(wls,digits=5),1)                        
                b = wls .> 0 
                phaploset[m,pp] = hls[b]
                phaploweight[m,pp] = wls[b]
            end
        end
    end
    # collect.() can speed up greatly!!
    collect.(eachcol(phaploset)), collect.(eachcol(phaploweight))    
end

function get_unphased2phased(g::AbstractString, nallele::Integer)
    if occursin("N", g)
        g2 = [i == "N" ? Int8.(0:nallele-1) : parse(Int8,i) for i in split(g, "/")]
        unique(reduce(vcat, [unique(permutations(i)) for i in Iterators.product(g2...)]))
    else
        get_unphased2phased(g)        
    end
end

function get_unphased2phased(g::AbstractString)
    g2 = parse.(Int8, split(g, "/"))
    unique(permutations(g2))     
end

function getchrparentgeno(chrparentgeno::AbstractMatrix)
    Matrix{Any}([eltype(i) <: AbstractString ? (in("N",i) ? missing : parse.(Int8, i)) : i for i in chrparentgeno])
end


function get_derivedgeno(phaploindex::AbstractVector,
    phaploset::AbstractVector,
    chrmultiallele::AbstractVector, 
    priorspace::AbstractDict; 
    popidls=keys(priorspace),
    parentindexls::AbstractVector=[])   
    phaplo = get_phaplo(phaploindex,phaploset)
    get_derivedgeno(phaplo, chrmultiallele, priorspace; popidls,parentindexls)
end

function get_derivedgeno(phaplo::AbstractMatrix,     
    chrmultiallele::AbstractVector, 
    priorspace::AbstractDict; 
    popidls=keys(priorspace),
    parentindexls::AbstractVector=[])       
    if isempty(parentindexls)
        popidls2 = popidls
    else
        popidls2 = collect(popidls)
        for (key, val) in priorspace
            pls = intersect(parentindexls,val["midparentindex"])
            isempty(pls) || push!(popidls2,key)    
        end
        unique!(popidls2)
    end
    res = Dict([begin
        states = priorspace[popid]["state"]
        pp = priorspace[popid]["parentindex"]
        subphaplo = reduce.(vcat,eachrow(phaplo[:,pp]))
        derivedgeno = [h[j] for h in subphaplo, j in states]
        # derived is a matrix with size rxc, r=#markers, c=#states
        # matrix element is an unphased/phased genotype code starting from 1 (==missing)   
        derived_unphased = unphased_genocode_2_index.(sort.(derivedgeno))
        derived_phased = similar(derived_unphased)        
        for na in unique(chrmultiallele)
            b = chrmultiallele .== na
            # dersub = view(derivedgeno,b,:)
            # ressub = view(derived_phased,b,:)
            # ressub .= phased_genocode_2_index.(dersub,na)
            map!(x->phased_genocode_2_index(x,na), view(derived_phased,b,:), view(derivedgeno,b,:))
        end        
        popid => (unphased = derived_unphased, phased = derived_phased)   
    end for popid = popidls2])
    res
end

function get_site_derivedgeno(snp::Integer, pindex::Integer,phaploindex::AbstractVector,
    phaploset::AbstractVector, priorspace::AbstractDict; site_nallele::Integer)
    popidls = popfrompindex(pindex,priorspace)
    sitephaplo = get_sitephaplo(snp,pindex,phaploindex, phaploset)   
    get_site_derivedgeno(sitephaplo, priorspace; site_nallele, popidls)
end


function get_site_derivedgeno(sitephaplo::AbstractMatrix, priorspace::AbstractDict; 
    site_nallele::Integer, popidls = keys(priorspace))    
    Dict([begin
        states = priorspace[popid]["state"]
        pp = priorspace[popid]["parentindex"]     
        subphaplo = reduce.(vcat,eachrow(sitephaplo[:,pp]))        
        hkey = hash((subphaplo, site_nallele,states))
        popid =>  get!(lru_site_derivedgeno, hkey) do
            get_site_subderivedgeno(subphaplo, site_nallele, states)    
        end                            
    end for popid in popidls])
end

function get_site_subderivedgeno(subphaplo::AbstractVector, site_nallele::Integer, 
    states::AbstractVector)
    derivedgeno = [h[j] for h in subphaplo, j in states]
    # derived is a matrix with size rxc, r=#possible founder haplotype combiinations, c=#states
    # matrix element is an unphased/phased genotype code starting from 1 (==missing)             
    derived_phased = phased_genocode_2_index.(derivedgeno,site_nallele)      
    sort!.(derivedgeno)
    derived_unphased = unphased_genocode_2_index.(derivedgeno)
    (unphased = derived_unphased,phased = derived_phased)   
end

const lru_site_derivedgeno = LRU{UInt, NamedTuple}(maxsize=10^3)

##########################################################################################

# gindex: unphased genotype index, tranformed from the observed dosage/unphased_genotypes
function genolike_GT_unphased(gindex::Union{Missing,Integer}, ngenotype::Integer, doseerror::Real)
    ismissing(gindex) && return ones(ngenotype)
    [gindex == i ? (1-doseerror) : doseerror/(ngenotype-1) for i in 1:ngenotype]
end

# gindex: phased genotype index, tranformed from the observed phased genoytpe
function genolike_GT_phased(gindex::Union{Missing,Integer}, ngenotype::Integer, phaseerror::Real)
    ismissing(gindex) && return ones(ngenotype)
    [gindex == i ? (1-phaseerror) : phaseerror/(ngenotype-1) for i in 1:ngenotype]
end

# allelicdepths: a list of numbers of reads for REF and ALt alleles. 
function genolike_AD(allelicdepths::Union{Missing,AbstractVector}, ploidy::Integer, 
    doseerror::Real,baseerror::Real,allelicbias::Real, overdispersion::Real)    
    (ismissing(allelicdepths) || sum(allelicdepths) == 0) && return missing        
    logpls = loglike_AD(allelicdepths, ploidy; baseerror, allelicbias,overdispersion)            
    maxlogp = maximum(logpls) # maxlogp depends on allelicdepths, but also on error rates. 
    pls = exp.(logpls .- maxlogp)    
    ngenotype = length(logpls)
    errormtx = get_genotype_errormtx(ngenotype,doseerror)
    pls2 = errormtx * pls 
    if maxlogp < - 700
        0
        # @info "test" maxlogp allelicdepths allelicbias doseerror baseerror overdispersion maxlog=20
    end
    pscale = exp(max(maxlogp,-700.0))  # upbound -700 to avoid underflow, exp(-750) == 0.0
    pls2 .* pscale
end


function loglike_AD(allelicdepths::Union{Missing,AbstractVector}, ploidy::Integer;  
    baseerror::Real=0.001,allelicbias::Real=1.0,overdispersion::Real=0.0)        
    nallele = length(allelicdepths)    
    genotypes = get_unphased_codeset(ploidy, nallele)    
    # Dirichlet-Multinomial
    # pooling alt alleles => Beta-binomial distributions
    n0 = allelicdepths[1]
    n1 = sum(allelicdepths[2:end])
    logfact = logfactorial(n0+n1) - logfactorial(n0) - logfactorial(n1)    
    loglls = zeros(length(genotypes))
    for g in eachindex(genotypes)
        genotype = genotypes[g]
        d = sum(genotype .> 0) # number of alt alleles
        f = (allelicbias * d) / (ploidy -d + allelicbias * d)
        q = f*(1-baseerror) + (1-f) * baseerror
        if overdispersion < 1e-5
            loglls[g]  = n0 * log(1-q) + n1 * log(q) + logfact        
        else
            alpha = q/overdispersion
            beta = (1-q)/overdispersion
            loglls[g] = logbeta(n1+alpha,n0+beta) - logbeta(alpha,beta) + logfact
        end
    end            
    # loglike for  dosage = 0, 1, ..., ploidy
    loglls
end

function postprob_AD(allelicdepths::Union{Missing,AbstractVector}, ploidy::Integer;  
    doseerror::Real = 0.0, 
    baseerror::Real=0.001,
    allelicbias::Real=1.0, 
    overdispersion::Real=0.0,
    digits::Integer=6)
    (ismissing(allelicdepths) || sum(allelicdepths) == 0) && return missing
    loglls = loglike_AD(allelicdepths, ploidy; baseerror, allelicbias,overdispersion)   
    pls = exp.(loglls .- logsumexp(loglls))
    if doseerror > 0 
        errormtx = get_genotype_errormtx(length(pls),doseerror)
        pls .= errormtx * pls 
    end
    round.(pls; digits)
end


function genocall_AD(allelicdepths::Union{Missing,AbstractVector}, ploidy::Integer;  
    doseerror::Real = 0.0, 
    baseerror::Real=0.001,
    allelicbias::Real=1.0, 
    overdispersion::Real=0.0,   
    callthreshold::Real=0.95)
    (ismissing(allelicdepths) || sum(allelicdepths) == 0) && return missing
    pls = postprob_AD(allelicdepths, ploidy; doseerror, baseerror, allelicbias,overdispersion)
    unphasedcalling(pls; callthreshold)    
end

# random genotype error model
function get_genotype_errormtx(ngenotype,doseerror)    
    [i==j ? 1-doseerror : doseerror/(ngenotype-1) for i in 1:ngenotype, j in 1:ngenotype]
end

# obsprob: a list of probabilties for each unphased genotype in order
function genolike_GPGL(obsprob::Union{Missing,AbstractVector}, doseerror::Real)    
    ismissing(obsprob) && return missing    
    ngenotype = length(obsprob)
    errormtx = get_genotype_errormtx(ngenotype,doseerror)
    if all(obsprob .>= 0) 
        # GP (genotype rawgenoprob)
        res = errormtx * obsprob
    elseif all(obsprob .<= 0) 
        # GL (log10 scaled likeliihood)        
        maxlogp = maximum(obsprob)
        pls = [10^(i - maxlogp) for i in obsprob]        
        pls2 = errormtx * pls
        pscale = exp(max(maxlogp,-700.0))  
        res = pls2 .* pscale
    elseif any(obsprob .> 0) || any(isnan.(obsprob))
        obsprob2 = copy(obsprob)
        b = isnan.(obsprob) 
        any(b) && (obsprob2[b] .= 1.0)
        b2 = isinf.(obsprob) 
        any(b2) && (obsprob2[b2] .= 100.0) 
        s = sum(obsprob2)
        if s ≈ 0.0 
            obsprob2 .= 1/length(obsprob2)
        else
            obsprob2 ./= s
        end
        res = errormtx * obsprob2
        if any(b) || any(b2) || any(obsprob .< 0)
            @warn string("unexpected obsprob = ", obsprob, " for GP or GL", ", new obsprob=",obsprob2)
        end
    else 
        @error string("unexpected obsprob = ", obsprob, " for GP or GL", ", new obsprob=",obsprob2)
    end  
    res    
end

# see update_remainfilels
function resolve_format(obsgeno,format::AbstractString)
    if format == "GT"
        format = isa(obsgeno,AbstractVector) ? "GT_phased" : "GT_unphased"
    elseif occursin("&",format) 
        formatls = sort(split(format, "&"))
        if formatls[2] == "GT"
            if typeof(obsgeno) <: Union{Integer,AbstractString}
                "GT_unphased"
            elseif isa(obsgeno,AbstractVector)
                if eltype(obsgeno) <: AbstractString
                    "GT_phased" 
                elseif eltype(obsgeno) <: Integer
                    formatls[1] == "AD" || @error string("unexpected mixed format = ",format," for obsgeno=",obsgeno)                
                    formatls[1]
                elseif eltype(obsgeno) <: AbstractFloat
                    in(formatls[1], ["GP","GL"]) || @error string("unexpected mixed format = ",format," for obsgeno=",obsgeno)                
                    formatls[1]
                else
                    @error string("unexpected mixed format = ",format," for obsgeno=",obsgeno)               
                end

            else
                @error string("unexpected mixed format = ",format," for obsgeno=",obsgeno)                
            end
        else
            @error string("unexpected mixed format = ",format," for obsgeno=",obsgeno)
        end
    else
        format 
    end    
end

function get_loglike_hidden(obsgeno, hidden_unphaseindex::Integer,ploidy::Integer; 
    doseerror::Real, baseerror::Real, allelicbias::Real,
    overdispersion::Real,  
    format::AbstractString, nallele::Integer)    
    ismissing(obsgeno) && return 0.0
    format = resolve_format(obsgeno,format)        
    if in(format, ["GT_unphased","GT_phased"])
        obsgeno2 = format == "GT_phased" ? join(sort(obsgeno),"/") : obsgeno
        gindex = unphased_genocode_2_index(obsgeno2)
        if ismissing(gindex)
            0.0
        else
            if gindex == hidden_unphaseindex
                log(1 - doseerror)
            else
                ngenotype = get_nunphased(ploidy, nallele)
                log(doseerror/(ngenotype - 1))
            end
        end  
    elseif format == "AD"                        
        logpls = loglike_AD(obsgeno, ploidy; baseerror, allelicbias,overdispersion)            
        nunphased = length(logpls)
        logpmax = maximum(logpls)
        logpls .- logpmax
        res = sum(exp(logpls[i]) * (i == hidden_unphaseindex ? (1-doseerror) : doseerror/(nunphased-1)) for i in eachindex(logpls))
        log(res) + logpmax
    elseif format in ["GP","GL"]            
        if format == "GP"
            log(obsgeno[hidden_unphaseindex]) # GP probabliiites are in order of unphased genotypes
        else
            obsgeno[hidden_unphaseindex] * log(10) # GL in log10 scale
        end
    else
        error(string("unknown genoformat=",format))
        0.0
    end    
end

function get_dataprob(obsgeno; 
    doseerror::Real, baseerror::Real, allelicbias::Real,
    overdispersion::Real,  
    format::AbstractString, nallele::Integer,ploidy::Integer)    
    format = resolve_format(obsgeno,format)        
    if format == "GT_unphased"
        ngenotype = get_nunphased(ploidy, nallele)
        # likels = genolike_GT_unphased(obsgeno, ngenotype, doseerror)                        
        lru_key = hash((obsgeno, ngenotype, doseerror))
        likels = get!(lru_genolike_GTunphased, lru_key) do
            genolike_GT_unphased(obsgeno, ngenotype, doseerror)                        
        end
    elseif format == "GT_phased"     
        if isa(obsgeno,AbstractVector)
            gtype = eltype(obsgeno)
            if gtype <: AbstractString
                if in("N",obsgeno)             
                    obsgeno2 = missing
                else
                    obsgeno2 = phased_genocode_2_index(parse.(Int8,obsgeno),nallele)                 
                end                
            elseif gtype <: Integer
                obsgeno2 = phased_genocode_2_index(obgeno,nallele) 
            else
                @error string("unexpected obsgeno = ", obsgeno, " for GT_phased")
            end            
        else
            obsgeno2 = obsgeno
        end
        ngenotype = get_nphased(ploidy, nallele)                       
        if ismissing(obsgeno2) 
            likels = missing
        else
            lru_key = hash((obsgeno2, ngenotype, doseerror))
            likels = get!(lru_genolike_GTphased, lru_key) do
                genolike_GT_phased(obsgeno2, ngenotype, doseerror)                           
            end
        end
    elseif format == "AD"                        
        # likels = genolike_AD(obsgeno, ploidy, doseerror, baseerror, allelicbias, overdispersion)                                
        lru_key = hash((obsgeno, ploidy, doseerror, baseerror, allelicbias, overdispersion))
        likels = get!(lru_genolike_AD, lru_key) do
            genolike_AD(obsgeno, ploidy, doseerror, baseerror, allelicbias, overdispersion)                                
        end

    elseif format in ["GP","GL"]            
        likels = genolike_GPGL(obsgeno, doseerror)        
    else
        error(string("unknown offspring genoformat=",format))
    end
    if ismissing(likels) 
        if format == "AD"                
            ngenotype = get_nphased(ploidy, nallele)     
        end
        ones(ngenotype)
    else
        likels 
    end
end

function get_dataprob(obsgeno, derivedgeno::AbstractVecOrMat;
    doseerror::Real, baseerror::Real, allelicbias::Real,
    overdispersion::Real,  
    format::AbstractString, nallele::Integer,ploidy::Integer)
    likels = get_dataprob(obsgeno; doseerror, baseerror, allelicbias, overdispersion, 
        format, nallele, ploidy)
    likels[derivedgeno]
end

const lru_genolike_GTunphased = LRU{UInt, Vector{Float64}}(maxsize=10^3)

const lru_genolike_GTphased = LRU{UInt, Vector{Float64}}(maxsize=10^3)

const lru_genolike_AD = LRU{UInt, Union{Missing,Vector{Float64}}}(maxsize=10^4)

function get_err_at_marker(err::Real, marker) # keep marker 
    err    
end 

function get_err_at_marker(err::AbstractVector, marker) 
    err[marker] 
end 

# indgeno is a list of observed genotypes for all markers in a linkage group for a given offspring/midparent
# popderivedgeno = derivedgeno[popid]
# dataprobls[i] has the fixed length given by the number of states for all possible bvpairprop
function cal_ind_dataprobls!(dataprobls, indgeno::AbstractVector,
    chrformatls::AbstractVector, 
    chrmultiallele::AbstractVector, 
    popderivedgeno::NamedTuple; 
    ploidy::Integer, 
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector})
    dgeno_unphased,dgeno_phased = popderivedgeno  # size of dgeno:  #snps x #states 
    # length(indgeno) == length(chrformatls) == length(chrmultiallele)  # the same length denotes the number of markers
    formatset = unique(chrformatls)
    for format in formatset
        markerls = findall(chrformatls .== format)                
        for m in markerls     
            if ismissing(indgeno[m]) 
                dataprobls[m] .= 1.0
                continue
            end
            doseerror_m = get_err_at_marker(doseerror, m)                                
            baseerror_m = get_err_at_marker(baseerror, m)                                
            allelicbias_m = get_err_at_marker(allelicbias, m)    
            overdispersion_m = get_err_at_marker(overdispersion, m)
            format2 = resolve_format(indgeno[m],format)                   
            dgeno_m =  format2 == "GT_phased" ? dgeno_phased[m,:] : dgeno_unphased[m,:]         
            dataprobls[m] .= get_dataprob(indgeno[m], dgeno_m; 
                doseerror = doseerror_m, baseerror = baseerror_m, 
                allelicbias = allelicbias_m, overdispersion = overdispersion_m,
                format = format2, nallele = chrmultiallele[m], ploidy
            )
        end
    end
    dataprobls
end

# dataprobls[i] has the same length as indspace given by the corresponding bvpair for the individual
function cal_ind_dataprobls(indgeno::AbstractVector, indspace::AbstractVector, 
    chrformatls::AbstractVector, 
    chrmultiallele::AbstractVector, 
    popderivedgeno::NamedTuple; 
    ploidy::Integer, 
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector}, 
    skeletonmarkers::Union{Nothing,AbstractVector}=nothing)    
    dgeno_unphased = popderivedgeno.unphased[:,indspace]
    dgeno_phased = popderivedgeno.phased[:,indspace]
    # length(indgeno) == length(chrformatls) == length(chrmultiallele)  # the same length denotes the number of markers
    dataprobls = [ones(length(indspace)) for _ in 1:length(indgeno)]
    formatset = unique(chrformatls)
    for format in formatset
        markerls = findall(chrformatls .== format)       
        if !isnothing(skeletonmarkers)
            intersect!(markerls, skeletonmarkers)
        end         
        for m in markerls     
            ismissing(indgeno[m]) && continue     #  default dataprob 1.0 in the case of missing
            doseerror_m = get_err_at_marker(doseerror, m)                                
            baseerror_m = get_err_at_marker(baseerror, m)                                
            allelicbias_m = get_err_at_marker(allelicbias, m)    
            overdispersion_m = get_err_at_marker(overdispersion, m)
            format2 = resolve_format(indgeno[m],format)        
            dgeno_m =  format2 == "GT_phased" ? dgeno_phased[m,:] : dgeno_unphased[m,:]         
            dataprobls[m] .= get_dataprob(indgeno[m], dgeno_m; 
                doseerror = doseerror_m, baseerror = baseerror_m, 
                allelicbias = allelicbias_m, overdispersion = overdispersion_m,
                format =format2, nallele = chrmultiallele[m], ploidy
            )
        end
    end
    dataprobls
end



##########################################################################################

# list of all subpopulations with pinex being one of the parents
function popfrompindex(pindex::Integer,priorspace::AbstractDict)
    popfrompindex([pindex],priorspace)
end

function popfrompindex(pindex::AbstractVector,priorspace::AbstractDict)
    res = String[]
    for (popid, val) in priorspace
        isempty(intersect(pindex,val["parentindex"])) || push!(res, popid)
    end
    res
end

# By precalculating sitedatalike, speed up in the case of  many repeats of calling cache_sitedataprob_multiphase
function cache_sitedataprob_multiphase(sitephaplo::AbstractMatrix,     
    sitederived::AbstractDict,     
    sitedatalike::AbstractVector,
    bvpair::NamedTuple,        
    priorspace::AbstractDict; 
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector},
    site_offformat::AbstractString, 
    site_nallele::Integer,
    popidls = keys(priorspace))    
    # dataprobls[ind] is a list of prob matrix for each individual
    # i-th row of the matrix corresonds to i-th  founder haplotype
    # the ij element is the prob of that offspring for  founder haplotype i and state j     
    d = setdiff(popidls, keys(sitederived))   
    isempty(d) || error(string("popid not in sitederived: ",d))   
    resls = []     
    for popid in popidls
        poppriorspace =  priorspace[popid]
        derived = site_offformat == "GT_phased" ? sitederived[popid].phased : sitederived[popid].unphased        
        ndstate = size(derived,2)
        ploidy = poppriorspace["offspringploidy"]  
        # leaf_offspring        
        offls = poppriorspace["offspringindex"]
        for o in offls
            likels = sitedatalike[o] 
            bv = bvpair.offspring[o]
            states = poppriorspace["condstate"][bv]
            der = states == 1:ndstate ? derived :  view(derived,:,states)
            res = likels[der]
            push!(resls,collect(eachrow(res)))
        end 
        # non-leaf offspring ( i.e., midparent parents) have the same ploidy as leaf offspring           
        pls = poppriorspace["midparentindex"]
        if !isempty(pls)
            phaseerror = isnothing(phaseerror_midparent) ? doseerror : phaseerror_midparent
            derived = sitederived[popid].phased
            ndstate = size(derived,2)
            for p in pls
                bv = bvpair.parent[p]
                states = poppriorspace["condstate"][bv]
                der = (states == 1:ndstate) ? derived :  view(derived,:,states)
                # statephalplo[1,p]:  only columns correspoind to the parents of pls  may have different phased genotypes; 
                pg = phased_genocode_2_index(sitephaplo[1,p], site_nallele)    
                res = get_dataprob(pg,der; 
                    doseerror=phaseerror, baseerror, allelicbias, overdispersion, 
                    format = "GT_phased", nallele = site_nallele, ploidy
                )                
                push!(resls,collect(eachrow(res)))
            end             
        end
    end
    nphase = size(sitephaplo, 1)
    # return res[h]: results for h-th founder haplotype; res[h][ind]]: a data probvec for individual ind for a given parental haplotype h
    [getindex.(resls,h) for h in 1:nphase] 
end



function cal_sitedataprob_multiphase(sitephaplo::AbstractMatrix,     
    sitederived::AbstractDict,     
    siteoffgeno::AbstractVector,
    bvpair::NamedTuple,        
    priorspace::AbstractDict; 
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector},
    exclparent::Bool=false,
    site_offformat::AbstractString, 
    site_nallele::Integer,
    popidls = keys(priorspace))    
    # dataprobls[ind] is a list of prob matrix for each individual
    # i-th row of the matrix corresonds to i-th  founder haplotype
    # the ij element is the prob of that offspring for  founder haplotype i and state j     
    d = setdiff(popidls, keys(sitederived))   
    isempty(d) || error(string("popid not in sitederived: ",d))       
    resls = [] 
    for popid in popidls
        poppriorspace =  priorspace[popid]
        derived = site_offformat == "GT_phased" ? sitederived[popid].phased : sitederived[popid].unphased        
        ndstate = size(derived,2)
        ploidy = poppriorspace["offspringploidy"]  
        # leaf_offspring        
        offls = poppriorspace["offspringindex"]
        for o in offls
            likels = get_dataprob(siteoffgeno[o];
                doseerror, baseerror, allelicbias, overdispersion, 
                format = site_offformat, nallele = site_nallele, ploidy,
            )
            bv = bvpair.offspring[o]
            states = poppriorspace["condstate"][bv]
            der = states == 1:ndstate ? derived :  view(derived,:,states)
            res = likels[der]
            push!(resls,collect(eachrow(res)))
        end 
        # non-leaf offspring ( i.e., midparent parents) have the same ploidy as leaf offspring    
        if !exclparent       
            pls = poppriorspace["midparentindex"]
            if !isempty(pls)
                phaseerror = isnothing(phaseerror_midparent) ? doseerror : phaseerror_midparent
                derived = sitederived[popid].phased
                ndstate = size(derived,2)                        
                for p in pls
                    bv = bvpair.parent[p]
                    states = poppriorspace["condstate"][bv]
                    der = (states == 1:ndstate) ? derived :  view(derived,:,states)
                    # statephalplo[1,p]:  only columns correspoind to the parents of pls  may have different phased genotypes; 
                    pg = phased_genocode_2_index(sitephaplo[1,p], site_nallele)    
                    res = get_dataprob(pg,der; 
                        doseerror=phaseerror, baseerror, allelicbias, overdispersion, 
                        format = "GT_phased", nallele = site_nallele, ploidy
                    )                
                    push!(resls,collect(eachrow(res)))
                end              
            end
        end
    end 
    nphase = size(sitephaplo, 1)    
    # return res[h]: results for h-th founder haplotype; res[h][ind]]: a data probvec for individual ind for a given parental haplotype h
    [getindex.(resls,h) for h in 1:nphase]
end


function cal_sitedataprob_singlephase(snp::Integer,        
    phaplo::AbstractMatrix, 
    derivedgeno::AbstractDict,     
    bvpair::NamedTuple,        
    priorspace::AbstractDict,
    chroffgeno::AbstractMatrix; 
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector},
    exclparent::Bool=false,
    site_offformat::AbstractString, 
    site_nallele::Integer)        
    sitephaplo = phaplo[[snp],:]    
    sitederived = Dict([popid => (unphased = val.unphased[[snp],:],phased = val.phased[[snp],:]) 
        for (popid, val) in derivedgeno])
    siteoffgeno = chroffgeno[snp,:]
    # return res[k]: a list of data  probvec for each individual in pouplation k
    first(cal_sitedataprob_multiphase(sitephaplo, sitederived, siteoffgeno, bvpair,priorspace; 
        phaseerror_midparent,doseerror, baseerror, allelicbias, overdispersion, exclparent,
        site_offformat, site_nallele
    ))        
end

##########################################################################################


function calmarglogl_bvpair(phaplo::AbstractMatrix, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector,
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    bvpair::NamedTuple,
    polyped::PolyPed;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector}, 
    popidls=keys(priorspace), 
    exclparent::Bool=false,
    skeletonsize::Union{Nothing,Integer}=nothing, 
    snporder::Union{Nothing,AbstractVector}=nothing)
    isnothing(snporder) && (snporder=1:size(chroffgeno,1))    
    nparent = size(phaplo,2)
    nmarker,noff = size(chroffgeno)        
    logl_off = Vector{Any}(missing,noff)                
    logl_parent = Vector{Any}(missing,nparent)
    parentformatls = repeat(["GT_phased"], nmarker)
    if isnothing(skeletonsize)
        skeletonmarkers = nothing
    else
        pri1=first(values(priorprocess))
        tseq = findall(pri1.markerincl)
        if length(tseq) > 1.5* skeletonsize
            thin = div(length(tseq),skeletonsize)
            tseq = tseq[1:thin:end]
        end
        skeletonmarkers = snporder[tseq]
    end
    # calculation    
    for popid = popidls
        poppriorspace = priorspace[popid]
        condstates = poppriorspace["condstate"]
        keyls = poppriorspace["valentkey"]
        offls = poppriorspace["offspringindex"]
        for off=offls            
            offgeno = chroffgeno[:,off]            
            ploidy = polyped.offspringinfo[off,:ploidy] 
            bvkey = keyls[bvpair.offspring[off]]
            offspace = condstates[bvpair.offspring[off]]                          
            dataprobls = cal_ind_dataprobls(offgeno,offspace, chroffformat, chrmultiallele, derivedgeno[popid]; 
                ploidy,doseerror, baseerror, allelicbias,overdispersion,skeletonmarkers)            
            logl_off[off] = caloffmarglogl_bvpair(dataprobls, bvkey, priorprocess,snporder; skeletonmarkers)
        end
        if !exclparent
            pls  = poppriorspace["midparentindex"]
            for p = pls        
                pgeno = phased_genocode_2_index.(phaplo[:, p], chrmultiallele)            
                ploidy = polyped.parentinfo[p,:ploidy]
                bvkey = keyls[bvpair.parent[p]]  
                pspace = condstates[bvpair.parent[p]]  
                phaseerror = isnothing(phaseerror_midparent) ? doseerror : phaseerror_midparent
                dataprobls = cal_ind_dataprobls(pgeno,pspace, parentformatls, chrmultiallele, derivedgeno[popid]; 
                    ploidy,doseerror=phaseerror, baseerror, allelicbias,overdispersion,skeletonmarkers)            
                logl_parent[p] = caloffmarglogl_bvpair(dataprobls, bvkey, priorprocess,snporder; skeletonmarkers)
            end
        end
    end
    (offspring=logl_off, parent=logl_parent)
end

function sum_margloglls(loglls::NamedTuple)
    sum(all(ismissing.(i)) ? 0.0 : sum(skipmissing(i)) for i in loglls)    
end

function cal_margloglls(phaploindex::AbstractVector,
    phaploset::AbstractVector,    
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector,
    chrmultiallele::AbstractVector,        
    priorspace::AbstractDict,
    priorprocess::AbstractDict,    
    bvpair::NamedTuple,
    polyped::PolyPed;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},     
    overdispersion::Union{Real,AbstractVector},     
    snporder::Union{Nothing,AbstractVector}=nothing,
    incldelmarker::Bool=true)
    phaplo = get_phaplo(phaploindex,phaploset)
    derivedgeno = get_derivedgeno(phaplo,chrmultiallele, priorspace)    
    siblogl = calmarglogl_bvpair(phaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
        priorspace, priorprocess, bvpair, polyped; phaseerror_midparent,doseerror, 
        baseerror, allelicbias, overdispersion, snporder)            
    if incldelmarker
        add_logl_delmarker!(siblogl, bvpair, phaplo, priorspace, priorprocess, 
            chroffgeno,chroffformat, chrmultiallele, polyped; 
            phaseerror_midparent,doseerror, baseerror, allelicbias, overdispersion, snporder)
    end
    siblogl    
end

function cal_marglogp(siblogl::NamedTuple, phaploindex::AbstractVector, phaploweight::AbstractVector)
    logl = sum_margloglls(siblogl)
    logl_parentdata = sum(cal_logl_parentdata(phaploindex, phaploweight))
    logl + logl_parentdata        
end

# for one offspring in one linkage group for a given bvpair
function caloffmarglogl_bvpair(dataprobls::AbstractVector,bvkey::AbstractString, 
    priorprocess::AbstractDict, snporder::AbstractVector; 
    skeletonmarkers::Union{Nothing,AbstractVector}=nothing)        
    markerincl, startprob, tranprobseq = strkey2prior(priorprocess, bvkey)
    tseq = findall(markerincl)        
    if !isnothing(skeletonmarkers) 
        b  = [in(snporder[i], skeletonmarkers) for i in tseq]
        keepat!(tseq, b)
    end
    dataprobseq = dataprobls[snporder[tseq]] 
    calloglike(startprob, tranprobseq, dataprobseq)
end

function calmarglogl_bvpairprop(phaplo::AbstractMatrix, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector,
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    bvpairprop::NamedTuple,
    polyped::PolyPed;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector},
    parentindexls::AbstractVector=[], 
    popidls=keys(priorspace), 
    snporder::Union{Nothing,AbstractVector}=nothing)
    isnothing(snporder) && (snporder=1:size(chroffgeno,1))    
    nparent = size(phaplo,2)
    nmarker,noff = size(chroffgeno)        
    logl_off = Vector{Any}(missing,noff)                
    logl_parent = Vector{Any}(missing,nparent)
    parentformatls = repeat(["GT_phased"], nmarker)
    # calculation    
    for popid = popidls
        poppriorspace = priorspace[popid]                
        offls = poppriorspace["offspringindex"]
        nstate = length(poppriorspace["state"])
        dataprobls = [ones(nstate) for _ in 1:nmarker]        
        for off=offls            
            offgeno = chroffgeno[:,off]            
            ploidy = polyped.offspringinfo[off,:ploidy]                         
            cal_ind_dataprobls!(dataprobls, offgeno,chroffformat, chrmultiallele, derivedgeno[popid]; 
                ploidy,doseerror, baseerror, allelicbias,overdispersion)            
            logl_off[off] = caloffmarglogl_bvpairprop(dataprobls, bvpairprop.offspring[off], poppriorspace,priorprocess,snporder)
        end
        pls  = poppriorspace["midparentindex"]        
        for p = pls        
            pgeno = phased_genocode_2_index.(phaplo[:, p], chrmultiallele)            
            ploidy = polyped.parentinfo[p,:ploidy]
            phaseerror = isnothing(phaseerror_midparent) ? doseerror : phaseerror_midparent            
            cal_ind_dataprobls!(dataprobls, pgeno,parentformatls, chrmultiallele, derivedgeno[popid]; 
                ploidy,doseerror=phaseerror, baseerror, allelicbias,overdispersion)            
            logl_parent[p] = caloffmarglogl_bvpairprop(dataprobls, bvpairprop.parent[p], poppriorspace,priorprocess,snporder)
        end
    end    
    parentinfo = polyped.parentinfo[parentindexls,:]
    ismidparentls = .!parentinfo[!,:isfounder]        
    if any(ismidparentls)
        midparentls = parentindexls[ismidparentls]
        familyls = parentinfo[ismidparentls,:family]          
        for popid = familyls
            poppriorspace = priorspace[popid]                 
            nstate = length(poppriorspace["state"])
            dataprobls = [ones(nstate) for _ in 1:nmarker]      
            pls  = intersect(poppriorspace["midparentindex"], midparentls)            
            for p = pls        
                pgeno = phased_genocode_2_index.(phaplo[:, p], chrmultiallele)            
                ploidy = polyped.parentinfo[p,:ploidy]  
                phaseerror = isnothing(phaseerror_midparent) ? doseerror : phaseerror_midparent                       
                cal_ind_dataprobls!(dataprobls, pgeno,parentformatls, chrmultiallele, derivedgeno[popid]; 
                    ploidy,doseerror=phaseerror, baseerror, allelicbias,overdispersion)            
                logl_parent[p] = caloffmarglogl_bvpairprop(dataprobls, bvpairprop.parent[p], poppriorspace,priorprocess,snporder)
            end
        end
    end
    (offspring=logl_off, parent=logl_parent)
end

# for one offspring in one linkage group for a list of proposal bvprop
function caloffmarglogl_bvpairprop(dataprobls::AbstractVector,bvprop::AbstractVector,
    poppriorspace::AbstractDict,priorprocess::AbstractDict,
    snporder::AbstractVector)
    condstates = poppriorspace["condstate"]
    keyls = poppriorspace["valentkey"]
    [begin
        markerincl, startprob, tranprobseq = strkey2prior(priorprocess, keyls[bv])
        tseq = findall(markerincl)        
        snpseq = snporder[tseq]        
        dataprobseq = [dataprobls[i][condstates[bv]] for i=snpseq]
        calloglike(startprob, tranprobseq, dataprobseq)
    end for bv=bvprop]
end

function upoffmarglogl!(loglgrid::AbstractMatrix,
    dataprobls::AbstractVector,popid::String,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    bvprop::AbstractVector,
    snporder::AbstractVector)
    condstates = priorspace[popid]["condstate"]
    keyls = priorspace[popid]["valentkey"]
    loglgrid[bvprop] = [if ismissing(loglgrid[bv])
        markerincl, startprob, tranprobseq = strkey2prior(priorprocess, keyls[bv])        
        tseq = findall(markerincl)
        snpseq = snporder[tseq]        
        dataprobseq = [dataprobls[i][condstates[bv]] for i=snpseq]
        calloglike(startprob, tranprobseq, dataprobseq)
    else
        loglgrid[bv]
    end for bv=bvprop]
end

function uploglgrid!(loglgrid::AbstractMatrix,bvnow::Integer,
    dataprob::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,popid::AbstractString,
    snporder::AbstractVector)
    # itmax=ceil(Int, sqrt(size(loglgrid,1)))
    itmax = size(loglgrid,1)
    nstuckmax=ceil(Int, itmax/5)
    nstuck=0
    loglnow = -Inf
    for it=1:itmax
        dims = rem(it,3)
        if dims == 0
            bvprop = get_prop_bvpair(priorspace,popid,loglgrid)
        else
            bvprop = priorspace[popid]["valentneighbor"][bvnow][dims]
        end
        loglls = upoffmarglogl!(loglgrid,dataprob, popid, priorspace,priorprocess,bvprop,snporder)
        logl = max(loglls...)
        bv = rand(bvprop[loglls .≈ logl])
        # println("   it=",it,",(bv,logl)=",(bv,logl))
        (bv == bvnow && logl ≈ loglnow) ? nstuck +=1 : nstuck = 0
        bvnow = bv
        loglnow = logl
        nstuck>=nstuckmax && break
    end
    bvnow
end

##########################################################################################
