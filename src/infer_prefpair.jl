

function get_prefpairing_pls(n::Integer, ppref::Real)
    # the first bivalent-only configuration is the prefered
    ppref2 = ppref ≈ 1.0 ? 1-1e-10 : ppref
    pls = ones(n) ./ n 
    pls[1] += (1-1/n) * ppref2
    pls[2:end] .-= ppref2/n
    pls
end

function get_prefparing_loglikels(bivalentcounts::AbstractVector, pprefls)
    yls = copy(bivalentcounts)
    sort!(yls, rev=true)
    # ylogfact = logfactorial(sum(yls)) - sum(logfactorial.(yls)) # works only for integer
    # ylogfact = loggamma(sum(yls)+1) - sum(loggamma.(yls .+ 1))
    n = length(yls)
    loglikels = [begin 
        pls = get_prefpairing_pls(n, ppref)
        sum(yls .* log.(pls)) 
    end for ppref in pprefls]
    loglikels
end


function get_pdf_cdf(xls::AbstractVector, loglls::AbstractVector)
    n = length(xls)
    n == length(loglls) || error("inconsistent length between xls and loglls")
    pdfls = exp.(loglls .- maximum(loglls))
    area = 0.0
    for i in 1:n-1
        area += (pdfls[i] + pdfls[i+1]) * (xls[i+1]-xls[i])/2
    end
    pdfls ./= area
    areals = [(pdfls[i] + pdfls[i+1]) * (xls[i+1]-xls[i])/2 for i in 1:n-1]
    cdfls = accumulate(+, areals)
    pushfirst!(cdfls, 0.0)
    pdfls, cdfls
end

function get_quantile(xls, cdfls, q)
    i = findfirst(>=(q), cdfls)
    if isnothing(i) 
        return last(xls)
    elseif i == 1
        return first(xls)
    elseif i == length(xls)
        return last(xls)
    else
        y1, y2 = cdfls[i-1], cdfls[i]
        x1, x2 = xls[i-1], xls[i]
        y1 ≈ y2 && return (x1+x2)/2
        x1 + (x2-x1)*(q-y1)/(y2-y1)
    end
end

function get_cdf_at(xls, cdfls, x)    
    i = findfirst(>=(x),xls)
    if isnothing(i) 
        return last(cdfls)
    elseif i == 1
        return first(cdfls)
    elseif i == length(xls)
        return last(cdfls)
    else
        y1, y2 = cdfls[i-1], cdfls[i]
        x1, x2 = xls[i-1], xls[i]        
        y1 + (x-x1)*(y2-y1)/(x2-x1)
    end
end


function infer_prefpair(bivalentdf::AbstractDataFrame; pstep::Real=0.001)
    qls = [0.025, 0.25, 0.5, 0.75,0.975]
    xstep = pstep
    xls = collect(0:xstep:1)    
    postdf = DataFrame(ones(0,length(xls)),string.("posteriorPDF",xls))
    col_1st_valent = 3
    namels = names(bivalentdf)[1:col_1st_valent-1]
    estdf = DataFrame(ones(0,length(qls) + col_1st_valent-1),vcat(["log10_Pr<=0.05","MAP"], string.("quantile",qls)))
    leftls = namels .=> [String[] for _ in eachindex(namels)]
    insertcols!(estdf, 1, leftls...)
    for row in eachrow(bivalentdf)
        bivalentcounts = Vector(row[col_1st_valent:end])
        loglls = get_prefparing_loglikels(bivalentcounts, xls)
        pdfls, cdfls = get_pdf_cdf(xls, loglls)
        estmle = xls[argmax(pdfls)]
        estls = [get_quantile(xls, cdfls,q) for q in qls]
        log10pvalue = log10(get_cdf_at(xls, cdfls, 0.05))
        push!(postdf, pdfls)        
        push!(estdf, vcat(Vector(row[1:col_1st_valent-1]),[log10pvalue, estmle],estls))
    end
    hcat(estdf, postdf)
end


function infer_prefpair_bivalentonly(valentcount::AbstractDataFrame; pstep::Real=0.001)
    col_1st_valent = 3
    valentnames = replace.(names(valentcount)[col_1st_valent:end],"valent"=>"")
    valentls = [[parse.(Int, j) for j in split.(split(i,"-"),":")] for i in valentnames]
    groupls = [maximum(length.(i)) for i in valentls]
    cols = vcat(1:col_1st_valent-1, findall(groupls .== 2) .+ (col_1st_valent - 1))
    bivalentdf = valentcount[!,cols]
    infer_prefpair(bivalentdf; pstep)
end
