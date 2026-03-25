

# (ma \kron mb)*x = vec(mb * mx * ma'), where x= vec(mx) 
function kronvec(ma::AbstractMatrix, mb::AbstractMatrix, vx::AbstractVector)
    mx = reshape(vx, size(mb,2),:)
    vec(mb * mx * ma')
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# fwprob[t] =  Pr(x[t]|y[1:t])
# fwscale[[t]] = Pr(y[t]|y[1:t-1])
function forward(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat
    fwprob= Vector{Vector{T}}()
    fwscale=Vector{T}()
    # prob=initprob .* dataprobseq[1]
    prob=kron(initprob...) .* dataprobseq[1]
    scale=sum(prob)
    prob /= scale
    push!(fwscale,scale)
    push!(fwprob,prob)
    for t=2:size(dataprobseq,1)
        # prob = (tranprobseq[t-1]' * fwprob[t-1]) .* dataprobseq[t]
        tranprob = [i[t-1] for i=tranprobseq]
        if any(isnothing.(tranprob))
            prob0 = fwprob[t-1]
        else
            prob0 = kronvec(tranprob[1]',tranprob[2]', fwprob[t-1])
        end
        prob = prob0 .* dataprobseq[t]
        scale=sum(prob)
        prob /= scale
        push!(fwscale,scale)
        push!(fwprob,prob)
    end
    (fwprob=fwprob,fwscale=fwscale)
end

function calloglike(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    calloglike(fwscale)
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# bwprob[t] =  Pr(y[t+1:T]|x[t]) / Pr(y[t+1:T]|y[1:t])
function backward(tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}},
    fwscale::Vector{T}) where T <: AbstractFloat
    bwprob=Vector{Vector{T}}()
    prob=ones(size(tranprobseq[end],2))
    pushfirst!(bwprob,prob)
    for t=size(dataprobseq,1)-1:-1:1
        # prob = tranprobseq[t] * (dataprobseq[t+1] .* bwprob[1])
        prob = (dataprobseq[t+1] .* bwprob[1])
        tranprob = [i[t] for i=tranprobseq]        
        if !any(isnothing.(tranprob))
            prob .= kronvec(tranprob[1],tranprob[2],prob)
        end
        prob /= fwscale[t+1]
        pushfirst!(bwprob,prob)
    end
    bwprob
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# logbwprob[t] = log(Pr(y[t+1:T],x[t]))
function logbackward(tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}},
    logbwinit::Vector{T}=zeros(T,length(dataprobseq[end]))) where T <: AbstractFloat
    nseq = length(dataprobseq)
    logbwprob = Vector{Vector{T}}(undef,nseq)
    logbwprob[end] = logbwinit
    for t=nseq-1:-1:1
        logmax = max(logbwprob[t+1]...)
        prob = exp.(logbwprob[t+1] .- logmax)        
        tranprob = [i[t] for i=tranprobseq]
        temp = (dataprobseq[t+1] .* prob)
        if !any(isnothing.(tranprob))
            temp .= kronvec(tranprob..., temp) 
        end
        logbwprob[t] = log.(temp) .+ logmax
    end
    logbwprob
end

function posteriordecode(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    loglike = calloglike(fwscale)
    bwprob=backward(tranprobseq,dataprobseq,fwscale)
    posteriorprob=[fwprob[i] .* bwprob[i] for i=1:size(dataprobseq,1)]
    (loglike=loglike,posteriorprob=posteriorprob)
end



function backwardsample(tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}},
    fwprob::Vector{Vector{T}}; 
    samplesize::Integer=1000) where T <: AbstractFloat
    tmax = length(fwprob)
    ns = length(first(dataprobseq))
    if !isnothing(tranprobseq[1][1]) && !isnothing(tranprobseq[2][1])
        ns1,ns2 = [size(tranprobseq[i][1],1) for i in 1:2]
        (ns == ns1 * ns2) || @error string("inconsistent dimensions!", " [ns, ns1,ns2]=",[ns,ns1,ns2])    
    end
    # statespace = [(i,j) for i in 1:ns1, j in 1:ns2]
    inttype = ns < typemax(Int8) ? Int8 : Int16
    bwsamples= zeros(inttype, samplesize, tmax)
    hidden = rand(Categorical(fwprob[end]),samplesize)
    nexthidden = similar(hidden)
    bwsamples[:,end] .= hidden
    for t = tmax-1:-1:1        
        sls = unique(hidden)
        tranprob1 = tranprobseq[1][t]
        tranprob2 = tranprobseq[2][t]
        if isnothing(tranprob1) || isnothing(tranprob2)
            nexthidden .= hidden
        else
            for s in sls
                # To speedup for column of kron
                # s1,s2 = statespace[s]                                                
                prob = kron(tranprob1,tranprob2)[:,s]                
                prob .*= fwprob[t]            
                prob ./= sum(prob)
                b = hidden .== s
                nexthidden[b] .= rand(Categorical(prob),sum(b))
            end
        end
        hidden .= nexthidden
        bwsamples[:,t] .= hidden
    end
    bwsamples
end
