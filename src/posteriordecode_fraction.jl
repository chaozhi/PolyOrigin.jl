
function kronvec(ma::AbstractMatrix, mb::AbstractMatrix, vx::AbstractVector)
    mx = reshape(vx, size(mb,2),:)
    vec(mb * mx * ma')
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# fwprob[t] =  Pr(x[t]|y[1:t])
# fwscale[[t]] = Pr(y[t]|y[1:t-1])
function forward(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{Matrix{T}}},
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
        tranprob = [i[t-1]' for i=tranprobseq]
        prob0 = kronvec(tranprob..., fwprob[t-1])
        prob = prob0 .* dataprobseq[t]
        scale=sum(prob)
        prob /= scale
        push!(fwscale,scale)
        push!(fwprob,prob)
    end
    (fwprob=fwprob,fwscale=fwscale)
end

function calloglike(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{Matrix{T}}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    calloglike(fwscale)
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# bwprob[t] =  Pr(y[t+1:T]|x[t]) / Pr(y[t+1:T]|y[1:t])
function backward(tranprobseq::Vector{Vector{Matrix{T}}},
    dataprobseq::Vector{Vector{T}},
    fwscale::Vector{T}) where T <: AbstractFloat
    bwprob=Vector{Vector{T}}()
    prob=ones(size(tranprobseq[end],2))
    pushfirst!(bwprob,prob)
    for t=size(dataprobseq,1)-1:-1:1
        # prob = tranprobseq[t] * (dataprobseq[t+1] .* bwprob[1])
        tranprob = [i[t]' for i=tranprobseq]
        prob = kronvec(tranprob..., (dataprobseq[t+1] .* bwprob[1]))
        prob /= fwscale[t+1]
        pushfirst!(bwprob,prob)
    end
    bwprob
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# logbwprob[t] = log(Pr(y[t+1:T],x[t]))
function logbackward(tranprobseq::Vector{Vector{Matrix{T}}},
    dataprobseq::Vector{Vector{T}},
    logbwinit::Vector{T}=zeros(T,length(dataprobseq[end]))) where T <: AbstractFloat
    nseq = length(dataprobseq)
    logbwprob = Vector{Vector{T}}(undef,nseq)
    logbwprob[end] = logbwinit
    for t=nseq-1:-1:1
        logmax = max(logbwprob[t+1]...)
        prob = exp.(logbwprob[t+1] .- logmax)
        #  log(0.0) = -inf
        tranprob = [i[t]' for i=tranprobseq]
        temp = kronvec(tranprob..., (dataprobseq[t+1] .* prob))
        logbwprob[t] = log.(temp) .+ logmax
    end
    logbwprob
end

function posteriorDecode(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{Matrix{T}}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    loglike = calloglike(fwscale)
    bwprob=backward(tranprobseq,dataprobseq,fwscale)
    posteriorprob=[fwprob[i] .* bwprob[i] for i=1:size(dataprobseq,1)]
    (loglike=loglike,posteriorprob=posteriorprob)
end
