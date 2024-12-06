
# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# fwprob[t] =  Pr(x[t]|y[1:t])
# fwscale[[t]] = Pr(y[t]|y[1:t-1])
function forward(initprob::Vector{Float64},tranprobseq::Vector{Matrix{Float64}},
    dataprobseq::Vector{Vector{Float64}})
    fwprob=Vector{Vector{Float64}}()
    fwscale=Vector{Float64}()
    prob=initprob .* dataprobseq[1]
    scale=sum(prob)
    prob /= scale
    push!(fwscale,scale)
    push!(fwprob,prob)
    for t=2:size(dataprobseq,1)
        prob = (tranprobseq[t-1]' * fwprob[t-1]) .* dataprobseq[t]
        scale=sum(prob)
        prob /= scale
        push!(fwscale,scale)
        push!(fwprob,prob)
    end
    (fwprob=fwprob,fwscale=fwscale)
end

calloglike(fwscale::Vector{Float64}) = sum(log.(fwscale))

function calloglike(initprob::Vector{Float64},tranprobseq::Vector{Matrix{Float64}},
    dataprobseq::Vector{Vector{Float64}})
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    calloglike(fwscale)
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# bwprob[t] =  Pr(y[t+1:T]|x[t]) / Pr(y[t+1:T]|y[1:t])
function backward(tranprobseq::Vector{Matrix{Float64}},
    dataprobseq::Vector{Vector{Float64}},fwscale::Vector{Float64})
    bwprob=Vector{Vector{Float64}}()
    prob=ones(size(tranprobseq[end],2))
    pushfirst!(bwprob,prob)
    for t=size(dataprobseq,1)-1:-1:1
        prob = tranprobseq[t] * (dataprobseq[t+1] .* bwprob[1])
        prob /= fwscale[t+1]
        pushfirst!(bwprob,prob)
    end
    bwprob
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# logbwprob[t] = log(Pr(y[t+1:T],x[t]))
function logbackward(tranprobseq::Vector{Matrix{Float64}},
    dataprobseq::Vector{Vector{Float64}},
    logbwinit::Vector{Float64}=zeros(length(dataprobseq[end])))
    nseq = length(dataprobseq)
    logbwprob = Vector{Vector{Float64}}(undef,nseq)
    logbwprob[end] = logbwinit
    for t=nseq-1:-1:1
        logmax = max(logbwprob[t+1]...)
        prob = exp.(logbwprob[t+1] .- logmax)
        #  log(0.0) = -inf
        logbwprob[t] = log.(tranprobseq[t] * (dataprobseq[t+1] .* prob)) .+ logmax
    end
    logbwprob
end

# """
#     posteriorDecode(initprob,tranprobseq, dataprobseq)
#
# calculate marginal posterior probablities using the HMM posterior decoding
# algorightm.
#
# # Arguments
#
# `initprob::Vector{Float64}`: initial probablity vector
#
# `tranprobseq::Vector{Matrix{Float64}}`: a vector of transition probaiblity
# matrix from one time to the next
#
# `dataprobseq::Vector{Vector{Float64}}`: a vector of emission probability
# vector at each time
#
# """
function posteriorDecode(initprob::Vector{Float64},tranprobseq::Vector{Matrix{Float64}},
    dataprobseq::Vector{Vector{Float64}})
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    loglike = calloglike(fwscale)
    bwprob=backward(tranprobseq,dataprobseq,fwscale)
    posteriorprob=[fwprob[i] .* bwprob[i] for i=1:size(dataprobseq,1)]
    (loglike=loglike,posteriorprob=posteriorprob)
end
