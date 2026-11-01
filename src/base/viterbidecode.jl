
# """
#     viterbidecode(initprob,tranprobseq, dataprobseq)

# calculate the optimized path using the HMM Viterbi decoding algorightm.

# # Arguments

# `initprob::AbstractVector`: the probablity vector at initial time t=1

# `tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}`: tranprobseq[t] is the transition probaiblity
# matrix from t to t+1. The value of nothing denotes an idenity matrix. If length of tranprobseq equals
# that of dataprobseq, the last element of tranprobseq is irrelevant. 

# `dataprobseq::AbstractVector{<: AbstractVector}`: dataprobseq[t] is the emission probability
# vector at time t. 

# """

function viterbidecode(initprob::Vector{T},
    tranprobseq::Vector{T2} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat
    tmax = length(dataprobseq)
    vprobls = [similar(i) for i in dataprobseq]
    vindexls = [zeros(Int, length(i)) for i in dataprobseq]
    vprobls[1] .= log.(initprob .* dataprobseq[1])
    for t = 2:tmax
        if isnothing(tranprobseq[t-1])        
            for j in eachindex(dataprobseq[t])                                
                val, vindexls[t][j] = vprobls[t-1][j], j
                vprobls[t][j] = val + log(dataprobseq[t][j])
            end
        else
            for j in eachindex(dataprobseq[t])
                tran2jcol = view(tranprobseq[t-1],:,j)            
                # temp = vprobls[t-1] .+ log.(tran2jcol)             
                temp = map((x,y)-> y<= 0.0 ? -Inf : x + log(y), vprobls[t-1], tran2jcol)              
                val, vindexls[t][j] = findmax(temp)
                vprobls[t][j] = val + log(dataprobseq[t][j])
            end
        end
    end
    vpath=zeros(Int,tmax)
    logprobpath,vpath[end]=findmax(vprobls[end])
    for t=tmax:-1:2
        vpath[t-1]= vindexls[t][vpath[t]]
    end
    (logprob=logprobpath, path=vpath, his_logprob=vprobls, his_index=vindexls)
end


function viterbidecode(initprob::Vector{Vector{T}},
    tranprobseq::Vector{Vector{T2}} where T2 <: Union{Nothing,Matrix{T}},
    dataprobseq::Vector{Vector{T}}) where T <: AbstractFloat    
    tmax = length(dataprobseq)
    vprobls = [similar(i) for i in dataprobseq]
    vindexls = [zeros(Int, length(i)) for i in dataprobseq]    
    vprobls[1] .= log.(kron(initprob...) .* dataprobseq[1])
    for t = 2:tmax
        tranprob = [i[t-1] for i=tranprobseq]
        if any(isnothing.(tranprob))            
            for j in eachindex(dataprobseq[t])                                
                val, vindexls[t][j] = vprobls[t-1][j], j
                vprobls[t][j] = val + log(dataprobseq[t][j])
            end
        else            
            tranprob_kron = kron(tranprob...)
            for j in eachindex(dataprobseq[t])
                tran2jcol = view(tranprob_kron,:,j)            
                # temp = vprobls[t-1] .+ log.(tran2jcol)             
                temp = map((x,y)-> y<= 0.0 ? -Inf : x + log(y), vprobls[t-1], tran2jcol)              
                val, vindexls[t][j] = findmax(temp)
                vprobls[t][j] = val + log(dataprobseq[t][j])
            end
        end
    end
    vpath=zeros(Int,tmax)
    logprobpath,vpath[end]=findmax(vprobls[end])
    for t=tmax:-1:2
        vpath[t-1]= vindexls[t][vpath[t]]
    end
    (logprob=logprobpath, path=vpath, his_logprob=vprobls, his_index=vindexls)
end


