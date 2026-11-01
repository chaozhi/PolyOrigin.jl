using Test
using PolyOrigin
exportall(PolyOrigin)

##################################################################

na = 3
ls = [[i,j,k] .- 1 for k in 1:na for j in 1:k for i in 1:j]
ils = genoindex_2_index.(ls; iszerobased=true)
@test eachindex(ls) == ils .+ 1
[ils ls]
ploidy = 3
parsedls = index_2_genoindex.(ils, ploidy, na; iszerobased=true)
@test ls == parsedls


##################################################################

ls = [[i,j,k] for k in 1:na for j in 1:k for i in 1:j]
ils = genoindex_2_index.(ls; iszerobased=false)
@test eachindex(ls) == ils
[ils ls]
ploidy = 3
parsedls = index_2_genoindex.(ils, ploidy, na; iszerobased=false)
@test ls == parsedls

##################################################################

na = 3
ls = [[i,j,k] .- 1 for k in 1:na for j in 1:na for i in 1:na]
ls2 = [[i,j,k] .- 1 for i in 1:na, j in 1:na, k in 1:na]
ls3 = vec(ls2)
@test ls == ls3
ils = diploindex_2_index.(ls,na; iszerobased=true)
@test eachindex(ls) == ils .+ 1
parsedls = index_2_diploindex.(ils,ploidy,na; iszerobased=true)
@test ls == parsedls

##################################################################

na = 3
ls = [[i,j,k]  for k in 1:na for j in 1:na for i in 1:na]
ls2 = [[i,j,k] for i in 1:na, j in 1:na, k in 1:na]
ls3 = vec(ls2)
@test ls == ls3
ils = diploindex_2_index.(ls,na; iszerobased=false)
[ils ls]
@test eachindex(ls) == ils 
parsedls = index_2_diploindex.(ils,ploidy,na; iszerobased=false)
@test ls == parsedls


##################################################################