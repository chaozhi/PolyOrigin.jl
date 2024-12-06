cd(@__DIR__)
dataid = "4x_diallel"
isgbs = false
isrefine = false
acc,tau, polyancestry = calacctau(dataid,isgbs,isrefine)

@test tau[1] >= 0.9
@test acc.ndoseerr <=1
@test acc.nphaseerr<=3
@test acc.assignerr<0.09
@test acc.callerr<0.09
