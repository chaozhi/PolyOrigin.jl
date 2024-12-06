cd(@__DIR__)
dataid = "4x_diallel"
isgbs = false
isrefine = true
acc,tau, polyancestry = calacctau(dataid,isgbs,isrefine)

# ancestryfile = string(dataid,"_polyancestry.csv")
# savePolyAncestry(ancestryfile,polyancestry)
# polyancestry2 = readPolyAncestry(ancestryfile)
# fn = fieldnames(PolyAncestry)
# fn = (fn[1:7]...,fn[9:end]...)
# for i = fn
#     @test getfield(polyancestry,i) == getfield(polyancestry2,i)
# end
# rm(ancestryfile)
@test tau[1] >= 0.9
@test acc.ndoseerr <=1
@test acc.nphaseerr<=3
@test acc.assignerr<0.09
@test acc.callerr<0.09
