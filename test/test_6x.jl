
cd(@__DIR__)
dataid = "6x_3waycross"
isrefine = false
acc = caltestacc(dataid,isrefine;verbose=true)


for key in setdiff(keys(acc),["markersum"])
    nunphasederr, ngeno, nheteroerr, npermerr, nopp = acc[key]
    @test nunphasederr/ngeno < 0.2
    @test nheteroerr/(nopp+1) < 0.2
    @test npermerr/nopp < 0.05
end
acc