
cd(@__DIR__)
dataid = "2x_3waycross"
isrefine = true
acc = caltestacc(dataid,isrefine; verbose=false)


for key in setdiff(keys(acc),["markersum"])
    nunphasederr, ngeno, nheteroerr, npermerr, nopp = acc[key]
    @test nunphasederr/ngeno < 0.05
    @test nheteroerr/(nopp+1) < 0.05
    @test npermerr/nopp < 0.05
end
