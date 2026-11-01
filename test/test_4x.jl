
cd(@__DIR__)
dataid = "4x_3waycross"
isrefine = true
acc = caltestacc(dataid,isrefine)



for key in setdiff(keys(acc),["markersum"])
    nunphasederr, ngeno, nheteroerr, npermerr, nopp = acc[key]
    @test nunphasederr/ngeno < 0.2
    @test nheteroerr/(nopp+1) < 0.2
    @test npermerr/nopp < 0.05
end
