using PolyOriginSim
cd(@__DIR__)

for ploidy = [2,4,6]
    outstem = string(ploidy, "x_3waycross")
    fhaplofile = PolyOriginSim.simfhaplo(;
        ploidy,
        chrlenls = 100, 
        nmarkerls = 150, 
        founderls = string.("P",1:3), 
        outstem, 
    )
    pedfile = PolyOriginSim.simnwaycross(;     
        ploidy,
        founderls = string.("P",1:3),
        famsizels = 5 * ploidy,
        outstem,    
    )
    polyoriginsim(fhaplofile,pedfile; 
        isphysmap = true, 
        recomrate  = 1, 
        isgbs = ploidy == 4, 
        outstem
    )
end

rm.(filter(x-> occursin(r"_polyoriginsim.log$",x),readdir()))
rm.(filter(x-> occursin(r"_fhaplo.vcf.gz$",x),readdir()))

