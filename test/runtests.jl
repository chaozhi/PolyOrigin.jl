
# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
using PolyOrigin
using Test

function calacctau(dataid::String, isgbs::Bool, isrefine::Bool)
    genofile =
        string(dataid, "_polyorigin_geno_", isgbs ? "GBS.csv" : "snparray.csv")
    pedfile = string(dataid, "_polyorigin_pedigree.csv")
    polyancestry = polyOrigin(
        genofile,
        pedfile,
        maxrun = 3,
        refinemap = isrefine,
        refineorder = isrefine,
        outstem = nothing,
        logfile = nothing,
        verbose = false,
    )
    truefile = string(dataid, "_polyorigin_truevalue_ancestral.csv")
    truegeno = readTruegeno!(truefile, polyancestry)
    acc = calAccuracy!(truegeno, polyancestry)
    tau = PolyOrigin.calmapkendall(truegeno.truemap, truegeno.estmap)
    acc, tau, polyancestry
end

@testset "PolyOrigin" begin
    @testset "4x" begin
        @time include("test_4x_diallel.jl")
    end
    @testset "2x" begin
        @time include("test_2x_diallel.jl")
    end
end
