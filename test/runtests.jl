using PolyOrigin
using Test
cd(@__DIR__)

function caltestacc(dataid::String, isrefine::Bool; verbose::Bool=false)    
    genofile =string(dataid, "_geno.vcf.gz")
    pedfile = string(dataid, "_ped.csv")
    truegenofile = string(dataid, "_truegeno.vcf.gz")
    trueorigfile = string(dataid, "_trueorig.vcf.gz")  
    joinedancestry = last(polyorigin(
        genofile,
        pedfile;
        refhapfile = truegenofile, 
        trueorigfile, 
        maxrun = 1,
        isrefine, 
        outstem = nothing,
        logfile = nothing,
        verbose,
    ))
    acc = calaccuracy!(truegenofile, trueorigfile, joinedancestry)
    # @info string(dataid) acc
    acc
end

@testset "PolyOrigin" begin
    @testset "coding" begin
        @time include("test_coding.jl")
    end
    @testset "2x" begin
        @time include("test_2x.jl")
    end
    @testset "4x" begin
        @time include("test_4x.jl")
    end
    # @testset "6x" begin
    #     @time include("test_6x.jl")
    # end    
end

