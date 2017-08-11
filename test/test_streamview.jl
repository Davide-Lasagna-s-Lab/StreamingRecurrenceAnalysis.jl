using BenchmarkTools
using Base.Test
using StreamingRecurrenceAnalysis

@testset "iteration protocol - width != 1        " begin
    g(x₀) = (x₀ .+= 1; x₀)
    output = [[1], [2], [3], [4], [5], [6], [7], [8], [9]]
    
    for w in 1:4
        sview = streamview(g, [1], w)
        @test sview == output[1:w]
        for i = 2:(9-w)
            @test step!(sview) == output[i:i+(w-1)]
        end
    end
end

@testset "error checking                         " begin
    @test_throws ArgumentError streamview((), [1], 0)
end

@testset "iteration                              " begin
    g(x₀) = (x₀ .+= 1; x₀)
    sview = streamview(g, [0], 3)
    @test length(sview) == 3
    @test sview == [[0], [1], [2]]
    @test sview[1] == [0]
    @test sview[2] == [1]
    @test sview[3] == [2]
end

@testset "example usage with allocation          " begin
    
    # initial condition and temporary
    x  = randn(100)
    u  = similar(x) 

    # increment state
    g(x₀) = (x₀ .+= 1.0; x₀)

    # make stream
    sview = streamview(g, x, 200)

    # sum
    function dowork!(sview::StreamView{X}, u::X) where {X}
        for i = 1:10000
            for snap in step!(sview)
                u .+= snap # do stuff
            end
        end
        u    
    end

    # warm up
    dowork!(sview, u)
        
    @test (@allocated dowork!(sview, u)) == 0
end

struct LogisticMap
    r::Float64
end
@inline (k::LogisticMap)(x) = x*k.r*(1-x)

@testset "example usage with allocation          " begin
    
    # make stream
    sview = streamview(LogisticMap(0.4), 0.5, 200)

    # sum
    function dowork!(sview::StreamView{X}) where {X}
        s = zero(X)
        for i = 1:10000
            step!(sview)
            # this enables vectorisation and obtain 4x speed up
            @simd for i in 1:length(sview)
                @inbounds s += sview[i]
            end
        end
        s
    end

    # warm up
    dowork!(sview)
            
    # still small allocation in this test
    @test (@allocated dowork!(sview)) == 16
end