using BenchmarkTools
using Base.Test
using StreamingRecurrenceAnalysis

@testset "iteration protocol - width != 1        " begin
    g(x₀) = (x₀ .+= 1; x₀)
    output = ([[1], [2], [3]],
              [[2], [3], [4]],
              [[3], [4], [5]],
              [[4], [5], [6]],
              [[5], [6], [7]],
              [[6], [7], [8]],
              [[7], [8], [9]])
    
    sview = streamview(g, [1], 3, 2)
    for (i, v) in enumerate(sview)
        @test v == output[i]
        @test i <= 1 # only one view is generated
    end

    sview = streamview(g, [1], 3, 3)
    for (i, v) in enumerate(sview)
        @test v == output[i]
        @test i <= 2 # only two views are generated
    end

    sview = streamview(g, [1], 3, 8)
    for (i, v) in enumerate(sview)
        @test v == output[i]
        @test i <= 7 # only seven views are generated
    end
end

@testset "iteration protocol - window width = 1  " begin
    g(x₀) = (x₀ .+= 1; x₀)
    output = ([[0]], [[1]], [[2]], [[3]], [[4]], [[5]], [[6]], [[7]])
    for N in 0:7
        sview = streamview(g, [0], 1, N)
        for (i, v) in enumerate(sview)
            @test v == output[i]
            @test i <= N+1 # N+1 views are generated
        end
    end
end

@testset "error checking                         " begin
    @test_throws ArgumentError streamview((), [1],  5, 1)
    @test_throws ArgumentError streamview((), [1],  5, 2)
    @test_throws ArgumentError streamview((), [1],  5, 3)
    @test_throws ArgumentError streamview((), [1],  0, 0)
    @test_throws ArgumentError streamview((), [1],  0, 1)
    @test_throws ArgumentError streamview((), [1], -1, 1)
end

@testset "length                                 " begin
    g(x₀) = (x₀ .+= 1; x₀)
    for width = 1:5
        for N = (width-1):10
            sview = streamview(g, [0], width, N)
            @test length(collect(sview)) == length(sview)
        end
    end
end

@testset "example usage with allocation          " begin
    
    # initial condition and temporary
    x  = randn(100)
    u  = similar(x) 

    # increment state
    g(x₀) = (x₀ .+= 1.0; x₀)

    # have a stream of 10 snapshots with a view of width 2
    sview = streamview(g, x, 20, 100)

    # sum
    function dowork!(sview::StreamView{X}, u::X) where {X}
        for v in sview     # the current view
            for snap in v  # each snapshot in the view 
                u .+= snap # do stuff
            end
        end
        u    
    end

    # warm up
    dowork!(sview, u)
        
    @test (@allocated dowork!(sview, u)) == 0
end

@testset "example usage with number              " begin
    
    # mapping
    g(x) = x + one(x)

    # have a stream of 10 snapshots with a view of width 2
    sview = streamview(g, 1, 2, 4)

    # this will have allocations, because we have a vector in the next! tuple
    function moving_average!(sview::StreamView{X}, out::Vector) where {X}
        for (i, v) in enumerate(sview)
            out[i] = mean(v)
        end
        out
    end

    out = moving_average!(sview, zeros(length(sview)))

    @test out[1] == 1.5
    @test out[2] == 2.5
    @test out[3] == 3.5
    @test out[4] == 4.5
end