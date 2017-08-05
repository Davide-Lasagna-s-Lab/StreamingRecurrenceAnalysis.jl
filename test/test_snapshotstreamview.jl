using Base.Test
using StreamingRecurrenceAnalysis

@testset "iteration protocol - width != 1        " begin
    g(x₀) = (x₀ .+= 1; x₀)
    output = ([[3], [2], [1]],
              [[4], [3], [2]],
              [[5], [4], [3]],
              [[6], [5], [4]],
              [[7], [6], [5]],
              [[8], [7], [6]],
              [[9], [8], [7]])
    
    sview = snapshot_stream_view(g, [1], 3, 2)
    for (i, v) in enumerate(sview)
        @test v == output[i]
        @test i <= 1 # only one view is generated
    end

    sview = snapshot_stream_view(g, [1], 3, 3)
    for (i, v) in enumerate(sview)
        @test v == output[i]
        @test i <= 2 # only two views are generated
    end

    sview = snapshot_stream_view(g, [1], 3, 8)
    for (i, v) in enumerate(sview)
        @test v == output[i]
        @test i <= 7 # only seven views are generated
    end
end

@testset "iteration protocol - window width = 1  " begin
    g(x₀) = (x₀ .+= 1; x₀)
    output = ([[0]], [[1]], [[2]], [[3]], [[4]], [[5]], [[6]], [[7]])
    for N in 0:7
        sview = snapshot_stream_view(g, [0], 1, N)
        for (i, v) in enumerate(sview)
            @test v == output[i]
            @test i <= N+1 # N+1 views are generated
        end
    end
end

@testset "error checking                         " begin
    @test_throws ArgumentError snapshot_stream_view((), [1],  5, 1)
    @test_throws ArgumentError snapshot_stream_view((), [1],  5, 2)
    @test_throws ArgumentError snapshot_stream_view((), [1],  5, 3)
    @test_throws ArgumentError snapshot_stream_view((), [1],  0, 0)
    @test_throws ArgumentError snapshot_stream_view((), [1],  0, 1)
    @test_throws ArgumentError snapshot_stream_view((), [1], -1, 1)
end

@testset "length                                 " begin
    g(x₀) = (x₀ .+= 1; x₀)
    for width = 1:5
        for N = (width-1):10
            sview = snapshot_stream_view(g, [0], width, N)
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
    sview = snapshot_stream_view(g, x, 20, 100)

    # sum
    function dowork!(sview::SnapshotStreamView{X}, u::X) where {X}
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