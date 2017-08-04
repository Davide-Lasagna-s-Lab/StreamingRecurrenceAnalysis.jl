using Base.Test
using StreamingRecurrenceAnalysis

@testset "iteration protocol                     " begin
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

@testset "error checking                         " begin
    @test_throws ArgumentError snapshot_stream_view((), [1], 5, 1)
    @test_throws ArgumentError snapshot_stream_view((), [1], 5, 2)
    @test_throws ArgumentError snapshot_stream_view((), [1], 5, 3)
end

@testset "length                                 " begin
    @test length(snapshot_stream_view((), [1], 3, 2)) == 1
    @test length(snapshot_stream_view((), [1], 3, 3)) == 2
    @test length(snapshot_stream_view((), [1], 3, 4)) == 3
end