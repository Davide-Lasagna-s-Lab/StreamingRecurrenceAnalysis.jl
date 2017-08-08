using BenchmarkTools
using Base.Test
import StreamingRecurrenceAnalysis: Tile, 
                                    isminimum, 
                                    TileIterator, 
                                    streamview,
                                    streamdistmat,
                                    tiles,
                                    recurrences,
                                    centre

@testset "Tile                                   " begin
    t = Tile(((1, 2, 3),
              (1, 0, 3),
              (1, 2, 3)))
    @test isminimum(t) == true

    t = Tile(((1, 0, 3),
              (1, 4, 3),
              (1, 2, 3)))
    @test isminimum(t) == false

    t = Tile(((1, 1, 1),
              (1, 1, 1),
              (1, 1, 1)))
    @test isminimum(t) == false

    t = Tile(((1, 2, 1),
              (1, 1, 1),
              (1, 1, 1)))
    @test isminimum(t) == false
end

@testset "TileIterator                           " begin
    dists = [[0, 2, 4, 0, 4, 1], 
             [0, 2, 2, 1, 8, 6],
             [0, 6, 1, 8, 3, 9]]

    # moving slice
    sl = TileIterator(dists, 1)

    # expected output
    output = [(1, Tile(((4, 2, 1),   # tile centred at index two
                        (2, 2, 6),   # for an offset equal to one
                        (0, 0, 0)))),
              (2, Tile(((0, 1, 8), 
                        (4, 2, 1), 
                        (2, 2, 6)))),
              (3, Tile(((4, 8, 3),
                        (0, 1, 8), 
                        (4, 2, 1)))),
              (4, Tile(((1, 6, 9),
                        (4, 8, 3),
                        (0, 1, 8))))]
    # typical usage             
    for (i, (offset, tile)) in enumerate(sl)
        @test output[i] == (offset, tile)
    end
end

@testset "streamdistmat                          " begin
    # generate a known stream of Ints
    counter = 1
    data = [1, 3, 5, 1, 4, 7, 9, 2, 8, 3, 1, 0]
    g(x₀) = data[counter+=1]
    dist(x, y) = abs(x-y)

    sview = streamview(g, data[1], 4, 10)

    # these will be the views generated
    views = Vector{Int}[[1, 3, 5, 1],
                        [3, 5, 1, 4],
                        [5, 1, 4, 7],
                        [1, 4, 7, 9],
                        [4, 7, 9, 2],
                        [7, 9, 2, 8],
                        [9, 2, 8, 3],
                        [2, 8, 3, 1]]
    # test they actually are                        
    @test [copy(view) for view in sview] == views             
             
    # streaming distance matrix
    R = streamdistmat(sview, dist, 1)

    # this will be the distances calculated             
    dists = [[0, 2, 4, 0], 
             [0, 2, 2, 1],
             [0, 4, 1, 2],
             [0, 3, 6, 8],
             [0, 3, 5, 2],
             [0, 2, 5, 1],
             [0, 7, 1, 6],
             [0, 6, 1, 1]]

    # these will be the triplets generated    
    output = [[[0, 2, 4, 0], 
               [0, 2, 2, 1],
               [0, 4, 1, 2]],
              [[0, 2, 2, 1],
               [0, 4, 1, 2],
               [0, 3, 6, 8]],
              [[0, 4, 1, 2],
               [0, 3, 6, 8],
               [0, 3, 5, 2]],
              [[0, 3, 6, 8],
               [0, 3, 5, 2],
               [0, 2, 5, 1]],
              [[0, 3, 5, 2],
               [0, 2, 5, 1],
               [0, 7, 1, 6]],
              [[0, 2, 5, 1],
               [0, 7, 1, 6],
               [0, 6, 1, 1]]]  

    # reset counter for data generation
    counter = 1     
    sview = streamview(g, data[1], 4, 10)
    for (i, r) in enumerate(R)
        @test output[i] == r
    end
end

@testset "tiles                                  " begin
    # generate a known stream of Ints
    counter = 1
    data = [1, 3, 5, 1, 4, 7, 9, 2, 8, 3, 1, 0]
    g(x₀) = data[counter+=1]
    dist(x, y) = abs(x-y)
    sview = streamview(g, data[1], 4, 10)
                    
    # streaming distance matrix
    R = streamdistmat(sview, dist, 1)

    # these will be the triplets generated    
    # output = [[[0, 2, 4, 0], 
    #            [0, 2, 2, 1],
    #            [0, 4, 1, 2]],
    #           [[0, 2, 2, 1],
    #            [0, 4, 1, 2],
    #            [0, 3, 6, 8]],
    #           [[0, 4, 1, 2],
    #            [0, 3, 6, 8],
    #            [0, 3, 5, 2]],
    #           [[0, 3, 6, 8],
    #            [0, 3, 5, 2],
    #            [0, 2, 5, 1]],
    #           [[0, 3, 5, 2],
    #            [0, 2, 5, 1],
    #            [0, 7, 1, 6]],
    #           [[0, 2, 5, 1],
    #            [0, 7, 1, 6],
    #            [0, 6, 1, 1]]]  
    # use this for generation
    # for o in output
    #     D = hcat(o...)
    #     display(flipdim(D[1:3, :], 1)); println()
    #     display(flipdim(D[2:4, :], 1)); println()
    # end

    ts = [(1, Tile(((4,  2,  1),
                    (2,  2,  4),
                    (0,  0,  0)))),
          (2, Tile(((0,  1,  2),
                    (4,  2,  1),
                    (2,  2,  4)))),
          (1, Tile(((2,  1,  6),
                    (2,  4,  3),
                    (0,  0,  0)))),
          (2, Tile(((1,  2,  8),
                    (2,  1,  6),
                    (2,  4,  3)))),
          (1, Tile(((1,  6,  5),
                    (4,  3,  3),
                    (0,  0,  0)))),
          (2, Tile(((2,  8,  2),
                    (1,  6,  5),
                    (4,  3,  3)))),
          (1, Tile(((6,  5,  5),
                    (3,  3,  2),
                    (0,  0,  0)))),
          (2, Tile(((8,  2,  1),
                    (6,  5,  5),
                    (3,  3,  2)))),
          (1, Tile(((5,  5,  1),
                    (3,  2,  7),
                    (0,  0,  0)))),
          (2, Tile(((2,  1,  6),
                    (5,  5,  1),
                    (3,  2,  7)))),
          (1, Tile(((5,  1,  1),
                    (2,  7,  6),
                    (0,  0,  0)))),
          (2, Tile(((1,  6,  1),
                    (5,  1,  1),
                    (2,  7,  6))))]

    for (i, (offset, tile)) in enumerate(tiles(R))
       @test ts[i] == (offset, tile)
    end
end

@testset "recurrences                            " begin
   # generate a known stream of Ints
    counter = 1
    data = [1, 2, 3, 1, 3, 4, 2, 6]
    g(x₀) = data[counter+=1]
    dist(x, y) = abs(x-y)

    # test views
    sview = streamview(g, data[1], 5, length(data)-1)

    output = [[1, 2, 3, 1, 3],
              [2, 3, 1, 3, 4],
              [3, 1, 3, 4, 2],
              [1, 3, 4, 2, 6]]

    for (i, s) in enumerate(sview)
        @test output[i] == s
    end

    # test streaming distance matrix
    counter = 1
    sview = streamview(g, data[1], 5, length(data)-1)
    R = streamdistmat(sview, dist, 1)

    output = [[[0, 1, 2, 0, 2], 
               [0, 1, 1, 1, 2], 
               [0, 2, 0, 1, 1]],
              [[0, 1, 1, 1, 2], 
               [0, 2, 0, 1, 1], 
               [0, 2, 3, 1, 5]]]

    for (i, s) in enumerate(R)
        @test output[i] == s
    end

    # test tiles
    output = [(1, Tile(((2,  1,  0),
                        (1,  1,  2),
                        (0,  0,  0)))),
              (2, Tile(((0,  1,  1),
                        (2,  1,  0),
                        (1,  1,  2)))),
              (3, Tile(((2,  2,  1),
                        (0,  1,  1),
                        (2,  1,  0)))),
              (1, Tile(((1,  0,  3),
                        (1,  2,  2),
                        (0,  0,  0)))),
              (2, Tile(((1,  1,  1),
                        (1,  0,  3),
                        (1,  2,  2)))),
              (3, Tile(((2,  1,  5),
                        (1,  1,  1),
                        (1,  0,  3))))]
    counter = 1
    sview = streamview(g, data[1], 5, length(data)-1)
    R = streamdistmat(sview, dist, 1)
    for (i, (offset, tile)) in enumerate(tiles(R))
        @test output[i] == (offset, tile)
    end

    # test recurrences
    output = [(2, Tile(((1,  1,  1),
                        (1,  0,  3),
                        (1,  2,  2))))]
    counter = 1
    sview = streamview(g, data[1], 5, length(data)-1)
    R = streamdistmat(sview, dist, 1)
    for (i, (offset, tile)) in enumerate(recurrences(R))
        @test output[i] == (x, offset)
    end
end