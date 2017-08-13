using BenchmarkTools
using Base.Test
using  StreamingRecurrenceAnalysis
import StreamingRecurrenceAnalysis: _isminimum, 
                                    _centre,
                                    unpack,
                                    step!


@testset "Tile                                   " begin
    t = ((1, 2, 3),
         (1, 0, 3),
         (1, 2, 3))
    @test _isminimum(t) == true
    @test _centre(t) == 0

    t = ((1, 0, 3),
         (1, 4, 3),
         (1, 2, 3))
    @test _isminimum(t) == false
    @test _centre(t) == 4

    t = ((1, 1, 1),
         (1, 1, 1),
         (1, 1, 1))
    @test _isminimum(t) == false
    @test _centre(t) == 1

    t = ((1, 2, 1),
         (1, 1, 1),
         (1, 1, 1))
    @test _isminimum(t) == false
    @test _centre(t) == 1
end

# truncated logistic map
struct LogisticMap
    r::Float64
end
(k::LogisticMap)(x) = round(x*k.r*(1-x), 5)
# # GENERATE DATA
# x = 0.1
# k = LogisticMap(4)
# for i = 1:25
#     println(x)
#     x = k(x)
# end
const DATA = [0.10000, 0.36000, 0.92160, 0.28901, 0.82193, 
              0.58544, 0.97080, 0.11339, 0.40213, 0.96169, 
              0.14737, 0.50261, 0.99997, 0.00012, 0.00048,
              0.00192, 0.00767, 0.03044, 0.11805, 0.41646, 
              0.97208, 0.10856, 0.38710, 0.94901, 0.19356]

@testset "DistMatrixView                         " begin
    # the distance function plus meta information
    dist(x, y) = (abs(x-y), x^2)

    @testset "streamviews                        " begin
        # test minimum shift
        @test_throws ArgumentError streamdistmat(LogisticMap(4), 0.1, dist, 0:5, 2)

        # test numbers in the views are consistent with what generated above
        R = streamdistmat(LogisticMap(4), 0.1, dist, 1:2, 2)
        @test R.distmatv.x      == DATA[1:3]
        @test R.distmatv.window == DATA[3:6]
        step!(R.distmatv)
        @test R.distmatv.x      == DATA[2:4]
        @test R.distmatv.window == DATA[4:7]

        # another case
        R = streamdistmat(LogisticMap(4), 0.1, dist, 2:4, 2)
        @test R.distmatv.x      == DATA[1:3]
        @test R.distmatv.window == DATA[4:8]
        step!(R.distmatv)
        @test R.distmatv.x      == DATA[2:4]
        @test R.distmatv.window == DATA[5:9]

        # test lengths of views
        R = streamdistmat(LogisticMap(4), 0.1, dist, 1:5, 2)
        @test length(R.distmatv.x) == 3
        @test length(R.distmatv.window) == 7

        R = streamdistmat(LogisticMap(4), 0.1, dist, 6:7, 2)
        @test length(R.distmatv.x) == 3
        @test length(R.distmatv.window) == 4

        # test consistency of the shift
        R = streamdistmat(LogisticMap(4), 0.1, dist, 1:5, 2)
        @test R.distmatv.x[3] == R.distmatv.window[1]
        step!(R.distmatv)
        @test R.distmatv.x[3] == R.distmatv.window[1]

        # another case, where x lags window by one 
        R = streamdistmat(LogisticMap(4), 0.1, dist, 2:5, 2)
        val = R.distmatv.window[1]
        step!(R.distmatv)
        @test R.distmatv.x[3] == val

        # a further case, where x lags window by three
        R = streamdistmat(LogisticMap(4), 0.1, dist, 4:5, 2)
        val = R.distmatv.window[1]
        step!(R.distmatv); step!(R.distmatv); step!(R.distmatv)
        @test R.distmatv.x[3] == val
    end

    @testset "distance function plus meta        " begin
        R = streamdistmat(LogisticMap(4), 0.1, dist, 1:2, 2)

        # fill distance matrices
        step!(R.distmatv)
        step!(R.distmatv)
        step!(R.distmatv)

        # test eltypes 
        @test eltype(R.distmatv.dist) == Vector{Float64}
        @test eltype(R.distmatv.meta) == Vector{Tuple{Float64}}

        # this is what should have been computed in the three steps
        # the j is the shift applied to the data. Extract the dist
        # part and discard the meta information
        d = [[ dist(DATA[i], DATA[i+j])[1]   for j = 0:3] for i = 4:6]
        m = [[(dist(DATA[i], DATA[i+j])[2],) for j = 0:3] for i = 4:6]
        @test d == R.distmatv.dist
        @test m == R.distmatv.meta

        # another example
        R = streamdistmat(LogisticMap(4), 0.1, dist, 2:4, 2)

        # fill distance matrices
        step!(R.distmatv)
        step!(R.distmatv)
        step!(R.distmatv)

        # this is what should have been computed in the three steps
        d = [[ dist(DATA[i], DATA[i+j])[1]   for j = 1:5] for i = 4:6]
        m = [[(dist(DATA[i], DATA[i+j])[2],) for j = 1:5] for i = 4:6]
        @test d == R.distmatv.dist
        @test m == R.distmatv.meta
    end

    @testset "iteration over the view            " begin
        for ΔminΔmax in [1:4, 2:4]
            R = streamdistmat(LogisticMap(4), 0.1, dist, ΔminΔmax, 1)
            # start up the view properly
            r = R.distmatv; step!(r); step!(r); step!(r)
            # this is the current state
            for c in [5, 6, 7, 8, 9] 
                for (i, el) in enumerate(r)
                    # at this point the current state is the fifth
                    @test snapshot(el) == DATA[c]
                    # the shift is what given via ΔminΔmax
                    @test shift(el) == ΔminΔmax[i]
                    # check distance and meta
                    @test distance(el) ==  dist(DATA[c], DATA[c+shift(el)])[1]
                    @test     meta(el) == (dist(DATA[c], DATA[c+shift(el)])[2], )
                end
                # shift and check again
                step!(r)
            end
        end

        # checkbounds
        ΔminΔmax = 2:4
        R = streamdistmat(LogisticMap(4), 0.1, dist, ΔminΔmax, 1)
        r = R.distmatv
        @test_throws BoundsError r[1]
        @test_throws BoundsError r[5]
    end
end

@testset "Distance Matrix Entries                " begin
    # the distance function plus meta information
    dist(x, y) = (abs(x-y), x^2)

    # allowed shifts
    for ΔminΔmax in [1:3, 5:6]

        # streaming object
        R = streamdistmat(LogisticMap(4), 0.1, dist, ΔminΔmax, 10)

        # fill with output
        out_d = Matrix{Float64}(       length(ΔminΔmax), 10)
        out_m = Matrix{Tuple{Float64}}(length(ΔminΔmax), 10)

        # indices are automatically shifted
        for (i, j, d, m) in entries(R, true)
            out_d[j, i] = d
            out_m[j, i] = m
        end

        # this is the expected output. We always start at five, because
        # we need three shifts from 2 to obtain the first full slice.
        # The streaming implementation loops over the entries of this
        # matrix, avoid the full storage of hyper-long simulations.
        # i runs from five to 14, because we asked for ten shifts
        d = [ dist(DATA[i], DATA[i+j])[1]   for j = ΔminΔmax, i = 5:14]
        m = [(dist(DATA[i], DATA[i+j])[2],) for j = ΔminΔmax, i = 5:14]
    
        @test out_d == d
        @test out_m == m
    end
end