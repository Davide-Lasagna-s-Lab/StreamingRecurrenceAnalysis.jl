using BenchmarkTools
using Base.Test
using  StreamingRecurrenceAnalysis
import StreamingRecurrenceAnalysis: _isrecurrence, 
                                    unpack,
                                    step!

@testset "Tile                                   " begin
    t = ((1, 2, 3),
         (1, 0, 3),
         (1, 2, 3))
    @test _isrecurrence(t) == true

    t = ((1, 0, 3),
         (1, 4, 3),
         (1, 2, 3))
    @test _isrecurrence(t) == false

    t = ((1, 1, 1),
         (1, 1, 1),
         (1, 1, 1))
    @test _isrecurrence(t) == false

    t = ((1, 2, 1),
         (1, 1, 1),
         (1, 1, 1))
    @test _isrecurrence(t) == false
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

@testset "Distance Matrix Entries                " begin
    # the distance function plus meta information
    dist(x, y) = (abs(x-y), x^2)

    # allowed shifts
    for ΔminΔmax in [1:2, 3:10, 2:5]
        for N = [1, 5]
            # streaming object
            R = streamdistmat(LogisticMap(4), 0.1, dist, ΔminΔmax, N)

            # fill with output
            out_d = Matrix{Float64}(       length(ΔminΔmax), N)
            out_m = Matrix{Tuple{Float64}}(length(ΔminΔmax), N)

            # indices are automatically shifted
            for (i, j, d, m) in entries(R, true)
                out_d[j, i] = d
                out_m[j, i] = m
            end

            # this is the expected output. We always start at four, because
            # we need three shifts from 2 to obtain the first full slice.
            # The streaming implementation loops over the entries of this
            # matrix, avoid the full storage of hyper-long simulations.
            d = [ dist(DATA[i], DATA[i+j])[1]   for j = ΔminΔmax, i = 4:(4+N-1)]
            m = [(dist(DATA[i], DATA[i+j])[2],) for j = ΔminΔmax, i = 4:(4+N-1)]
        
            @test out_d == d
            @test out_m == m
        end
    end
end