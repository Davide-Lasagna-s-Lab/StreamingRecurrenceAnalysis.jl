@testset "3 by 3 tuple                           " begin
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

# the distance function plus meta information
dist(x, y) = (round(abs(x-y); digits=5), round(x^2, digits=5))

@testset "Distance Matrix Entries                " begin
    # allowed shifts
    for ΔminΔmax in [1:2, 3:10, 2:5]
        for N = [2, 5]
            # streaming object
            R = streamdistmat(LogisticMap(4), 0.1, dist, ΔminΔmax, N)
           
            # collect data
            out = full(R)

            # this is the expected output. We always start at four, because
            # we need three shifts from 2 to obtain the first full slice.
            # The streaming implementation loops over the entries of this
            # matrix, avoid the full storage of hyper-long simulations.
            # Note that the meta information is a matrix of tuples!
            d = [dist(DATA[i], DATA[i+j]) for j = ΔminΔmax, i = 4:(4+N-1)]
        
            @test out == d
        end
    end
end

@testset "recurrences                            " begin
    # real distance matrix, with additional rows/columns for visual inspection
    # d = [dist(DATA[i], DATA[i+j])[1] for j = 0:4, i = 3:9]
    # d = [dist(DATA[i], DATA[i+j])[2] for j = 0:4, i = 3:9]
    # THIS IS THE DISTANCE VALUE
    #       3         4        5        6        7        8         9
    # 0.63259 | 0.53292  0.23649  0.38536  0.85741  0.28874 | 0.55956 1
    # --------------------------------------------------------------- 
    # 0.09967 | 0.29643  0.14887  0.47205  0.56867  0.84830 | 0.25476 2 
    # 0.33616 | 0.68179  0.70854  0.18331  0.00911  0.03398 | 0.10048 3 
    # ---------------------------------------------------------------
    # 0.0492  | 0.17562  0.41980  0.37625  0.82343  0.38922 | 0.59784 4
    
    # THIS IS THE META INFORMATION
    # 0.84935 | 0.08353  0.67557  0.34274  0.94245  0.01286  | 0.16171
    # ---------------------------------------------------------------
    # 0.84935 | 0.08353  0.67557  0.34274  0.94245  0.01286  | 0.16171
    # 0.84935 | 0.08353  0.67557  0.34274  0.94245  0.01286  | 0.16171
    # ---------------------------------------------------------------
    # 0.84935 | 0.08353  0.67557  0.34274  0.94245  0.01286  | 0.16171

    # streaming object
    R = streamdistmat(LogisticMap(4), 0.1, dist, 2:3, 5)

    # expected value from table above
    expected = [(DATA[5], 4, 2, (0.14887, 0.67557)), 
                (DATA[7], 6, 3, (0.00911, 0.94245))]

    # check
    for (i, rec) in enumerate(recurrences(R))
        @test rec == expected[i]
    end

    # test with custom predicate
    predicate(dist) = dist[1] < 0.1

    # streaming object
    R = streamdistmat(LogisticMap(4), 0.1, dist, 2:3, 5)

    # expected value from table above
    expected = [(DATA[7], 6, 3, (0.00911, 0.94245))]
    
    # check
    for (i, rec) in enumerate(recurrences(R, predicate))
        @test rec == expected[i]
    end
end