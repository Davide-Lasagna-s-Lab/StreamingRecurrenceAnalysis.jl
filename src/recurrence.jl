export streamdistmat, recurrences, entries, distance, meta, StreamDistMatrix

# ~~~ DISTANCE INFO BETWEEN TWO STATES ~~~
struct DistInfo{T, D}
    dist::NTuple{3,NTuple{3,T}}
    meta::D
end

isrecurrence(d::DistInfo) = _isrecurrence(d.dist)
        dist(d::DistInfo) = d.dist[2][2]
        meta(d::DistInfo) = d.meta
        pack(d::DistInfo) = (dist(d), meta(d)...)

# helper functions on 3 by 3 tuples
@inline function _isrecurrence(tup::NTuple{3,NTuple{3}})
    (a, b, c), (d, e, f), (g, h, i) = tup
    e < min(min(a, b, c, d), min(f, g, h, i))
end

# ~~~ VIEW OVER ENTRIES OF dist MATRIX ~~~
struct DistMatrixView{T,D,F}
    distfun::F
       dist::Vector{Vector{T}}
       meta::Vector{Vector{D}}
end

# Update dist data with new current state and window
function update!(dmv::DistMatrixView, x::X, window::StreamView{X}) where {X}
    d, m = _kernel(shift!(dmv.dist), shift!(dmv.meta),
                   x, window, dmv.distfun)
    push!(dmv.dist, d)
    push!(dmv.meta, m)
    dmv
end

# update dist information with new snapshots
function _kernel(dist, meta, x::X, window::StreamView{X}, distfun) where {X}
    # use threads here
    for i in eachindex(dist)
        dist[i], meta[i] = unpack(distfun(x, window[i])...)
    end
    dist, meta
end

unpack(x, rest...) = (x, rest)

# Iterator interface 
Base.start(dmv::DistMatrixView) = 2
Base.done(dmv::DistMatrixView, j) = j == length(dmv.dist[1])
function Base.next(dmv::DistMatrixView, j)
    @inbounds begin
        dist = ((dmv.dist[1][j+1], dmv.dist[2][j+1], dmv.dist[3][j+1]),
                (dmv.dist[1][j  ], dmv.dist[2][j  ], dmv.dist[3][j  ]),
                (dmv.dist[1][j-1], dmv.dist[2][j-1], dmv.dist[3][j-1]))
    end
    DistInfo(dist, dmv.meta[2][j]), j+1
end

# ~~~ ITERATION OVER SLICES OF dist MATRIX ~~~
struct StreamDistMatrix{T,D,X,F,DMV<:DistMatrixView{T,D,F},S1<:StreamView{X},S2<:StreamView{X}}
    distmatv::DMV            # current view over dist matrix
    ΔminΔmax::UnitRange{Int} # range of shifts
      window::S1             # view over future snapshots
           x::S2             # view over current snapshots
           N::Int            # number of views to generate
end

function streamdistmat(g, x₀::X, distfun, ΔminΔmax::UnitRange, N::Int) where {X}
      Δmax = last(ΔminΔmax)
      Δmin = first(ΔminΔmax)
    Δmin > 0 || throw(ArgumentError("Δmin must be positive, got $Δmin"))
    # make view large enough so that we can find minima with
    # shifts between Δmin and Δmax, included
     width = Δmax-Δmin+3
         x = streamview(g, copy(x₀), 3)
    # shift forward window ahead such that the first minimum
    # can be found for a discrete shift equal to Δmin.
    window = streamview(g, copy(x₀), width); for i = 1:Δmin+1; step!(window); end
    # obtain the type of the dist and meta information, then allocate.
      d, m = unpack(distfun(x[1], window[1])...)
      T, D = typeof(d), typeof(m)
      dist = Vector{T}[Vector{T}(width) for i = 1:3]
      meta = Vector{D}[Vector{D}(width) for i = 1:3]
    # instantiate the view and then the StreamDistMatrix object
       dmv = DistMatrixView(distfun, dist, meta)
    StreamDistMatrix{T, D, X, typeof(distfun),
                     typeof(dmv), typeof(window), 
                     typeof(x)}(dmv, ΔminΔmax, window, x, N)
end

# Iterator interface
function Base.start(sdm::StreamDistMatrix)
    update!(sdm.distmatv, last(      sdm.x),        sdm.window)
    update!(sdm.distmatv, last(step!(sdm.x)), step!(sdm.window))
    # this is the index of the first state 
    # for which we can determine recurrences
    return 4 
end

function Base.next(sdm::StreamDistMatrix, Δi)
    update!(sdm.distmatv, last(step!(sdm.x)), step!(sdm.window))
    # elements of the flattened iterator will be (x, Δi, Δj, dinfo)
    return zip(Iterators.repeated(sdm.x[end-1]), 
               Iterators.repeated(Δi), 
               sdm.ΔminΔmax, 
               sdm.distmatv), Δi+1
end

# generate N views in total
Base.done(sdm::StreamDistMatrix, Δi) = Δi == sdm.N+4

# Fill dist matrix (used mainly for plotting?)
function Base.full(R::StreamDistMatrix{T, D}) where {T, D}
    Δmax, Δmin = last(R.ΔminΔmax), first(R.ΔminΔmax)
    shape = (Δmax-Δmin+1, R.N)
    d, m = Matrix{T}(shape), Matrix{D}(shape)
    for (i, r) in enumerate(R)
        for (j, (x, Δi, Δj, dinfo)) in enumerate(r)
            d[j, i] = dist(dinfo)
            m[j, i] =     meta(dinfo)
        end
    end
    d, m
end

# ~~~ Iterator over recurrences  ~~~
function recurrences(R::StreamDistMatrix, predicate::Function=x->true)
    itr = Iterators.filter(Iterators.flatten(R)) do args
        x, Δi, Δj, dinfo = args
        isrecurrence(dinfo) && predicate(pack(dinfo)) 
    end
    # transform
    ((rec[1], rec[2], rec[3], pack(rec[4])) for rec in itr)
end