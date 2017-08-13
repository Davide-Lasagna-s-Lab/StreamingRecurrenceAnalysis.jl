export isminimum,
       snapshot,
       ishift,
       jshift,
       distance,
       meta,
       streamdistmat,
       recurrences,
       entries

# ~~~ MAIN OUTPUT OF ITERATION ~~~
struct DistInfo{X, T, D}
       x::X
      Δj::Int
      Δi::Int
    dist::NTuple{3,NTuple{3,T}}
    meta::NTuple{3,NTuple{3,D}}
end

isminimum(d::DistInfo) = _isminimum(d.dist)
 snapshot(d::DistInfo) = d.x
   jshift(d::DistInfo) = d.Δj
   ishift(d::DistInfo) = d.Δi
 distance(d::DistInfo) = _centre(d.dist)
     meta(d::DistInfo) = _centre(d.meta)

# helper functions on 3 by 3 tuples
@inline function _isminimum(tup::NTuple{3,NTuple{3}})
    (a, b, c), (d, e, f), (g, h, i) = tup
    e < min(min(a, b, c, d), min(f, g, h, i))
end
@inline _centre(tup::NTuple{3,NTuple{3}}) = tup[2][2]

# ~~~ View over entries of distance matrix ~~~
mutable struct DistMatrixView{F,T,D,X,S<:StreamView{X}} <: AbstractVector{DistInfo{X, T, D}}
     distfun::F
        dist::Vector{Vector{T}}
        meta::Vector{Vector{D}}
      window::S
           x::S
    ΔminΔmax::UnitRange{Int}
          Δi::Int
end

# outer constructor
DistMatrixView(distfun::F,
               dist::Vector{Vector{T}},
               meta::Vector{Vector{D}},
               window::S,
               x::S,
               ΔminΔmax::UnitRange,
               Δi::Int) where {F, T, D, S<:StreamView{X}} where {X} =
    DistMatrixView{F, T, D, X, S}(distfun, dist, meta, window, x, ΔminΔmax, Δi)

unpack(x, rest...) = (x, rest)

@inline function _kernel(dist, meta, x, window, distfun)
    # use threads here
    for i in eachindex(dist)
        dist[i], meta[i] = unpack(distfun(x, window[i])...)
    end
    dist, meta
end

function step!(dmv::DistMatrixView)
    d, m = _kernel(shift!(dmv.dist),
                   shift!(dmv.meta),
                   last(step!(dmv.x)),
                   step!(dmv.window),
                   dmv.distfun)
    push!(dmv.dist, d)
    push!(dmv.meta, m)
    # update current index
    dmv.Δi += 1
    dmv
end

# ~~~ Iterator interface ~~~
Base.indices(dmv::DistMatrixView) = (dmv.ΔminΔmax, )
@inline function Base.getindex(dmv::DistMatrixView, Δj::Int)
    checkbounds(dmv, Δj)
    # Δ is the shift, so we have to calculate the 
    # appropriate index in the dist and meta vectors
    j = Δj - first(dmv.ΔminΔmax)+2
    @inbounds begin
        dist = ((dmv.dist[1][j+1], dmv.dist[2][j+1], dmv.dist[3][j+1]),
                (dmv.dist[1][j  ], dmv.dist[2][j  ], dmv.dist[3][j  ]),
                (dmv.dist[1][j-1], dmv.dist[2][j-1], dmv.dist[3][j-1]))
        meta = ((dmv.meta[1][j+1], dmv.meta[2][j+1], dmv.meta[3][j+1]),
                (dmv.meta[1][j  ], dmv.meta[2][j  ], dmv.meta[3][j  ]),
                (dmv.meta[1][j-1], dmv.meta[2][j-1], dmv.meta[3][j-1]))
    end
    # return current state
    DistInfo(dmv.x[end-1], Δj, dmv.Δi, dist, meta)
end

# ~~~ Iteration over slices of distance matrix ~~~
struct StreamDistMatrix{DMV<:DistMatrixView}
    distmatv::DMV # view over the distance matrix
           N::Int # number of views to generate
end

function streamdistmat(g, x₀, distfun, ΔminΔmax::UnitRange, N::Int)
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
    # obtain the type of the distance and meta information by computing
    # the distance between to random snapshots. Then allocate.
      d, m = unpack(distfun(x[1], window[1])...)
      T, D = typeof(d), typeof(m)
      dist = Vector{T}[Vector{T}(width) for i = 1:3]
      meta = Vector{D}[Vector{D}(width) for i = 1:3]
    # instantiate the view and then the StreamDistMatrix object
    StreamDistMatrix(DistMatrixView(distfun, dist, meta, window, x, ΔminΔmax, 2), N)
end

@inline Base.start(sdm::StreamDistMatrix) = (step!(sdm.distmatv);
                                             step!(sdm.distmatv); 1)

@inline Base.next(sdm::StreamDistMatrix, state) = (step!(sdm.distmatv), state+1)
@inline Base.done(sdm::StreamDistMatrix, state) = state == sdm.N+1


# ~~~ Iteration over the entries of the distance matrix ~~~
struct StreamDistMatrixEntries{I}
    itr::I
    normalise::Bool
    shifts::Tuple{Int, Int}
end

entries(R::StreamDistMatrix, normalise::Bool=true) = 
    StreamDistMatrixEntries(Iterators.flatten(R), 
                            normalise, 
                            (4, first(R.distmatv.ΔminΔmax)-1))

@inline Base.start(s::StreamDistMatrixEntries) = start(s.itr)
@inline function Base.next(s::StreamDistMatrixEntries, state) 
    # update iterator
    val, state = next(s.itr, state)
    # unpack DistInfo and return items of interest
    Δi = s.normalise ? ishift(val) - s.shifts[1] : ishift(val)
    Δj = s.normalise ? jshift(val) - s.shifts[2] : jshift(val)
    (Δi, Δj, distance(val), meta(val)), state
end
@inline Base.done(s::StreamDistMatrixEntries, state) = done(s.itr, state)