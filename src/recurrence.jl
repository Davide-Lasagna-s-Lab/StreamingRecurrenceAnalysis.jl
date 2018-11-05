export streamdistmat, recurrences, StreamDistMatrix, repeatf

# repeat function function f n times on x using recursion
repeatf(f, x, n::Int) = n > 1 ? f(repeatf(f, x, n-1)) : f(x)

# ~~~ VIEW OVER ENTRIES OF dist MATRIX ~~~
struct DistMatrixView{D, F}
     distfun::F
    dist::Vector{Vector{D}}
    DistMatrixView(distfun::F, ::Type{D}, width::Int) where {F, D} =
        new{D, F}(distfun, Vector{D}[Vector{D}(width) for i = 1:3])
end

# Update dist data with new current state and window
update!(dmv::DistMatrixView, x::X, window::StreamView{X}) where {X} =
    push!(dmv.dist, _kernel(shift!(dmv.dist), x, window, dmv.distfun))

# update dist information with new snapshots
function _kernel(dist, x::X, window::StreamView{X}, distfun) where {X}
    # use threads here
    for i in eachindex(dist)
        dist[i] = distfun(x, window[i])
    end
    dist
end

# Iterator interface
Base.start(dmv::DistMatrixView) = 2
Base.done(dmv::DistMatrixView, j) = j == length(dmv.dist[1])
function Base.next(dmv::DistMatrixView, j)
    @inbounds begin
        # pre check we have a minimum
        dist = dmv.dist
        isrec = _isrecurrence(((dist[1][j+1][1], dist[2][j+1][1], dist[3][j+1][1]),
                               (dist[1][j  ][1], dist[2][j  ][1], dist[3][j  ][1]),
                               (dist[1][j-1][1], dist[2][j-1][1], dist[3][j-1][1])))
    end
    (dmv.dist[2][j], isrec), j+1
end

# helper functions on 3 by 3 tuples
function _isrecurrence(tup::NTuple{3,NTuple{3}})
    (a, b, c), (d, e, f), (g, h, i) = tup
    e < minimum((a, b, c, d, f, g, h, i))
end

# ~~~ ITERATION OVER SLICES OF dist MATRIX ~~~
struct StreamDistMatrix{D,X,DMV<:DistMatrixView{D},S1<:StreamView{X},S2<:StreamView{X}}
         dmv::DMV            # current view over dist matrix
    ΔminΔmax::UnitRange{Int} # range of shifts
      window::S1             # view over future snapshots
           x::S2             # view over current snapshots
           N::Int            # number of views to generate
end

function streamdistmat(g, x₀::X, distfun, ΔminΔmax::UnitRange, N::Int) where {X}
    Δmin, Δmax = extrema(ΔminΔmax)
    Δmin > 0 || throw(ArgumentError("Δmin must be positive, got $Δmin"))
    # this is the sliding view of width 3, i.e. the minimum required
    # to determine whether we have a local minimum. Essentially, this
    # view is centered around the snapshot at time t, so that we will 
    # find recurrences with snapshots at time t+T
    x = streamview(g, copy(x₀), 3)
    # make the width of the second streamview large enough so that we
    # can find minima with shifts between Δmin and Δmax, included. This
    # view contains "future" snapshots, that are compare to snapshots in
    # the view x. Note we shift forward the window  such that the first 
    # recurrence can be found for a discrete shift equal to Δmin.
    window = streamview(g, copy(x₀), Δmax-Δmin+3)
    for i = 1:Δmin+1; step!(window); end
    # obtain the type of the distance information, then allocate. 
    # D could be just a float, i.e. just a distance, or there could 
    # be other meta information, e.g. the required translation for
    # problems with continuous symmetries
    D = typeof(distfun(x[1], window[1]))
    # instantiate the view and then the StreamDistMatrix object
    dmv = DistMatrixView(distfun, D, length(window))
    StreamDistMatrix{D, X, typeof(dmv), typeof(window),
                       typeof(x)}(dmv, ΔminΔmax, window, x, N)
end

# Iterator interface
function Base.start(sdm::StreamDistMatrix)
    update!(sdm.dmv, last(      sdm.x),        sdm.window)
    update!(sdm.dmv, last(step!(sdm.x)), step!(sdm.window))
    # this is the index of the first state
    # for which we can determine recurrences
    return 3
end

function Base.next(sdm::StreamDistMatrix, Δi)
    update!(sdm.dmv, last(step!(sdm.x)), step!(sdm.window))
    # elements of the flattened iterator will be (x, Δi, Δj, (d, isrec))
    return zip(Iterators.repeated(copy(sdm.x[end-1])),
               Iterators.repeated(Δi),
               sdm.ΔminΔmax,
               sdm.dmv), Δi+1
end

# generate N views in total
Base.done(sdm::StreamDistMatrix, Δi) = Δi == sdm.N+3

# Fill dist matrix (used mainly for plotting?)
function Base.full(R::StreamDistMatrix{D}) where {D}
    Δmax, Δmin = last(R.ΔminΔmax), first(R.ΔminΔmax)
    dist = Matrix{D}(Δmax-Δmin+1, R.N)
    for (i, r) in enumerate(R)
        for (j, (x, Δi, Δj, (d, isrec))) in enumerate(r)
            dist[j, i] = d
        end
    end
    dist
end

# ~~~ Iterator over recurrences  ~~~
struct Recurrences{D,X,I}
    itr::I
end

function recurrences(R::StreamDistMatrix{D,X}, predicate::Function=x->true) where {D,X}
    itr = Iterators.filter(Iterators.flatten(R)) do args
        x, Δi, Δj, (d, isrec) = args
        isrec && predicate(d)
    end
    Recurrences{D,X,typeof(itr)}(itr)
end

Base.start(recs::Recurrences) = start(recs.itr)
function Base.next(recs::Recurrences{D,X}, state) where {D,X}
    rec, state = next(recs.itr, state)
    (rec[1], rec[2], rec[3], rec[4][1])::Tuple{X,Int,Int,D}, state
end
Base.done(recs::Recurrences, state) = done(recs.itr, state)