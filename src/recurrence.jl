# ~~~ TILE ~~~
struct Tile{T} <: AbstractMatrix{T}
    data::NTuple{3, NTuple{3, T}}
end
Base.size(t::Tile) = (3, 3)
Base.getindex(t::Tile, i, j) = t.data[i][j] 

centre(t::Tile) = t[2, 2] 

@inline function isminimum(t::Tile) 
    (a, b, c), (d, e, f), (g, h, i) = t.data
    e < min(min(a, b, c, d), min(f, g, h, i))
end

# ~~~ Tile iterator over a slice ~~~
struct TileIterator{D, X}
              x::X
              R::Vector{Vector{D}}
    startoffset::Int
end

@inline Base.start(s::TileIterator)        = s.startoffset
@inline Base.done( s::TileIterator, state) = state == length(s.R[1])-1
@inline function Base.next(s::TileIterator, state)
    @inbounds begin
        a = s.R[1][state+2]; b = s.R[2][state+2]; c = s.R[3][state+2]
        d = s.R[1][state+1]; e = s.R[2][state+1]; f = s.R[3][state+1]
        g = s.R[1][state  ]; h = s.R[2][state  ]; i = s.R[3][state  ]
    end
    # return the offset of the current tile and the tile
    (s.x, state, Tile(((a, b, c),
                       (d, e, f),
                       (g, h, i)))), state+1
end   

# ~~~ THE STREAMING RECURRENCE MATRIX ~~~
struct StreamDistMatrix{S<:StreamView, F, D}
          sview::S
           dist::F
    startoffset::Int
              R::Vector{Vector{D}} # three vectors of eltype D
end

function streamdistmat(sview::StreamView, dist, startoffset::Int)
    # get one element of the stream, compute distance and 
    # initialise vectors R with appropriate type
    x = sview.buffer[1]
    res = dist(x, x)
    R = Vector{typeof(res)}[[zero(res) for i = 1:width(sview)] for j = 1:3]
    StreamDistMatrix(sview, dist, startoffset, R)
end

# ~~~ Iteration over the slices ~~~
@inline function Base.start(srm::StreamDistMatrix)
    _, state = _advance(srm, start(srm.sview))
    _, state = _advance(srm, state)
    state
end

@inline Base.next(srm::StreamDistMatrix, state) = _advance(srm, state)

@inline function _advance(srm::StreamDistMatrix, state)
    # get new snapshots
    window, state = next(srm.sview, state)
    # compute distances and shift
    insert!(srm.R, 
            length(srm.R), 
            _kernel(shift!(srm.R), srm.startoffset, window, srm.dist))
    (srm.R, window[2]), state
end

@inline function _kernel(r, startoffset, window, dist)
    for i in startoffset:length(r)
        r[i] = dist(window[1], window[i])
    end
    r
end

@inline Base.done(srm::StreamDistMatrix, state) = done(srm.sview, state)

# ~~~ Iteration over all tiles ~~~
struct Tiles{D<:StreamDistMatrix}
    d::D
end

@inline Base.start(t::Tiles) = start(t.d)
@inline function Base.next(t::Tiles, state) 
    (R, x), state = next(t.d, state)
    TileIterator(x, R, t.d.startoffset), state
end
@inline Base.done(t::Tiles, state) = done(t.d, state)

tiles(d::StreamDistMatrix) = Iterators.flatten(Tiles(d))

# ~~~ Iterator over Recurrences  ~~~
function recurrences(R::StreamDistMatrix, predicate::Function=(x->true))
    # filter based on minima and custom predicate
    g(data) = isminimum(data[2]) && predicate(centre(data[3]))
    f = Iterators.filter(g, tiles(R))
    # unpack and return: current state, offset and cargo
    ((data[1], data[2]) for data in f)
end