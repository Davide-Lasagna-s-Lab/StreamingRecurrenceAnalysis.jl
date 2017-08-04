export snapshot_stream_view, _step!, SnapshotStreamView

# This cannot be a subtype of AbstractVector, because its elements are
# inherently sequential.
struct SnapshotStreamView{X, G}
         g::G         # the forward map
    buffer::Vector{X} # buffer with snapshots
         N::Int       # number of snapshots produced
end

"""
    snapshot_stream_view(g, x₀::X, width::Int, N::Int)

Construct a view over a stream of snapshots, starting from an initial condition
`x₀` and generated by the forward map operator `g`. The map `g` is supposed to
work in-place. The view is a vector with elements of type `X` and length `width`.
The snapshot stream terminates after `N` elements have been generated.
"""
function snapshot_stream_view(g, x₀::X, width::Int, N::Int) where {X}
    N ≥ width - 1 ||
        throw(ArgumentError("length of snapshot stream must be higher than buffer width"))
    width > 0 ||
        throw(ArgumentError("width must be positive, got $width"))    
    buffer = X[similar(x₀) for i = 1:width]
    buffer[1] .= x₀
    SnapshotStreamView(g, buffer, N)
end


# ~~~ Iteration Protocol ~~~
function Base.start(s::SnapshotStreamView)
    # fill all buffer except one element
    for i = 1:length(s.buffer) - 2
        _step!(s)
    end
    return length(s.buffer) - 2
end

# handle the width = 1 case
Base.next(s::SnapshotStreamView, state) = (state < 0 ? s.buffer : _step!(s), state += 1)
Base.done(s::SnapshotStreamView, state) = state == s.N

# advance time and return buffer
function _step!(s::SnapshotStreamView)
    # copy current state `buffer[1]` to storage that will be overwritten
    s.buffer[end] .= s.buffer[1]

    # remove last, advance in time, then insert at the beginning
    insert!(s.buffer, 1, s.g(pop!(s.buffer)))
end

# number of views returned
Base.length(s::SnapshotStreamView) = s.N - 1