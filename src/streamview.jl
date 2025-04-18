export StreamView, streamview

struct StreamView{X,G} <: AbstractVector{X}
    g::G        
    buf::Vector{X}
end


"""
    sview = streamview(g, x₀::X, width::Int)

Construct a view over a stream of snapshots, starting from an initial condition
`x₀` and generated by the forward map operator `g : X → X` . The map `g` works
in-place. The view is a vector with elements of type `X` and length `width`.
The snapshot stream is updated using `step!(sview)`.

Examples
========

julia> g(x₀) = (x₀ .+= 1.0; x₀)

julia> sview = streamview(g, [1.0], 2)
2-element StreamingRecurrenceAnalysis.StreamView{Array{Float64,1},#g}:
 [1.0]
 [2.0]

julia> step!(sview)
2-element StreamingRecurrenceAnalysis.StreamView{Array{Float64,1},#g}:
 [2.0]
 [3.0]

"""
function streamview(g, x₀::X, width::Int) where {X}
    width > 0 || throw(ArgumentError("width must be positive, got $(length(width))"))
    # allocate memory and initialise
    buf = push!(sizehint!(X[], width), x₀)
    for i = 1:width-1
        push!(buf, g(copy(buf[end])))
    end
    return StreamView(g, buf)
end

# ~~~ Array interface ~~~
@inline Base.size(s::StreamView) = (length(s.buf), )
@inline Base.getindex(s::StreamView, i::Int) = s.buf[i]
@inline Base.eltype(::Type{StreamView{X}}) where {X} = X

@inline function step!(s::StreamView{X}) where {X}
    @inbounds s.buf[1] .= s.buf[end]
    push!(s.buf, s.g(popfirst!(s.buf)))
    return s
end

# hack to work with numbers or tuple, add if needed
_Immutable = Union{Number, Tuple{Number, Vararg{Number}}}
@inline function step!(s::StreamView{<:_Immutable})
    popfirst!(push!(s.buf, s.g(s.buf[end])))
    return s
end