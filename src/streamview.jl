export StreamView, streamview, step!

struct StreamView{X,G} <: AbstractVector{X}
    g::G        
    buf::Vector{X}
end

function streamview(g, x₀::X, width::Int) where {X}
    width > 0 ||
        throw(ArgumentError("width must be positive, got $(length(width))"))
    # allocate memory and initialise
    buf = push!(sizehint!(X[], width), x₀)
    for i = 1:width-1
        push!(buf, g(deepcopy(buf[end])))
    end
    StreamView(g, buf)
end

# ~~~ Array interface ~~~
@inline Base.size(s::StreamView) = (length(s.buf), )
@inline Base.getindex(s::StreamView, i::Int) = s.buf[i]
@inline Base.eltype(::Type{StreamView{X}}) where {X} = X

@inline function step!(s::StreamView{X}) where {X}
    @inbounds s.buf[1] .= s.buf[end]
    push!(s.buf, s.g(shift!(s.buf)))
    s
end