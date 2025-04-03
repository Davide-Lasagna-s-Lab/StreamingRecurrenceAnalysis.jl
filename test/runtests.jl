using Test

using  StreamingRecurrenceAnalysis
import StreamingRecurrenceAnalysis: _isrecurrence, step!

struct LogisticMap
    r::Float64
end
@inline (k::LogisticMap)(x) = round(x*k.r*(1-x),; digits=5)

include("test_streamview.jl")
include("test_recurrence.jl")