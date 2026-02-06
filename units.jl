const _UNITS_SYNONYMS = [
    ["mm", "kg m-2"],
    ["mm h-1", "kg m-2 h-1"],
    ["mm s-1", "kg m-2 s-1"],
    ["degC", "°C", "℃", "degree_Celsius"],
]

# Conversions accept two parameters: scalar value and timestep width (seconds).
const _UNITS_CONVERSIONS = Dict{
    Tuple{String, String},
    Function
}(
    ("degC", "K") => (x, _) -> x + 273.15,
    ("K", "degC") => (x, _) -> x - 273.15,
    ("Pa", "kPa") => (x, _) -> x / 1000,
    ("kPa", "Pa") => (x, _) -> x * 1000,
    ("mm h-1", "mm") => (x, t) -> x * t / 3600,
    ("mm s-1", "mm") => (x, t) -> x * t,
)

const _CANONICAL = Dict{String, String}(
    u => group[1] for group in _UNITS_SYNONYMS for u in group
)

canonical(u::String)::String = get(_CANONICAL, u, u)

#
# Convert units of an array of data.
#
# @param data Array of data to convert.
# @param from Original units of the data.
# @param to Desired units of the data.
# @param timestep Timestep width (seconds) for time-varying units.
function convert_units(data::AbstractArray{<:AbstractFloat},
                       from::String, to::String,
                       timestep::Int)::AbstractArray{<:AbstractFloat}
    f = canonical(from)
    t = canonical(to)
    if f == t
        return data
    end

    if haskey(_UNITS_CONVERSIONS, (f, t))
        @debug "Converting $(length(data)) values from $from to $to"
        return _UNITS_CONVERSIONS[(f, t)].(data, timestep)
    end
    error("No conversion function found for units $from -> $to")
end
