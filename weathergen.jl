#!/usr/bin/env julia

using ArgParse
using Logging
using NCDatasets
using Dates

################################################################################
# Constants
################################################################################

# Standard name of latitude axes (as per CF spec).
const STD_LAT = "latitude"

# Standard name of longitude axes (as per CF spec).
const STD_LON = "longitude"

# Standard name of time axes (as per CF spec).
const STD_TIME = "time"

struct Options
    seed::Int
    in_tmin::String
    in_tmax::String
    in_rs::String
    in_pr::String
    in_ps::String
    in_ws::String
    out_temp::String
    out_pr::String
    out_ps::String
    out_rs::String
    out_ws::String
    out_vpd::String
    name_tmax::String
    name_tmin::String
    out_name_temp::String
    out_name_vpd::String
    show_progress::Bool
    log_level::Logging.LogLevel
end

struct DimensionIndices
    index_lon::Int
    index_lat::Int
    index_time::Int
end

function parse_log_level(level::Int)::Logging.LogLevel
    if level == 0
        return Logging.Error
    elseif level == 1
        return Logging.Warn
    elseif level == 2
        return Logging.Info
    elseif level == 3
        return Logging.Debug
    else
        error("Invalid log level: $level")
    end
end

function parse_cli()::Options
    parser = ArgParseSettings(description="Generate hourly meteorology from daily inputs")
    @add_arg_table! parser begin
        "--seed", "-s"
            arg_type=Int
            default=0
            help="Seed for deterministic PRNG"
        "--file-tmin"
            arg_type=String
            required=true
            help="Input file with daily minimum temperature (℃)"
        "--file-tmax"
            arg_type=String
            required=true
            help="Input file with daily maximum temperature (℃)"
        "--file-rs"
            arg_type=String
            required=true
            help="Input file with daily shortwave radiation (W m-2)"
        "--file-pr"
            arg_type=String
            required=true
            help="Input file with daily precipitation (mm)"
        "--file-ps"
            arg_type=String
            required=true
            help="Input file with daily air pressure (Pa)"
        "--file-ws"
            arg_type=String
            required=true
            help="Input file with daily wind speed (m s-1)"
        "--out-temp"
            arg_type=String
            required=true
            help="Path to hourly air temperature output file (℃)"
        "--out-pr"
            arg_type=String
            required=true
            help="Path to hourly precipitation output file (mm)"
        "--out-ps"
            arg_type=String
            required=true
            help="Path to hourly air pressure output file (Pa)"
        "--out-rs"
            arg_type=String
            required=true
            help="Path to hourly shortwave radiation output file (W m-2)"
        "--out-ws"
            arg_type=String
            required=true
            help="Path to hourly wind speed output file (m s-1)"
        "--out-vpd"
            arg_type=String
            required=true
            help="Path to hourly vapour pressure deficit output file (kPa)"
        "--name-tmax"
            arg_type=String
            default="tasmax"
            help="Name of daily maximum temperature variable"
        "--name-tmin"
            arg_type=String
            default="tasmin"
            help="Name of daily minimum temperature variable"
        "--out-name-temp"
            arg_type=String
            default="tas"
            help="Name of output hourly air temperature variable"
        "--out-name-vpd"
            arg_type=String
            default="vpd"
            help="Name of output hourly vapour pressure deficit variable"
        "--show-progress"
            action = :store_true
            help="Show progress bar"
        "--verbosity", "-v"
            arg_type=Int
            default=2
            help="Verbosity level (0: errors, 1: warnings, 2: info, 3: debug)"
    end
    parsed = parse_args(parser)
    log_level = parse_log_level(parsed["verbosity"])
    return Options(parsed["seed"], parsed["file-tmin"], parsed["file-tmax"],
                   parsed["file-rs"], parsed["file-pr"], parsed["file-ps"],
                   parsed["file-ws"], parsed["out-temp"], parsed["out-pr"],
                   parsed["out-ps"], parsed["out-rs"], parsed["out-ws"],
                   parsed["out-vpd"], parsed["name-tmax"], parsed["name-tmin"],
                   parsed["out-name-temp"], parsed["out-name-vpd"],
                   parsed["show-progress"], log_level)
end

################################################################################
# Native wrappers
################################################################################
const libwg = joinpath(@__DIR__, "libweathergen.so")

# Seed wrapper
function wg_seed(seed::Int)::Cint
    return ccall((:wg_seed, libwg), Cint, (Int64,), Int64(seed))
end

# One-day generation wrapper (in-place outputs)
function wg_generate_day!(;
    idate::Int,
    alat::Float32, dayl::Float32, dec::Float32,
    deltat::NTuple{12,Float32},   # or use Vector{Float32} length 12
    tmin::Float32, tmax::Float32, sw_mean_wm2::Float32, precip_mm::Float32,
    wind_ms::Float32, press_pa::Float32, ca_umol_mol::Float32,
    tair::Vector{Float32}, tsoil::Vector{Float32}, rh::Vector{Float32},
    vpd::Vector{Float32}, vmfd::Vector{Float32},
    radabv::Vector{Float32}, fbeam::Vector{Float32},
    ppt::Vector{Float32}, winda::Vector{Float32},
    press::Vector{Float32}, ca::Vector{Float32},
)::Cint
    nhrs = Cint(length(tair))
    @assert length(tsoil) == nhrs
    @assert length(rh)    == nhrs
    @assert length(vpd)   == nhrs
    @assert length(vmfd)  == nhrs
    @assert length(ppt)   == nhrs
    @assert length(winda) == nhrs
    @assert length(press) == nhrs
    @assert length(ca)    == nhrs
    @assert length(radabv) == nhrs*3
    @assert length(fbeam)  == nhrs*3

    return ccall((:wg_generate_day, libwg), Cint,
        (Cint, Cfloat, Cfloat, Cfloat, Ptr{Cfloat}, Cfloat, Cfloat, Cfloat, Cfloat,
         Cfloat, Cfloat, Cfloat, Cint,
         Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat},
         Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}),
        Cint(idate), Cfloat(alat), Cfloat(dayl), Cfloat(dec),
        Base.unsafe_convert(Ptr{Cfloat}, pointer_from_objref(Ref{NTuple{12,Cfloat}}(Cfloat.(deltat)))),
        Cfloat(tmin), Cfloat(tmax), Cfloat(sw_mean_wm2), Cfloat(precip_mm),
        Cfloat(wind_ms), Cfloat(press_pa), Cfloat(ca_umol_mol), nhrs,
        pointer(tair), pointer(tsoil), pointer(rh), pointer(vpd), pointer(vmfd),
        pointer(radabv), pointer(fbeam), pointer(ppt), pointer(winda), pointer(press), pointer(ca))
end

function unique_by_path(ds::AbstractVector{<:NCDataset})::Vector{<:NCDataset}
    seen = Set{String}()
    out = NCDataset[]
    for d in ds
        p = NCDatasets.path(d)
        if !(p in seen)
            push!(seen, p)
            push!(out, d)
        end
    end
    return out
end

function var_from_std_name(nc::NCDataset, std_name::String)
    matching = varbyattrib(nc, standard_name = std_name)
    if length(matching) == 1
        return matching[1]
    else
        error("No variable found with standard name $std_name")
    end
end

function validate_time_axes(datasets::AbstractVector{<:NCDataset})
    # Group datasets by path/URI.
    datasets = unique_by_path(datasets)

    # Get time variable from the first file.
    time = var_from_std_name(datasets[1], STD_TIME)

    # Read timestamps from the first file.
    t0 = time[:]

    # Compare to timestamps in other files.
    for i in 2:length(datasets)
        t1 = var_from_std_name(datasets[i], STD_TIME)[:]
        if t0 != t1
            error("Time axes are not equal")
        end
    end

    @debug "Time axes are consistent between input files"

    # Ensure daily spacing.
    # TODO: allow for the absence of Feb 29.
    d = diff(t0)
    MS_PER_DAY = Dates.Millisecond(86400 * 1000)
    if any(d .!= MS_PER_DAY)
        # Iterate to discover error, so we can produce a targeted message.
        for i in 1:length(d)
            if d[i] != MS_PER_DAY
                @warn "Time axis is not daily at index $i: $(t0[i]) to $(t0[i+1]); expected $MS_PER_DAY but got $(d[i])"
                break
            end
        end
        error("Time axis is not daily")
    end

    @debug "Time axes are daily with no gaps"
end

function validate_spatial_axes(datasets::AbstractVector{<:NCDataset}, axis::String)
    # Group datasets by path/URI.
    datasets = unique_by_path(datasets)

    if length(datasets) < 2
        return
    end

    ax0 = var_from_std_name(datasets[1], axis)[:]

    for i in 2:length(datasets)
        ax1 = var_from_std_name(datasets[i], axis)[:]
        if ax0 != ax1
            error("$axis axes are not equal in files $(NCDatasets.path(datasets[1])) and $(NCDatasets.path(datasets[i]))")
        end
    end
end

function validate_variable(nc::NCDataset, var::NCDatasets.CFVariable, units::String)::DimensionIndices
    if get(var, "units", "") != units
        error("Variable $(var.name) has units $(get(var, "units", "")) but expected $units")
    end

    # Ensure that the variable is 3-dimensional.
    if ndims(var) != 3
        error("Variable $(var.name) has $(ndims(var)) dimensions but expected 3")
    end

    # Ensure that the three dimensions are latitude, longitude, and time.
    var_lon = var_from_std_name(nc, STD_LON)
    var_lat = var_from_std_name(nc, STD_LAT)
    var_time = var_from_std_name(nc, STD_TIME)

    if ndims(var_lon) != 1
        error("Longitude variable has $(ndims(var_lon)) dimensions but expected 1")
    end

    if ndims(var_lat) != 1
        error("Latitude variable has $(ndims(var_lat)) dimensions but expected 1")
    end

    if ndims(var_time) != 1
        error("Time variable has $(ndims(var_time)) dimensions but expected 1")
    end

    dim_lon = dimnames(var_lon)[1]
    dim_lat = dimnames(var_lat)[1]
    dim_time = dimnames(var_time)[1]

    index_lon = findfirst(==(dim_lon), dimnames(var))
    index_lat = findfirst(==(dim_lat), dimnames(var))
    index_time = findfirst(==(dim_time), dimnames(var))

    if index_lon < 0
        error("Variable $(var.name) does not have dimension $dim_lon")
    end
    if index_lat < 0
        error("Variable $(var.name) does not have dimension $dim_lat")
    end
    if index_time < 0
        error("Variable $(var.name) does not have dimension $dim_time")
    end

    return DimensionIndices(index_lon, index_lat, index_time)
end

function validate_variable_from_name(nc::NCDataset, name::String, units::String)::DimensionIndices
    if !haskey(nc, name)
        error("Variable $name not found in dataset")
    end
    var = nc[name]
    return validate_variable(nc, var, units)
end

function validate_variable_from_std_name(nc::NCDataset, std_name::String, units::String)::DimensionIndices
    var = var_from_std_name(nc, std_name)
    return validate_variable(nc, var, units)
end

function process_data(opts::Options, tmin::NCDataset, tmax::NCDataset, rs::NCDataset, pr::NCDataset, ps::NCDataset, ws::NCDataset)
    # Validate dimensions.
    validate_time_axes([tmin, tmax, rs, pr, ps, ws])
    validate_spatial_axes([tmin, tmax, rs, pr, ps, ws], STD_LON)
    validate_spatial_axes([tmin, tmax, rs, pr, ps, ws], STD_LAT)

    # Validate variables and get dimension indices.
    idx_tmin = validate_variable_from_name(tmin, opts.name_tmin, "degC")
    idx_tmax = validate_variable_from_name(tmax, opts.name_tmax, "degC")
    idx_rs = validate_variable_from_std_name(rs, "surface_downwelling_shortwave_flux_in_air", "W m-2")
    idx_pr = validate_variable_from_std_name(pr, "precipitation_amount", "mm")
    idx_ps = validate_variable_from_std_name(ps, "air_pressure", "Pa")
    idx_ws = validate_variable_from_std_name(ws, "wind_speed", "m s-1")

    # Initialise PRNG seed.
    wg_seed(opts.seed)

    # Iterate through gridcells. (All input files use the same grid.)
    var_lon = var_from_std_name(tmin, STD_LON)
    var_lat = var_from_std_name(tmin, STD_LAT)
    var_time = get_time_variable(tmin)
end

function main(opts::Options)
    @info "Running weather generator..."

    # Open input files for reading.
    NCDataset(opts.in_tmin) do tmin
        NCDataset(opts.in_tmax) do tmax
            NCDataset(opts.in_rs) do rs
                NCDataset(opts.in_pr) do pr
                    NCDataset(opts.in_ps) do ps
                        NCDataset(opts.in_ws) do ws
                            process_data(opts, tmin, tmax, rs, pr, ps, ws)
                        end
                    end
                end
            end
        end
    end
end

opts = parse_cli()

logger = ConsoleLogger(stdout, opts.log_level)
global_logger(logger)

main(opts)
