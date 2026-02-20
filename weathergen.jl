#!/usr/bin/env -S julia --project=@.
#
# Usage:
# ./weathergen.jl -i infile.nc -o outfile.nc
#
# Parallel usage:
# mpiexecjl -n 4 ./weathergen.jl -i infile.nc -o outfile.nc
#
# Note: mpiexecjl must first be installed (typically to ~/.julia/bin):
# using Pkg
# Pkg.add("MPIPreferences")
# using MPIPreferences
# MPIPreferences.install_mpiexecjl()
#

module Weathergen

using ArgParse
using Logging
using NCDatasets
using Dates
using DataStructures
using MPI

################################################################################
# Includes
################################################################################
include(joinpath(@__DIR__, "units.jl"))

################################################################################
# Constants
################################################################################

# Standard name of latitude axes (as per CF spec).
const STD_LAT = "latitude"

# Standard name of longitude axes (as per CF spec).
const STD_LON = "longitude"

# Standard name of time axes (as per CF spec).
const STD_TIME = "time"

# Standard name attribute name (as per CF spec).
const ATTR_STD_NAME = "standard_name"

# Units attribute name (as per CF spec).
const ATTR_UNITS = "units"

# Long name attribute name (as per CF spec).
const ATTR_LONG_NAME = "long_name"

# Standard name of surface downwelling shortwave flux in air (as per CF spec).
const STD_RS = "surface_downwelling_shortwave_flux_in_air"

# Standard name of precipitation amount (as per CF spec).
const STD_PR = "precipitation_amount"

# Standard name of air pressure (as per CF spec).
const STD_PS = "air_pressure"

# Number of timesteps per day.
const DAY_LENGTH = 24

# Name of the air pressure variable created in the output file when using fixed
# air pressure.
const NAME_PS = "ps"

# Units of the air pressure variable in the output file.
const UNITS_PS = "Pa"

# Long name of the air pressure variable in the output file.
const LONG_PS = "Air pressure"

# Maximum allowed size of a single chunk in bytes: 4GiB.
const MAX_CHUNK_SIZE = 4 * 1024^3

# Standard name of air temperature (as per CF spec).
const STD_TEMP = "air_temperature"

# Long name of air temperature (as per CF spec).
const LONG_TEMP = "Air temperature"

# Units of the air temperature variable in the output file.
const UNITS_TEMP = "degC"

# Units of the shortwave radiation variable in the output file.
const UNITS_RS = "W m-2"

# Units of the VPD variable in the output file.
const UNITS_VPD = "kPa"

# Units of the precipitation variable in the output file.
const UNITS_PR = "mm"

# Default constant air pressure if not using dynamic air pressure (Pa).
const DEFAULT_PS = 101300

# Standard name of VPD (as per CF spec).
const STD_VPD = "vapour_pressure_deficit"

# Long name of VPD (as per CF spec).
const LONG_VPD = "Vapour pressure deficit"

# Required units of air temperature.
const REQ_UNITS_TEMP = "degC"

# Required units of radiation.
const REQ_UNITS_RS = "W m-2"

# Required units of precipitation.
const REQ_UNITS_PR = "mm"

# Required units of air pressure.
const REQ_UNITS_PS = "Pa"

# Constants for splitmix64.
const C1 = 0x9e3779b97f4a7c15
const C2 = 0xbf58476d1ce4e5b9
const C3 = 0x94d049bb133111eb
const WG_SEED_DEFAULT = 88172645463393265

################################################################################
# Types
################################################################################

# Parsed CLI options.
struct Options
    seed::Int
    in_tmin::String
    in_tmax::String
    in_rs::String
    in_pr::String
    in_ps::String
    out_temp::String
    out_pr::String
    out_ps::String
    out_rs::String
    out_vpd::String
    name_tmax::String
    name_tmin::String
    out_name_temp::String
    out_name_vpd::String
    log_level::Logging.LogLevel
    compression_level::Int
    chunk_size_lon::Int
    chunk_size_lat::Int
    chunk_size_time::Int

    # Constant air pressure value (Pa) to be used if input data doesn't include
    # air pressure.
    default_ps::Float32
    parallel::Bool
    show_progress::Bool
end

# Struct to hold a variable along with the indices of its dimensions.
struct DimensionIndices
    var::NCDatasets.CFVariable
    index_lon::Int
    index_lat::Int
    index_time::Int
end

# Struct to hold input and output file paths for a given variable, along with
# the corresponding CLI argument names for error reporting.
struct PerVariablePaths
    input_file::String
    output_file::String
    input_arg_name::String
    output_arg_name::String
end

# Dimension orders for input variables.
struct InputDimensionOrders
    idx_tmin::DimensionIndices
    idx_tmax::DimensionIndices
    idx_rs::DimensionIndices
    idx_pr::DimensionIndices
    idx_ps::Union{DimensionIndices, Nothing}
end

# Dimension orders for output variables.
struct OutputDimensionOrders
    idx_rs::DimensionIndices
    idx_pr::DimensionIndices
    idx_temp::DimensionIndices
    idx_vpd::DimensionIndices
    idx_ps::DimensionIndices
end

# Chunk sizes for each dimension.
struct ChunkSizes
    lon::Int
    lat::Int
    time::Int
end

################################################################################
# CLI parsing.
################################################################################

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
        "-c", "--compression-level"
            arg_type=Int
            default=5
            help="Compression level for output NetCDF file (0-9; 0 = no compression, 9 = maximum compression)."
        "--chunk-lon"
            arg_type=Int
            default=1
            help="Chunk size to use on the longitude dimension."
        "--chunk-lat"
            arg_type=Int
            default=1
            help="Chunk size to use on the latitude dimension."
        "--chunk-time"
            arg_type=Int
            default=8760
            help="Chunk size to use on the time dimension."
        "--default-ps"
            arg_type=Float32
            default=DEFAULT_PS
            help="Default air pressure (Pa) to use if no input file is provided."
        "--verbosity", "-v"
            arg_type=Int
            default=2
            help="Verbosity level (0: errors, 1: warnings, 2: info, 3: debug)"
        "-P", "--parallel"
            action=:store_true
            help="Use MPI for parallel processing. All MPI-specific behaviour is hidden behind this flag."
        "-i", "--input-file"
            arg_type=String
            help="Input file with daily meteorology. Use this if all variables are in a single file."
        "-o", "--output-file"
            arg_type=String
            help="Desired path to the hourly output file. Use this if all variables are in a single file."
        "--file-tmin"
            arg_type=String
            help="Input file with daily minimum temperature (℃). Mutually exclusive with --input-file."
        "--file-tmax"
            arg_type=String
            help="Input file with daily maximum temperature (℃). Mutually exclusive with --input-file."
        "--file-rs"
            arg_type=String
            help="Input file with daily shortwave radiation (W m-2). Mutually exclusive with --input-file."
        "--file-pr"
            arg_type=String
            help="Input file with daily precipitation (mm). Mutually exclusive with --input-file."
        "--file-ps"
            arg_type=String
            help="Input file with daily air pressure (Pa). Mutually exclusive with --input-file."
        "--out-temp"
            arg_type=String
            help="Path to hourly air temperature output file (℃). Mutually exclusive with --output-file."
        "--out-rs"
            arg_type=String
            help="Path to hourly shortwave radiation output file (W m-2). Mutually exclusive with --output-file."
        "--out-pr"
            arg_type=String
            help="Path to hourly precipitation output file (mm). Mutually exclusive with --output-file."
        "--out-ps"
            arg_type=String
            help="Path to hourly air pressure output file (Pa). Mutually exclusive with --output-file."
        "--out-vpd"
            arg_type=String
            help="Path to hourly vapour pressure deficit output file (kPa). Mutually exclusive with --output-file."
        "--show-progress"
            action=:store_true
            help="Show overall progress during processing."
    end

    parsed = parse_args(parser)
    log_level = parse_log_level(parsed["verbosity"])

    # Validate file paths.
    file_tmin = parsed["file-tmin"]
    file_tmax = parsed["file-tmax"]
    file_rs = parsed["file-rs"]
    file_pr = parsed["file-pr"]
    file_ps = parsed["file-ps"]
    out_temp = parsed["out-temp"]
    out_pr = parsed["out-pr"]
    out_ps = parsed["out-ps"]
    out_rs = parsed["out-rs"]
    out_vpd = parsed["out-vpd"]

    if parsed["input-file"] !== nothing || parsed["output-file"] !== nothing
        if parsed["input-file"] === nothing || parsed["output-file"] === nothing
            error("Must specify both --input-file and --output-file if using a single file.")
        end
        file_tmin = parsed["input-file"]
        file_tmax = parsed["input-file"]
        file_rs = parsed["input-file"]
        file_pr = parsed["input-file"]
        file_ps = parsed["input-file"]
        out_temp = parsed["output-file"]
        out_pr = parsed["output-file"]
        out_ps = parsed["output-file"]
        out_rs = parsed["output-file"]
        out_vpd = parsed["output-file"]

        validate_file_path(parsed["input-file"], "--input-file")

        if parsed["input-file"] == parsed["output-file"]
            error("Input and output files must be different.")
        end

        ensure_not_set(parsed["file-tmin"], "--file-tmin")
        ensure_not_set(parsed["file-tmax"], "--file-tmax")
        ensure_not_set(parsed["file-rs"], "--file-rs")
        ensure_not_set(parsed["file-pr"], "--file-pr")
        ensure_not_set(parsed["file-ps"], "--file-ps")

        ensure_not_set(parsed["out-temp"], "--out-temp")
        ensure_not_set(parsed["out-pr"], "--out-pr")
        ensure_not_set(parsed["out-ps"], "--out-ps")
        ensure_not_set(parsed["out-rs"], "--out-rs")
        ensure_not_set(parsed["out-vpd"], "--out-vpd")
    else
        validate_file_path(file_tmin, "--file-tmin")
        validate_file_path(file_tmax, "--file-tmax")
        validate_file_path(file_rs, "--file-rs")
        validate_file_path(file_pr, "--file-pr")
        if file_ps !== nothing
            validate_file_path(file_ps, "--file-ps")
        end

        ensure_set(out_temp, "--out-temp")
        ensure_set(out_pr, "--out-pr")
        ensure_set(out_ps, "--out-ps")
        ensure_set(out_rs, "--out-rs")
        ensure_set(out_vpd, "--out-vpd")

        paths = [
            PerVariablePaths(file_tmin, out_temp, "--file-tmin", "--out-temp"),
            PerVariablePaths(file_tmax, out_temp, "--file-tmax", "--out-temp"),
            PerVariablePaths(file_rs, out_rs, "--file-rs", "--out-rs"),
            PerVariablePaths(file_pr, out_pr, "--file-pr", "--out-pr")
        ]
        if file_ps !== nothing
            push!(paths, PerVariablePaths(file_ps, out_ps, "--file-ps", "--out-ps"))
        else
            # FIXME: this is not ideal.
            # Dummy path to get past validation. This won't contain ps data,
            # so constant air pressure will be used instead.
            file_ps = file_tmin
        end
        validate_per_variable_paths(paths)
    end

    return Options(parsed["seed"], file_tmin, file_tmax, file_rs, file_pr,
                   file_ps, out_temp, out_pr, out_ps, out_rs, out_vpd,
                   parsed["name-tmax"], parsed["name-tmin"],
                   parsed["out-name-temp"], parsed["out-name-vpd"],
                   log_level, parsed["compression-level"],
                   parsed["chunk-lon"], parsed["chunk-lat"],
                   parsed["chunk-time"], parsed["default-ps"],
                   parsed["parallel"], parsed["show-progress"])
end

function validate_per_variable_paths(paths::Vector{PerVariablePaths})
    for variable in paths
        for var2 in paths
            if var2.input_arg_name == variable.input_arg_name ||
               var2.output_arg_name == variable.output_arg_name
                continue
            end

            if variable.input_file == var2.input_file
                error("$(variable.input_arg_name) and $(var2.input_arg_name) cannot have the same input file (currently: $(variable.input_file)).")
            end

            if variable.output_file == var2.output_file
                error("$(variable.output_arg_name) and $(var2.output_arg_name) cannot have the same output file (currently: $(variable.output_file)).")
            end
        end
    end
end

function ensure_not_set(path::Union{String, Nothing}, arg_name::String)
    if path !== nothing
        error("Invalid $arg_name value: cannot specify both $arg_name and --input-file/--output-file")
    end
end

function ensure_set(path::Union{String, Nothing}, arg_name::String)
    if path === nothing
        error("Missing required argument: $arg_name")
    end
end

function validate_file_path(path::Union{String, Nothing}, arg_name::String)
    ensure_set(path, arg_name)

    if !isfile(path)
        error("Invalid $arg_name value: file $(path) does not exist.")
    end
end

################################################################################
# Native wrappers
################################################################################
const libwg = joinpath(@__DIR__, "libweathergen.so")

function wg_init()::Nothing
    ccall((:wg_init, libwg), Cvoid, ())
    return nothing
end

# Seed wrapper
function wg_seed(seed::Int64)::Cint
    # Seed value of 0 causes RNG to depend on prior state, which breaks
    # reproducibility. In that case, use default seed value instead.
    if seed == 0
        @warn "Seed value of 0 is not allowed; using default seed value of $WG_SEED_DEFAULT instead."
        seed = WG_SEED_DEFAULT
    end
    return ccall((:wg_seed, libwg), Cint, (Int64,), seed)
end

function wg_seed(seed::UInt64)::Cint
    return wg_seed(reinterpret(Int64, seed))
end

"""Generate hourly meteorology for a single day.

This is a thin wrapper over the Fortran `wg_generate_day` C ABI.
All scalars are passed by value. Output buffers must be preallocated.
"""
function wg_generate_day(
    idate::Int,
    alat::Real,
    tmin::Real,
    tmax::Real,
    sw_mean_wm2::Real,
    precip_mm::Real,
    press_pa::Real,
    nhrs::Int,
    tair::Ptr{Cfloat},
    tsoil::Ptr{Cfloat},
    rh::Ptr{Cfloat},
    vpd::Ptr{Cfloat},
    vmfd::Ptr{Cfloat},
    radabv::Ptr{Cfloat},
    fbeam::Ptr{Cfloat},
    ppt::Ptr{Cfloat},
    press::Ptr{Cfloat},
)::Cint
    return ccall(
        (:wg_generate_day, libwg),
        Cint,
        (Cint, Cfloat, Cfloat, Cfloat, Cfloat, Cfloat, Cfloat, Cint,
            Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat},
            Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}),
        Cint(idate),
        Cfloat(alat),
        Cfloat(tmin),
        Cfloat(tmax),
        Cfloat(sw_mean_wm2),
        Cfloat(precip_mm),
        Cfloat(press_pa),
        Cint(nhrs),
        tair,
        tsoil,
        rh,
        vpd,
        vmfd,
        radabv,
        fbeam,
        ppt,
        press,
    )
end

################################################################################
# Main script.
################################################################################

function get_rank()::Int
    if @isdefined MPI
        return MPI.Comm_rank(MPI.COMM_WORLD)
    else
        return 0
    end
end

function get_world_size()::Int
    if @isdefined MPI
        return MPI.Comm_size(MPI.COMM_WORLD)
    else
        return 1
    end
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

function has_var_with_std_name(nc::NCDataset, std_name::String)::Bool
    return length(varbyattrib(nc, standard_name = std_name)) > 0
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

function validate_spatial_axes(datasets::AbstractVector{<:NCDataset},
                               axis::String)
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

function validate_variable(nc::NCDataset,
                           var::NCDatasets.CFVariable)::DimensionIndices
    if !haskey(var.attrib, ATTR_UNITS)
        error("Variable $(name(var)) has no units attribute")
    end

    # Ensure that the variable is 3-dimensional.
    if ndims(var) != 3
        error("Variable $(name(var)) has $(ndims(var)) dimensions but expected 3")
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
        error("Variable $(name(var)) does not have dimension $dim_lon")
    end
    if index_lat < 0
        error("Variable $(name(var)) does not have dimension $dim_lat")
    end
    if index_time < 0
        error("Variable $(name(var)) does not have dimension $dim_time")
    end

    return DimensionIndices(var, index_lon, index_lat, index_time)
end

function validate_variable_from_name(nc::NCDataset,
                                     name::String)::DimensionIndices
    if !haskey(nc, name)
        error("Variable $name not found in dataset")
    end
    var = nc[name]
    return validate_variable(nc, var)
end

function validate_variable_from_std_name(nc::NCDataset,
                                         std_name::String)::DimensionIndices
    var = var_from_std_name(nc, std_name)
    return validate_variable(nc, var)
end

function hyperslab(indices::DimensionIndices, i::Int, j::Int, k::Int)
    ndim = length(indices.var.var.dimids)
    @assert ndim == 3

    # C-style indexing: start at 0.
    start = zeros(Int, ndim)
    count = ones(Int, ndim)
    stride = ones(Int, ndim)

    # i-1, j-1 for C-style indexing.
    start[indices.index_lat] = i - 1
    start[indices.index_lon] = j - 1
    start[indices.index_time] = 0
    count[indices.index_time] = k

    # C-style ordering for the netcdf API.
    return reverse(start), reverse(count), reverse(stride)
end

function read_raw(indices::DimensionIndices, i::Int, j::Int, ntime::Int)
    v = indices.var.var
    std_name = get(v.attrib, ATTR_STD_NAME, "unknown")
    @debug "Reading $ntime values of variable $(v.varid) ($std_name) at gridcell ($i, $j)"

    start, count, stride = hyperslab(indices, i, j, ntime)
    data = Vector{eltype(v)}(undef, ntime)

    NCDatasets.nc_get_vars!(v.ds.ncid, v.varid, start, count, stride, data)
    @debug "Successfully read $(length(data)) values for variable $(v.varid) ($std_name) at gridcell ($i, $j)"
    return data
end

function read(indices::DimensionIndices, i::Int, j::Int, units::String, dt::Int, ntime::Int)
    data = read_raw(indices, i, j, ntime)

    # Error if any data is missing.
    if any(ismissing, data)
        error("Missing data in variable $(name(indices.var)) at gridcell ($i, $j)")
    end

    # Convert from Vector{Union{Float32, Missing}} to Vector{Float32}.
    data = convert(Vector{Float32}, data)

    return convert_units(data, indices.var.attrib[ATTR_UNITS], units, dt)
end

"""Convert a date-like value to Maespa's `idate` (days since 1950-01-01).

This matches the Fortran calendar math used by Maespa (Julian-style leap
years: every 4 years, with no century exception), so `JDATE(idate)` returns
the correct day-of-year.
"""
function date_to_idate(date)
    # Note: can't use Dates.value() because the underlying fortran code uses a
    # Julian calendar, with leap days exactly every 4 years.

    yearValue = Dates.year(date)
    monthValue = Dates.month(date)
    dayValue = Dates.day(date)

    ifd = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

    dayOfYear = ifd[monthValue] + dayValue
    isLeapYear = (4 * (yearValue / 4) == yearValue)
    if isLeapYear && monthValue >= 3
        dayOfYear += 1
    end

    yearsSince1950 = yearValue - 1950
    daysBeforeYear = 365 * yearsSince1950 + div(yearsSince1950 - 1, 4)

    return daysBeforeYear + dayOfYear - 1
end

function write(indices::DimensionIndices, data::Vector{Float32}, i::Int, j::Int)
    v = indices.var.var
    start, count, stride = hyperslab(indices, i, j, length(data))

    std_name = get(v.attrib, ATTR_STD_NAME, "unknown")
    @debug "Writing $(length(data)) values of variable $(v.varid) ($std_name) at gridcell ($i, $j)"

    NCDatasets.nc_put_vars(v.ds.ncid, v.varid, start, count, stride, data)
end

function barrier(message::AbstractString)
    @debug "Entering MPI barrier: $message"
    MPI.Barrier(MPI.COMM_WORLD)
    @debug "Exiting MPI barrier: $message"
end

const NcAttribContainer = Union{NCDataset, NCDatasets.CFVariable}
function copy_attributes(in::NcAttribContainer, out::NcAttribContainer)
    for attr in keys(in.attrib)
        # Attributes starting with _ are for netcdf internal use only.
        if attr[1] != '_'
            out.attrib[attr] = in.attrib[attr]
        end
    end
end

# Create a variable in the output file using the specified metadata.
function create_var(path::String, name::String, units::String,
                    std_name::String, long_name::String,
                    dims::Vector{String}, compression_level::Int,
                    chunk_sizes::Vector{Int})

    # No need to use the open_netcdf() MPI-mode wrapper here, because variable
    # creation only occurs on the master node, so we can just use independent
    # access mode (the default).
    NCDataset(path, "a") do nc_out
        var = defVar(nc_out, name, Float32, dims,
                     deflatelevel=compression_level,
                     shuffle=compression_level > 0,
                     chunksizes=chunk_sizes,
                     fillvalue=fillvalue(Float32))
        var.attrib[ATTR_UNITS] = units
        var.attrib[ATTR_STD_NAME] = std_name
        var.attrib[ATTR_LONG_NAME] = long_name
    end
end

# Convenience function for when the output variable name is the same as the
# input variable name.
function create_var_from_existing(path::String, idx_in::DimensionIndices,
                                  compression_level::Int, units::String,
                                  sizes::ChunkSizes)
    var = idx_in.var
    std_name = var.attrib[ATTR_STD_NAME]
    long_name = var.attrib[ATTR_LONG_NAME]
    dims = [dimnames(var)...]
    chunk_sizes = get_chunk_size(idx_in, sizes)

    create_var(path, name(var), units, std_name, long_name, dims,
               compression_level, chunk_sizes)
end

function create_outfile(nc_in::NCDataset, nc_out::NCDataset, opts::Options)
    compression = opts.compression_level
    shuffle = compression > 0

    # Copy all dimensions from the input file.
    for dim in keys(nc_in.dim)
        if dim != "time"
            size = nc_in.dim[dim]
            defDim(nc_out, dim, size)
        end
    end

    # Copy all global attributes from the input file.
    copy_attributes(nc_in, nc_out)

    # Create and populate coordinate variables in output file.
    in_lon = var_from_std_name(nc_in, STD_LON)
    in_lat = var_from_std_name(nc_in, STD_LAT)
    in_time = var_from_std_name(nc_in, STD_TIME)

    # Read coordinate values from input file.
    lons = in_lon[:]
    lats = in_lat[:]

    # Get the size of each longitude value in bytes.
    lon_size = sizeof(eltype(lons))

    # Calculate the maximum chunk size for the longitude dimension to
    # ensure it does not exceed MAX_CHUNK_SIZE.
    chunk_size_lon = min(Int(MAX_CHUNK_SIZE / lon_size), length(lons))
    @debug "Longitude axis contains $(length(lons)) values, each of size $lon_size bytes"
    @debug "Using longitude chunk size of $chunk_size_lon"

    # Create longitude variable in the output file.
    out_lon = defVar(nc_out, name(in_lon), lons, dimnames(in_lon),
                     deflatelevel=compression, shuffle=shuffle,
                     chunksizes=[chunk_size_lon],
                     fillvalue=fillvalue(typeof(lons[1])))
    copy_attributes(in_lon, out_lon)

    # Get the size of each latitude value in bytes.
    lat_size = sizeof(eltype(lats))

    # Calculate the maximum chunk size for the latitude dimension to
    # ensure it does not exceed MAX_CHUNK_SIZE.
    chunk_size_lat = min(Int(MAX_CHUNK_SIZE / lat_size), length(lats))
    @debug "Latitude axis contains $(length(lats)) values, each of size $lat_size bytes"
    @debug "Using latitude chunk size of $chunk_size_lat"

    # Create latitude variable in the output file.
    out_lat = defVar(nc_out, name(in_lat), lats, dimnames(in_lat),
                        deflatelevel=compression, shuffle=shuffle,
                        chunksizes=[chunk_size_lat],
                        fillvalue=fillvalue(typeof(lats[1])))
    copy_attributes(in_lat, out_lat)

    # Construct hourly timeseries from each day in the input time
    # variable.
    times = in_time[:]
    hours = [t + Hour(h) for t in times for h in 0:(DAY_LENGTH - 1)]

    # Not writing a fill value attribute for time, since there should be no
    # missing values.
    out_time = defVar(nc_out, name(in_time), hours, dimnames(in_time),
                        attrib = OrderedDict(
                            ATTR_UNITS => in_time.attrib[ATTR_UNITS],
                            "calendar" => in_time.attrib["calendar"],
                        ),
                        deflatelevel=compression, shuffle=shuffle,
                        chunksizes=[opts.chunk_size_time])
    copy_attributes(in_time, out_time)
end

function create_output_files(opts::Options, nc_in::NCDataset)
    paths = unique([opts.out_temp, opts.out_rs, opts.out_pr, opts.out_ps,
                    opts.out_vpd])
    for path in paths
        # Create directory if it doesn't already exist.
        dir = dirname(path)
        if dir != "" && !isdir(dir)
            mkdir(dir)
        end

        # Note: we never enable MPI parallel access when creating output files,
        # because this is only done by the master node.
        NCDataset(path, "c") do nc_out
            create_outfile(nc_in, nc_out, opts)
        end
    end
end

function get_chunk_size(order::DimensionIndices, sizes::ChunkSizes)::Vector{Int}
    chunks = [1, 1, 1]
    chunks[order.index_time] = sizes.time
    chunks[order.index_lat] = sizes.lat
    chunks[order.index_lon] = sizes.lon
    return chunks
end

# Deterministic 64-bit mix (SplitMix64-style constants).
function mix64(base_seed::UInt64, i::UInt64, j::UInt64)::UInt64
    u = base_seed
    u ⊻= i * C1
    u ⊻= j * C2

    u += C1
    u = (u ⊻ (u >> 30)) * C2
    u = (u ⊻ (u >> 27)) * C3
    u = u ⊻ (u >> 31)
    return u
end

function mix64(seed::Int, i::Int, j::Int)::UInt64
    return mix64(reinterpret(UInt64, Int64(seed)),
                 reinterpret(UInt64, Int64(i)),
                 reinterpret(UInt64, Int64(j)))
end

function mix64(seed::Int, lat::Float64, lon::Float64)::UInt64
    return mix64(reinterpret(UInt64, Int64(seed)),
                 reinterpret(UInt64, lat),
                 reinterpret(UInt64, lon))
end

function mix64(seed::Int, lat::Float32, lon::Float32)::UInt64
    return mix64(reinterpret(UInt64, Int64(seed)),
                 reinterpret(UInt64, Float64(lat)),
                 reinterpret(UInt64, Float64(lon)))
end

function get_partition(total::Int, world_size::Int, rank::Int)::UnitRange{Int}
    # Partition a 1D index space into contiguous blocks, distributed as evenly
    # as possible across ranks.
    base = total ÷ world_size
    extra = total % world_size

    count = base + (rank < extra ? 1 : 0)
    start = rank * base + min(rank, extra) + 1

    if count == 0
        return 1:0
    end

    return start:(start + count - 1)
end

function get_workload(opts::Options,
                      nlat::Int,
                      nlon::Int)::AbstractUnitRange{Int}
    ncells = nlat * nlon

    # If running in serial mode, return all gridcells.
    if !opts.parallel
        return 1:ncells
    end

    return get_partition(ncells, get_world_size(), get_rank())
end

function get_max_workload_size(opts::Options,
                               nlat::Int,
                               nlon::Int)
    ncells = nlat * nlon

    if !opts.parallel
        return ncells
    end

    return cld(ncells, get_world_size())
end

function generate_weather(opts::Options, indices_in::InputDimensionOrders,
                          indices_out::OutputDimensionOrders)
    # Iterate through gridcells. (All input files use the same grid.)
    var_lon = var_from_std_name(indices_in.idx_tmin.var.var.ds, STD_LON)
    var_lat = var_from_std_name(indices_in.idx_tmin.var.var.ds, STD_LAT)
    var_time = var_from_std_name(indices_in.idx_tmin.var.var.ds, STD_TIME)

    lons = var_lon[:]
    lats = var_lat[:]
    times = var_time[:]

    # Get timestep width in seconds.
    # TODO: more robust time delta handling.
    dt = Second(Dates.value(times[2] - times[1]) / 1000).value

    workitems = get_workload(opts, length(lats), length(lons))
    workload_size = length(workitems)
    if opts.parallel
        if workload_size == 0
            @info "Workload contains 0 gridcells"
        else
            nlon = length(lons)
            first_cell = first(workitems)
            last_cell = last(workitems)
            first_i = ((first_cell - 1) ÷ nlon) + 1
            first_j = ((first_cell - 1) % nlon) + 1
            last_i = ((last_cell - 1) ÷ nlon) + 1
            last_j = ((last_cell - 1) % nlon) + 1
            @info "Workload contains $workload_size gridcells (cell range $first_cell:$last_cell, first ($first_i,$first_j), last ($last_i,$last_j))"
        end
    else
        @info "Processing $workload_size gridcells"
    end

    # Record start time for progress reporting.
    start_time = time()

    # Iterate through gridcells. Generate climate one gridcell at a time.
    nlon = length(lons)
    ntime = length(times)
    for cell in workitems
        i = ((cell - 1) ÷ nlon) + 1
        j = ((cell - 1) % nlon) + 1

        lat = lats[i]
        lon = lons[j]

        # Initialise PRNG seed.
        wg_seed(mix64(opts.seed, lat, lon))

        @info "Processing gridcell $i, $j ($lat, $lon)"

        # Read timeseries for this gridcell.
        tmin_data = read(indices_in.idx_tmin, i, j, UNITS_TEMP, dt, ntime)
        tmax_data = read(indices_in.idx_tmax, i, j, UNITS_TEMP, dt, ntime)
        rs_data = read(indices_in.idx_rs, i, j, UNITS_RS, dt, ntime)
        pr_data = read(indices_in.idx_pr, i, j, UNITS_PR, dt, ntime)
        if indices_in.idx_ps === nothing
            ps_data = fill(opts.default_ps, length(times))
        else
            ps_data = read(indices_in.idx_ps, i, j, UNITS_PS, dt, ntime)
        end

        tair_out = Vector{Float32}(undef, DAY_LENGTH * length(times))
        vpd_out = Vector{Float32}(undef, DAY_LENGTH * length(times))
        rs_out = Vector{Float32}(undef, DAY_LENGTH * length(times))
        pr_out = Vector{Float32}(undef, DAY_LENGTH * length(times))
        ps_out = Vector{Float32}(undef, DAY_LENGTH * length(times))

        # Define daily arrays for the outputs we don't care about. These
        # will just be overwritten each day.
        tsoil_out = Vector{Float32}(undef, DAY_LENGTH)
        rh_out = Vector{Float32}(undef, DAY_LENGTH)
        vmfd_out = Vector{Float32}(undef, DAY_LENGTH)
        radabv_out = Vector{Float32}(undef, DAY_LENGTH * 3)
        fbeam_out = Vector{Float32}(undef, DAY_LENGTH * 3)

        # wg_generate_day() operates at the day level, so we need to iterate
        # over the days in the input file.
        for k in eachindex(times)
            # @debug "Generating climate for day: $(times[k])"

            # idate: days since 1950
            idate = date_to_idate(times[k])
            # alat: latitude in radians
            alat = deg2rad(lat)

            # Compute offset into output arrays.
            start = (k - 1) * DAY_LENGTH + 1

            # Get pointers to today's data in the long output arrays.
            GC.@preserve tair_out vpd_out rs_out pr_out ps_out tsoil_out rh_out vmfd_out radabv_out fbeam_out begin

                tair_day = pointer(tair_out, start)
                vpd_day = pointer(vpd_out, start)
                pr_day = pointer(pr_out, start)
                ps_day = pointer(ps_out, start)

                tsoil_day = pointer(tsoil_out, 1)
                rh_day = pointer(rh_out, 1)
                vmfd_day = pointer(vmfd_out, 1)
                radabv_day = pointer(radabv_out, 1)
                fbeam_day = pointer(fbeam_out, 1)

                # Call the weather generator.
                wg_generate_day(idate, alat, tmin_data[k], tmax_data[k],
                                rs_data[k], pr_data[k],
                                ps_data[k], DAY_LENGTH, tair_day,
                                tsoil_day, rh_day, vpd_day, vmfd_day,
                                radabv_day, fbeam_day, pr_day,
                                ps_day)

                # rs output is the sum of the PAR and NIR components of
                # radabv_day. radabv is a 3-column matrix flattened to a
                # column-major buffer.
                for ihr in 1:DAY_LENGTH
                    rs_out[start + ihr - 1] = radabv_out[ihr] +
                                              radabv_out[ihr + DAY_LENGTH]
                end
            end
        end # iteration through times
        @debug "Successfully generated entire timeseries for gridcell ($i, $j)"

        # Convert VPD from Pa to kPa.
        vpd_out /= 1000

        # Write data for this gridcell to the output files.
        write(indices_out.idx_temp, tair_out, i, j)
        write(indices_out.idx_vpd, vpd_out, i, j)
        write(indices_out.idx_rs, rs_out, i, j)
        write(indices_out.idx_pr, pr_out, i, j)
        write(indices_out.idx_ps, ps_out, i, j)

        # Write progress message after processing each gridcell.
        # In MPI mode, effectively all IO is collective, so we can assume that
        # all workers are making progress at the same rate.
        if opts.show_progress && (!opts.parallel || get_rank() == 0)
            progress = (cell - first(workitems) + 1) / workload_size
            percent = progress * 100
            elapsed = time() - start_time
            total = elapsed / progress
            remaining = total - elapsed
            @info "Progress: $(round(percent, digits=2))% (Elapsed: $(round(elapsed, digits=2))s, Remaining: $(round(remaining, digits=2))s)"
        end
    end # iteration through assigned gridcells

    if !opts.parallel
        return
    end

    # Perform busy wait consisting of zero-length reads and writes that mimics
    # the access patterns of the remaining workers.
    # Number of iterations is max workload size - workload size of this worker.
    niter = get_max_workload_size(opts, length(lats), length(lons)) - workload_size
    @info "Performing busy wait for $niter iterations to allow other workers to finish"
    for _ in 1:niter
        # For each gridcell, a real worker reads the entire timeseries of
        # tmin, tmax, rs, pr, and conditionally, ps. Use full-shape reads
        # here so all workers participate in matching collective operations.

        read_raw(indices_in.idx_tmin, 1, 1, 0)
        read_raw(indices_in.idx_tmax, 1, 1, 0)
        read_raw(indices_in.idx_rs, 1, 1, 0)
        read_raw(indices_in.idx_pr, 1, 1, 0)
        if indices_in.idx_ps !== nothing
            read_raw(indices_in.idx_ps, 1, 1, 0)
        end

        write(indices_out.idx_temp, Float32[], 1, 1)
        write(indices_out.idx_vpd, Float32[], 1, 1)
        write(indices_out.idx_rs, Float32[], 1, 1)
        write(indices_out.idx_pr, Float32[], 1, 1)
        write(indices_out.idx_ps, Float32[], 1, 1)
    end
end

function initialise_output_files(opts::Options)
    if opts.parallel && get_rank() != 0
        # Only the master node should create output files, to avoid conflicts.
        return
    end

    # User-specified output chunk sizes.
    sizes = ChunkSizes(opts.chunk_size_lon, opts.chunk_size_lat,
                       opts.chunk_size_time)

    # Get a reference to the tmin dataset. This will be useful as a template.
    dim_order_default, chunk_size_default = NCDataset(opts.in_tmin) do tmin
        # Create output files with coordinate variables.
        create_output_files(opts, tmin)

        # Default dimension order for created-from-scratch variables can be taken
        # from tmin input file.
        idx_tmin = validate_variable_from_name(tmin, opts.name_tmin)
        ([dimnames(tmin[opts.name_tmin])...], get_chunk_size(idx_tmin, sizes))
    end

    # Initialise data variables in output files.
    NCDataset(opts.in_rs) do nc_rs
        idx_rs = validate_variable_from_std_name(nc_rs, STD_RS)
        create_var_from_existing(opts.out_rs, idx_rs, opts.compression_level,
                                 UNITS_RS, sizes)
    end

    NCDataset(opts.in_pr) do nc_pr
        idx_pr = validate_variable_from_std_name(nc_pr, STD_PR)
        create_var_from_existing(opts.out_pr, idx_pr, opts.compression_level,
                                 UNITS_PR, sizes)
    end

    # Create air temperature variable.
    create_var(opts.out_temp, opts.out_name_temp, UNITS_TEMP, STD_TEMP,
               LONG_TEMP, dim_order_default, opts.compression_level,
               chunk_size_default)

    # Create VPD variable.
    create_var(opts.out_vpd, opts.out_name_vpd, UNITS_VPD, STD_VPD, LONG_VPD,
               dim_order_default, opts.compression_level,
               chunk_size_default)

    # Air pressure can use same dimension order and suitable chunk sizes as
    # input file (if one is provided). Otherwise, use defaults based on tmin
    # input file.
    dims_ps, chunks_ps = NCDataset(opts.in_ps) do nc_ps
        if has_var_with_std_name(nc_ps, STD_PS)
            idx_ps = validate_variable_from_std_name(nc_ps, STD_PS)
            ([dimnames(idx_ps.var)...], get_chunk_size(idx_ps, sizes))
        else
            @info "Using fixed air pressure = $(opts.default_ps) $(UNITS_PS)"
            (dim_order_default, chunk_size_default)
        end
    end

    # Create air pressure variable.
    create_var(opts.out_ps, NAME_PS, UNITS_PS, STD_PS, LONG_PS,
               dims_ps, opts.compression_level, chunks_ps)
end

function open_netcdf(f::Function, path::AbstractString, mode::AbstractString,
                     opts::Options)
    open() = begin
        if opts.parallel
            return NCDataset(MPI.COMM_WORLD, path, mode)
        else
            return NCDataset(path, mode)
        end
    end
    nc = open()
    if opts.parallel
        @info "Setting collective access mode for file $(NCDatasets.path(nc))"
        NCDatasets.paraccess(nc, :collective)
    end

    @debug "Successfully opened path $(path): ncid=$(nc.ncid)"

    try
        return f(nc)
    finally
        close(nc)
    end
end

function with_open_netcdfs(f::Function, paths::Vector{String},
                           mode::AbstractString, opts::Options)
    unique_paths = unique(paths)
    datasets = Dict{String, NCDataset}()

    function open_next(i::Int)
        if i > length(unique_paths)
            return f(datasets)
        end

        path = unique_paths[i]
        open_netcdf(path, mode, opts) do nc
            datasets[path] = nc
            return open_next(i + 1)
        end
    end

    return open_next(1)
end

function process_data(opts::Options, tmin::NCDataset, tmax::NCDataset,
                      rs::NCDataset, pr::NCDataset, ps::NCDataset)
    # Validate dimensions.
    validate_time_axes([tmin, tmax, rs, pr, ps])
    validate_spatial_axes([tmin, tmax, rs, pr, ps], STD_LON)
    validate_spatial_axes([tmin, tmax, rs, pr, ps], STD_LAT)

    # Validate variables and get dimension indices.
    idx_tmin = validate_variable_from_name(tmin, opts.name_tmin)
    idx_tmax = validate_variable_from_name(tmax, opts.name_tmax)
    idx_rs = validate_variable_from_std_name(rs, STD_RS)
    idx_pr = validate_variable_from_std_name(pr, STD_PR)
    dynamic_ps = has_var_with_std_name(ps, STD_PS)

    idx_ps = nothing
    if dynamic_ps
        idx_ps = validate_variable_from_std_name(ps, STD_PS)
    end

    idx_in = InputDimensionOrders(idx_tmin, idx_tmax, idx_rs, idx_pr, idx_ps)

    with_open_netcdfs([opts.out_temp, opts.out_vpd, opts.out_rs, opts.out_pr,
                       opts.out_ps], "a", opts) do out_datasets
        nc_out_temp = out_datasets[opts.out_temp]
        nc_out_vpd = out_datasets[opts.out_vpd]
        nc_out_rs = out_datasets[opts.out_rs]
        nc_out_pr = out_datasets[opts.out_pr]
        nc_out_ps = out_datasets[opts.out_ps]

        idx_temp = validate_variable_from_std_name(nc_out_temp, STD_TEMP)
        idx_vpd = validate_variable_from_std_name(nc_out_vpd, STD_VPD)
        idx_rs_out = validate_variable_from_std_name(nc_out_rs, STD_RS)
        idx_pr_out = validate_variable_from_std_name(nc_out_pr, STD_PR)
        idx_ps_out = validate_variable_from_std_name(nc_out_ps, STD_PS)

        idx_out = OutputDimensionOrders(idx_rs_out, idx_pr_out, idx_temp,
                                        idx_vpd, idx_ps_out)

        generate_weather(opts, idx_in, idx_out)
    end
end

function make_mpi_logger(level::Logging.LogLevel, rank::Int)
    return ConsoleLogger(stdout, level; meta_formatter = (lvl, _module, group, id, file, line) -> begin
        color, prefix, suffix = Logging.default_metafmt(lvl, _module, group, id, file, line)
        return color, "[rank=$rank] $prefix", suffix
    end)
end

function make_std_logger(level::Logging.LogLevel)
    return ConsoleLogger(stdout, level)
end

function make_logger(opts::Options)
    if opts.parallel
        return make_mpi_logger(opts.log_level, get_rank())
    else
        return make_std_logger(opts.log_level)
    end
end

function main(opts::Options)
    wg_init()

    initialise_output_files(opts)
    if opts.parallel
        barrier("post-init")
    end
    @info "Output files initialised; starting weather generation"

    # Open input files for reading.
    with_open_netcdfs([opts.in_tmin, opts.in_tmax, opts.in_rs, opts.in_pr,
                       opts.in_ps], "r", opts) do in_datasets
        tmin = in_datasets[opts.in_tmin]
        tmax = in_datasets[opts.in_tmax]
        rs = in_datasets[opts.in_rs]
        pr = in_datasets[opts.in_pr]
        ps = in_datasets[opts.in_ps]

        process_data(opts, tmin, tmax, rs, pr, ps)
    end
end

# Main CLI entrypoint function. Initialises logging and MPI, and runs the
# generator. Does not swallow exceptions.
function cli_main(opts::Options)
    if opts.parallel
        @eval using MPI
        MPI.Init()
    end

    # Logging initialisation must happen after MPI, because in parallel mode,
    # the logger needs to know the MPI rank to include in log messages.
    logger = make_logger(opts)
    global_logger(logger)

    # Can't emit this warning until after logging has been initialised.
    if opts.parallel
        @info "Running in MPI mode. World size is $(get_world_size()), rank is $(get_rank())"
        if get_world_size() == 1
            @warn "MPI parallelism enabled but only one process detected; running in serial. This is almost certainly not what you want. To fix, run with multiple processes (e.g. using mpirun)."
        end
    end

    try
        main(opts)
        if opts.parallel && MPI.Initialized() && !MPI.Finalized()
            MPI.Finalize()
        end
    catch err
        if opts.parallel
            @error "Error on rank $(get_rank()): $err"
            Base.display_error(err, catch_backtrace())
            MPI.Abort(MPI.COMM_WORLD, 1)
        end
        rethrow()
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    opts = Weathergen.parse_cli()
    Weathergen.cli_main(opts)
end
