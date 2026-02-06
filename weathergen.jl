#!/usr/bin/env julia

using ArgParse
using Logging
using NCDatasets
using Dates
using DataStructures

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

# Standard name of the air pressure variable in the output file.
const STD_PS = "air_pressure"

# Long name of the air pressure variable in the output file.
const LONG_PS = "Air pressure"

################################################################################
# Types
################################################################################

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
end

struct DimensionIndices
    var::NCDatasets.CFVariable
    index_lon::Int
    index_lat::Int
    index_time::Int
end

struct PerVariablePaths
    input_file::String
    output_file::String
    input_arg_name::String
    output_arg_name::String
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
            default=101300.0f0
            help="Default air pressure (Pa) to use if no input file is provided."
        "--verbosity", "-v"
            arg_type=Int
            default=2
            help="Verbosity level (0: errors, 1: warnings, 2: info, 3: debug)"
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
        validate_file_path(file_ps, "--file-ps")

        ensure_set(out_temp, "--out-temp")
        ensure_set(out_pr, "--out-pr")
        ensure_set(out_ps, "--out-ps")
        ensure_set(out_rs, "--out-rs")
        ensure_set(out_vpd, "--out-vpd")

        validate_per_variable_paths([
            PerVariablePaths(file_tmin, out_temp, "--file-tmin", "--out-temp"),
            PerVariablePaths(file_tmax, out_temp, "--file-tmax", "--out-temp"),
            PerVariablePaths(file_rs, out_rs, "--file-rs", "--out-rs"),
            PerVariablePaths(file_pr, out_pr, "--file-pr", "--out-pr"),
            PerVariablePaths(file_ps, out_ps, "--file-ps", "--out-ps"),
        ])
    end

    return Options(parsed["seed"], file_tmin, file_tmax, file_rs, file_pr,
                   file_ps, out_temp, out_pr, out_ps, out_rs, out_vpd,
                   parsed["name-tmax"], parsed["name-tmin"],
                   parsed["out-name-temp"], parsed["out-name-vpd"],
                   log_level, parsed["compression-level"],
                   parsed["chunk-lon"], parsed["chunk-lat"],
                   parsed["chunk-time"], parsed["default-ps"])
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

# Seed wrapper
function wg_seed(seed::Int)::Cint
    return ccall((:wg_seed, libwg), Cint, (Int64,), Int64(seed))
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

function read_variable(var::NCDatasets.CFVariable, idx::DimensionIndices,
                       i::Int, j::Int, units::String)
    hyperslab = selectdim(var, idx.index_lat, i)
    hyperslab = selectdim(hyperslab, idx.index_lon, j)
    data = hyperslab[:]

    # Error if any data is missing.
    if any(ismissing, data)
        error("Missing data in variable $(name(var)) at gridcell ($i, $j)")
    end

    # Convert from Vector{Union{Float32, Missing}} to Vector{Float32}.
    data = convert(Vector{Float32}, data)

    # Get timestep width in seconds.
    time = var_from_std_name(var.var.ds, STD_TIME)
    dt = Second(Dates.value(time[2] - time[1]))

    return convert_units(data, var.attrib[ATTR_UNITS], units, dt.value)
end

"""Convert a `DateTime` to Maespa's `idate` (days since 1950-01-01).

This matches the Fortran calendar math used by Maespa (Julian-style leap
years: every 4 years, with no century exception), so `JDATE(idate)` returns
the correct day-of-year.
"""
function date_to_idate(date::DateTime)
    # Note: can't use Dates.value() because the underlying fortran code uses a
    # Julian calendar, with leap days exactly every 4 years.

    yearValue = year(date)
    monthValue = month(date)
    dayValue = day(date)

    ifd = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

    dayOfYear = ifd[monthValue] + dayValue
    isLeapYear = (4 * (yearValue ÷ 4) == yearValue)
    if isLeapYear && monthValue >= 3
        dayOfYear += 1
    end

    yearsSince1950 = yearValue - 1950
    daysBeforeYear = 365 * yearsSince1950 + div(yearsSince1950 - 1, 4)

    return daysBeforeYear + dayOfYear - 1
end

function write_outputs(out_file::String, name::String, data::Vector{Float32},
                       ilat::Int, ilon::Int, indices::DimensionIndices)
    NCDataset(out_file, "a") do nc
        var = nc[name]

        hyperslab = selectdim(var, indices.index_lat, ilat)
        hyperslab = selectdim(hyperslab, indices.index_lon, ilon)
        hyperslab[:] = data
    end
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

function init_outfile(path::String, name::String, dims::Vector{String},
                      units::String, std_name::String, long_name::String,
                      compression_level::Int, chunk_sizes::Vector{Int})
    NCDataset(path, "a") do nc_out
        var = defVar(nc_out, name, Float32, dims,
                     deflatelevel=compression_level,
                     shuffle=compression_level > 0,
                     chunksizes=chunk_sizes)
        var.attrib[ATTR_UNITS] = units
        var.attrib[ATTR_STD_NAME] = std_name
        var.attrib[ATTR_LONG_NAME] = long_name
    end
end

function init_outfile(path::String, nc_in::NCDataset, out_var_name::String,
                      in_var_name::String, compression_level::Int,
                      chunk_sizes::Vector{Int})
    # Open the file and create the required output variable.
    NCDataset(path, "a") do nc_out
        # Get the input variable.
        in_var = nc_in[in_var_name]

        # Create the output variable.
        # var = add_variable(nc_out, out_var_name, Float32, dimnames(in_var))
        # defVar(nc_out, name(in_lon), in_lon[:], dimnames(in_lon))
        var = defVar(nc_out, out_var_name, Float32, dimnames(in_var),
                     deflatelevel=compression_level,
                     shuffle=compression_level > 0,
                     chunksizes=chunk_sizes)
        copy_attributes(in_var, var)
    end
end

# Convenience function for when the output variable name is the same as the
# input variable name.
function init_outfile(path::String, nc_in::NCDataset, var_name::String,
                      compression_level::Int, chunk_sizes::Vector{Int})
    init_outfile(path, nc_in, var_name, var_name, compression_level,
                 chunk_sizes)
end

function init_outfile(path::String, var_name::String, units::String,
                      std_name::String, long_name::String,
                      dims::Vector{String},
                      compression_level::Int,
                      chunk_sizes::Vector{Int})
    NCDataset(path, "a") do nc
        var = defVar(nc, var_name, Float32, dims,
               deflatelevel=compression_level,
               shuffle=compression_level > 0,
               chunksizes=chunk_sizes)

        var.attrib[ATTR_UNITS] = units
        var.attrib[ATTR_STD_NAME] = std_name
        var.attrib[ATTR_LONG_NAME] = long_name
    end
end

function init_outfiles(opts::Options, nc_in::NCDataset)
    paths = unique([opts.out_temp, opts.out_rs, opts.out_pr, opts.out_ps,
                    opts.out_vpd])

    compression = opts.compression_level
    shuffle = compression > 0

    for path in paths
        # Create directory if it doesn't already exist.
        dir = dirname(path)
        if !isdir(dir)
            mkdir(dir)
        end

        NCDataset(path, "c") do nc_out
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

            # TODO: prevent chunk size from exceeding 4GiB. Same for lat.
            lons = in_lon[:]
            out_lon = defVar(nc_out, name(in_lon), lons, dimnames(in_lon),
                             deflatelevel=compression, shuffle=shuffle,
                             chunksizes=[length(lons)])
            copy_attributes(in_lon, out_lon)

            lats = in_lat[:]
            out_lat = defVar(nc_out, name(in_lat), lats, dimnames(in_lat),
                             deflatelevel=compression, shuffle=shuffle,
                             chunksizes=[length(lats)])
            copy_attributes(in_lat, out_lat)

            # Construct hourly timeseries from each day in the input time
            # variable.
            # TODO: rethink time chunking. Is input (typically 8760) too small?
            times = in_time[:]
            hours = [t + Hour(h) for t in times for h in 0:(DAY_LENGTH - 1)]
            out_time = defVar(nc_out, name(in_time), hours, dimnames(in_time),
                              attrib = OrderedDict(
                                  ATTR_UNITS => in_time.attrib[ATTR_UNITS],
                                  "calendar" => in_time.attrib["calendar"],
                              ),
                              deflatelevel=compression, shuffle=shuffle,
                              chunksizes=[opts.chunk_size_time])
            copy_attributes(in_time, out_time)
        end
    end
end

function get_chunk_size(indices::DimensionIndices, opts::Options)
    chunks = [1, 1, 1]
    chunks[indices.index_time] = opts.chunk_size_time
    chunks[indices.index_lat] = opts.chunk_size_lat
    chunks[indices.index_lon] = opts.chunk_size_lon
    return chunks
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

    if dynamic_ps
        idx_ps = validate_variable_from_std_name(ps, STD_PS)
    else
        idx_ps = idx_tmin
    end

    # Initialise PRNG seed.
    wg_seed(opts.seed)

    # Iterate through gridcells. (All input files use the same grid.)
    var_lon = var_from_std_name(tmin, STD_LON)
    var_lat = var_from_std_name(tmin, STD_LAT)
    var_time = var_from_std_name(tmin, STD_TIME)

    lons = var_lon[:]
    lats = var_lat[:]
    times = var_time[:]

    # Create output files with coordinate variables.
    init_outfiles(opts, tmin)

    # Initialise data variables in output files.
    # path::String, nc_in::NCDataset, out_var_name::String, in_var_name::String
    init_outfile(opts.out_rs, rs, name(idx_rs.var), opts.compression_level,
                 get_chunk_size(idx_rs, opts))
    init_outfile(opts.out_pr, pr, name(idx_pr.var), opts.compression_level,
                 get_chunk_size(idx_pr, opts))

    if dynamic_ps
        init_outfile(opts.out_ps, ps, name(idx_ps.var), opts.compression_level,
                    get_chunk_size(idx_ps, opts))
        ps_var = name(idx_ps.var)
    else
        @info "Using fixed air pressure = $(opts.default_ps) $(UNITS_PS)"
        init_outfile(opts.out_ps, NAME_PS, UNITS_PS, STD_PS, LONG_PS,
                     [dimnames(tmin[opts.name_tmin])...],
                     opts.compression_level, get_chunk_size(idx_tmin, opts))
        ps_var = NAME_PS
    end

    # Temperature can be created by copying metadata from tmin input file.
    init_outfile(opts.out_temp, tmin, opts.out_name_temp, opts.name_tmin,
                 opts.compression_level, get_chunk_size(idx_tmin, opts))
    idx_temp = idx_tmin # same dimension order as tmin

    # VPD must be created from scratch. We can use same dimension order as tmin
    # input file.
    init_outfile(opts.out_vpd, opts.out_name_vpd,
                 [dimnames(tmin[opts.name_tmin])...],
                 "kPa", "vapour_pressure_deficit", "Vapour pressure deficit",
                 opts.compression_level, get_chunk_size(idx_tmin, opts))
    idx_vpd = idx_tmin # same dimension order as tmin

    # Iterate through gridcells. Generate climate one gridcell at a time.
    for i in eachindex(lats)
        for j in eachindex(lons)
            lat = lats[i]
            lon = lons[j]

            @info "Processing gridcell $i, $j ($lat, $lon)"

            # Read timeseries for this gridcell.
            tmin_data = read_variable(idx_tmin.var, idx_tmin, i, j, "degC")
            tmax_data = read_variable(idx_tmax.var, idx_tmax, i, j, "degC")
            rs_data = read_variable(idx_rs.var, idx_rs, i, j, "W m-2")
            pr_data = read_variable(idx_pr.var, idx_pr, i, j, "mm")
            if dynamic_ps
                ps_data = read_variable(idx_ps.var, idx_ps, i, j, "Pa")
            else
                ps_data = fill(opts.default_ps, length(times))
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
                @debug "Generating climate for day: $(times[k])"

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
                    rs_day = pointer(rs_out, start)
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

            # Write data for this gridcell to the output files.
            name_temp = opts.out_name_temp
            name_vpd = opts.out_name_vpd
            write_outputs(opts.out_temp, name_temp, tair_out, i, j, idx_temp)
            write_outputs(opts.out_vpd, name_vpd, vpd_out, i, j, idx_vpd)
            write_outputs(opts.out_rs, name(idx_rs.var), rs_out, i, j, idx_rs)
            write_outputs(opts.out_pr, name(idx_pr.var), pr_out, i, j, idx_pr)
            write_outputs(opts.out_ps, ps_var, ps_out, i, j, idx_ps)
        end # iteration through lons
    end # iteration through lats
end

function main(opts::Options)
    # Open input files for reading.
    NCDataset(opts.in_tmin) do tmin
        NCDataset(opts.in_tmax) do tmax
            NCDataset(opts.in_rs) do rs
                NCDataset(opts.in_pr) do pr
                    NCDataset(opts.in_ps) do ps
                        process_data(opts, tmin, tmax, rs, pr, ps)
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
