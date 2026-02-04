#!/usr/bin/env julia

using ArgParse
using Logging
using NCDatasets
using Dates
using DataStructures

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

# Number of timesteps per day.
const DAY_LENGTH = 24

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
end

struct DimensionIndices
    var::NCDatasets.CFVariable
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
        "--verbosity", "-v"
            arg_type=Int
            default=2
            help="Verbosity level (0: errors, 1: warnings, 2: info, 3: debug)"
    end
    parsed = parse_args(parser)
    log_level = parse_log_level(parsed["verbosity"])
    return Options(parsed["seed"], parsed["file-tmin"], parsed["file-tmax"],
                   parsed["file-rs"], parsed["file-pr"], parsed["file-ps"],
                   parsed["out-temp"], parsed["out-pr"],
                   parsed["out-ps"], parsed["out-rs"],
                   parsed["out-vpd"], parsed["name-tmax"], parsed["name-tmin"],
                   parsed["out-name-temp"], parsed["out-name-vpd"],
                   log_level)
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
    if !haskey(var.attrib, ATTR_UNITS)
        error("Variable $(name(var)) has no units attribute")
    end

    if var.attrib[ATTR_UNITS] != units
        error("Variable $(name(var)) has units $(var.attrib[ATTR_UNITS]) but expected $units")
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

function read_variable(var::NCDatasets.CFVariable, idx::DimensionIndices, i::Int, j::Int)
    hyperslab = selectdim(var, idx.index_lat, i)
    hyperslab = selectdim(hyperslab, idx.index_lon, j)
    data = hyperslab[:]

    # Error if any data is missing.
    if any(ismissing, data)
        error("Missing data in variable $(name(var)) at gridcell ($i, $j)")
    end

    return data
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

function init_outfile(path::String, name::String, dims::Vector{String}, units::String, std_name::String, long_name::String)
    NCDataset(path, "a") do nc_out
        var = defVar(nc_out, name, Float32, dims)
        var.attrib[ATTR_UNITS] = units
        var.attrib[ATTR_STD_NAME] = std_name
        var.attrib[ATTR_LONG_NAME] = long_name
    end
end

function init_outfile(path::String, nc_in::NCDataset, out_var_name::String, in_var_name::String)
    # Open the file and create the required output variable.
    NCDataset(path, "a") do nc_out
        # Get the input variable.
        in_var = nc_in[in_var_name]

        # Create the output variable.
        # var = add_variable(nc_out, out_var_name, Float32, dimnames(in_var))
        # defVar(nc_out, name(in_lon), in_lon[:], dimnames(in_lon))
        var = defVar(nc_out, out_var_name, Float32, dimnames(in_var))
        copy_attributes(in_var, var)
    end
end

# Convenience function for when the output variable name is the same as the
# input variable name.
function init_outfile(path::String, nc_in::NCDataset, var_name::String)
    init_outfile(path, nc_in, var_name, var_name)
end

function init_outfiles(opts::Options, nc_in::NCDataset)
    paths = unique([opts.out_temp, opts.out_rs, opts.out_pr, opts.out_ps,
                    opts.out_vpd])
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

            out_lon = defVar(nc_out, name(in_lon), in_lon[:], dimnames(in_lon))
            copy_attributes(in_lon, out_lon)

            out_lat = defVar(nc_out, name(in_lat), in_lat[:], dimnames(in_lat))
            copy_attributes(in_lat, out_lat)

            # Construct hourly timeseries from each day in the input time
            # variable.
            times = in_time[:]
            hours = [t + Hour(h) for t in times for h in 0:(DAY_LENGTH - 1)]
            out_time = defVar(nc_out, name(in_time), hours, dimnames(in_time),
                              attrib = OrderedDict(
                ATTR_UNITS => in_time.attrib[ATTR_UNITS],
                "calendar" => in_time.attrib["calendar"],
            ))
            copy_attributes(in_time, out_time)
        end
    end
end

function process_data(opts::Options, tmin::NCDataset, tmax::NCDataset, rs::NCDataset, pr::NCDataset, ps::NCDataset)
    # Validate dimensions.
    validate_time_axes([tmin, tmax, rs, pr, ps])
    validate_spatial_axes([tmin, tmax, rs, pr, ps], STD_LON)
    validate_spatial_axes([tmin, tmax, rs, pr, ps], STD_LAT)

    # Validate variables and get dimension indices.
    idx_tmin = validate_variable_from_name(tmin, opts.name_tmin, "degC")
    idx_tmax = validate_variable_from_name(tmax, opts.name_tmax, "degC")
    idx_rs = validate_variable_from_std_name(rs, "surface_downwelling_shortwave_flux_in_air", "W m-2")
    idx_pr = validate_variable_from_std_name(pr, "precipitation_amount", "mm")
    idx_ps = validate_variable_from_std_name(ps, "air_pressure", "Pa")

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
    init_outfile(opts.out_rs, rs, name(idx_rs.var))
    init_outfile(opts.out_pr, pr, name(idx_pr.var))
    init_outfile(opts.out_ps, ps, name(idx_ps.var))

    # Temperature can be created by copying metadata from tmin input file.
    init_outfile(opts.out_temp, tmin, opts.out_name_temp, opts.name_tmin)
    idx_temp = idx_tmin # same dimension order as tmin

    # VPD must be created from scratch. We can use same dimension order as tmin
    # input file.
    init_outfile(opts.out_vpd, opts.out_name_vpd,
                 [dimnames(tmin[opts.name_tmin])...],
                 "kPa", "vapour_pressure_deficit", "Vapour pressure deficit")
    idx_vpd = idx_tmin # same dimension order as tmin

    # Iterate through gridcells. Generate climate one gridcell at a time.
    for i in eachindex(lats)
        for j in eachindex(lons)
            lon = lons[j]
            lat = lats[i]

            @info "Processing gridcell $i, $j ($lon, $lat)"

            # Read timeseries for this gridcell.
            tmin_data = read_variable(idx_tmin.var, idx_tmin, i, j)
            tmax_data = read_variable(idx_tmax.var, idx_tmax, i, j)
            rs_data = read_variable(idx_rs.var, idx_rs, i, j)
            pr_data = read_variable(idx_pr.var, idx_pr, i, j)
            ps_data = read_variable(idx_ps.var, idx_ps, i, j)

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
                        rs_out[start + ihr - 1] = radabv_out[ihr] + radabv_out[ihr + DAY_LENGTH]
                    end
                end
            end # iteration through times

            # Write data for this gridcell to the output files.
            write_outputs(opts.out_temp, opts.out_name_temp, tair_out, i, j, idx_temp)
            write_outputs(opts.out_vpd, opts.out_name_vpd, vpd_out, i, j, idx_vpd)
            write_outputs(opts.out_rs, name(idx_rs.var), rs_out, i, j, idx_rs)
            write_outputs(opts.out_pr, name(idx_pr.var), pr_out, i, j, idx_pr)
            write_outputs(opts.out_ps, name(idx_ps.var), ps_out, i, j, idx_ps)
        end # iteration through lons
    end # iteration through lats
end

function main(opts::Options)
    @info "Running weather generator..."

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
