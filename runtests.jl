#!/usr/bin/env -S julia --project=@.
using Test
using Dates
using Statistics
using NCDatasets
using Logging

# If Revise is available, use it to include the weathergen module in a way that
# allows for interactive development. Otherwise, just include the file directly.
if Base.find_package("Revise") !== nothing
    using Revise
    includet("weathergen.jl")
else
    include("weathergen.jl")
end

using .Weathergen

function write_synthetic_daily_input(path::AbstractString)
    lats = Float32[-35.0f0, -34.0f0]
    lons = Float32[149.0f0, 150.0f0]
    times = collect(DateTime(2001, 1, 1):Day(1):DateTime(2001, 1, 6))

    NCDataset(path, "c") do ds
        defDim(ds, "lat", length(lats))
        defDim(ds, "lon", length(lons))
        defDim(ds, "time", length(times))

        vlat = defVar(ds, "lat", lats, ("lat",))
        vlat.attrib["standard_name"] = "latitude"
        vlat.attrib["units"] = "degrees_north"
        vlat.attrib["long_name"] = "latitude"

        vlon = defVar(ds, "lon", lons, ("lon",))
        vlon.attrib["standard_name"] = "longitude"
        vlon.attrib["units"] = "degrees_east"
        vlon.attrib["long_name"] = "longitude"

        defVar(ds, "time", times, ("time",), attrib = Dict(
            "standard_name" => "time",
            "units" => "days since 2000-01-01 00:00:00",
            "calendar" => "standard",
            "long_name" => "time",
        ))

        dims = ("lon", "lat", "time")
        nl = length(lons)
        na = length(lats)
        nt = length(times)

        tmin = Array{Float32}(undef, nl, na, nt)
        tmax = similar(tmin)
        rs = similar(tmin)
        pr = similar(tmin)

        for i in 1:nl, j in 1:na, k in 1:nt
            base = Float32(8 + i + j + 0.4f0 * k)
            tmin[i, j, k] = base
            tmax[i, j, k] = base + 8f0
            rs[i, j, k] = Float32(180 + 5 * i + 3 * j + k)
            pr[i, j, k] = Float32((k % 3 == 0) ? 4.5 : 0.8)
        end

        vtmin = defVar(ds, "tasmin", tmin, dims)
        vtmin.attrib["units"] = "degC"
        vtmin.attrib["standard_name"] = "air_temperature"
        vtmin.attrib["long_name"] = "Daily minimum near-surface air temperature"

        vtmax = defVar(ds, "tasmax", tmax, dims)
        vtmax.attrib["units"] = "degC"
        vtmax.attrib["standard_name"] = "air_temperature"
        vtmax.attrib["long_name"] = "Daily maximum near-surface air temperature"

        vrs = defVar(ds, "rsds", rs, dims)
        vrs.attrib["units"] = "W m-2"
        vrs.attrib["standard_name"] = "surface_downwelling_shortwave_flux_in_air"
        vrs.attrib["long_name"] = "Surface downwelling shortwave radiation"

        vpr = defVar(ds, "pr", pr, dims)
        vpr.attrib["units"] = "mm"
        vpr.attrib["standard_name"] = "precipitation_amount"
        vpr.attrib["long_name"] = "Daily precipitation"
    end

    return (times = times, lons = lons, lats = lats)
end

function mkopts(inpath::AbstractString, outpath::AbstractString; seed::Int=123,
                name_tmax::String="tasmax", name_tmin::String="tasmin",
                out_name_temp::String="tas", out_name_vpd::String="vpd",
                log_level::Logging.LogLevel=Logging.Error,
                compression_level::Int=5,
                chunk_size_lon::Int=1, chunk_size_lat::Int=1, chunk_size_time::Int=24,
                default_ps::Float32=101300f0, parallel::Bool=false)
    return Weathergen.Options(
        seed,
        inpath, inpath, inpath, inpath, inpath,
        outpath, outpath, outpath, outpath, outpath,
        name_tmax, name_tmin, out_name_temp, out_name_vpd,
        log_level, compression_level,
        chunk_size_lon, chunk_size_lat, chunk_size_time,
        default_ps, parallel,
    )
end

function run_weathergen(opts::Weathergen.Options)
    Weathergen.cli_main(opts)
end

function log_level_to_verbosity(level::Logging.LogLevel)
    if level === Logging.Error
        return 0
    elseif level === Logging.Warn
        return 1
    elseif level === Logging.Info
        return 2
    elseif level === Logging.Debug
        return 3
    else
        return 2
    end
end

# Run the weather generator via the CLI in a subprocess, passing options as
# command-line arguments. This will use mpiexecjl for parallel execution if
# opts.parallel is true, using the specified number of processors. If parallel
# is false, nprocs will be ignored.
#
# Per-variable file paths are ignored. in_tmin and out_temp are used as the
# input and output file paths.
#
# This is very inefficient. Don't use this for serial execution.
function run_weathergen_cli(opts::Weathergen.Options; nprocs::Int=2)
    # Map Logging level back to verbosity integer used by the CLI parser.
    verbosity = log_level_to_verbosity(opts.log_level)

    args = ["julia", "--project=@.", "weathergen.jl"]
    push!(args, "-s"); push!(args, string(opts.seed))
    push!(args, "--chunk-lon"); push!(args, string(opts.chunk_size_lon))
    push!(args, "--chunk-lat"); push!(args, string(opts.chunk_size_lat))
    push!(args, "--chunk-time"); push!(args, string(opts.chunk_size_time))
    push!(args, "--verbosity"); push!(args, string(verbosity))
    push!(args, "--default-ps"); push!(args, string(opts.default_ps))
    if opts.parallel
        push!(args, "--parallel")
    end

    # Assume single-file mode: provide input-file and output-file
    push!(args, "-i"); push!(args, opts.in_tmin)
    push!(args, "-o"); push!(args, opts.out_temp)

    if opts.parallel
        args = vcat(["mpiexecjl", "-n", string(nprocs)], args)
    end

    # Build command string and run via shell to avoid splicing issues
    run(`sh -c $(join(args, " "))`)
end

function read_output(path::AbstractString)
    NCDataset(path) do ds
        return (
            tas = Float32.(coalesce.(ds["tas"][:], NaN32)),
            vpd = Float32.(coalesce.(ds["vpd"][:], NaN32)),
            rs = Float32.(coalesce.(ds["rsds"][:], NaN32)),
            pr = Float32.(coalesce.(ds["pr"][:], NaN32)),
            ps = Float32.(coalesce.(ds["ps"][:], NaN32)),
            pr3d = Float32.(coalesce.(ds["pr"][:, :, :], NaN32)),
            dims_tas = size(ds["tas"]),
            dims_time = length(ds["time"][:]),
        )
    end
end

@testset "calendar-preserving hourly timestamps" begin
    calendar_types = (
        NCDatasets.CFTime.DateTimeStandard,
        NCDatasets.CFTime.DateTimeJulian,
        NCDatasets.CFTime.DateTimeProlepticGregorian,
        NCDatasets.CFTime.DateTimeAllLeap,
        NCDatasets.CFTime.DateTimeNoLeap,
        NCDatasets.CFTime.DateTime360Day,
    )

    for T in calendar_types
        daily_time = T(2016, 2, 28, 12)
        midnight = Weathergen.start_of_day(daily_time)
        hourly_times = [midnight + Hour(h) for h in 0:23]

        @test typeof(midnight) == typeof(daily_time)
        @test (year(midnight), month(midnight), day(midnight)) == (2016, 2, 28)
        @test (hour(first(hourly_times)), hour(last(hourly_times))) == (0, 23)
        @test all(typeof(t) == typeof(daily_time) for t in hourly_times)
    end

    daily_time = DateTime(2016, 2, 28, 12)
    midnight = Weathergen.start_of_day(daily_time)
    @test midnight == DateTime(2016, 2, 28)
    @test typeof(midnight) == typeof(daily_time)
end

@testset "weathergen integration" begin
    mktempdir() do tmp
        input_path = joinpath(tmp, "input_daily.nc")
        output_path = joinpath(tmp, "output_hourly.nc")
        input = write_synthetic_daily_input(input_path)

        # Run serial generation via the API (same as existing tests)
        opts = mkopts(input_path, output_path)
        run_weathergen(opts)
        out = read_output(output_path)

        expected_time = length(input.times) * 24
        @test out.dims_tas == (length(input.lons), length(input.lats), expected_time)
        @test out.dims_time == expected_time

        @test all(isfinite, out.tas)
        @test all(isfinite, out.vpd)
        @test all(isfinite, out.rs)
        @test all(isfinite, out.pr)
        @test all(isfinite, out.ps)

        @test minimum(out.rs) >= 0
        @test minimum(out.pr) >= 0
        @test minimum(out.vpd) >= 0

        # Daily precipitation should be conserved within numerical tolerance.
        expected_daily_pr = [0.8f0, 0.8f0, 4.5f0, 0.8f0, 0.8f0, 4.5f0]
        for k in eachindex(input.times)
            t0 = (k - 1) * 24 + 1
            t1 = k * 24
            @test out.pr3d[1, 1, t0:t1] |> sum ≈ expected_daily_pr[k] atol=1e-5
        end

        @test minimum(out.ps) ≈ 101300f0 atol=1f-3
        @test maximum(out.ps) ≈ 101300f0 atol=1f-3

        # Compact golden signatures to detect behavior regressions.
        @test mean(out.tas) ≈ 16.108633f0 atol=1f-5
        @test sum(out.tas) ≈ 9278.572f0 atol=1f-2

        @test mean(out.vpd) ≈ 0.41560256f0 atol=1f-6
        @test sum(out.vpd) ≈ 239.38707f0 atol=1f-3

        @test mean(out.rs) ≈ 195.5f0 atol=1f-4
        @test sum(out.rs) ≈ 112608.0f0 atol=1f-1

        @test mean(out.pr) ≈ 0.08472223f0 atol=1f-7
        @test sum(out.pr) ≈ 48.800003f0 atol=1f-4

        @test out.tas[1:8] ≈ Float32[
            12.190443, 13.190443, 13.192387, 14.192387,
            11.771329, 12.771329, 12.777161, 13.777161,
        ] atol=1f-5

        @test out.vpd[1:8] ≈ Float32[
            0.159512, 0.169053, 0.169246, 0.179282,
            0.120668, 0.127898, 0.128464, 0.136096,
        ] atol=1f-5

        @test out.pr[1:12] ≈ Float32[
            0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0, 0, 0,
        ] atol=1f-6
    end
end

@testset "weathergen parallel parity" begin
    mktempdir() do tmp
        input_path = joinpath(tmp, "input_daily.nc")
        out_serial = joinpath(tmp, "output_serial.nc")
        out_parallel = joinpath(tmp, "output_parallel.nc")
        input = write_synthetic_daily_input(input_path)

        # Run serial generation via the API (same as existing tests)
        opts_serial = mkopts(input_path, out_serial)
        run_weathergen(opts_serial)

        # If mpiexecjl is not available, skip the parallel parity test.
        if Sys.which("mpiexecjl") === nothing
            @info "mpiexecjl not found; skipping parallel parity test"
            return
        end

        # Run the CLI in parallel using mpiexecjl. Use 2 processes.
        # Pass the same deterministic seed and small chunk sizes so NetCDF
        # chunking is valid for the small test dataset.
        opts_parallel = mkopts(input_path, out_parallel; parallel=true)
        try
            run_weathergen_cli(opts_parallel; nprocs=2)
        catch err
            @error "Parallel weathergen CLI failed" err
            rethrow()
        end

        # Helper to open dataset and extract info for comparison.
        function ds_summary(path)
            NCDataset(path) do ds
                vars = Dict{String,Any}()
                for vn in ["tas", "vpd", "rsds", "pr", "ps"]
                    if haskey(ds, vn)
                        arr = Array(ds[vn][:])
                        attrs = Dict{String,Any}(collect(ds[vn].attrib))
                        vars[vn] = (data = arr, attrib = attrs, dims = size(ds[vn]))
                    end
                end

                coords = Dict(
                    "lat" => Array(ds["lat"][:]),
                    "lon" => Array(ds["lon"][:]),
                    "time" => Array(ds["time"][:]),
                )

                global_attribs = Dict{String,Any}(collect(ds.attrib))

                return (vars = vars, coords = coords, global_attribs = global_attribs)
            end
        end

        s = ds_summary(out_serial)
        p = ds_summary(out_parallel)

        # Compare coordinates
        @test keys(s.coords) == keys(p.coords)
        for k in keys(s.coords)
            @test s.coords[k] == p.coords[k]
        end

        # Compare global attributes
        @test s.global_attribs == p.global_attribs

        # Compare variables: existence, attributes, dimensions and data
        @test keys(s.vars) == keys(p.vars)
        for vn in keys(s.vars)
            sv = s.vars[vn]
            pv = p.vars[vn]
            @test sv.dims == pv.dims
            @test sv.attrib == pv.attrib
            @test sv.data == pv.data
        end
    end
end
