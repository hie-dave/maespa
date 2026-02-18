#!/usr/bin/env julia
using Test
using Dates
using Statistics
using NCDatasets
using Logging

include("weathergen.jl")

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

function run_weathergen(input_path::AbstractString, output_path::AbstractString)
    seed = 123
    opts = Options(
        seed,          # seed::Int
        input_path,    # in_tmin::String
        input_path,    # in_tmax::String
        input_path,    # in_rs::String
        input_path,    # in_pr::String
        input_path,    # in_ps::String
        output_path,   # out_temp::String
        output_path,   # out_pr::String
        output_path,   # out_ps::String
        output_path,   # out_rs::String
        output_path,   # out_vpd::String
        "tasmax",      # name_tmax::String
        "tasmin",      # name_tmin::String
        "tas",         # out_name_temp::String
        "vpd",         # out_name_vpd::String
        Logging.Error, # log_level::Logging.LogLevel
        5,             # compression_level::Int
        1,             # chunk_size_lon::Int
        1,             # chunk_size_lat::Int
        24,            # chunk_size_time::Int
        101300,        # default_ps::Float32
        false,         # parallel::Bool
    )
    cli_main(opts)
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

@testset "weathergen integration" begin
    mktempdir() do tmp
        input_path = joinpath(tmp, "input_daily.nc")
        output_path = joinpath(tmp, "output_hourly.nc")
        input = write_synthetic_daily_input(input_path)

        run_weathergen(input_path, output_path)
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
