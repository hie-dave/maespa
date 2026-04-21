#!/usr/bin/env -S julia --project=@.

using Dates
using NCDatasets

const DRY_RUN = true

function usage_and_exit()
    println("Usage: ./fix_time_axis.jl <path/to/file.nc>")
    exit(1)
end

length(ARGS) == 1 || usage_and_exit()
path = ARGS[1]
mode = DRY_RUN ? "r" : "a"

NCDataset(path, mode) do ds
    haskey(ds, "time") || error("Variable 'time' not found in file: $path")
    tvar = ds["time"]

    times = tvar[:]
    isempty(times) && error("Variable 'time' is empty in file: $path")

    t0 = first(times)
    offset = Hour(hour(t0)) + Minute(minute(t0)) + Second(second(t0))

    if offset == Hour(0)
        println("No change needed: first timestamp already at midnight ($t0)")
        exit(0)
    end

    if DRY_RUN
        println("DRY_RUN enabled: would shift time axis by -$offset")
        println("First timestamp would be: $t0 -> $(t0 - offset)")
        exit(0)
    end

    tvar[:] = [t - offset for t in times]
    println("Shifted time axis by -$offset")
    println("First timestamp: $t0 -> $(tvar[1])")
end
