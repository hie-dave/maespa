#!/usr/bin/env julia

using NCDatasets
using PlotlyJS

site = "CumberlandPlain"
daily_file = "data/silo/processed/$site.nc"
hourly_file = "data/silo/hourly/$site.nc"

NCDataset(daily_file) do nc_daily
    time_daily = nc_daily["time"][:]
    tmax = nc_daily["tasmax"][:]
    tmin = nc_daily["tasmin"][:]

    NCDataset(hourly_file) do nc_hourly
        time_hourly = nc_hourly["time"][:]
        tas = nc_hourly["tas"][:]

        # Plot daily min/max temperature alongside hourly temperature on the
        # same plot.
        title = "$site air temperature"
        xlab = "Date"
        ylab = "Air Temperature (°C)"
        plot([
            scatter(x=time_daily, y=tmax, name="Daily max"),
            scatter(x=time_daily, y=tmin, name="Daily min"),
            scatter(x=time_hourly, y=tas, name="Hourly"),
        ], title=title, xlab=xlab, ylab=ylab)
    end
end
