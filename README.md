### averaging_rastermaps
This repository contains all the codes prepared to get the average value of the raster maps.

* forcing_avg_series_basin.py - The grid values are averaged over the basin and if we want long-term monthly, the monthly sum is calculated for each year and then averaged over the long term.
* average_series_basin_usingfunction.py - Same thing using the function to make it easy.
* rastermaps_average_monthly.py and annual.py take the monthly and annual maps generated by SPHY, which are each year's monthly and annual sum. These maps are then summed and divided by the number of maps to get an average


