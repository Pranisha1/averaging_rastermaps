# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:14:28 2024

@author: Pokhr002

############ It basically reads the maps, calculates the average from all grid cells and store in df with the date ###########
##### This can be later used to manipulate and get averages for the month and year and season so on ############
"""

import os
import rasterio
import pandas as pd
from datetime import datetime
from rasterstats import zonal_stats
import geopandas as gpd  # Import geopandas
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# Directory containing raster files
directory_forcing = 'C:/SPHY_input/forcing/'
dir_model = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/03Model/04_final_calib/results/'
basin_map = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/03Model/01_SPHY_Input001(500x500)/Data/Catchment_area/'
output_dir = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/06Programming/01Python/02droughts/processed_data/for_SPEI/'
fig_dir = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/06Programming/01Python/02droughts/figures/'

basin = 'Karnali'

#%% 

start_date = datetime(1991, 1, 1)
end_date = datetime(2022, 12, 31)
# Create date range with custom date format (Y-m-d)
date_range = pd.date_range(start=start_date, end=end_date, freq='D').strftime('%Y-%m-%d')

# Read the shapefile using geopandas
basin_gdf = gpd.read_file(basin_map + basin + "_basin.shp")

# List all files in the directory
files = os.listdir(directory_forcing)
model_files = os.listdir(dir_model)


#%%
######################################Calculating the Forcing average ############################################################


parameter = "prec"

# Filter files that start with "prec"
parameter_files = [f for f in files if f.startswith(parameter)]

print (parameter_files)

# Create lists to store data
average_values = []


for file in parameter_files :
    # Open the raster file
    with rasterio.open(os.path.join(directory_forcing, file)) as src:
        
        affine_transform = src.transform
        
        # Read the raster data as a numpy array
        raster_data = src.read(1)  # Assuming a single band raster
        
        # Compute zonal statistics based on the raster that is within basin geometry
        stats = zonal_stats(basin_gdf, raster_data, affine=affine_transform, nodata=-999, stats=["mean"])
        
        # Extract the mean value from the zonal stats dictionary
        mean_value = round(stats[0]["mean"], 2)

        # Append the zonal statistics to the list along with the filename
        average_values.append({f"Mean_{parameter}": mean_value, 'Filename': file})
        
# Create a DataFrame from the list of dictionaries
df_average = pd.DataFrame(average_values)
df_dates = pd.DataFrame({'Dates': date_range})

# Merge the date range DataFrame with the average precipitation DataFrame
average_series = pd.merge(df_dates, df_average, left_index=True, right_index=True)

# Save the merged DataFrame to a CSV file
average_series.to_csv(output_dir + f"{basin}_average_{parameter}.csv", encoding="utf-8", sep='\t', index=False)

print(f"The average {parameter} series for the {basin} basin area is produced.")


#%%

######################################## Model result averages #####################################################################

# Filter files that 
m_parameter = "ETp"
not_parameter = ["SnowS", "SnowDTS", "SnowR_GLAC", "SnowWatStore_GLAC", "Snow_GLAC", "ETpDTS"]

# Filter files that start with "Snow" but not with any prefix in not_parameter
model_parameter_files = [f for f in model_files if f.startswith(m_parameter) and not any(f.startswith(prefix) for prefix in not_parameter)]

print(model_parameter_files)


# Create lists to store data
average_values = []


for file in model_parameter_files:
    # Open the raster file
    with rasterio.open(os.path.join(dir_model, file)) as src:
        
        # Assuming `src` is the rasterio dataset object opened with your raster file
        affine_transform = src.transform
        
        # Read the raster data as a numpy array
        raster_data = src.read(1)  # Assuming a single band raster
        
        # Compute zonal statistics based on the cropped raster and basin geometry
        stats = zonal_stats(basin_gdf, raster_data, affine=affine_transform, nodata=-999, stats=["mean"])
        
        # Extract the mean value from the zonal stats dictionary
        mean_value = round(stats[0]["mean"], 2)

        # Append the zonal statistics to the list along with the filename
        average_values.append({f"Mean_{m_parameter}": mean_value, 'Filename': file})
        
# Create a DataFrame from the list of dictionaries
df_average = pd.DataFrame(average_values)
df_dates = pd.DataFrame({'Dates': date_range})

# Merge the date range DataFrame with the average precipitation DataFrame
average_series = pd.merge(df_dates, df_average, left_index=True, right_index=True)

# Save the merged DataFrame to a CSV file
average_series.to_csv(output_dir + f"{basin}_average_{m_parameter}.csv", encoding="utf-8", sep='\t', index=False)

print(f"The average {m_parameter} series for the {basin} basin area is produced.")



#%%   ## this is a code to plot snow map for the 


snow_df = pd.read_csv(f"{output_dir}karnali/Karnali_average_snow.csv", delimiter='\t', parse_dates=[0], usecols=lambda x: x != 'Filename')

# Extract year and month from the 'date' column
snow_df['year'] = snow_df['Dates'].dt.year
snow_df['month'] = snow_df['Dates'].dt.month

# Group by year and month, and calculate the sum of snowfall for each group
monthly_sum = snow_df.groupby(['year', 'month'])['Mean_Snow'].sum().reset_index()

# Merge year and month as a datetime object
monthly_sum['date'] = pd.to_datetime(monthly_sum['year'].astype(str) + '-' + monthly_sum['month'].astype(str), format='%Y-%m')

# Calculate the annual average snowfall
longterm_monthly_average = monthly_sum.groupby('month')['Mean_Snow'].mean().reset_index()
longterm_monthly_average.rename(columns={'Mean_Snow': 'Long-term monthly average snowfall'}, inplace=True)

# Set up subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# Plot the monthly snowfall sum
ax1.plot(monthly_sum['date'], monthly_sum['Mean_Snow'], linestyle='-', color='b')
ax1.set_xlabel('Date (m-y)')
ax1.set_ylabel('Snowfall (mm)')
ax1.set_title('Monthly Snowfall Sum')
ax1.xaxis.set_major_locator(mdates.YearLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%y'))
plt.setp(ax1.get_xticklabels(), rotation=90)  # Rotate x-axis tick labels by 45 degrees


# Plot the long-term monthly average snowfall
ax2.plot(longterm_monthly_average ['month'], longterm_monthly_average ['Long-term monthly average snowfall'], marker='o', linestyle='-', color='r')
ax2.set_xlabel('Year')
ax2.set_ylabel('Snowfall')
ax2.set_title('Long-term Monthly Average Snowfall')

# Adjust layout
plt.tight_layout()

# Save the plot as an image in the output directory
output_path = os.path.join(fig_dir, 'snowfall_plot.png')
plt.savefig(output_path)

# Show the plot
plt.show()



#%%

# ----------------------------------------------------------------------------------------------------------------------

# import numpy as np
# from rasterstats.io import Raster
# from rasterio.mask import mask
# import matplotlib.pyplot as plt
# from rasterio.plot import show

# idea to use mask and view in python

# # Iterate over the filtered files
# for file in prec_files:
#     # Open the raster file
#     with rasterio.open(os.path.join(directory, file)) as src:
        
#         # Assuming `src` is the rasterio dataset object opened with your raster file
#         affine_transform = src.transform
        
#         # Read the raster data as a numpy array
#         raster_data = src.read(1)  # Assuming a single band raster
        
#         # # Crop the raster using the basin geometry
#         # raster_data_crop, out_transform = mask(src, basin_gdf.geometry, crop=True)
        
#         # # Display the cropped raster
#         # plt.figure(figsize=(8, 6))
#         # show(raster_data_crop[0], transform=out_transform)
#         # plt.title('Cropped Raster Data')
#         # plt.xlabel('X')
#         # plt.ylabel('Y')
#         # # plt.colorbar(label='Pixel Value')

#         # Calculate the average of the cropped raster values
#         # raster_average = round(np.nanmean(raster_data_crop), 2)  # NaN values are ignored in the calculation
        
#         # print(raster_average)

#         # Compute zonal statistics based on the cropped raster and basin geometry
#         stats = zonal_stats(basin_gdf, raster_data, stats=["mean"])

#         # Append the average precipitation value and zonal statistics to the list
#         average_values.append({'Average_precip': stats, 'Filename': file})

