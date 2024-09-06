# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 18:01:21 2024

@author: Pokhr002
"""

import os
import rasterio
import pandas as pd
from rasterstats import zonal_stats
import geopandas as gpd  # Import geopandas

directory_forcing = 'C:/SPHY_input/forcing/'
directory_basin = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/03Model/01_SPHY_Input001(500x500)/Data/Catchment_area/'
# List all files in the directory
files = os.listdir(directory_forcing)



# basins = ['Karnali_basin.shp', 'Seti_basin.shp', 'Bheri_Basin.shp']

# # Path to the subbasin map
# basin_map = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/03Model/01_SPHY_Input001(500x500)/Data/Catchment_area/Karnali_basin.shp'

# # Read the shapefile using geopandas
# basin_gdf = gpd.read_file(basin_map)


#%%

def process_raster_files(directory, basin_dir, files, basin="", prefix=""):
    # Create lists to store data
    average_values = []
    # Read the shapefile using geopandas
    basin_geometry = gpd.read_file(os.path.join(basin_dir, basin + '_basin.shp'))
    
    prefix_files = []

    for file in files:
# Check if the file starts with the specified prefix
        if file.startswith(prefix):
# Add the file to the list of files with the prefix
            prefix_files.append(file)
            print(prefix_files)
            with rasterio.open(os.path.join(directory, file)) as src:
                affine_transform = src.transform                
                # Read the raster data as a numpy array
                raster_data = src.read(1)  # Assuming a single band raster                
                # Compute zonal statistics based on the cropped raster and basin geometry
                stats = zonal_stats(basin_geometry, raster_data, affine=affine_transform, nodata=-999, stats=["mean"])                
                # Extract the mean value from the zonal stats dictionary
                mean_value = round(stats[0]["mean"], 2)
                # Append the zonal statistics to the list along with the filename
                average_values.append({f'{prefix}_mean_value': mean_value, 'Filename': file})

    # Create a DataFrame from the list of dictionaries
    df_average = pd.DataFrame(average_values)

    # Save the DataFrame to a CSV file
    output_filename = os.path.join(directory, f'{basin}_{prefix}_averages.csv')
    df_average.to_csv(output_filename, encoding="utf-8", sep='\t', index=False)

    print(f"The average values for the files starting with 'prec', 'tmin', 'tmax', and 'tavg' are saved to: {output_filename}")

process_raster_files(directory_forcing, directory_basin, files, basin="Karnali", prefix="tavg" )


