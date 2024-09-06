# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 15:17:05 2024

@author: Pokhr002
"""
import os
import shutil
import numpy as np
import rasterio
import pandas as pd
import geopandas as gpd
from rasterio.mask import mask

#%%

## I have output from SPHY for Prec, ETa and runoff components in Y maps that is yearly sum. I am using this code to get long term averages.

input_dir = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/03Model/04_final_calib/results/'
output_dir = 'C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\output_data\\'
output_dir_maps = 'C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\output_data\\maps\\yearly\\'

catchment_shapefile = "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\01GIS\\Main_map\\Karnali_basin.shp"  
catchment_data = gpd.read_file(catchment_shapefile)
basin_area_m2 = catchment_data.geometry.area.sum()
print(basin_area_m2)


basin_area_km2 = catchment_data.geometry.area.sum() / 1e6
print(f"Basin area: {basin_area_km2} square kilometers")

# Create a date range from '1991-01-01' to '2022-12-31'
dates = pd.date_range(start='1991-01-01', end='2022-12-31')

prefix = "BaserY"

# Generate the filenames in the specified format
filenames = []
for i in range(12):
    for j in range(1, 1001):
        if j < 1000:
            filenames.append(f'{prefix}{i:02d}.{j:03d}')
        else:
            i += 1
            filenames.append(f'{prefix}{i:02d}.000')

# Combine dates and filenames as per the length of dates
combined_data = list(zip(dates, filenames[:len(dates)]))

# Create a DataFrame from the combined data
df = pd.DataFrame(combined_data, columns=['Dates', 'Filenames'])
df['year'] = df['Dates'].dt.year


df_yearly = df[df['Dates'].dt.is_year_end].reset_index(drop=True)
df_yearly['Year'] = df_yearly['Dates'].dt.year


#%%

## This code exists because I changes the file name name to Y when I printed in the 
# reporting table so I changed back to prec to avoid future confusion

# for index, row in df_yearly.iterrows():
#     original_filename = row['Filenames']
    
#     # Replace the first 'Y' with 'Prec' in the filename
#     new_filename = original_filename.replace('Y', 'Prec', 1)
    
#     # Construct the original file path
#     original_filepath = os.path.join(input_dir, original_filename)
    
#     # Construct the new file path
#     new_filepath = os.path.join(output_dir_maps, new_filename)
    
#     # Check if the file exists in the input directory
#     if os.path.exists(original_filepath):
#         # Copy the file to the new location with the new name
#         shutil.copy(original_filepath, new_filepath)
#         print(f"Renamed {original_filename} to {new_filename}")
#     else:
#         print(f"File {original_filepath} does not exist. Skipping.")

#%%

long_term_avg_sum = None
year_count = 0

# Group by the 'year' and process each group (where each file is already a yearly sum map)
  
for filename in df_yearly['Filenames']:
    raster_path = os.path.join(input_dir, filename)
 
    if not os.path.exists(raster_path):
        print(f"File {raster_path} does not exist. Skipping.")
        continue
    
    with rasterio.open(raster_path) as src:
        data = src.read(1)  # Read the raster data (yearly sum map)
        
        if long_term_avg_sum is None:
            # Initialize the long_term_avg_sum array with zeros
            long_term_avg_sum = np.zeros_like(data, dtype=np.float64)
            print(f"Initialized long_term_avg_sum array for year with shape: {long_term_avg_sum.shape}")
        
        # Accumulate the yearly sum map data
        long_term_avg_sum += data
        year_count += 1
        print(f"Processed file: {filename}, Current file count: {year_count}")
    
# After processing all yearly maps, calculate the long-term average
if year_count > 0:
    long_term_avg = long_term_avg_sum / year_count  # Average over the number of years
    print("Calculated long-term average.")

    # Clip the long-term average raster to the catchment boundary
    output_file_clipped = os.path.join(output_dir_maps, f'{prefix}_long_term_average_clipped.tif')
    
    with rasterio.open(raster_path) as src:
        # Create a temporary raster to clip
        meta = src.meta.copy()
        meta.update({
            "driver": "GTiff",
            "height": long_term_avg.shape[0],
            "width": long_term_avg.shape[1],
            "count": 1,
            "dtype": 'float32'
        })
        
        with rasterio.MemoryFile() as memfile:
            with memfile.open(**meta) as temp_dst:
                temp_dst.write(long_term_avg.astype(rasterio.float32), 1)
                
                # Now clip the temporary raster
                out_image, out_transform = mask(temp_dst, catchment_data.geometry, crop=True, filled=True, nodata=src.nodata)
                
                # Convert to 2D if necessary and to float64 for accurate calculations
                out_image = out_image[0].astype(np.float64)
                
                # Update the metadata for the clipped raster
                out_meta = src.meta.copy()
                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[0],
                    "width": out_image.shape[1],
                    "transform": out_transform,
                    "count": 1,
                    "dtype": 'float32'
                })
                
                # Write the clipped raster to a new file
                with rasterio.open(output_file_clipped, "w", **out_meta) as dst:
                    dst.write(out_image, 1)
    
    print(f"Saved clipped long-term average to {output_file_clipped}")
else:
    print("No data processed. No output generated.")
