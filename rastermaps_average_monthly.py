# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 15:17:05 2024

@author: Pokhr002
"""
import os
import numpy as np
import rasterio
import pandas as pd
import geopandas as gpd
from rasterio.mask import mask

#%%

input_dir = 'C:/Users/pokhr002/OneDrive - Universiteit Utrecht/03Model/04_final_calib/results/'
forcing_dir = 'C:\\SPHY_input\\forcing\\'
output_dir = 'C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\output_data\\'
output_dir_maps = 'C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\output_data\\maps\\'
catchment_shapefile = "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\01GIS\\Main_map\\Karnali_basin.shp"  
catchment_data = gpd.read_file(catchment_shapefile)
basin_area_m2 = catchment_data.geometry.area.sum()
print(basin_area_m2)


basin_area_km2 = catchment_data.geometry.area.sum() / 1e6
print(f"Basin area: {basin_area_km2} square kilometers")

# Create a date range from '1991-01-01' to '2022-12-31'
dates = pd.date_range(start='1991-01-01', end='2022-12-31')

prefix = "PrecM"

# Generate the filenames in the specified format
filenames = []
for i in range(12):
    for j in range(1, 1001):
        if j < 1000:
            filenames.append(f'{prefix}{i:03d}.{j:03d}')
        else:
            i += 1
            filenames.append(f'{prefix}{i:03d}.000')

# Combine dates and filenames as per the length of dates
combined_data = list(zip(dates, filenames[:len(dates)]))

# Create a DataFrame from the combined data
df = pd.DataFrame(combined_data, columns=['Dates', 'Filenames'])

# Filter the DataFrame for the last day of each month
df_monthly = df[df['Dates'].dt.is_month_end].reset_index(drop=True)
df_monthly['Month'] = df_monthly['Dates'].dt.month



#%%

for month, group in df_monthly.groupby('Month'):
    print(f"Processing month: {month}")
    print("Group DataFrame:")
    print(group)
    
    monthly_sum = None
    file_count = 0
    
    for filename in group['Filenames']:
        raster_path = os.path.join(input_dir, filename)
        
        if not os.path.exists(raster_path):
            print(f"File {raster_path} does not exist. Skipping.")
            continue
        
        with rasterio.open(raster_path) as src:
            data = src.read(1)
            if monthly_sum is None:
                # Initialize the monthly_sum array with zeros
                monthly_sum = np.zeros_like(data, dtype=np.float64)
                print(f"Initialized monthly_sum array for month {month} with shape: {monthly_sum.shape}")
                
                
            monthly_sum += data
            file_count += 1
            print(f"Processed file: {filename}, Current file count: {file_count}")
    
    if file_count > 0:
        # Calculate the average by dividing the sum by the number of files
        monthly_avg = monthly_sum / file_count
        
        # Define the output file path for the unclipped average raster
        output_file_unclipped = os.path.join(output_dir_maps, f'{prefix}_long_term_average_M{month:02d}.tif')
        
        # Save the unclipped monthly average raster
        with rasterio.open(raster_path) as src:
            meta = src.meta.copy()
            meta.update({
                "driver": "GTiff",
                "height": monthly_avg.shape[0],
                "width": monthly_avg.shape[1],
                "count": 1,
                "dtype": 'float32'
            })
            
            with rasterio.open(output_file_unclipped, "w", **meta) as dst:
                dst.write(monthly_avg.astype(rasterio.float32), 1)
        
        print(f"Saved unclipped monthly average for month {month:02d} to {output_file_unclipped}")
        
        # Now, clip the raster using the catchment boundary
        output_file_clipped = os.path.join(output_dir_maps, f'_{prefix}_{os.path.basename(catchment_shapefile).split(".")[0]}_long_term_average_M{month:02d}.tif')
        
        with rasterio.open(output_file_unclipped) as src:
            out_image, out_transform = mask(src, catchment_data.geometry, crop=True)
            
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
        
        print(f"Saved clipped monthly average for month {month:02d} to {output_file_clipped}")
        
#%%

# Open the clipped raster file, calculate the prec value for the basin and store that in the df

basins = ['Bheri', 'Seti', 'Karnali']
months = [f'{i:02d}' for i in range(1, 13)]
prefixes = ['PrecM', 'EtaM']

results_list = []
# Loop through each basin and month
for prefix in prefixes:
    for basin in basins:
        for month in months:
            # Construct the filename for the current basin and month
            filename = f'_{prefix}_{basin}_basin_long_term_average_M{month}.tif'
            filepath = os.path.join(output_dir_maps, filename)
                        
            if os.path.exists(filepath):
                # Open the clipped raster file
                with rasterio.open(filepath) as src:
                    # Read the data
                    out_image = src.read(1).astype(np.float64)
                    affine = src.transform
                    
                    # Handle nodata values (if any)
                    nodata_value = src.nodata
                    if nodata_value is not None:
                        valid_data = out_image[out_image != nodata_value]
                    else:
                        valid_data = out_image
                    
                    # Calculate the average precipitation in mm
                    value_mm = np.sum(valid_data)
                    
                    # cell_area = affine[0] * -affine[4]
                    total_precipitation = np.sum(out_image[out_image > 0])  # * cell_area #Check if y ou need to multiply with the cell area
                    print(total_precipitation)
                    
                    # Append the result to the list
                    results_list.append({
                        'Prefix': prefix,
                        'Basin': basin,
                        'Month': month,
                        'value_mm': value_mm
                    })
                    
            #         print(f"Processed {filename}: value = {value_mm} mm")
            # else:
            #     print(f"File {filename} does not exist. Skipping.")

# Convert the results list to a DataFrame
results_df = pd.DataFrame(results_list)

final_df = results_df.pivot(index='Month', columns=['Basin','Prefix'], values='value_mm')

# Rename the columns to have 'prec_' as a prefix
final_df.columns = [f'{col[0]}_{col[1]}' for col in final_df.columns]

# Reset the index to make 'Month' a column again
final_df = final_df.reset_index()

# Display the final DataFrame
# print(final_df)
final_df.to_csv(os.path.join(output_dir, 'final_PE_sum_data.csv'), index=False)



