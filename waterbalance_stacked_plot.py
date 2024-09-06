# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 23:30:37 2024

@author: Pokhr002
"""

#### THIS IS THE CODE TO HAVE STACKED PLOR FOR THE WATER BALANCE COMPONENTS OF THE MODEL 

##### UPDATE THIS WITH MELT COMPONENT AND THE NON-MELT COMPONENTS


import os
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
from calendar import monthrange

#%% 


data_dir = "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\input_data\\"
output_dir = "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\output_data\\"
fig_dir = "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\06Programming\\01Python\\06_waterbalance\\figures\\"

# Plot properties

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.titlesize': 10,
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 15
})

plot_width_cm = 21 - 2 * 1.6
plot_height_in = 12
plot_width_in = plot_width_cm / 2.54
print(plot_width_in)



#Area for each basins 

basin_shapefiles = {
    "Karnali": "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\01GIS\\Main_map\\Karnali_basin.shp",
    "Seti": "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\01GIS\\Main_map\\Seti_basin.shp",
    "Bheri": "C:\\Users\\pokhr002\\OneDrive - Universiteit Utrecht\\01GIS\\Main_map\\Bheri_basin.shp"
}

area_data = []
# Loop through each basin and calculate the area
for basin_name, shapefile_path in basin_shapefiles.items():
    catchment_data = gpd.read_file(shapefile_path)
    basin_area_m2 = catchment_data.geometry.area.sum()  # Calculate area in square meters
    area_data.append({"Basin": basin_name, "Area_m2": basin_area_m2})
    print(f"{basin_name} area: {basin_area_m2} m²")
area_df = pd.DataFrame(area_data)

basin_areas = {
    "West Seti_260": area_df.loc[area_df['Basin'] == 'Seti', 'Area_m2'].values[0],
    "Bheri_270": area_df.loc[area_df['Basin'] == 'Bheri', 'Area_m2'].values[0],
    "Karnali_280": area_df.loc[area_df['Basin'] == 'Karnali', 'Area_m2'].values[0]
}
#%% 

#### READING DATA ( RUNOFF COMPONENTS, PRECIPITATION AND EVAPOTRANSPIRATION ACTUAL)

# Station number for modeled QTot (adjust as needed)
stn_number = 6

# List of desired output dataframes (modify names as needed)
stations = [
"Seti_259.2",
"West Seti_260",
"Bheri_270",
"Kar_215",
"Kar_240",
"Karnali_280"
]

# Function to read data from file and create dataframe with monthly averages
def read_data(file_prefix, columns):
    file_name = [f for f in os.listdir(data_dir) if f.startswith(file_prefix)][0]  # Assuming only one file matches the pattern
    data = pd.read_table(os.path.join(data_dir, file_name), skiprows=stn_number+3, header=None, delim_whitespace=True)
    
    start_date = datetime(1991, 1, 1)
    end_date = datetime(2022, 12, 31)
    modeltime = pd.date_range(start=start_date, end=end_date, freq='D')
    
    # Rename the first column to "Dates" and assign the modeltime
    data.rename(columns={data.columns[0]: "Dates"}, inplace=True)
    data.iloc[:, 0] = modeltime
    data['Dates'] = pd.to_datetime(data['Dates'])
      
    # Add "Month" and "Year" columns to the data
    data['Month'] = data['Dates'].dt.month
    data = data.drop(columns= ['Dates'])  
    
    long_term_monthly_data = data.groupby('Month').mean().reset_index()
    long_term_monthly_data.columns = ['Month'] + columns
        
    return long_term_monthly_data
    
# Reading data and renaming columns based on stations
runoff_long_term_monthly = read_data("QAllDTS", stations)
snowr_long_term_monthly = read_data("STotDTS", stations)
baser_long_term_monthly = read_data("BTotDTS", stations)
rainr_long_term_monthly = read_data("RTotDTS", stations)
glacr_long_term_monthly = read_data("GTotDTS", stations)

# Display first few rows of the long-term monthly dataframes
print(baser_long_term_monthly.head())
print(rainr_long_term_monthly.head())
print(runoff_long_term_monthly.head())


#%% 

### THIS IS THE CODE TO PLOT THE WATER BALANCE DATA 

days_in_month = {month: monthrange(2022, month)[1] for month in range(1, 13)}
print(days_in_month)
hours_per_day = 24

# Reading the values for the Prec and Eta that is calculated from the raster file
df_Eta_Prec = pd.read_csv(os.path.join(output_dir,"final_PE_data.csv"))
df_Eta_Prec = df_Eta_Prec.rename(columns=lambda col: col.replace('Karnali_', 'Karnali_280_')
                                 .replace('Seti_', 'West Seti_260_')
                                 .replace('Bheri_', 'Bheri_270_'))

station_plot = ['West Seti_260', 'Bheri_270', 'Karnali_280'] 
# fig, axes = plt.subplots(len(station_plot), 1, figsize=(plot_width_in, plot_height_in), sharex=True)
fig, axes = plt.subplots(len(station_plot), 1, figsize=(5.5, 8), sharex=True)

for idx, station in enumerate(station_plot):
    ax = axes[idx]
    
    # Since these df contains data for all the station, selecting for each staton using loop
    runoff_station = runoff_long_term_monthly[['Month', station]].copy()
    baser_station = baser_long_term_monthly[['Month', station]].copy()
    snowr_station = snowr_long_term_monthly[['Month', station]].copy()
    glacr_station = glacr_long_term_monthly[['Month', station]].copy()
    rainr_station = rainr_long_term_monthly[['Month', station]].copy()
    eta_prec_station = df_Eta_Prec.filter(regex=f'^{station}.*(_EtaM|_PrecM)$')
    eta_prec_station = pd.concat([df_Eta_Prec[['Month']], eta_prec_station], axis=1)
    basin_area = basin_areas[station]
    
    # Merging dataframes for the selected station
    combined_data = baser_station.merge(snowr_station, on='Month', suffixes=('_Base', '_Snow'), how='inner')
    combined_data = combined_data.merge(rainr_station, on='Month', suffixes=('', '_Rain'), how='inner')   
    combined_data = combined_data.merge(glacr_station, on='Month', suffixes=('', '_Glacier'), how='inner') 
    combined_data = combined_data.merge(runoff_station, on='Month', suffixes=('', '_Runoff'), how='inner') 
    combined_data = combined_data.merge(eta_prec_station, on='Month', how='inner') 

    # Renaming columns
    combined_data.rename(columns={
    f'{station}_Base': 'Baser',
    f'{station}_Snow': 'Snowr',
    station: 'Rainr',
    f'{station}_Glacier': 'Glacr',
    f'{station}_Runoff': 'Runoff',
    f'{station}_EtaM': 'Eta_mm',
    f'{station}_PrecM': 'Prec_mm', 
}, inplace=True)
    
  
    # Invert the sign of the selected columns because they are negative in the water balance equation
    columns_to_convert_sign = ['Baser', 'Snowr', 'Rainr', 'Glacr', 'Runoff', 'Eta_mm']
    for col in columns_to_convert_sign:
        if col in combined_data.columns:
            combined_data[col] = combined_data[col] * -1
        else:
            print(f"Warning: Column {col} not found in combined_data for {station}")
  
    # The followong columns are then changed to unit mm so that they can be comparable with prec and eta    
    columns_to_convert_mm = ['Baser', 'Snowr', 'Rainr', 'Glacr', 'Runoff']    
    for col in columns_to_convert_mm:
        combined_data[f'{col}_mm'] = combined_data.apply(
            lambda row: row[col] / basin_area * days_in_month[row['Month']] * hours_per_day * 3600 * 1000, axis=1
        )              
        
        # Calculate storage to check if the storage is pos or neg in each month
    if 'Prec_mm' in combined_data.columns:
        combined_data['Storage'] = combined_data['Prec_mm'] + combined_data['Eta_mm'] + combined_data['Runoff_mm']
    else:
        print(f"ERROR: 'Prec_mm' is missing when calculating storage for {station}!")
     
    print(combined_data)
   
    # Define pastel colors
    set_colors = sns.color_palette("colorblind")
    set_colors1 = sns.color_palette("Paired")
    
    # Calculate the bottom for storage
    bottom_storage = combined_data.apply(
        lambda row: row['Prec_mm'] if row['Storage'] >= 0 else row['Eta_mm'] + row['Rainr_mm'] + row['Baser_mm'] + row['Snowr_mm'],
        axis=1
    )

    # Plotting the stacked bar plot
    ax.bar(combined_data['Month'], combined_data['Prec_mm'], label='Precipitation', color=set_colors[9], alpha = 0.8)
    ax.bar(combined_data['Month'], combined_data['Eta_mm'], label='Evapotranspiration', color=set_colors1[2])
    ax.bar(combined_data['Month'], combined_data['Rainr_mm'], bottom=combined_data['Eta_mm']
           , label='Rain runoff', color=set_colors1[8])
    ax.bar(combined_data['Month'], combined_data['Baser_mm'], bottom=combined_data['Eta_mm'] + combined_data['Rainr_mm']
           , label='Baseflow runoff', color=set_colors[5])
    ax.bar(combined_data['Month'], combined_data['Snowr_mm'], bottom=combined_data['Eta_mm'] + combined_data['Rainr_mm'] + combined_data['Baser_mm']
           , label='Snow runoff', color=set_colors[7])
    ax.axhline(0, color='black', linewidth=1, alpha = 0.75) # this is for pos and neg side to clearly see the balance

    # Plotting the storage as a line plot    
    ax.plot(combined_data['Month'], combined_data['Storage'], color=set_colors1[9], marker='o', markersize= 3 ,  linewidth=1 , label ='Storage change')

    # Set title for each sub-plot        
    subplot_title = f"({chr(97+ idx)})"  # chr(97) is 'a', chr(98) is 'b', and so on.
    ax.set_title(subplot_title,  loc='right')
                     
    # ax.set_title(station.split('_')[0], fontsize=12, fontweight='bold', backgroundcolor=set_colors[8], loc='right', 
    #               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # this was to have the name of the sation on the right of the figure
    
    ax.set_ylabel('Value [mm]')
    # ax.set_xticks(combined_data['Month'])
    # ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize=10)
    
    if idx == len(station_plot) - 1:
        # Set x-ticks and labels only for the last subplot
        ax.set_xticks(combined_data['Month'])
        ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    else:
        # Remove x-ticks and labels for all other subplots
        ax.tick_params(labelbottom=False, bottom=False)  
        
    ax.set_ylim(-600, 600)
    ax.grid(True, linestyle='--', linewidth=0.5)


# Adjust space between sub-plots
plt.subplots_adjust(hspace=0.0)

# Add legend at the bottom
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True)
        
# n_stations = len(station_plot)
# spacing = 1 / n_stations
# for idx, station in enumerate(station_plot):
#     # Calculate the vertical position to ensure the text spans the entire figure length
#     y_position = (1 - spacing / 2) - idx * spacing

#     fig.text(0.98, y_position, station.split('_')[0], fontsize=12, fontweight='bold', 
#              rotation=90, va='center', ha='center',
#              bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

# Adjust layout to make space for the text and legend
plt.tight_layout(rect=[0, 0.06, 0.97, 0.9])
plt.savefig(os.path.join(fig_dir, 'combined_water_balance.png'), dpi=300, bbox_inches='tight')


#bbbox to anchor and set the coordinates