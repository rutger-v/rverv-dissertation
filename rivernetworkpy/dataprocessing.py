from .analysisvisualisation import remove_directionality, get_strikes_and_dip_orientations_from_dip_directions, find_diff_segm_or_and_mean_strike_or_dip_or, fix_list

import pandas as pd
import numpy as np
import geopandas as gpd
import utm
from pyproj import CRS
from shapely.geometry import Point, Polygon, LineString, shape
from shapely.ops import unary_union
import json
from scipy.stats import circmean


#### CHI DATA MAP MERGING #### 

def merge_segment_bearings_with_chi_data_map_for_basin(nickname, part):
    """
    MERGE SEGMENT BEARINGS WITH CHI DATA MAPS
    - to get rid of segments outside the selected basin
    - to add the basin_key and source_key columns to the segments_bearings.csv
    """
    Dataset_prefix = nickname
    DataDirectory = "./Basins/"+part+"/"+nickname+"/"
    
    # load the csv's in as dataframes
    chi_data_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_chi_data_map.csv")
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings.csv")
    
    # use pandas' merge function
    merged = pd.merge(chi_data_df, segment_bearings_df, how = "right", on = ["latitude", "longitude"])
    
    from IPython.display import display
    display(merged.head())
 
    # drop the rows with no basin key
    merged['basin_key'].replace('', np.nan, inplace=True)
    merged.dropna(subset=['basin_key'], inplace=True)

    # drop the information from the chi data map that we don't want
    merged.drop(columns=['chi', 'elevation', 'flow_distance', 'drainage_area'], inplace=True)
    
    from IPython.display import display
    display(merged.head())

    # overwrite the segment_bearings csv
    merged.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings.csv", index = False)  
    
def add_basin_key_to_junction_angles_csv(nickname, part):
    """
    Similar to segment bearings merging above; goal is to add basin key to junction angles csv
    """
    Dataset_prefix = nickname
    DataDirectory = "./Basins/"+part+"/"+nickname+"/"
    
    # load the csv's in as dataframes
    chi_data_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_chi_data_map.csv")
    junction_angles_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_FULL_Basins_JAngles.csv")
    
    # use panda's merge function
    merged = pd.merge(chi_data_df, junction_angles_df, how = "right", on = ["latitude", "longitude"])
    
    # drop the information from the chi data map that we don't want
    merged.drop(columns=['chi', 'elevation', 'flow_distance', 'drainage_area'], inplace=True)
    
    from IPython.display import display
    display(merged.head())
    
    # overwrite the junction angles csv
    merged.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_FULL_Basins_JAngles.csv", index = False)      
    
#### GETTING RID OF SEGMENT BEARINGS / JUNCTION ANGLES IN BASINS WITH A ORDER LOWER THAN A CERTAIN THRESHOLD ########

def remove_low_order_basins_from_segment_bearings_csv(nickname, part, threshold = 5):
    """
    This function removes segment bearings from basins whose maximum stream order is lower than a certain threshold.
    
    10/01/2024  
    
    Returns: 
        A NEW segment bearings csv
    """
    Dataset_prefix = nickname
    DataDirectory = "./Basins/"+part+"/"+nickname+"/"
    
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings.csv")
    
    # group by basin key, and compute the maximum junction stream order for each basin key
    temp_df = segment_bearings_df.groupby("basin_key")["junction_stream_order"].max().reset_index()
    
    # merge temp_df to the the segment_bearings_df
    segment_bearings_df = pd.merge(segment_bearings_df, temp_df, on="basin_key", suffixes = ("", "_max"))

    # filter out rows with a maximum stream order smaller than the threshold (default 5)
    segment_bearings_big_basins_df = segment_bearings_df[segment_bearings_df["junction_stream_order_max"] >= threshold]
    
    # write csv
    segment_bearings_big_basins_df.to_csv(DataDirectory+Dataset_prefix+"_COP30_UTM_segment_bearings_filtered.csv", index = False)

def remove_low_order_basins_from_junction_angles_csv(nickname, part, threshold = 5):
    """
    This function removes junction angles from basins whose maximum stream order is lower than a certain threshold.
    
    10/01/2024  
    
    Returns: 
        A NEW junction angles csv
    """
    Dataset_prefix = nickname
    DataDirectory = "./Basins/"+part+"/"+nickname+"/"
    
    junction_angles_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_FULL_Basins_JAngles.csv")
    
    # group by basin key, and compute the maximum junction stream order for each basin key
    temp_df = junction_angles_df.groupby("basin_key")["junction_stream_order"].max().reset_index()
    
    # merge temp_df to the the segment_bearings_df
    junction_angles_df = pd.merge(junction_angles_df, temp_df, on="basin_key", suffixes = ("", "_max"))

    # filter out rows with a maximum stream order smaller than the threshold (default 5)
    junction_angles_big_basins_df = junction_angles_df[junction_angles_df["junction_stream_order_max"] >= threshold]
    
    # write csv
    junction_angles_big_basins_df.to_csv(DataDirectory+Dataset_prefix+"_COP30_UTM_FULL_Basins_JAngles_filtered.csv", index = False)

    
    
#### REMOVING LAKE SEGMENT BEARINGS AND JUNCTION ANGLES (PART 2 ONLY, INCORPORATED IN BEDDING MEASUREMENTS FUNCTION [see below] IN PART 1) ######   
def remove_lake_segment_bearings_and_junction_angles(nickname, filepath_to_hydrolakes = "./hydrolakes/HydroLAKES_VRP_AP/HydroLAKES_VRP_AP.shp"):
    """
    This function creates new segment bearings and junction angle csvs with only the junctions and segment bearings with segments fully outside lakes.
    
    10/01/2024
    """
    # Open chi data map and segment bearings csv and put in geodataframes
    print("Starting to remove lake segment bearings and junction angles of  " + nickname)
    Dataset_prefix = nickname
    DataDirectory = "./Basins/PART2/"+nickname+"/"
    
    # load the chi data map and segment bearings in as dataframes. we need to merge the two of them because chi data map has the entire river network and segment bearings has the junctions
    chi_data_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_chi_data_map.csv")
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings_filtered.csv")
    
    # add a Geometry column containing shapely points. Generate from longitude and latitude. (The construction df.apply(lambda row: f()) allows for applying one function to every row in a pd Dataframe)
    chi_data_df['geometry_points'] = chi_data_df.apply(lambda row: Point(row.longitude, row.latitude), axis = 1)
    segment_bearings_df['geometry_points'] = segment_bearings_df.apply(lambda row: Point(row.longitude, row.latitude), axis = 1)
    
    # turn both dataframes into a geodataframe 
    chi_data_gdf = gpd.GeoDataFrame(chi_data_df, geometry = 'geometry_points')
    segment_bearings_gdf = gpd.GeoDataFrame(segment_bearings_df, geometry = 'geometry_points')
    
    # get segment lines, store as separate gdf
    segment_lines_gdf = turn_streams_into_lines(chi_data_gdf, segment_bearings_gdf)
    
    # write the lines to GeoJSON so we can manually check whether turning the nodes into lines worked well
    segment_lines_gdf['lines'].to_file(DataDirectory+ Dataset_prefix + "_segment_lines.geojson", driver = 'GeoJSON', crs = "EPSG:4326")
    print("I saved you a GeoJSON file with the river network as lines instead of points. I deserve a sticker.")
    
    # remove all lines that are (partially) inside a HydroLAKE
    segment_lines_gdf = find_lines_fully_outside_lakes(segment_lines_gdf, filepath_to_hydrolakes)
    
    # write the lines to GeoJSON so we can manually check whether we successfully got rid of (partially) inside a hydroLAKE
    segment_lines_gdf['lines'].to_file(DataDirectory+ Dataset_prefix + "_segment_lines_nolakes.geojson", driver = 'GeoJSON', crs = "EPSG:4326")
    
    # all segment bearings and junction angle measurements have a unique junction number. 
    # so, load in the segment bearings and junction angle csv's as dataframes, and remove all datapoints with a junction number that is NOT in segment lines_gdf
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings_filtered.csv")
    junction_angles_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_FULL_Basins_JAngles_filtered.csv")
    
    valid_junction_numbers = segment_lines_gdf["junction_number"].unique()
    
    segment_bearings_df = segment_bearings_df[segment_bearings_df["junction_number"].isin(valid_junction_numbers)]
    junction_angles_df = junction_angles_df[junction_angles_df["junction_number"].isin(valid_junction_numbers)]    
    
    # save result
    segment_bearings_df.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings_filtered_nolakes.csv", index = False)
    junction_angles_df.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_FULL_Basins_JAngles_filtered_nolakes.csv", index = False)
    
#### ADDING BEDDING MEASUREMENTS TO SEGMENT BEARINGS DATASET #### 
def add_bedding_measurements_to_segments(nickname, buffer_radius = 1000, filepath_to_beddingdata = "./wvbed/WV_BED_VRP.shp", filepath_to_beddingquads = "./wvbed/WV_QUADS_VRP.shp", filepath_to_hydrolakes = "./hydrolakes/HydroLAKES_VRP_AP/HydroLAKES_VRP_AP.shp"):
    """
    This function adds a column of bedding measurements to the segment bearings csv.
    
    Arguments:
        nickname (str): the nickname of the basin
        buffer_radius (int): the radius in meters defining the buffer around the river segment in which to find bedding datapoints
        filepath_to_beddingdata (str): is obvious
    
    Returns:
        a NEW segment_bearings file, with the name _segment_bearings_bedding_measurements.csv
    """
    # Open chi data map and segment bearings csv and put in geodataframes
    print("Starting the bedding measurements procedure of " + nickname + ". Radius is set to " + str(buffer_radius) + "m")
    Dataset_prefix = nickname
    DataDirectory = "./Basins/PART1/"+nickname+"/"
    
    # load the chi data map and segment bearings in as dataframes. we need to merge the two of them because chi data map has the entire river network and segment bearings has the junctions
    chi_data_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_chi_data_map.csv")
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings.csv")
    
    # add a Geometry column containing shapely points. Generate from longitude and latitude. (The construction df.apply(lambda row: f()) allows for applying one function to every row in a pd Dataframe)
    chi_data_df['geometry_points'] = chi_data_df.apply(lambda row: Point(row.longitude, row.latitude), axis = 1)
    segment_bearings_df['geometry_points'] = segment_bearings_df.apply(lambda row: Point(row.longitude, row.latitude), axis = 1)
    
    # turn both dataframes into a geodataframe 
    chi_data_gdf = gpd.GeoDataFrame(chi_data_df, geometry = 'geometry_points')
    segment_bearings_gdf = gpd.GeoDataFrame(segment_bearings_df, geometry = 'geometry_points')
    
    # the chi_data_gdf includes all the basins originally extracted. I only want those that occur in the segment_bearings_gdf
    valid_basin_keys = segment_bearings_gdf['basin_key'].unique()
    chi_data_gdf = chi_data_gdf[chi_data_gdf['basin_key'].isin(valid_basin_keys)]
    
    # get segment lines, store as separate gdf
    segment_lines_gdf = turn_streams_into_lines(chi_data_gdf, segment_bearings_gdf)
    
    # write the lines to GeoJSON so we can manually check whether turning the nodes into lines worked well
    segment_lines_gdf['lines'].to_file(DataDirectory+ Dataset_prefix + "_segment_lines.geojson", driver = 'GeoJSON', crs = "EPSG:4326")
    print("I saved you a GeoJSON file with the river network as lines instead of points. I deserve a sticker.")
        
    # remove all lines that are not fully contained wihtin the WV_QUADS_VRP.shp polygon. 
    segment_lines_gdf = find_lines_fully_inside_quads(segment_lines_gdf, filepath_to_beddingquads)    
    
    # remove all lines that are (partially) inside a HydroLAKE
    segment_lines_gdf = find_lines_fully_outside_lakes(segment_lines_gdf, filepath_to_hydrolakes)
    
    # write the lines to GeoJSON so we can manually check whether we successfully got rid of the segments (partially) outside the quads and (partially) inside a hydroLAKE
    segment_lines_gdf['lines'].to_file(DataDirectory+ Dataset_prefix + "_segment_lines_inquads_nolakes.geojson", driver = 'GeoJSON', crs = "EPSG:4326")
    
    # buffer the remaining lines
    segment_lines_buffered_gdf = buffer_lines(segment_lines_gdf, chi_data_gdf, buffer_radius)
    
    # find the bedding measurements within that buffer polygon (this results in the gdf being extended with columns for dips, azimuths, dip_directions and strikes containing any measurements found)
    gdf, df = find_bedding_measurements_within_polygons(segment_lines_buffered_gdf, filepath_to_beddingdata)
    
    # write gdf and df to CSV - if something goes wrong in the following lines we can just import the csv w/o having to rerun the previous function
    gdf.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_bedding_measurements.csv", index = False)
    df.to_csv(DataDirectory + Dataset_prefix + "_unique_bedding_measurements.csv", index = False)
    
def add_bedding_measurements_to_segment_bearings(nickname):
    """
    This function calculates some stats based on bedding data and adds these to the segment_bearings dataset. The function `add_bedding_measurements_to_segments` should be executed first.
    
    Arguments:
        nickname (str): the name of the basin
    
    Returns:
        a NEW segment_bearings file, with the name _segment_bearings_bedding_stats.csv
    """

    # open file
    Dataset_prefix = nickname
    DataDirectory = "./Basins/PART1/"+nickname+"/"
    gdf1 = gpd.read_file(DataDirectory + Dataset_prefix + "_COP30_UTM_bedding_measurements.csv")
    
    # fix the lists (i.e. after opening the csv the list-like strings should become actual python-lists)
    gdf1['azimuths'] = gdf1.apply(lambda row: fix_list(row.azimuths), axis = 1)
    gdf1['dips'] = gdf1.apply(lambda row: fix_list(row.dips), axis = 1)
    gdf1['dip_directions'] = gdf1.apply(lambda row: fix_list(row.dip_directions), axis = 1)
    
    # remove the columns with the line and the buffer
    gdf1.drop(columns = ["lines", "lines_buffer"], inplace = True)      
                    
    # remove the rows with NO bedding measurements (i.e. rows where dips/azimuths/dip_directions lists are empty)
    gdf = gdf1[gdf1.dips.map(tuple) != tuple([])]
    print("I removed all rows without bedding measurements")
        
    # derive strike from dip_direction column
    print("Deriving strikes from dip directions now")
    gdf[['strikes', 'dip_orientations']] = gdf.apply(lambda row: get_strikes_and_dip_orientations_from_dip_directions(row.dip_directions), axis = 1, result_type = "expand")
    
    """
    # calculate n, mean direction and circular variance of strikes and dip orienations using the calculate_directional_stats function and put it in new columns
    print("Calculating some circular statistics on your strikes")
    gdf[['n_strike', 'mean_strike', 'circ_var_strike']] = gdf.apply(lambda row: calculate_directional_stats(row.strikes), axis = 1, result_type = "expand")
    gdf[['n_dip_orientation', 'mean_dip_orientation', 'circ_var_dip_orientations']] = gdf.apply(lambda row: calculate_directional_stats(row.dip_orientations), axis = 1, result_type = "expand")
    """
    gdf["mean_strike"] = gdf.apply(lambda row: np.rad2deg(circmean(np.deg2rad(np.array(row.strikes)), high = np.pi, low = 0)), axis=1, result_type="expand")
    gdf["mean_dip_orientation"] = gdf.apply(lambda row: np.rad2deg(circmean(np.deg2rad(np.array(row.dip_orientations)), high = np.pi, low = 0)), axis=1, result_type="expand")
    
    
    # merge the segment_bearings dataset with the gdf based on the column "junction_number"
    print("Now merging your segment_bearings with the segments and bedding dataset we just made")
    
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings.csv")
    gdf['junction_number']=gdf['junction_number'].astype(float) # somehow junction_number in gdf is an object and not an int. To merge both datasets should be ints.
    gdf['junction_number']=gdf['junction_number'].astype(int)  
    merged = pd.merge(segment_bearings_df, gdf, on = ["junction_number"], how = "right")
    
    # drop the rows with NO bedding measurements
    print("I am now removing all the segment that are not associated with bedding measurements")
    merged['mean_strike'].replace('', np.nan, inplace=True)
    merged.dropna(subset=['mean_strike'], inplace=True)
    
    # add column with segment_orientation
    merged["segment_orientation"] = merged.apply(lambda row: remove_directionality(row.segment_bearing_whole_segment), axis = 1)
    
    # add column with difference between segment orientation and mean strike (result is on [-90,90])
    merged["difference_segment_orientation_and_mean_strike"] = merged.apply(lambda row: find_diff_segm_or_and_mean_strike_or_dip_or(row.segment_orientation, row.mean_strike), axis = 1)
    
    # and column with difference between segment orientation and mean dip orientation (result also on [-90,90])
    merged["difference_segment_orientation_and_mean_dip_orientation"] = merged.apply(lambda row: find_diff_segm_or_and_mean_strike_or_dip_or(row.segment_orientation, row.mean_dip_orientation), axis = 1)
                
    # make a new segment_bearings csv with beddings
    print("I am now saving a new segment_bearings csv for segments with bedding stats")
    merged.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings_bedding_stats.csv", index = False) 

    
##### FUNCTIONS TO MERGE CSVs FROM INDIVIDUAL BASIN FOLDERS TO ONE CSV IN PARENT FOLDER ####    
######## for part 1 ############
def merge_segment_bearings_bedding_stats(basins_dict):
    """
    Function to merge the segment_bearings_bedding_stats csv files from all basins in the hybas_dict to one csv, and sve it in the PART1 directory.
    """
    basins_dfs = []
    for hybas_id in basins_dict:
        basin_nickname=basins_dict[hybas_id]
        df = pd.read_csv("./basins/PART1/"+basin_nickname+"/"+basin_nickname+"_COP30_UTM_segment_bearings_bedding_stats.csv")
        basins_dfs.append(df)
    total_df = pd.concat(basins_dfs)
    # Save CSV
    total_df.to_csv("./basins/PART1/all_segment_bearings_bedding_stats.csv")
    print("Done")
    
def merge_unique_bedding_measurements(basins_dict):
    """
    Function to merge the unique_bedding_measurements csv files from all basins in the hybas_dict to one csv, and sve it in the PART1 directory.
    """
    basins_dfs = []
    for hybas_id in basins_dict:
        basin_nickname=basins_dict[hybas_id]
        df = pd.read_csv("./basins/PART1/"+basin_nickname+"/"+basin_nickname+"_unique_bedding_measurements.csv")
        basins_dfs.append(df)
    total_df = pd.concat(basins_dfs)
    
    # Drop rows with duplicate values in the 'ID' column (this could happen when two basins are next to each other)
    total_df = total_df.drop_duplicates(subset='ID')
    # Save CSV
    total_df.to_csv("./basins/PART1/all_unique_bedding_measurements.csv")
    print("Done")
######## for part 2 ####################

def create_all_segment_bearings_and_junction_angles_csvs(subdivisions):
    """
    Function to put together each segment bearings csv together into 1 big one, and also for junction angles.
    On the way, NEW unique Basin IDs are generated such that basins can be separated. Also a new column will be introduced stating whether it is a VRP or AP basin.
    """
    segment_bearings_dfs = []
    junction_angles_dfs = []
    count = 1
    
    for subdivision in subdivisions:
        basin_nickname = subdivision
        # open dfs
        segment_bearings_df = pd.read_csv("./basins/PART2/"+basin_nickname+"/"+basin_nickname+"_COP30_UTM_segment_bearings_filtered_nolakes.csv")
        junction_angles_df =  pd.read_csv("./basins/PART2/"+basin_nickname+"/"+basin_nickname+"_COP30_UTM_FULL_Basins_JAngles_filtered_nolakes.csv")
        
        # add new basin key such that all basins have a unique id when merged
        segment_bearings_df["new_basin_key"] = str(count) + segment_bearings_df["basin_key"].astype(str)
        junction_angles_df["new_basin_key"] = str(count) + junction_angles_df["basin_key"].astype(str)
        count+=1 
            
        # append df to list
        segment_bearings_dfs.append(segment_bearings_df)
        junction_angles_dfs.append(junction_angles_df)
       
    # Concatenate the modified dataframes into a single dataframe
    segment_bearings_df = pd.concat(segment_bearings_dfs, ignore_index=True)
    
    junction_angles_df = pd.concat(junction_angles_dfs, ignore_index=True)
        
    # Extract the prefix ("AP" or "VRP") from the keys of the subdivisions dictionary
    prefix_list = [key.split('_')[0] for key in subdivisions.keys()]
    
    # Add a new column "key_prefix" to result_df based on the prefix_list
    segment_bearings_df['Province'] = pd.concat([pd.Series([prefix] * len(df)) for prefix, df in zip(prefix_list, segment_bearings_dfs)], ignore_index=True)
    junction_angles_df['Province'] = pd.concat([pd.Series([prefix] * len(df)) for prefix, df in zip(prefix_list, junction_angles_dfs)], ignore_index=True)
    
    # Save the result
    segment_bearings_df.to_csv("./basins/PART2/all_segment_bearings.csv", index = False)
    junction_angles_df.to_csv("./basins/PART2/all_junction_angles.csv", index = False)
    
        
    
##################### INTERNAL FUNCTIONS #####################

def turn_streams_into_lines(chi_data_gdf, segment_bearings_gdf):
    """
    This function returns a geodataframe containing stream segments (along-stream) as lines.
    Requires chi_data_df and corresponding segment_bearings_df
    
    Author: RV
    
    Date: 11/10/2023  
    """
    
    # put all the source keys of segment bearings in list (this will be used later)
    source_keys = segment_bearings_gdf["source_key"].tolist()

    chi_data_gdf_merged = pd.merge(chi_data_gdf, segment_bearings_gdf, on=['geometry_points'], how = "left")
    print("I have now mentally merged your chi data map and segment bearings, so I will now loop through this merging result to produce line shapes of river segments")
    
    # now loop through the merged dataset, collecting points in a points_list which is used to create
    # a line when hitting a junction. Then points list is cleared. This is repeated until the end of the merged dataset is reached.
    lines_list = []
    junction_numbers_list = []
    previous_source_key = 0
    
    points_list = []
    chi_data_gdf_merged = chi_data_gdf_merged.reset_index()
    for index, row in chi_data_gdf_merged.iterrows():
        if pd.isna(row['junction_number']) == True:
            # checks if this node is a junction: this is no junction
            if row["source_key_x"] in source_keys:
                # this 'source key check' is there because without it, there would be occurences of streams consisting of only 1 node
                # that are not in the segment_bearings file, so it would not be recognised as a junction by this loop.
                # it resulted in this 1 node stream being added to the stream directly above it in the chi data file which is wrong
                # this error can be filtered out with introducing the source_key check
                
                # we will also need to check for a rare case where ther is a one-noded stream that is NOT in the segment_bearings file but DOES HAVE
                # a source_key that occurs in the segment_bearings file because of the stream immediately downstream. 
                # we solve this by checking if the source_key is different from the source_key one above
                if row["source_key_x"] == previous_source_key: 
                    points_list.append(row['geometry_points']) #only append geometry points if source key is the same (otherwise it will be different stream)
    
        elif pd.isna(row['junction_number']) == False:
            # checks if this node is a junciton: this is a junction
            if len(points_list) == 0:
                # this happens at the very first data point, which is technically a junction, but we don't want to 'cut off' a river segment before because there is none
                # so add this point to the points list, and add to junction_number_list
                points_list = [row['geometry_points']]
                junction_numbers_list.append(row['junction_number'])
                previous_source_key = row['source_key_x']
                
            elif len(points_list) == 1:
                # this is a rare occasion where the river segment only consist of 1 node: a junction. (i.e., two junctions must be next to each other)
                # we cannot make a LineString from just the 1 point. It's segment bearing won't tell anything anyways. So lets just remove the node and remove it from the junction_numbers_list
                # and report that this happened:
                print("Just had a 1-node long segment. I removed it. It's junction number (upstream) was", str(junction_numbers_list[-1]))
                del junction_numbers_list[-1]
                # start a new points_list with the current junction, and add the new junction number to the junction_number_list.
                points_list = [row['geometry_points']]
                junction_numbers_list.append(row['junction_number']) 
                previous_source_key = row['source_key_x']
                
            elif len(points_list) > 1:
                # hurray, we can no make a LineString of river segments with at least 2 nodes.
                lines_list.append(LineString(points_list))
                # and then do the following two things to make a new river segment starting at the junction
                points_list = [row['geometry_points']]
                junction_numbers_list.append(row['junction_number'])
                previous_source_key = row['source_key_x']
    
    # make a dataframe with the lines and junction numbers that we just made 
    segment_lines_df = pd.DataFrame(list(zip(junction_numbers_list, lines_list)), columns = ["junction_number", "lines"])
  
    # turn that into a geodataframe 
    segment_lines_gdf = gpd.GeoDataFrame(segment_lines_df, geometry = "lines")
    
    # return the segment_lines_gdf
    return segment_lines_gdf

def find_lines_fully_inside_quads(segment_lines_gdf, filepath_to_beddingquads):
    """
    This function finds which lines in segment_lines_gdf are fully contained within the bedding quads and gets rid of those that don't and returns segment_lines_gdf
    """
    print("Now starting to remove segment lines that are (partially) outside the quadrangles with bedding data.")
    # read the beddingquads file 
    quads_gdf = gpd.read_file(filepath_to_beddingquads)
    
    # Put the beddinquads in EPSG:4326 (that is, the coordinate reference system of the lines)
    quads_gdf = quads_gdf.to_crs("EPSG:4326")
    
    # Dissolve all polygons into a single polygon
    dissolved_polygon = quads_gdf.unary_union
    
    # Check if each line is fully within the dissolved polygon and filter the GeoDataFrame
    segment_lines_gdf = segment_lines_gdf[segment_lines_gdf.geometry.within(dissolved_polygon)]
    
    # return segment_lines_gdf
    return segment_lines_gdf

def find_lines_fully_outside_lakes(segment_lines_gdf, filepath_to_hydrolakes):
    """
    This function finds which lines in segment_lines_gdf are fully outside the hydrolakes and deletes the segment_lines that are not.
    """
    print("Now starting to remove segment lines that are (partially) inside lakes.")
    
    # open hydrolakes
    hydrolakes_gdf = gpd.read_file(filepath_to_hydrolakes)
    
    # Put the hydrolakes in EPSG:4326 (that is, the coordinate reference system of the lines)
    hydrolakes_gdf = hydrolakes_gdf.to_crs("EPSG:4326")
 
    # Dissolve all hydrolake polygons into a single polygon
    dissolved_polygon = hydrolakes_gdf.unary_union

    # Check if each line is fully outside the dissolved polygon and filter the GeoDataFrame
    segment_lines_gdf = segment_lines_gdf[~segment_lines_gdf.geometry.intersects(dissolved_polygon)]
    
    # return segment_lines_gdf
    return segment_lines_gdf

def buffer_lines(segment_lines_gdf, chi_data_gdf, radius = 1000):
    """
    This function buffers lines of a geodataframe by a given radius.
    Returns the same geodataframe with an extra column (buffer polygons)
    
    Args:
        segment_lines_gdf: the geodataframe with the lines
        chi_data_gdf: important! the geodataframe of chi data map, needed to check for the right EPSG before buffering
        radius: default 1000m    
    """
    print("I am so excited. I am now buffering lines!!")    
    # assign crs to the lines geoseries 
    segment_lines_gdf["lines"].crs = "EPSG:4326"
    
    # determine ESPG of basin and convert to the right EPSG (need to do before buffering). Code snippets from lsdmap_otgrabber.py's to_UTM function and from https://gis.stackexchange.com/questions/365584/convert-utm-zone-into-epsg-code
    lat = chi_data_gdf["latitude"].iloc[0] #this is just the first datapoint in the chi data map. assume that the dataset does not span multiple WGS zones.
    lon = chi_data_gdf["longitude"].iloc[0]
    temp_info = utm.from_latlon(lat,lon)
    
    if(temp_info[3] in ['X','W','V','U','T','S','R','Q','P','N']):
        south = False
    else:
        south = True
    crs = CRS.from_dict({'proj':'utm', 'zone':temp_info[2],'south': south})
    print("Oooh wait I need to transform this first to: ")
    print(crs.to_authority())
    print("in order to buffer correctly!")
    segment_lines_gdf['lines'] = segment_lines_gdf["lines"].to_crs(crs.to_authority())
        
    # buffer all the lines and put the resulting polygons in another column
    segment_lines_gdf['lines_buffer'] = segment_lines_gdf.geometry.buffer(radius)
    segment_lines_gdf['lines_buffer'].crs = crs.to_authority()
    
    # go back to EPSG:4326 (lat/lons)
    print("Finished buffering, returning all data to EPSG:4326 (lat/lons)")
    segment_lines_gdf['lines'] = segment_lines_gdf['lines'].to_crs(4326)
    segment_lines_gdf['lines_buffer'] = segment_lines_gdf['lines_buffer'].to_crs(4326)    

    # we now have a gdf with junction number, gds of (downstream) river line path, and gds of (downstream) river line path buffer
    # the junction number can be used to merge with segment bearings dataset    
    return segment_lines_gdf
    
def find_bedding_measurements_within_polygons(geodataframe, filepath_beddingdata):
    """
    This function takes a geodataframe with at least 1 column that is called 'lines_buffer', and then adds three columns (dip direction, azimuth, dip) with lists of all those measurements found in that buffer. 
    Requires a filepath to the bedding dataset to be set.
    Returns the new geodataframe
    
    Author: RV
    
    Date: 11/10/2023
    """
    print("now opening bedding data file at " + filepath_beddingdata)
    
    # load the bedding data into a geodataframe 
    bedding_gdf =  gpd.read_file(filepath_beddingdata)
    
    # loop through the geodataframe with segment buffers, and for every buffer polygon ask the bedding dataset "are your points within this polygon?
    print("Going to connect bedding measurements to segments. This is going to take a while.")
    number_of_segments = len(geodataframe["junction_number"])
    
    bed_ids = set([]) # keep track of the unique ids of bed measurements connected to segments. will use later 
    
    azimuths = []
    dips = []
    dip_directions = []
    count = 1 #for progress
    for index, row in geodataframe.iterrows():
        # print percentage of progress to terminal
        progress = count/number_of_segments * 100
        print(str(progress)+"%")
        count+=1
        
        # extract the polygon of this segment
        polygon = row['lines_buffer']
        
        # copy bedding gdf into new df
        this_bedding_gdf = bedding_gdf
        
        # check if bedding is in polygon. results in a true/false series that can be appended as a row to this_bedding_gdf
        this_bedding_gdf["inbuffer"] = bedding_gdf.within(polygon)
        
        # then we need to loop through this_bedding_gdf, putting azimuths, dips, dip_directions in a list 
        az = []
        di = []
        dd = []
        for ind, rw in this_bedding_gdf.iterrows():
            if rw["inbuffer"] == True: 
                if (int(rw["Azimuth"]), int(rw["Dip_Angle"]), int(rw["Dip_Direct"])) != (0,0,0): # only append bedding measurements if not all of azimuth, dip angle and dip direction are 0. It is very, very, very highly likely these (0,0,0) bed measurements are not representative of the real situation so do not include.
                    az.append(int(rw["Azimuth"]))
                    di.append(int(rw["Dip_Angle"]))
                    dd.append(int(rw["Dip_Direct"]))
                    bed_ids.add(int(rw["ID"])) 
                
        # then the lists az, di, dd need to be appended to the azimuths, dips, dip_directions lists respectively
        azimuths.append(az)
        dips.append(di)
        dip_directions.append(dd)
            
    # append the azimuth, dips and dip_directions lists with lists to a new dataframe containing the segments. first check if they are all equal length (should be)
    print("I have "+ str(len(geodataframe["junction_number"])) + "segments, and your azimuths, dips, and dip direction lists have lengths of " + str(len(azimuths)) + ", " + str(len(dips)) + ", " + str(len(dip_directions)) + "respectively")
    
    geodataframe["azimuths"] = azimuths
    geodataframe["dips"] = dips
    geodataframe["dip_directions"] = dip_directions
    
    # drop the columns with the lines and lines_buffer because we don't need it anymore
    geodataframe.drop(['lines', 'lines_buffer'], axis = 1)
    
    # remove all bedding data points that are not in bed_ids from the bedding_gdf
    bedding_gdf = bedding_gdf[bedding_gdf["ID"].isin(bed_ids)]
    bedding_df = bedding_gdf.drop(columns='geometry')
    
    
    # now return the geodataframe (which is kind of not a gdf anymore since it doesn't have geoseries contained in it now we dropped lines and lines buffer)
    return geodataframe, bedding_df


#### THIS FOLLOWING NOT USED IN THE CURRENT VERSION OF THE DATA ANALYSIS NOTEBOOK BUT INCLUDED FOR REFERENCE #### 

def add_junction_angles_to_segment_bearings(hybas_id):
    """
    Merges junction_angles dataset to segment_bearings_bedding_stats
    
    Args:
        hybas_id (int)
        
    Returns:
        saves a NEW csv called with the name segment_bearings_bedding_stats_junction_angles.csv
    
    Date:
        13/10/23
    """
    # open junction angles file, put in dataframe
    Dataset_prefix = str(hybas_id)
    DataDirectory = "./Basins/"+str(hybas_id)+"/"
    junction_angles_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_FULL_Basins_JAngles.csv")
                    
    # make a drop all columns except junction_number, donors_junction_angle, donor1_receiver_junction_angle and donor2_receiver_junction_angle
    junction_angles_df = junction_angles_df[["junction_number", "donors_junction_angle", "donor1_receiver_junction_angle", "donor2_receiver_junction_angle"]]
            
    # rename columns to include the fact that it concerns the upstream junction angle
    junction_angles_df.rename(columns={"donors_junction_angle":"upstream_donors_junction_angle", "donor1_receiver_junction_angle":"upstream_donor1_receiver_junction_angle", "donor2_receiver_junction_angle":"upstream_donor2_receiver_junction_angle"})
    
    # open segment bearings with beddings dataset
    segment_bearings_df = pd.read_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings_bedding_stats.csv")
                    
    # merge the junction angles and segment bearings datasets based on junction angle
    # note that this merge result will have NaN values at first order streams since they obviously don't have an upstream junction angle
    merged = pd.merge(segment_bearings_df, junction_angles_df, on = ["junction_number"], how = "left")
                    
    # write result out to csv
    merged.to_csv(DataDirectory + Dataset_prefix + "_COP30_UTM_segment_bearings_bedding_stats_junction_angles.csv", index = False) 
                    
                    