import utm
from pyproj import CRS
import pandas as pd
import geopandas as gpd
import numpy as np
import fiona 
import rasterio
import rasterio.mask as mask
import json
import lsdviztools.lsdbasemaptools as bmt
from lsdviztools.lsdplottingtools import lsdmap_gdalio as gio
import os
from rasterio.warp import reproject, Resampling, calculate_default_transform

def get_buffered_basin(hybas_id, filepath = './hydrobasins'):
    """
    This function returns a geopandas geoseries of a buffered basin from a given hybas_id and filepath to hydrobasins dataset.
    
    Author: RV
    
    Date: 15/08/2023  
    """
    # get filepath right
    pfaf_lvl = str(hybas_id)[1:3]
    filepath = filepath + '/hybas_na_lev'+pfaf_lvl+'_v1c.shp'
    print("Reading the Hydrobasins dataset at: "+filepath)
    
    # open file, get selected basin
    dataset = gpd.read_file(filepath)
    selected_basin = dataset.loc[dataset["HYBAS_ID"]==hybas_id]
        
    # find buffer radius (according to Strong and Mudd, 2022)
    basin_area = selected_basin.iloc[0,6]
    if basin_area < 100:
        buffer_radius = 2*1000
    elif basin_area > 15625:
        buffer_radius = 25*1000
    else: 
        buffer_radius = 0.25*np.sqrt(basin_area)*1000
    
    # get geometry of selected basin
    selected_basin_geometry = selected_basin.geometry
    
    # determine ESPG of basin and convert to the right EPSG (need to do before buffering). Code snippets from lsdmap_otgrabber.py's to_UTM function and from https://gis.stackexchange.com/questions/365584/convert-utm-zone-into-epsg-code
    lat = selected_basin_geometry.centroid.y.iloc[0]
    lon = selected_basin_geometry.centroid.x.iloc[0]
    temp_info = utm.from_latlon(lat,lon)
    
    if(temp_info[3] in ['X','W','V','U','T','S','R','Q','P','N']):
        south = False
    else:
        south = True
    crs = CRS.from_dict({'proj':'utm', 'zone':temp_info[2],'south': south})        
    selected_basin_geometry = selected_basin_geometry.to_crs(crs.to_authority())
    
    # buffer basin     
    selected_basin_buffered = selected_basin_geometry.buffer(buffer_radius)
    
    # return buffered basin geoseries
    print("Finished buffering basin "+str(hybas_id))
    return selected_basin_buffered    
    
def get_ll_ur_coords_of_geoseries(geoseries):
    """
    This function returns a list of the lower left and upper right coordinates of a rectangle enveloping a geopandas geoseries
    
    Author: RV
    
    Date: 15/08/2023
    """
    geoseries = geoseries.to_crs(4326)
    envelope_coords = geoseries.envelope.get_coordinates()
    ur = [max(envelope_coords["y"].tolist()),max(envelope_coords["x"].tolist())]
    ll = [min(envelope_coords["y"].tolist()),min(envelope_coords["x"].tolist())]
    return [ll, ur]

def clip_raster_to_geoseries(geoseries, in_rasterpath, out_rasterpath):
    """
    This function clips any raster to any geoseries
    
    Author: RV
    
    Date: 16/08/2023
    """
    coords = [json.loads(geoseries.to_json())['features'][0]['geometry']]
    with rasterio.open(in_rasterpath) as src:
        out_image, out_transform = mask.mask(src, shapes=coords, crop=True)
        out_meta = src.meta

    out_meta.update({"driver":"ENVI",
                     "height":out_image.shape[1],
                     "width":out_image.shape[2],
                     "transform":out_transform})

    with rasterio.open(out_rasterpath, "w", **out_meta) as dest:
        dest.write(out_image)
    print("Finished clipping "+in_rasterpath+" to geoseries, new file is: "+out_rasterpath)
    
def setup_directory_with_DEM(hybas_id, nickname, part):
    """
    This function sets up a directory with DEM for one basin based on hybas_id.
    
    Author: RV
    
    Date: 16/08/2023  
    """
    # set dataset prefix and data directory
    Dataset_prefix = nickname
    DataDirectory = "./Basins/"+part+"/"+nickname+"/"
    
    # check if there is already a data directory with the same hybas id. Will not continue if it already exists.
    if os.path.exists(DataDirectory) == True:
        print("The directory " + DataDirectory + " already exists. I am not going to download you data for this basin: " + str(hybas_id) + ":" + nickname)
        return

    # get basin buffered
    bb = get_buffered_basin(hybas_id)

    # get upper right and lower left coordinates of buffered basin
    coords = get_ll_ur_coords_of_geoseries(bb)

    # read the OpenTopography API key
    your_OT_api_key_file = "my_OT_api_key.txt"

    with open(your_OT_api_key_file, 'r') as file:
        print("I am reading you OT API key from the file "+your_OT_api_key_file)
        api_key = file.read().rstrip()
        print("Your api key starts with: "+api_key[0:4])

    # set source to COP30
    source_name = "COP30"

    # download the DEM and convert for LSDTT
    DEM = bmt.ot_scraper(source = source_name,
                         lower_left_coordinates = coords[0],
                         upper_right_coordinates = coords[1],
                         prefix = Dataset_prefix,
                         api_key_file = your_OT_api_key_file,
                         path = DataDirectory)
    DEM.print_parameters()
    DEM.download_pythonic()
    Fname = Dataset_prefix+"_"+source_name+".tif"
    gio.convert4lsdtt(DataDirectory,Fname)   

    # clip / cookie-cut download to buffered basin
    clip_raster_to_geoseries(geoseries = bb, in_rasterpath = DataDirectory + Dataset_prefix+"_"+source_name + "_UTM.bil", out_rasterpath= DataDirectory + Dataset_prefix+"_"+source_name + "_UTM_clipped.bil")
    
    # now remove the tif files and the non-clipped bil/hdr files as we don't need them anymore. Saving pc storage.
    os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.bil")
    os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.bil.aux.xml")
    os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.hdr")
    os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.tif")
    os.remove(DataDirectory+Dataset_prefix+"_"+source_name+".tif")
    
def setup_directories_with_subdivision_outlines(subdivisions, buffer_radius = 2000, filepath_to_hydrobasins = "./hydrobasins"):
    """
    This function sets up directories and makes outlines for PART 2
    
    Author: RV
    
    Date: 10/01/2023  
    """
    for subdivision in subdivisions:
        # set dataset prefix and data directory
        Dataset_prefix = subdivision
        DataDirectory = "./Basins/PART2/"+subdivision+"/"

        # check if there is already a data directory with the same subdivision name. Will not continue if it already exists.
        if os.path.exists(DataDirectory) == True: 
            print("The directory " + DataDirectory + " already exists. I am not going to set up a folder for: "+ subdivision)
        else:   
            # get list of hydrobasins
            hydrobasins_list = subdivisions[subdivision]
            
            # create a geodataframe with one polygon describing all selected hydrobasins
            basins_gdfs = []
            for hybas_id in hydrobasins_list: 
                # get the right hydrobasins
                pfaf_lvl = str(hybas_id)[1:3]
                filepath = filepath_to_hydrobasins + '/hybas_na_lev'+pfaf_lvl+'_v1c.shp'
                print("Reading the Hydrobasins dataset at: "+filepath)

                # open file, get selected basin
                dataset = gpd.read_file(filepath)
                selected_basin = dataset.loc[dataset["HYBAS_ID"]==hybas_id]
                
                # append selected basin gdf to basins_gdfs list
                basins_gdfs.append(selected_basin)
            combined_basins_gdf = gpd.GeoDataFrame(pd.concat(basins_gdfs, ignore_index=True))
            total_basins_gs = gpd.GeoSeries([combined_basins_gdf.unary_union])
            total_basins_gs.crs = "EPSG:4326"
            
            # determine ESPG of total basins outline and convert to the right EPSG (need to do before buffering). Code snippets from lsdmap_otgrabber.py's to_UTM function and from https://gis.stackexchange.com/questions/365584/convert-utm-zone-into-epsg-code
            lat = total_basins_gs.centroid.y.iloc[0]
            lon = total_basins_gs.centroid.x.iloc[0]
            temp_info = utm.from_latlon(lat,lon)

            if(temp_info[3] in ['X','W','V','U','T','S','R','Q','P','N']):
                south = False
            else:
                south = True
            crs = CRS.from_dict({'proj':'utm', 'zone':temp_info[2],'south': south})        
            total_basins_gs = total_basins_gs.to_crs(crs.to_authority())

            # buffer basin     
            total_basins_gs_buffered = total_basins_gs.buffer(buffer_radius)

            # load the outline of either the VRP or the AP
            if subdivision[:2] == "VR":
                outline_gdf = gpd.read_file("./vrpapoutlines/VRP_OUTLINE_NONGLAC.shp")
            elif subdivision[:2] == "AP":
                outline_gdf = gpd.read_file("./vrpapoutlines/AP_OUTLINE_NONGLAC.shp")
           
            
            # get geometry of outline and convert to the right crs
            outline_gs = gpd.GeoSeries([outline_gdf.unary_union])
            outline_gs.crs = "EPSG:4326"
            outline_gs = outline_gs.to_crs(crs.to_authority())
            
            # now find the intersection of the outline_gs and the total_basins_gs_buffered
            this_subdivision_gs = gpd.GeoSeries(outline_gs.geometry.intersection(total_basins_gs_buffered.geometry))
                                                
            # save the subdivision outline as a GeoJSON for inspection
            os.makedirs(DataDirectory, exist_ok=True)
            this_subdivision_gs.to_file(DataDirectory+Dataset_prefix+"_subdivision_outline.geojson", driver = "GeoJSON")
            print("I saved the outline of subdivision " + subdivision + " as a geojson. Now go read the Jupyter Notebook instruction.")
        
def prepare_DEMS_for_subdivision_outlines(subdivisions):
    """
    This function prepares DEMs for PART 2
    
    Author: RV
    
    Date: 10/01/2023  
    """
    
    for subdivision in subdivisions:
        # set dataset prefix and data directory
        Dataset_prefix = subdivision
        DataDirectory = "./Basins/PART2/"+Dataset_prefix+"/"

        # check if there is already a data directory with a clipped DEM. Will not continue if it already exists.
        if os.path.exists(DataDirectory+Dataset_prefix+"_COP30_UTM_clipped.bil") == True: 
            print(DataDirectory+Dataset_prefix+"_COP30_UTM_clipped.bil already exists. I am not going to download you a DEM again!")
        else:
            # load geojson
            this_subdivision_gs = gpd.read_file(DataDirectory+Dataset_prefix+"_subdivision_outline.geojson")
            
            # request the upper right and lower left coordinates of subdivision geometry
            coords = get_ll_ur_coords_of_geoseries(this_subdivision_gs)
           
            # read the OpenTopography API key
            your_OT_api_key_file = "my_OT_api_key.txt"

            with open(your_OT_api_key_file, 'r') as file:
                print("I am reading you OT API key from the file "+your_OT_api_key_file)
                api_key = file.read().rstrip()
                print("Your api key starts with: "+api_key[0:4])

            # set source to COP30
            source_name = "COP30"

            # LSDTT OT grabber script wants to make the DataDirectory. If it already exists (which is the case) it will give an error.
            # Thus first remove the DataDirectory. 
            os.remove(DataDirectory+Dataset_prefix+"_subdivision_outline.geojson")
            os.rmdir(DataDirectory)
            
            # download the DEM and convert for LSDTT
            DEM = bmt.ot_scraper(source = source_name,
                                 lower_left_coordinates = coords[0],
                                 upper_right_coordinates = coords[1],
                                 prefix = Dataset_prefix,
                                 api_key_file = your_OT_api_key_file,
                                 path = DataDirectory)
            DEM.print_parameters()
            DEM.download_pythonic()
            Fname = Dataset_prefix+"_"+source_name+".tif"
            gio.convert4lsdtt(DataDirectory,Fname) 
        
            # clip / cookie-cut download to outline
            clip_raster_to_geoseries(geoseries = this_subdivision_gs, in_rasterpath = DataDirectory + Dataset_prefix+"_"+source_name + "_UTM.bil", out_rasterpath= DataDirectory + Dataset_prefix+"_"+source_name + "_UTM_clipped.bil")

            # now remove the tif files and the non-clipped bil/hdr files as we don't need them anymore. Saving pc storage.
            os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.bil")
            os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.bil.aux.xml")
            os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.hdr")
            os.remove(DataDirectory+Dataset_prefix+"_"+source_name+"_UTM.tif")
            os.remove(DataDirectory+Dataset_prefix+"_"+source_name+".tif")
            
            # becasue we removed the direcotry including the GeoJSON of the outline, we want to save the this_subdivision_gs again
            this_subdivision_gs.to_file(DataDirectory+Dataset_prefix+"_subdivision_outline.geojson", driver = "GeoJSON")
        
def reproj_match(infile, match, outfile):
    # 21/11/2023: THIS FUNCTION HAS BECOME REDUNDANT IN NEWER VERSIONS OF THE DATA ANALYSIS NOTEBOOK - BUT IS STILL INCLUDED FOR FUTURE REFERENCE!
    
    """Reproject a file to match the shape and projection of existing raster. FROM https://pygis.io/docs/e_raster_resample.html
    
    Parameters
    ----------
    infile : (string) path to input file to reproject
    match : (string) path to raster with desired shape and projection 
    outfile : (string) path to output file tif
    """
    # open input
    with rasterio.open(infile) as src:
        src_transform = src.transform
        
        # open input to match
        with rasterio.open(match) as match:
            dst_crs = match.crs
            
            # calculate the output transform matrix
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src.crs,     # input CRS
                dst_crs,     # output CRS
                match.width,   # input width
                match.height,  # input height 
                *match.bounds,  # unpacks input outer boundaries (left, bottom, right, top)
            )

        # set properties for output
        dst_kwargs = src.meta.copy()
        dst_kwargs.update({"crs": dst_crs,
                           "transform": dst_transform,
                           "width": dst_width,
                           "height": dst_height,
                           "nodata": 0})
        print("Coregistered to shape:", dst_height,dst_width,'\n Affine',dst_transform)
        # open output
        with rasterio.open(outfile, "w", **dst_kwargs) as dst:
            # iterate through bands and write using reproject function
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)