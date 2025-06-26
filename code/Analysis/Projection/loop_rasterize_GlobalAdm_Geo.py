import os
from os import path
import click
import glob
import time
import yaml
import fnmatch
from utils import rasterize_shapefile

@click.command()
@click.option("--dir_raster", "-r", default='/shared/zhulab/Tian/VIIRS_NTL/Tiff_ref/', type=str, help="The directory of the NTL change maps")
@click.option("--dir_shp", "-s", default='/shared/zhulab/Tian/Analysis/AdmBoundary/Country/', type=str, help="The directory of shapefile/destination folder")
@click.option("--fname_shp", "-i", default='WB_countries_Admin0_10m.shp', type=str, help="The file name of the input shapefile")

def main(dir_raster, dir_shp, fname_shp):
    # msg of core
    print('\n*****************************************************************************************************')
    print('* function: create the map of World Bank administration rasters for NTL tiles')
    print('* core: {:04d}/{:04d}'.format(ci, cn))
    print('* config file: {}'.format(fname_shp))

    # Get the tile list
    tilelist = os.listdir(dir_raster)
    tilelist = [os.path.splitext(fname)[0] for fname in tilelist]
    # Get the polygon path
    path_shp = dir_shp + fname_shp

    ## Loop converting
    # Create task record object
    for tile_name in tilelist:
        if fnmatch.fnmatch(tile_name, 'h*v*'):
            fname_NTL = "{}.tif".format(tile_name)
            path_NTL = dir_raster + fname_NTL

            fname_save = "Country_{}.tif".format(tile_name)
            path_save = dir_shp + fname_save

            print("Converting the {} to the tile {} like raster".format(fname_shp, tile_name))
            rasterize_shapefile(path_shp, path_NTL, path_save, field='IDs', dtype='uint16', fill=65535)

# Main function
if __name__ == "__main__":
    main()
