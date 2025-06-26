import os
import click
import glob
from utils import rasterize_shapefile

@click.command()
@click.option("--ci",  "-i", default=1, type=int, help="The core's id")
@click.option("--cn",  "-n", default=1, type=int, help="The number of cores")
@click.option("--dir_raster", "-r", default='/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Result_ChangeMap_new/Mosaic_Sin/', type=str, help="The directory of the NTL change maps")
@click.option("--dir_shp", "-s", default='/shared/zhulab/Tian/Analysis/AdmBoundary/', type=str, help="The directory of shapefile/destination folder")
@click.option("--fname_shp", "-f", default='WB_countries_Admin0_10m.shp', type=str, help="The file name of the input shapefile")

def main(ci, cn, dir_raster, dir_shp, fname_shp):
    # msg of core
    print('\n*****************************************************************************************************')
    print('* function: create the map of World Bank administration rasters for sinusoidal tiles')
    print('* core: {:04d}/{:04d}'.format(ci, cn))

    # Get the tile list
    tile_refs = sorted(glob.glob(os.path.join(dir_raster, 'h*v*')))
    # select tiles by the core id and the number of cores
    tilelist = [os.path.join(dir_raster, tile_refs[i]) for i in range(ci-1, len(tile_refs), cn) ]
    # Get the polygon path
    path_shp = dir_shp + fname_shp

    ## Loop converting
    # Create task record object
    for tile in tilelist:
        tile_name = os.path.basename(tile)
        fname_ref = tile_name + '_AbruptChangeFrequency_20142022.tif'
        path_ref = os.path.join(tile, fname_ref)

        fname_save_continent = "Continent_{}_sin.tif".format(tile_name)
        path_save_continent = os.path.join(dir_shp, 'Continent_sin', fname_save_continent)
        fname_save_country = "Country_{}_sin.tif".format(tile_name)
        path_save_country = os.path.join(dir_shp, 'Country_sin', fname_save_country)

        print("Converting the {} to the tile {} like raster".format(fname_shp, tile_name))
        rasterize_shapefile(path_shp, path_ref, path_save_continent, field='CNTNT_ID')
        rasterize_shapefile(path_shp, path_ref, path_save_country, field='CNTRY_ID')

# Main function
if __name__ == "__main__":
    main()
