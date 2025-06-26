# This script aims to conduct the reprojection from geographic coordinate system to sinusoidal coordinate system
# The input is the folder containing the raster files in geographic coordinate system
# The output is the folder containing the raster files in sinusoidal coordinate system
# The script will iterate through the input folder and reproject each raster file to the output folder

# import the necessary libraries
import os
import glob
import rasterio
from rasterio import warp
import click

# define the main function
@click.command()
@click.option("--ci",  "-i", default=1, type=int, help="The core's id")
@click.option("--cn",  "-n", default=1, type=int, help="The number of cores")
@click.option("--ref", "-r", default="/shared/zhulab/Tian/Analysis/PixelArea/Sinusoidal/", type=str, help="The filepath of the data's location")
@click.option("--src", "-s", default="/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Result_ChangeMap_new/Mosaic", type=str, help="The filepath of the data's location")
@click.option("--des", "-d", default="/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Result_ChangeMap_new/Mosaic_Sin", type=str, help="The filepath of the data's location")
@click.option("--fn", "-f", default="Global_*.tif", type=str, help="File name of the data")

#@click.option("--src", "-s", default="/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Analysis/DarkPixel/DarkMask", type=str, help="The filepath of the data's location")
#@click.option("--des", "-d", default="/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Analysis/DarkPixel/DarkMask_Sin", type=str, help="The filepath of the data's location")
#@click.option("--fn", "-f", default="DarkMask_Global.tif", type=str, help="File name of the data")

#@click.option("--src", "-s", default="/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/BlackMarble_BRDFcorrected_new/Mosaic_IniRad/", type=str, help="The filepath of the data's location")
#@click.option("--des", "-d", default="/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Result_ChangeMap_new/Mosaic_Sin", type=str, help="The filepath of the data's location")
#@click.option("--fn", "-f", default="Global_initialRadiance_VZAint4.tif", type=str, help="File name of the data")

def main(ci, cn, src, ref, des, fn):
    # print the input parameters
    print('ci:', ci)
    print('cn:', cn)
    print('src:', src)
    print('ref:', ref)
    print('des:', des)
    print('fn:', fn)

    # create des folder if it does not exist
    if not os.path.exists(des):
        os.makedirs(des)

    # get reference layers
    file_refs = sorted(glob.glob(os.path.join(ref, 'M*.A2020*.hdf')))
   
    # find the file in the folder from src
    file_srcs = sorted(glob.glob(os.path.join(src, fn)))

    # assign the file_srcs to different cores according to ic and in by parallel computing
    file_refs = file_refs[ci-1::cn]

    for file_ref in file_refs:
        # get the tile name
        tile_name = os.path.basename(file_ref).split('.')[2]
        des_tile = os.path.join(des, tile_name)

        # get the target image profile 
        with rasterio.open(file_ref) as ref:
            with rasterio.open(ref.subdatasets[0]) as subdataset:
                like_profile = subdataset.profile.copy()
                   
        for file_src in file_srcs:
            fname = os.path.basename(file_src).split('.')[0]
            # save as local file # warp as the destination geotif
            des_file = os.path.join(des_tile, '{:}_{:}.tif'.format(tile_name, fname[7:]))
            # Skip if the file already exists
            if os.path.exists(des_file):
                # TRY TO LOAD THe FILE as 2d array to make sure it is correct, if not, delete the file
                try:
                    with rasterio.open(des_file) as des_file_mask:
                        des_array = des_file_mask.read(1)
                        #print('Exist:', des_file)
                        del des_array
                    continue
                except:
                    os.remove(des_file)
                    print('Removed:', des_file)

            if not os.path.exists(des_tile):
                    os.makedirs(des_tile)

            # warp the src image into the target image
            try: 
                with rasterio.open(file_src, 'r') as src_mask: # obtain the resource dataset
                    des_profile = src_mask.meta.copy()
    
                    # update the profile to the image-like
                    des_profile.update({
                    'crs':          like_profile['crs'],
                    'transform':    like_profile['transform'],
                    'width':        like_profile['width'],
                    'height':       like_profile['height'],
                    # Ensure proper data type to avoid overflow
                    'dtype': src_mask.dtypes[0],  # Ensures matching data type with source
                    'nodata': 0  # Ensure NoData is explicitly set
                    })
                    with rasterio.open(des_file, 'w', **des_profile) as dst_mask:
                        for i in range(1, src_mask.count + 1):
                            warp.reproject(
                                source          =   rasterio.band(src_mask, i), 
                                destination     =   rasterio.band(dst_mask, i),
                                src_transform   =   src_mask.transform,
                                src_crs         =   src_mask.crs,
                                dst_transform   =   des_profile['transform'],
                                dst_crs         =   des_profile['crs'],
                                resampling      =   warp.Resampling.nearest,
                                src_nodata = src_mask.nodata,  # Handle NoData from the source
                                dst_nodata      =   0
                                )
                        print('Successfully saved:', des_file)
            except:
                print('Failed to save:', file_src)
                try:
                    os.remove(des_file)
                    print('Removed:', des_file)
                except:
                    print('Failed to remove:', des_file)
                continue

# run the main function
if __name__ == "__main__":
    main()