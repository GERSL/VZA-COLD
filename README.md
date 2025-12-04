# VZA-COLD: Viewing Zenith Angle (VZA) stratified COntinuous monitoring of Land Disturbance (COLD)
VZA-COLD detects artificial light at night (ALAN) change at 15 arc-second spatial resolution with daily updating capability based on [NASA's Black Marble products](https://blackmarble.gsfc.nasa.gov/).

This page has been created for peer review purposes.

## Code Explanation
We provide the code resources (package name associated with the folder provided in the code) for our publication, which are programmed in MATLAB (2022b) and Python (3.10):
- **VZACOLD**: Nighttime light change detection package (MATLAB)
- **Analysis**: Code for data analysis and visualization (MATLAB + Python)

## NTL Change Product
Our global NTL change product provides artificial light at night change information, including change time (Day-of-year and Year) and change intensity for abrupt change and gradual change, for each year at the VIIRS 15-arc-second pixel level. 

### Product Preview
You can view the dataset through [Google Earth Engine-based application](https://ee-downloading.projects.earthengine.app/view/alan-change).

### Product Access
You can download Collection 1 of the global NTL product dataset (2014–2022), described in Li et al. (2026), from [TBD].

### Product Format
The products are provided in the Suomi-NPP VIIRS linear latitude/longitude (or geographic) (see Figure 2 in [this document](https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_Collection2.0.pdf)).
Annual change maps are delivered in GeoTIFF format, including:
 1) Abrupt change time in Day-of-Year<br>
    Filename: AbruptChangeTime_YYYY.tif (where YYYY means year)<br>
    Pixel Value: Day-of-Year<br>
   
3) Abrupt change intensity:<br>
    Filename: AbruptChangeIntensity_YYYY.tif (where YYYY means year)<br>
    Pixel value: Intensity of the change in units nW·cm⁻²·sr⁻¹<br>
   
4) Gradual change intensity: <br>
    Filename: GradualChangeIntensity_YYYY.tif (where YYYY means year)<br>
    Pixel value: Intensity of the change in units nW·cm⁻²·sr⁻¹<br>

## User Discussion and Feedback
For user discussions or comments, please feel free to contact us by email or visit [the project’s issue page] (https://github.com/GERSL/VZA-COLD/issues)

## Contact US
Tian Li (tianli@uconn.edu) and Zhe Zhu (zhe@uconn.edu) at the Department of Natural Resources and the Environment, University of Connecticut.

## Reference

Li, T., Zhu, Z., Wang, Z., Kyba, C. C. M., Seto, K. C., Yang, Y., Qiu, S., Kuester, T., Fragkias, M., Chen, X., Meyer, T. H., Rittenhouse, C. D., Tai, X., Cullerton, M., Hong, F., Grinstead, A., Song, K., Suh, J. W., Yang, X., Kalb, V. L., Deng, C., & Román, M. O. Increasing volatility in human nighttime activity revealed by daily and high-resolution satellite imagery. 

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.
