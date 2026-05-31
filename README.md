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
Visualization of the dataset can be found through [Google Earth Engine](https://ee-downloading.projects.earthengine.app/view/alan-change)-based application.

### Product Access
Collection 1 (beta version) of the global ALAN change product dataset (2014–2022) and the related codes, described in Li et al. (2026), can be downloaded from [this Zenodo repository](https://doi.org/10.5281/zenodo.18264642).

### Product Format
The products are provided in the Suomi-NPP VIIRS 15-arc-second linear latitude/longitude grids (or geographic) (see Figure 2 in [this document](https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_Collection2.0.pdf)).
Annual change maps are delivered in GeoTIFF format, including:
 1) Abrupt change time in Day-of-Year<br>
    Filename: AbruptChangeTime_YYYY.tif (where YYYY means year)<br>
    Pixel Value: Day-of-year of the detected abrupt ALAN change events<br>
   
3) Abrupt change intensity:<br>
    Filename: AbruptChangeIntensity_YYYY.tif (where YYYY means year)<br>
    Pixel value: Intensity of the abrupt ALAN change in units nW·cm⁻²·sr⁻¹<br>
   
4) Gradual change intensity: <br>
    Filename: GradualChangeIntensity_YYYY.tif (where YYYY means year)<br>
    Pixel value: Intensity of the gradual ALAN change in units nW·cm⁻²·sr⁻¹<br>

## User Discussion and Feedback
For user discussions or comments, please feel free to contact us by email or visit [the project’s issue page](https://github.com/GERSL/VZA-COLD/issues)

## Contact US
Tian Li (tianli@uconn.edu) and Zhe Zhu (zhe@uconn.edu) at the Department of Natural Resources and the Environment, University of Connecticut.

## Reference

Li, T., Wang, Z., Kyba, C., Román, M. O., Seto, K. C., Yang, Y., ... & Zhu, Z. (2026). Satellite imagery reveals increasing volatility in human night-time activity. Nature, 652(8109), 379-386.

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.
