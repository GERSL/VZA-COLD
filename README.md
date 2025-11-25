# VZA-COLD: Viewing Zenith Angle (VZA) stratified COntinuous monitoring of Land Disturbance (COLD)
VZA-COLD detects nighttime light (NTL) change at 15 arc-second spatial resolution with daily updating capability based on [NASA's Black Marble products](https://blackmarble.gsfc.nasa.gov/).

This page has been created for peer review purposes.

## Code Explanation
We provide the code resources (package name associated with the folder provided in the code) for our publication, which are programmed in MATLAB (2022b) and Python (3.10):
- **VZACOLD**: Nighttime light change detection package (MATLAB)
- **Analysis**: Code for data analysis and visualization (MATLAB + Python)

## NTL Change Product
Our global NTL change product provides nighttime light change information, including change time (Day-of-year and Year), change type (abrupt or gradual), and change intensity, for each year at the VIIRS 15-arc-second pixel level. 

### Product Preview
You can view the dataset through [Google Earth Engine](https://ee-downloading.projects.earthengine.app/view/alan-change)-based application.

### Product Access
You can download Version 1 of the global NTL product dataset (2014–2022), described in Tian et al. (2026), from [this link] (xx GB compressed; 80 GB uncompressed).

### Product Format
The products are provided in the Suomi-NPP VIIRS linear latitude/longitude (or geographic) grid with 10°×10° (see Figure 2 in [this document](https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_Collection2.0.pdf)).
For each grid tile, annual change maps are delivered in GeoTIFF format, including:
 1) Change time in Day-of-Year: Please insert the definition here
    e.g., h08v05_change_time_2020.tif, where h08v05 means the grid tile name and 2020 indicates the year.
    Pixel value indicates the Day-of-Year.
    
 3) Change type (abrupt or gradual): Please insert the definition here
    e.g., h08v05_change_type_2020.tif, where h08v05 means the grid tile name and 2020 indicates the year.
    | Pixel value | Change type        |
    |-------------|----------------|
    | 1           | Abrupt  |
    | 2           | Gradual |
    
3) Change intensity: Please insert the definition here
   e.g., h08v05_change_intensity_2020.tif, where h08v05 means the grid tile name and 2020 indicates the year.
   Pixel value is the intensity of the change in units: xxxx.

## User Discussion and Feedback
For user discussions or comments, please feel free to contact us by email or visit [the project’s issue page] (https://github.com/GERSL/VZA-COLD/issues)

## Contact US
Tian Li(tianli@uconn.edu) and Zhe Zhu (zhe@uconn.edu) at the Department of Natural Resources and the Environment, University of Connecticut.

## Reference

Li, T., Zhu, Z., Wang, Z., Kyba, C. C. M., Seto, K. C., Yang, Y., Qiu, S., Kuester, T., Fragkias, M., Chen, X., Meyer, T. H., Rittenhouse, C. D., Tai, X., Cullerton, M., Hong, F., Grinstead, A., Song, K., Suh, J. W., Yang, X., Kalb, V. L., Deng, C., & Román, M. O. Increasing volatility in human nighttime activity revealed by daily and high-resolution satellite imagery. 

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.
