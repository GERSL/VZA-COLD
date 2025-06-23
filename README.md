# VZA-COLD: Viewing Zenith Angle (VZA) stratified COntinuous monitoring of Land Disturbance (COLD)
VZA-COLD detects nighttime light (NTL) change at 15 arc-second spatial resolution with daily updating capability based on [NASA's Black Marble products](https://blackmarble.gsfc.nasa.gov/).


## Code Explanation
We provide the code resources (package name associated with the folder provided in the code) for our publication, which are programmed in MATLAB (2022b) and Python (3.10):
- **LandsatData**: Data processing and density analysis (Python)
- **COLD**: Disturbance detection package (MATLAB)
- **ODACA**: Disturbance agent classification package (MATLAB + Python)
- **Analysis**: Code for data analysis and visualization (MATLAB + Python (Jupyter Notebook))
- **Release**: Code for organizing the disturbance dataset  (MATLAB + Python)

## NTL change Product
This product provides disturbance time, agent, agent classification confidence, magnitude, and severity for each year at the Landsat 30m pixel level. For detailed descriptions of each layer, see [Land Disturbance Dataset page](https://github.com/GERSL/usdist/wiki/Land-Disturbance-Dataset).

You can view the dataset through [Google Earth Engine](https://ee-gers.projects.earthengine.app/view/us-disturbance)-based application.

You can download the Collection 1 dataset (1988–2022; 28 GB compressed, 80 GB uncompressed) from [this link](https://uconn-my.sharepoint.com/:u:/g/personal/shi_qiu_uconn_edu/EfcrNnvj2jpElWNygovkbcQBWRaBvFnQuvPCQfHujNSP-Q?e=9zj3qd).


## Contact US
Tian Li(tianli@uconn.edu) and Zhe Zhu (zhe@uconn.edu) at the Department of Natural Resources and the Environment, University of Connecticut.

## Reference

Li, T., Zhu, Z., Wang, Z., Kyba, C. C. M., Seto, K. C., Yang, Y., Qiu, S., Kuester, T., Fragkias, M., Chen, X., Meyer, T. H., Rittenhouse, C. D., Tai, X., Cullerton, M., Hong, F., Grinstead, A., Song, K., Suh, J. W., Yang, X., Kalb, V. L., Deng, C., & Román, M. O. Increasing volatility in human nighttime activity revealed by daily and high-resolution satellite imagery. Submitted to Nature.

Li, T., Zhu, Z., Wang, Z., Román, M. O., Kalb, V. L., & Zhao, Y. (2022). Continuous monitoring of nighttime light changes based on daily NASA's Black Marble product suite. Remote Sensing of Environment, 282, 113269.
