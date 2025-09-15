# water-mass-clustering
This repository contains the R script used for a comprehensive analysis of oceanographic data. The script provides a full pipeline for data processing, dimensionality reduction via Principal Component Analysis (PCA), and subsequent K-means clustering to identify distinct water mass assemblages.

Features
Data Pre-processing: Cleans, interpolates, and standardizes raw oceanographic profiles.

Dimensionality Reduction: Performs Principal Component Analysis (PCA) to reduce the data's complexity.

Clustering Analysis: Applies K-means clustering on the PCA results to group similar water profiles.

Visualization: Generates high-quality plots to visualize cluster distributions in the PCA space.

Prerequisites
To run this script, you need R installed. The following R packages are required and can be installed from CRAN:

FactoMineR

factoextra

tidyverse (includes dplyr, ggplot2, tidyr, etc.)

zoo

readr

future

furrr

janitor

gsw

RColorBrewer

Usage
Download the Data: Obtain the updated_dataframe.csv file from the Zenodo repository linked in the Data Availability section below.

Place the Data: Save the data file in the specified path in the script.

Run the Script: Open and execute the water-mass-clustering.R script in an R environment.

Data Availability
The nutrient data for this study are available in a Zenodo repository (DOI: https://doi.org/10.5281/zenodo.17106712). Seawater measurements of temperature, salinity, and DO presented in this document were made available by the Amundsen Science program through their data access portal. Seawater measurements of temperature, salinity, and DO presented in this document were made available by the Amundsen Science program through their data access portal (https://amundsenscience.com/data/data-access/). The views expressed in this publication do not necessarily represent the views of Amundsen Science or its partners. 


License
This project is licensed under the MIT License - see the LICENSE.md file for details.
