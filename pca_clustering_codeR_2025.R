# --- R Script for Oceanographic Data Analysis (PCA and Clustering) ---
# This script performs Principal Component Analysis (PCA) followed by K-means clustering
# on oceanographic data to identify distinct water mass assemblages and their characteristics.

# --- 1. Load Libraries and Setup Parallel Processing ---
# Loads the necessary R packages for data manipulation, analysis, visualization, and parallel computing.
library("FactoMineR") # For Principal Component Analysis (PCA)
library("factoextra") # For visualizing PCA results
library("zoo") # For time series analysis, specifically interpolation (na.approx)
library("readr") # For efficient data import (read_csv)
library("tidyr") # For data tidying, e.g., pivot_wider
library("tidyverse") # A collection of R packages (includes dplyr, ggplot2, tidyr, purrr)
library("future") # For setting up parallel processing
library("furrr") # For parallel processing with purrr-style functions
library("janitor") # For data cleaning, e.g., clean_names
library("gsw") # For Gibbs Seawater Oceanographic Toolbox functions
library('ggplot2')
library('RColorBrewer')
# Configures parallel processing to speed up computations.
# Uses all but one available CPU core for parallel tasks.
no_cores <- availableCores() - 1
future::plan(multisession, workers = no_cores)

# --- 2. Data Preparation and Pre-processing ---
# The goal is to produce a clean, standardized matrix ('ori_pca') where rows represent
# individual oceanographic profiles and columns represent scaled variables at specific depths.

# Loads your data file.
# Note: A portable file path is used here. For your specific path, replace the string below.
nutrient_data_df <- read.csv('path_to_folder/updated_dataframe.csv')

# Selects key variables for PCA and removes rows with missing values.
nut_prepca <- nutrient_data_df %>%
  dplyr::select(station, years, julian_day, depth_m, sa_g.kg, ct_degC, aou_µmol.kg, nstar_µmol.kg, sistar_µmol.kg) %>%
  drop_na(sa_g.kg, ct_degC, aou_µmol.kg, nstar_µmol.kg, sistar_µmol.kg) %>%
  # Renames columns for easier use in the script.
  rename(depth = depth_m, sa = sa_g.kg, ct = ct_degC, aou = aou_µmol.kg, nstar = nstar_µmol.kg, sistar = sistar_µmol.kg)

# Interpolates missing depth values within profiles and filters the data.
# The 'group_by' is now applied once to streamline the code.
prepca_df1 <- nut_prepca %>%
  group_by(years, julian_day, station) %>%
  arrange(depth) %>%
  mutate(depth = round(depth)) %>%
  # Ensures a continuous depth series for each profile before interpolation.
  complete(depth = min(depth):max(depth)) %>%
  # Uses linear interpolation to fill in any NA values.
  mutate(across(c(sa, ct, nstar, aou, sistar), ~ zoo::na.approx(., na.rm = TRUE))) %>%
  ungroup() %>%
  # Filters the data to the specific depth range of 10 to 250 meters.
  filter(depth <= 250, depth >= 10)

# Scales the selected variables (mean=0, std dev=1) for PCA to ensure equal weighting.
cols_to_scale <- c("sa", "ct", "nstar", "sistar", "aou")
prepca_df1_scaled <- prepca_df1 %>%
  mutate(across(all_of(cols_to_scale), scale)) %>%
  group_by(station, years, julian_day) %>%
  arrange(depth)

# Reshapes the data from a 'long' format to a 'wide' format for PCA.
# Each row becomes a single profile, and each column is a variable at a specific depth.
prepca_df2 <- prepca_df1_scaled %>%
  group_by(station, years, julian_day) %>%
  pivot_wider(
    names_from = depth,
    names_sep = "_",
    values_from = all_of(cols_to_scale),
    values_fn = mean
  )

# Final cleaning for PCA input: removes any remaining NAs and formats the data matrix.
ori_pca <- prepca_df2 %>%
  filter_all(all_vars(!is.na(.))) %>% # Removes rows with any remaining NAs
  mutate(profile_id = paste(station, years, julian_day, sep = "XXXX")) %>% # Creates a unique profile ID
  column_to_rownames(., var = "profile_id") %>% # Sets the profile ID as row names
  dplyr::select(-station, -years, -julian_day) # Removes the metadata columns

# 'ori_pca' is now the final prepared data matrix for PCA.

# --- 3. Principal Component Analysis (PCA) ---
# This section executes the PCA and extracts key results for interpretation.
res.pca <- PCA(ori_pca, scale.unit = TRUE, ncp = 10, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
print("Eigenvalues (Variance Explained by PCs):")
print(eig.val)
eigen_plot <- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100))
print(eigen_plot)
# Extracts and processes results for variables (loadings, contributions).
var_pca_results <- get_pca_var(res.pca)
cor_var_df <- as.data.frame(var_pca_results$cor) %>%
  rownames_to_column(., var = "var_depth") %>%
  separate(col = "var_depth", into = c("variable", "depth"), sep = "\\_") %>%
  janitor::clean_names() %>%
  mutate(depth = as.numeric(depth)) %>%
  pivot_wider(id_cols = "depth", names_from = variable, values_from = c(dim_1, dim_2, dim_3, dim_4, dim_5))
contrib_var_df <- as.data.frame(var_pca_results$contrib) %>%
  rownames_to_column(., var = "var_depth") %>%
  separate(col = "var_depth", into = c("variable", "depth"), sep = "\\_") %>%
  janitor::clean_names() %>%
  mutate(depth = as.numeric(depth)) %>%
  pivot_wider(id_cols = "depth", names_from = variable, values_from = c(dim_1, dim_2, dim_3, dim_4, dim_5))
# Extracts results for individual profiles in PCA space.
ind_pca_results <- get_pca_ind(res.pca)
coord_ind_df <- as.data.frame(ind_pca_results$coord) 

coord_ind_df2 <- as.data.frame(ind_pca_results$coord) %>%
  rownames_to_column(., var = "profile_id") %>%
  separate(profile_id, into = c("station", "years", "julian_day"), sep = "XXXX") %>%
  mutate(years = as.numeric(years))
# Describes the dimensions by identifying the most correlated variables, aiding in interpretation.
dim_descriptions <- dimdesc(res.pca, axes = c(1, 2, 3), proba = 0.05)
print("Description of Dimension 1 (Top Correlated Variables):")
print(as.data.frame(dim_descriptions$Dim.1) %>% top_n(20, wt = "correlation"))
print("Description of Dimension 2 (Top Correlated Variables):")
print(as.data.frame(dim_descriptions$Dim.2) %>% top_n(20, wt = "correlation"))
print("Description of Dimension 3 (Top Correlated Variables):")
print(as.data.frame(dim_descriptions$Dim.3) %>% top_n(20, wt = "correlation"))

# --- 4. Clustering Analysis on PCA Results ---
# This section applies K-means clustering to the PCA individual coordinates
# to identify distinct water mass groups or assemblages.
# Determines the optimal number of clusters using the 'Within Sum of Squares' (WSS) method.
# This plot helps identify the "elbow" point for a suitable 'k'.
fviz_nbclust(coord_ind_df[, c("Dim.1", "Dim.2", "Dim.3")],
             FUNcluster = kmeans, method = "wss", verbose = TRUE, k.max = 10, nboot = 1000)
# Performs K-means clustering with a chosen number of centers.
# The number of clusters (centers = 3) was likely determined by the fviz_nbclust step.
res.km <- kmeans(coord_ind_df[, c("Dim.1", "Dim.2", "Dim.3")], centers = 3, nstart = 100)
# Extracts cluster assignments and merges with profile metadata.
cluster_assignments_df <- as.data.frame(res.km$cluster) %>%
  rownames_to_column(., var = "profile_id") %>%
  separate(profile_id, into = c("station", "years", "julian_day"), sep = "XXXX") %>%
  janitor::clean_names() %>%
  mutate(grp = as.factor(res_km_cluster),
         years = as.numeric(years)) %>%
  dplyr::select(-res_km_cluster)
# Reorders cluster labels if a specific order is desired.
cluster_reordered_df <- cluster_assignments_df %>%
  mutate(grp = if_else(grp == "1", "2", if_else(grp == "2", "1", grp)))
# Joins reordered cluster assignments with PCA coordinates for visualization.
final_clustered_pca_data <- cluster_reordered_df %>%
  left_join(coord_ind_df2, by = c("station", "years", "julian_day"))
# Visualizes the PCA individuals colored by their assigned cluster.
# Create the scatter plot

set1_colors <- brewer.pal(n = 3, name = "Set1")

color_mapping <- c(
  "1" = set1_colors[2],  
  "2" = set1_colors[1],  
  "3" = set1_colors[3]  
)

pca_cluster_plot_ggplot <- ggplot(final_clustered_pca_data, aes(x = Dim.1, y = Dim.2)) +
  # Add data points, colored by the cluster group (grp)
  geom_point(aes(color = grp), size = 3) +
  # Add confidence ellipses around each cluster
  stat_ellipse(aes(color = grp), geom = "polygon", alpha = 0.2, linetype = "dashed") +
  # Customize titles and labels
  labs(
    title = "PCA of Oceanographic Data Profiles",
    subtitle = "Individuals colored by K-means cluster",
    x = "PC1", # Adjust the axis title if you know the explained variance
    y = "PC2", # Adjust the axis title if you know the explained variance
    color = "Cluster"
  ) +
  # Use a more readable color palette
  scale_color_manual(values = color_mapping) +
  # Apply a minimalist theme for better readability
  theme_minimal() +
  # Adjust legend position
  theme(legend.position = "bottom")

# Display the plot
print(pca_cluster_plot_ggplot)

