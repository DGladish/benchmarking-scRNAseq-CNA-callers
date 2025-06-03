# install dependencies
if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")

# import dependencies

# import all of the rds files
## obtain all of the paths
mwa_out_dir <- "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/mwa_script_output"
mwa_out_rds_paths <- list.files(mwa_out_dir, full.names=TRUE)

## do the importing
mwa_outs <- list()

for (i in 1:length(mwa_out_rds_paths)) {
rds <- readRDS(mwa_out_rds_paths[i])
mwa_outs[[i]] <- rds
basename <- basename(mwa_out_rds_paths[i])
element_name <- sub("_mwa_outs.rds", "", basename)

names(mwa_outs)[i] <- element_name
}

# check the success of the import / integrity of the list
if (!any(duplicated(mwa_outs))) {
    print("None of the RDS files are duplicates of one another :)")
} else {
    stop("Stopping execution, at least 1 rds file has been duplicated.")
}


# summarizing
## merge element-wise
mwa_outs_trans <- transpose(mwa_outs) # NTS: this creates a dataset duplicate

ascat_windows <- unlist(mwa_outs_trans[[1]])
scevan_windows <- unlist(mwa_outs_trans[[2]])
distances <- unlist(mwa_outs_trans[[3]])
differences <- unlist(mwa_outs_trans[[4]])
# add distances summary
# add differences summary



## clean
### remove NA values
ascat_windows <- ascat_windows[!is.na(ascat_windows)]
scevan_windows <- scevan_windows[!is.na(scevan_windows)]
distances <- distances[!is.na(distances)]
differences <- differences[!is.na(differences)]

### eliminate sample identifiers from the vector names
#### obtain element names with and without the sample identifiers if I should need them
ascat_windows_names_nosampID <- sub(".*(chr.*)", "\\1", names(ascat_windows))
scevan_windows_names_nosampID <- sub(".*(chr.*)", "\\1", names(scevan_windows))
distances_names_nosampID <- sub(".*(chr.*)", "\\1", names(distances))
differences_names_nosampID <- sub(".*(chr.*)", "\\1", names(differences))

ascat_windows_names_withsampID <- names(ascat_windows)
scevan_windows_names_withsampID <- names(scevan_windows)
distances_names_withsampID <- names(distances)
differences_names_withsampID <- names(differences)

## now analysis!
### summarize similarity metrics
#### convert the distance and difference vectors to data frames (required by ggplot)
dists_df <- data.frame(distances = distances)
diffs_df <- data.frame(differences = differences)


#### distances
##### histogram
dist_h <- ggplot(dists_df, aes(x = distances)) +
    geom_histogram(binwidth=1) + 
    #coord_cartesian(xlim = c(0, 35)) + 
    ggtitle("Distribution of element-wise distances between ASCAT & SCEVAN 1MB Windows")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_histogram_distances.png", plot = dist_h)

dist_hz <- ggplot(dists_df, aes(x = distances)) +
    geom_histogram(binwidth=1) + 
    coord_cartesian(xlim = c(0, 35)) + 
    ggtitle("Distribution of element-wise distances between ASCAT & SCEVAN 1MB Windows (zoomed)")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_histogram_distances_zoomed.png", plot = dist_hz)


##### boxplot
dist_bp <- ggplot(dists_df, aes(y=distances)) + 
  geom_boxplot(fill="tomato") + 
  ggtitle("Boxplot of distances between copy number values")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_boxplot_distances.png", plot = dist_bp)


dist_bpz <- ggplot(dists_df, aes(y=distances)) + 
  geom_boxplot(fill="tomato") + 
  coord_cartesian(ylim = c(0, 35)) + 
  ggtitle("Boxplot of distances between copy number values (Zoomed)")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_boxplot_distances_zoomed.png", plot = dist_bpz)


#### differences
##### histogram
diff_h <- ggplot(diffs_df, aes(x = differences)) +
    geom_histogram(binwidth=1) + 
    #coord_cartesian(xlim = c(0, 35)) + 
    ggtitle("Distribution of element-wise differences between ASCAT & SCEVAN 1MB Windows")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_histogram_differences.png", plot = diff_h)

diff_hz <- ggplot(diffs_df, aes(x = differences)) +
    geom_histogram(binwidth=1) + 
    coord_cartesian(xlim = c(-10, 10)) + 
    ggtitle("Distribution of element-wise differences between ASCAT & SCEVAN 1MB Windows (zoomed)")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_histogram_differences_zoomed.png", plot = diff_hz)


##### box plot
diff_bp <- ggplot(diffs_df, aes(y=differences)) + 
  geom_boxplot(fill="blue") + 
  ggtitle("Distribution of element-wise differences between ASCAT & SCEVAN 1MB Windows")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_boxplot_differences.png", plot = diff_bp)


diff_bpz <- ggplot(diffs_df, aes(y=differences)) + 
  geom_boxplot(fill="blue") + 
  coord_cartesian(ylim = c(-10, 10)) + 
  ggtitle("Distribution of element-wise differences between ASCAT & SCEVAN 1MB Windows (Zoomed)")

ggsave(filename = "/scratch/dgladish/benchmarking_spatial_scale/analysis_outs/single_parameter/dist_summary/plots/single_param_boxplot_differences_zoomed.png", plot = diff_bpz)
