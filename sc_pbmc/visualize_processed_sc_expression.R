library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}

make_number_of_cells_per_individual_bar_plot <- function(filtered_covariate_data) {
  unique_indis <- as.character(unique(filtered_covariate_data$ind_cov))
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(filtered_covariate_data$ind_cov == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
  ordering <- order(cells_arr)
  df <- data.frame(num_cells=cells_arr[ordering], indi=factor(unique_indis[ordering], levels=unique_indis[ordering]))

  p<-ggplot(data=df, aes(x=indi, y=num_cells)) +
    geom_bar(stat="identity", fill="steelblue")+
    figure_theme() +
    labs(x = "Individual", y = "Number of cells") +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  return(p) 
}

make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label) {
	# Put everthing in compact data frame
	unique_categories = as.character(sort(unique(categorical_variable)))

	df <- data.frame(categorical_variable=factor(categorical_variable, levels=unique_categories), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	

	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				guides(colour = guide_legend(override.aes = list(size=3)))
  	
  	
  	return(scatter)
}

make_pc_variance_explained_line_plot <- function(variance_explained, num_pcs) {
	variance_explained <- variance_explained[1:num_pcs]
	df <- data.frame(variance_explained = variance_explained, pc_num = 1:num_pcs)

	# PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + .002) + 
                scale_x_continuous(breaks=seq(0,(num_pcs-1),5)) +
                labs(x = "PC number", y = "Variance Explained") + 
                figure_theme() 

    return(line_plot)
}

boxplot_of_number_of_cells_per_individual_in_each_cluster <- function(cluster_assignments, filtered_covariate_data, cluster_assignment_name) {
	unique_cluster_assignments = as.character(sort(unique(cluster_assignments)))
	unique_indis <- as.character(unique(filtered_covariate_data$ind_cov))

	cluster_assignment_arr <- c()
	cells_per_indi_arr <- c()

	for (cluster_iter in 1:length(unique_cluster_assignments)) {
		cluster_assignment = as.character(unique_cluster_assignments[cluster_iter])
		cluster_assignment_covariate_data = filtered_covariate_data[as.character(cluster_assignments)==cluster_assignment,]
		for (indi_iter in 1:length(unique_indis)) {
			indi_name = unique_indis[indi_iter]
			num_cells = sum(cluster_assignment_covariate_data$ind_cov == indi_name)
			cluster_assignment_arr <- c(cluster_assignment_arr, cluster_assignment)
			cells_per_indi_arr <- c(cells_per_indi_arr, num_cells)
		}

	}

	df <- data.frame(cells_per_individual=cells_per_indi_arr, cluster_assignment=factor(cluster_assignment_arr, levels=unique_cluster_assignments))
	p <- ggplot(df, aes(x=cluster_assignment, y=cells_per_individual, fill=cluster_assignment)) + 
  		geom_boxplot() +
  		figure_theme() + 
  		theme(legend.position="none") +
      labs(x="Cluster assignment", y="Cells / individual") +
      geom_hline(yintercept = 5)
  	return(p)
}

#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_expression_dir <- args[1]  # Input dir
visualize_processed_expression_dir <- args[2]  # Output Dir

processed_expression_dir <- paste0(processed_expression_dir, "scran_normalization_hvg_2000_regress_batch_True_")




##########################
# Load in Data
##########################
# Load in cluster assignments
cluster_assignment_file <- paste0(processed_expression_dir, "cluster_assignments_at_various_resolutions.txt")
#cluster_assignment_data <- read.table(cluster_assignment_file, header=TRUE, sep="\t")
#saveRDS(cluster_assignment_data, "cluster_assignments.rds")
cluster_assignment_data <- readRDS("cluster_assignments.rds")


# Load in Covariates
filtered_covariate_file <- paste0(processed_expression_dir, "cell_covariates.txt")
#filtered_covariate_data <- read.table(filtered_covariate_file, header=TRUE, sep="\t")
#saveRDS(filtered_covariate_data, "cov.rds")
filtered_covariate_data <- readRDS("cov.rds")

# Load in PCS
pc_file <- paste0(processed_expression_dir, "sc_expression_pcs.txt")
#pcs <- read.table(pc_file, header=FALSE, sep="\t")
#saveRDS(pcs, "pcs.rds")
pcs <- readRDS("pcs.rds")

# Load in PC PVE
pc_pve_file <- paste0(processed_expression_dir, "sc_expression_pcs_percent_variance_explained.txt")
pc_pve <- read.table(pc_pve_file, header=FALSE, sep="\t")

# Load in umap_loadings
umap_file <- paste0(processed_expression_dir, "sc_expression_umaps.txt")
#umap_loadings <- read.table(umap_file, header=FALSE, sep="\t")
#saveRDS(umap_loadings, "umap.rds")
umap_loadings <- readRDS("umap.rds")




##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(visualize_processed_expression_dir, "pca_variance_explained_", num_pcs, "_pcs_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pc_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")



##########################
# Make bar plot showing number of cells per individual
##########################
num_cells_bar_plot <- make_number_of_cells_per_individual_bar_plot(filtered_covariate_data)
output_file <- paste0(visualize_processed_expression_dir, "number_of_cells_per_individual_bar_plot.pdf")
ggsave(num_cells_bar_plot, file=output_file, width=7.2, height=5, units="in")


##########################
# Make UMAP Plot colored by cluster assignments
##########################
umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(factor(cluster_assignment_data$leiden_0.001), umap_loadings[,1], umap_loadings[,2], "Leiden cluster 0.001", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_leiden_cluster_.001.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")
# boxplot of cells per individual in each cluster
boxplot_of_cells_per_indi_in_this_cluster <- boxplot_of_number_of_cells_per_individual_in_each_cluster(factor(cluster_assignment_data$leiden_0.001), filtered_covariate_data, "Leiden cluster 0.001")
output_file <- paste0(visualize_processed_expression_dir, "boxplot_of_number_of_cells_per_indi_in_cluster_leiden_cluster_.001.pdf")
ggsave(boxplot_of_cells_per_indi_in_this_cluster, file=output_file, width=7.2, height=5, units="in")
# Joint cowplot
joint <- plot_grid(umap_scatter_colored_by_cell_type + theme(legend.position="bottom"), boxplot_of_cells_per_indi_in_this_cluster, ncol=1, rel_heights=c(1,.45))
output_file <- paste0(visualize_processed_expression_dir, "joint_leiden_cluster_.001.pdf")
ggsave(joint, file=output_file, width=7.2, height=6, units="in")


umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(factor(cluster_assignment_data$leiden_0.01), umap_loadings[,1], umap_loadings[,2], "Leiden cluster 0.01", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_leiden_cluster_.01.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")
# boxplot of cells per individual in each cluster
boxplot_of_cells_per_indi_in_this_cluster <- boxplot_of_number_of_cells_per_individual_in_each_cluster(factor(cluster_assignment_data$leiden_0.01), filtered_covariate_data, "Leiden cluster 0.01")
output_file <- paste0(visualize_processed_expression_dir, "boxplot_of_number_of_cells_per_indi_in_cluster_leiden_cluster_.01.pdf")
ggsave(boxplot_of_cells_per_indi_in_this_cluster, file=output_file, width=7.2, height=5, units="in")
# Joint cowplot
joint <- plot_grid(umap_scatter_colored_by_cell_type + theme(legend.position="bottom"), boxplot_of_cells_per_indi_in_this_cluster, ncol=1, rel_heights=c(1,.45))
output_file <- paste0(visualize_processed_expression_dir, "joint_leiden_cluster_.01.pdf")
ggsave(joint, file=output_file, width=7.2, height=6, units="in")

umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(factor(cluster_assignment_data$leiden_0.1), umap_loadings[,1], umap_loadings[,2], "Leiden cluster 0.1", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_leiden_cluster_.1.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")
# boxplot of cells per individual in each cluster
boxplot_of_cells_per_indi_in_this_cluster <- boxplot_of_number_of_cells_per_individual_in_each_cluster(factor(cluster_assignment_data$leiden_0.1), filtered_covariate_data, "Leiden cluster 0.1")
output_file <- paste0(visualize_processed_expression_dir, "boxplot_of_number_of_cells_per_indi_in_cluster_leiden_cluster_.1.pdf")
ggsave(boxplot_of_cells_per_indi_in_this_cluster, file=output_file, width=7.2, height=5, units="in")
# Joint cowplot
joint <- plot_grid(umap_scatter_colored_by_cell_type + theme(legend.position="bottom"), boxplot_of_cells_per_indi_in_this_cluster, ncol=1, rel_heights=c(1,.45))
output_file <- paste0(visualize_processed_expression_dir, "joint_leiden_cluster_.1.pdf")
ggsave(joint, file=output_file, width=7.2, height=6, units="in")

umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(factor(cluster_assignment_data$leiden_1), umap_loadings[,1], umap_loadings[,2], "Leiden cluster 1", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_leiden_cluster_1.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")
# boxplot of cells per individual in each cluster
boxplot_of_cells_per_indi_in_this_cluster <- boxplot_of_number_of_cells_per_individual_in_each_cluster(factor(cluster_assignment_data$leiden_1), filtered_covariate_data, "Leiden cluster 1")
output_file <- paste0(visualize_processed_expression_dir, "boxplot_of_number_of_cells_per_indi_in_cluster_leiden_cluster_1.pdf")
ggsave(boxplot_of_cells_per_indi_in_this_cluster, file=output_file, width=7.2, height=5, units="in")
# Joint cowplot
joint <- plot_grid(umap_scatter_colored_by_cell_type + theme(legend.position="bottom"), boxplot_of_cells_per_indi_in_this_cluster, ncol=1, rel_heights=c(1,.45))
output_file <- paste0(visualize_processed_expression_dir, "joint_leiden_cluster_1.pdf")
ggsave(joint, file=output_file, width=7.2, height=6, units="in")

umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(factor(cluster_assignment_data$leiden_3), umap_loadings[,1], umap_loadings[,2], "Leiden cluster 3", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_leiden_cluster_3.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")
# boxplot of cells per individual in each cluster
boxplot_of_cells_per_indi_in_this_cluster <- boxplot_of_number_of_cells_per_individual_in_each_cluster(factor(cluster_assignment_data$leiden_3), filtered_covariate_data, "Leiden cluster 3")
output_file <- paste0(visualize_processed_expression_dir, "boxplot_of_number_of_cells_per_indi_in_cluster_leiden_cluster_3.pdf")
ggsave(boxplot_of_cells_per_indi_in_this_cluster, file=output_file, width=7.2, height=5, units="in")
# Joint cowplot
joint <- plot_grid(umap_scatter_colored_by_cell_type + theme(legend.position="bottom"), boxplot_of_cells_per_indi_in_this_cluster, ncol=1, rel_heights=c(1,.45))
output_file <- paste0(visualize_processed_expression_dir, "joint_leiden_cluster_3.pdf")
ggsave(joint, file=output_file, width=7.2, height=6, units="in")

##########################
# Make UMAP Plot colored by cell type
##########################
umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$cg_cov, umap_loadings[,1], umap_loadings[,2], "Cell Type", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_cell_type.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

##########################
# Make UMAP Plot colored by status
##########################
umap_scatter_colored_by_batch <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$Status, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_status.pdf")
ggsave(umap_scatter_colored_by_batch, file=output_file, width=7.2, height=5, units="in")

##########################
# Make UMAP Plot colored by Lupus
##########################
umap_scatter_colored_by_lupus <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$SLE_status, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(visualize_processed_expression_dir, "umap_1_2_scatter_colored_by_disease_cov.pdf")
ggsave(umap_scatter_colored_by_lupus, file=output_file, width=7.2, height=5, units="in")





