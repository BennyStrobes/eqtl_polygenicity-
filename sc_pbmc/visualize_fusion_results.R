args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')



figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}




load_in_fusion_results <- function(data_set_names, num_clusters, fusion_processed_output_dir) {
	num_data_sets <- length(data_set_names)

	gene_name_arr <- c()
	h2_arr <- c()
	h2_err_arr <-c()
	h2_p_arr <- c()
	fraction_non_zero_weights_arr <- c()
	data_set_names_arr <- c()
	data_set_cluster_size_arr <- c()

	for (data_set_num in 1:num_data_sets) {
		data_set_name <- data_set_names[data_set_num]
		num_cluster <- num_clusters[data_set_num]

		# load in data set info
		data_set_file_name <- paste0(fusion_processed_output_dir, data_set_name, "_1KG_only_fusion_results_summary.txt")
		data_set_info <- read.table(data_set_file_name, header=TRUE, sep="\t")

		gene_name_arr <- c(gene_name_arr, as.character(data_set_info$gene_name))
		h2_arr <- c(h2_arr, data_set_info$h2)
		h2_err_arr <- c(h2_err_arr, data_set_info$h2_err)
		h2_p_arr <- c(h2_p_arr, data_set_info$h2_p)
		fraction_non_zero_weights_arr <- c(fraction_non_zero_weights_arr, data_set_info$fraction_non_zero_weights)
		data_set_names_arr <- c(data_set_names_arr, rep(data_set_name, length(data_set_info$h2_err)))
		data_set_cluster_size_arr <- c(data_set_cluster_size_arr, rep(num_cluster, length(data_set_info$h2_err)))
	}
	df <- data.frame(h2=h2_arr, h2_err=h2_err_arr, h2_p=h2_p_arr, gene_name=gene_name_arr, fraction_non_zero_weights=fraction_non_zero_weights_arr, data_set_name=factor(data_set_names_arr, levels=data_set_names), data_set_cluster_resolution=data_set_cluster_size_arr)
	return(df)
}

create_boxplot_comparing_fraction_of_non_zero_weights <- function(df) {

  boxplot <- ggplot(df, aes(x=data_set_name, y=fraction_non_zero_weights, fill=data_set_cluster_resolution)) + geom_boxplot() +
            figure_theme() + 
            labs(x="Data Set name", y = paste0("Fraction non-zero lasso weights/gene"), fill="cluster resolution") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  return(boxplot)

}

create_boxplot_comparing_h2_g <- function(df) {
	  boxplot <- ggplot(df, aes(x=data_set_name, y=h2, fill=data_set_cluster_resolution)) + geom_boxplot() +
            figure_theme() + 
            labs(x="Data Set name", y = paste0("h2/gene"), fill="cluster resolution") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  return(boxplot)
}

create_scatterplot_comparing_h2_and_fraction_non_zero_snps <- function(df) {
	corry = cor(df$h2, df$fraction_non_zero_weights, method="spearman")
	scatter <- ggplot(df) + geom_point(aes(x=h2,y=fraction_non_zero_weights,color=data_set_cluster_resolution), size=.01) +
		figure_theme() +
		labs(x="h2", y="Fraction non-zero lasso weights", color="cluster resolution", title=paste0("spearman correlation: ", corry))

	return(scatter)

}

create_scatterplot_comparing_h2_err_and_fraction_non_zero_snps <- function(df) {
	corry = cor(df$h2_err, df$fraction_non_zero_weights, method="spearman")
	scatter <- ggplot(df) + geom_point(aes(x=h2_err,y=fraction_non_zero_weights,color=data_set_cluster_resolution), size=.01) +
		figure_theme() +
		labs(x="h2 standard error", y="Fraction non-zero lasso weights", color="cluster resolution", title=paste0("spearman correlation: ", corry))

	return(scatter)

}

create_scatter_comparing_mean_fraction_non_zero_and_avg_cells_per_individual <- function(df) {
	scatter <- ggplot(df) + geom_point(aes(x=avg_cells_per_individual,y=mean_fraction_non_zero_weights,color=data_set_cluster_resolution)) +
		figure_theme() +
		labs(x="Mean(cells/individual)", y="Mean(fraction non-zero weights)", color="cluster resolution")

	return(scatter)	
}

create_scatter_comparing_mean_fraction_non_zero_and_sd_cells_per_individual <- function(df) {
	scatter <- ggplot(df) + geom_point(aes(x=sd_cells_per_individual,y=mean_fraction_non_zero_weights,color=data_set_cluster_resolution)) +
		figure_theme() +
		labs(x="stdev(cells/individual)", y="Mean(fraction non-zero weights)", color="cluster resolution")

	return(scatter)	
}
create_scatter_comparing_mean_fraction_non_zero_and_mean_h2 <- function(df) {
	scatter <- ggplot(df) + geom_point(aes(x=mean_h2,y=mean_fraction_non_zero_weights,color=data_set_cluster_resolution)) +
		figure_theme() +
		labs(x="Mean(h^2)", y="Mean(fraction non-zero weights)", color="cluster resolution")

	return(scatter)	
}


get_df_for_sem <- function(df, full_df, data_set_names, num_clusters, avg_cells_per_individual, sd_cells_per_individual) {
	mean_h2_arr <- c()
	mean_h2_se_arr <- c()
	mean_nonzero_arr <-c()
	mean_nonzero_se_arr <- c()
	num_genes_arr <- c()

	for (data_set_iter in 1:length(data_set_names)) {
		data_set_name <- data_set_names[data_set_iter]
		data_set_indices <- as.character(df$data_set_name) == data_set_name
		data_set_df <- df[data_set_indices,]

		full_data_set_indices <- as.character(full_df$data_set_name) == data_set_name
		full_data_set_df <- full_df[full_data_set_indices,]

		mean_h2 <- mean(data_set_df$h2)
		mean_h2_se <- sd(data_set_df$h2)/sqrt(length(data_set_df$h2))

		mean_fraction_nonzero <- mean(data_set_df$fraction_non_zero_weights)
		mean_fraction_nonzero_se <- sd(data_set_df$fraction_non_zero_weights)/sqrt(length(data_set_df$fraction_non_zero_weights))


		mean_h2_arr <- c(mean_h2_arr, mean_h2)
		mean_h2_se_arr <- c(mean_h2_se_arr, mean_h2_se)

		mean_nonzero_arr <- c(mean_nonzero_arr, mean_fraction_nonzero)
		mean_nonzero_se_arr <- c(mean_nonzero_se_arr, mean_fraction_nonzero_se)

		num_genes_arr <- c(num_genes_arr, length(full_data_set_df$h2))
	}
	df <- data.frame(num_genes=num_genes_arr,data_set_name=factor(data_set_names,levels=data_set_names),data_set_cluster_resolution=factor(num_clusters), avg_cells_per_individual=avg_cells_per_individual, sd_cells_per_individual=sd_cells_per_individual, mean_h2=mean_h2_arr, mean_h2_se=mean_h2_se_arr, mean_fraction_non_zero_weights=mean_nonzero_arr, mean_fraction_non_zero_weights_se=mean_nonzero_se_arr)
	return(df)
}



create_sem_barplot_comparing_fraction_of_non_zero_weights <- function(df) {
	p <- ggplot(data=df, aes(x=data_set_name, y=mean_fraction_non_zero_weights, fill=data_set_cluster_resolution)) +
		geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
		figure_theme() +
		theme(legend.position="top") +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
		labs(y="Mean(Fraction non-zero gene weights)", fill="Cluster resolution") +
		geom_errorbar(aes(ymin=mean_fraction_non_zero_weights-(1.96*mean_fraction_non_zero_weights_se), ymax=mean_fraction_non_zero_weights+(1.96*mean_fraction_non_zero_weights_se)), position = position_dodge(), width = .75, size=.2)

	return(p)
}

create_sem_barplot_comparing_h2_g <- function(df) {
	p <- ggplot(data=df, aes(x=data_set_name, y=mean_h2, fill=data_set_cluster_resolution)) +
		geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
		figure_theme() +
		theme(legend.position="top") +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
		labs(y="Mean(h2)", fill="Cluster resolution") +
		geom_errorbar(aes(ymin=mean_h2-(1.96*mean_h2_se), ymax=mean_h2+(1.96*mean_h2_se)), position = position_dodge(), width = .75, size=.2)

	return(p)
}

get_average_cells_per_individual_in_each_data_set <- function(data_set_names, pseudobulk_data_dir) {
	avg_cells_per_indi_arr <- c()
	for (data_set_iter in 1:length(data_set_names)) {
		data_set_name <- data_set_names[data_set_iter]

		cells_per_indi_data_set_file <- paste0(pseudobulk_data_dir, data_set_name, "_num_cells_per_individual.txt")
		cells_per_indi_df <- read.table(cells_per_indi_data_set_file, header=TRUE, sep="\t")

		avg_cells_per_indi_in_data_set <- mean(cells_per_indi_df$Number_of_cells)
		avg_cells_per_indi_arr <- c(avg_cells_per_indi_arr, avg_cells_per_indi_in_data_set)
	}
	return(avg_cells_per_indi_arr)
}

get_stdev_cells_per_individual_in_each_data_set <- function(data_set_names, pseudobulk_data_dir) {
	avg_cells_per_indi_arr <- c()
	for (data_set_iter in 1:length(data_set_names)) {
		data_set_name <- data_set_names[data_set_iter]

		cells_per_indi_data_set_file <- paste0(pseudobulk_data_dir, data_set_name, "_num_cells_per_individual.txt")
		cells_per_indi_df <- read.table(cells_per_indi_data_set_file, header=TRUE, sep="\t")

		avg_cells_per_indi_in_data_set <- sd(cells_per_indi_df$Number_of_cells)
		avg_cells_per_indi_arr <- c(avg_cells_per_indi_arr, avg_cells_per_indi_in_data_set)
	}
	return(avg_cells_per_indi_arr)
}


############################
# Command line args
############################
data_set_summary_file <- args[1]
fusion_processed_output_dir <- args[2]
pseudobulk_data_dir <- args[3]
fusion_visualization_dir <- args[4]

# load in data set names
ordered_data_sets <- c("leiden_0_0", "leiden_0.001_0", "leiden_0.001_1", "leiden_0.001_2", "perez_et_al__T4", "perez_et_al__cM", "perez_et_al__T8", "perez_et_al__B", "perez_et_al__NK", "perez_et_al__ncM", "perez_et_al__cDC", "perez_et_al__pDC", "perez_et_al__Prolif")
data_set_df <- read.table(data_set_summary_file, header=TRUE, sep="\t")
data_set_df$data_set_name = factor(data_set_df$data_set_name, levels=ordered_data_sets)
data_set_names <- as.character(data_set_df$data_set_name)
num_clusters <- as.character(c(1, 3, 3, 3, 9, 9, 9, 9,9,9,9,9,9))

avg_cells_per_individual <- get_average_cells_per_individual_in_each_data_set(data_set_names, pseudobulk_data_dir)
sd_cells_per_individual <- get_stdev_cells_per_individual_in_each_data_set(data_set_names, pseudobulk_data_dir)

# Load in fusion results
fusion_results_full_df <- load_in_fusion_results(data_set_names, num_clusters, fusion_processed_output_dir)
fusion_results_df <- fusion_results_full_df[!is.na(fusion_results_full_df$fraction_non_zero_weights),]


fusion_results_sem_df <- get_df_for_sem(fusion_results_df,fusion_results_full_df, data_set_names, num_clusters, avg_cells_per_individual, sd_cells_per_individual)
fusion_results_sem_df$data_set_name = factor(fusion_results_sem_df$data_set_name, levels=ordered_data_sets)

print(fusion_results_sem_df)

# Create SEM barplot comparing distributions of mean fraction of non-zero weights
output_file <- paste0(fusion_visualization_dir, "sem_barplot_comparing_fraction_of_non_zero_weights.pdf")
sem_barplot_non_zero <- create_sem_barplot_comparing_fraction_of_non_zero_weights(fusion_results_sem_df)
ggsave(sem_barplot_non_zero, file=output_file, width=7.2, height=6.0, units="in")

# Create SEM barplot comparing distributions of h2
output_file <- paste0(fusion_visualization_dir, "sem_barplot_comparing_h2_g_for_sig_genes.pdf")
sem_barplot_h2_g <- create_sem_barplot_comparing_h2_g(fusion_results_sem_df)
ggsave(sem_barplot_h2_g, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(fusion_visualization_dir, "joint_sem_barplot.pdf")
joint_sem_barplot <- plot_grid(sem_barplot_h2_g + theme(legend.position="none"), sem_barplot_non_zero+ theme(legend.position="none"), ncol=1)
ggsave(joint_sem_barplot, file=output_file, width=7.2, height=8.5, units="in")


# Create scatter plot comparing mean_fraction_non_zero_weights with avg_cells_per_individual
output_file <- paste0(fusion_visualization_dir, "scatter_mean_fraction_non_zero_and_avg_cells_per_individual.pdf")
scatter <- create_scatter_comparing_mean_fraction_non_zero_and_avg_cells_per_individual(fusion_results_sem_df)
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

# Create scatter plot comparing mean_fraction_non_zero_weights with avg_cells_per_individual
output_file <- paste0(fusion_visualization_dir, "scatter_mean_fraction_non_zero_and_sd_cells_per_individual.pdf")
scatter <- create_scatter_comparing_mean_fraction_non_zero_and_sd_cells_per_individual(fusion_results_sem_df)
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

# Create scatter plot comparing mean_fraction_non_zero_weights with avg_cells_per_individual
output_file <- paste0(fusion_visualization_dir, "scatter_mean_fraction_non_zero_and_mean_h2.pdf")
scatter <- create_scatter_comparing_mean_fraction_non_zero_and_mean_h2(fusion_results_sem_df)
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

# Create boxplot comparing distributions of fraction of non-zero weights
output_file <- paste0(fusion_visualization_dir, "boxplot_comparing_fraction_of_non_zero_weights.pdf")
boxplot <- create_boxplot_comparing_fraction_of_non_zero_weights(fusion_results_df)
ggsave(boxplot, file=output_file, width=7.2, height=6.0, units="in")

# Create boxplot comparing distributions of h2_g
output_file <- paste0(fusion_visualization_dir, "boxplot_comparing_h2_g_for_all_genes.pdf")
boxplot <- create_boxplot_comparing_h2_g(fusion_results_full_df)
ggsave(boxplot, file=output_file, width=7.2, height=6.0, units="in")
output_file <- paste0(fusion_visualization_dir, "boxplot_comparing_h2_g_for_sig_genes.pdf")
boxplot <- create_boxplot_comparing_h2_g(fusion_results_df)
ggsave(boxplot, file=output_file, width=7.2, height=6.0, units="in")

# Scatter plot comparing h2 and fraction_of_non_zero_snps
output_file <- paste0(fusion_visualization_dir, "scatterplot_comparing_h2_and_fraction_non_zero_snps.pdf")
scatter <- create_scatterplot_comparing_h2_and_fraction_non_zero_snps(fusion_results_df)
ggsave(scatter, file=output_file, width=7.2, height=3.0, units="in")

# Scatter plot comparing h2 and fraction_of_non_zero_snps
output_file <- paste0(fusion_visualization_dir, "scatterplot_comparing_h2_err_and_fraction_non_zero_snps.pdf")
scatter <- create_scatterplot_comparing_h2_err_and_fraction_non_zero_snps(fusion_results_df)
ggsave(scatter, file=output_file, width=7.2, height=3.0, units="in")
