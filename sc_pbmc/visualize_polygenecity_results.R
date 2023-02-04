args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


extract_kurtosis_df <- function(gene_heritability_normalization, cell_types, input_dir) {
	iteration="600"
	kurtosis_arr <- c()
	cell_type_arr <- c()
	subset_arr <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type = cell_types[cell_type_iter]
		real_file_name <- paste0(input_dir, "perez_et_al__", cell_type, "_real_polygenecity_results_non_negative_h2_gene_h2_norm_", gene_heritability_normalization, "_temp_", iteration, "_expected_kurtosis.txt")

		temp_df <- read.table(real_file_name, header=FALSE)
		kurtosis <- temp_df$V1[2]

		kurtosis_arr <- c(kurtosis_arr, kurtosis)
		cell_type_arr <- c(cell_type_arr, cell_type)
		subset_arr <- c(subset_arr, "real")

		for (subsample_iter in 0:9) {
			subsample_file_name <- paste0(input_dir, "perez_et_al__", cell_type, "_null_", subsample_iter, "_polygenecity_results_non_negative_h2_gene_h2_norm_", gene_heritability_normalization, "_temp_", iteration, "_expected_kurtosis.txt")
			temp_df <- read.table(subsample_file_name, header=FALSE)
			kurtosis <- temp_df$V1[2]

			kurtosis_arr <- c(kurtosis_arr, kurtosis)
			cell_type_arr <- c(cell_type_arr, cell_type)
			subset_arr <- c(subset_arr, "subsample_full")
		}
	}

	df <- data.frame(cell_type=factor(cell_type_arr, levels=cell_types), kurtosis=kurtosis_arr, subset=factor(subset_arr, levels=c("real", "subsample_full")))
}

make_cell_type_specific_boxplot <- function(ct_kurtosis_df, cell_type) {
	df2 = ct_kurtosis_df[ct_kurtosis_df$subset=="subsample_full",]
	real_kurtosis = ct_kurtosis_df[ct_kurtosis_df$subset!="subsample_full",]$kurtosis[1]
	p<-ggplot(df2, aes(x=cell_type, y=kurtosis)) +
 	 geom_boxplot() +
 	 figure_theme() +
 	 theme(legend.position="bottom") +
 	 geom_hline(yintercept=real_kurtosis, linetype="dashed", color = "red") +
 	 labs(x="", title=cell_type)
	return(p)
}

make_cell_type_subsampled_boxplot <- function(kurtosis_df, cell_types) {
	plot_list <- list()
	for (cell_type_iter in 1:length(cell_types)) {

		cell_type <- cell_types[cell_type_iter]

		ct_kurtosis_df <- kurtosis_df[kurtosis_df$cell_type==cell_type,]

		plot_list[[cell_type_iter]] <- make_cell_type_specific_boxplot(ct_kurtosis_df, cell_type)
	}

	p <- plot_grid(plotlist=plot_list, ncol=4)

	return(p)
}





visualization_dir <- args[1]
input_dir <- args[2]

cell_types <- c("B", "cDC", "cM", "ncM", "NK", "Prolif", "T4", "T8")



gene_heritability_normalization="True"
kurtosis_df <- extract_kurtosis_df(gene_heritability_normalization, cell_types, input_dir)

output_file <- paste0(visualization_dir, "cell_type_subsampled_kurtosis_estimates_normalize_genes_", gene_heritability_normalization, ".pdf")
boxplot <- make_cell_type_subsampled_boxplot(kurtosis_df, cell_types)
ggsave(boxplot, file=output_file, width=7.2, height=4.5, units="in")


gene_heritability_normalization="False"
kurtosis_df <- extract_kurtosis_df(gene_heritability_normalization, cell_types, input_dir)

output_file <- paste0(visualization_dir, "cell_type_subsampled_kurtosis_estimates_normalize_genes_", gene_heritability_normalization, ".pdf")
boxplot <- make_cell_type_subsampled_boxplot(kurtosis_df, cell_types)
ggsave(boxplot, file=output_file, width=7.2, height=4.5, units="in")


