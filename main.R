#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  data <- read.delim(filename)
  return (tibble(data))
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  nonzero_genes <- rowSums(verse_counts[-1])!=0
  return (verse_counts[nonzero_genes, ])
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(x) {
  return (substring(x, 2, 3))
}

#' Grab sample replicate number from sample name
#'
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(x) {
  return (substring(x, 5))
}


#' Generate sample-level metadata from sample names.
#'
#' Will include columns named "sample", "timepoint", and "replicate" that store
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample",
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint"
#' stores sample time points; and "replicate" stores sample replicate
#'
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  timepoints <- unlist(lapply(sample_names, timepoint_from_sample))
  replicates <- unlist(lapply(sample_names, sample_replicate))
  metadata <- tibble(
    sample = sample_names, 
    timepoint = timepoints, 
    replicate = replicates
  )
  return (metadata)
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble or named vector of read totals from each sample. Vectors must
#' be length `_S_ `, a tibble can be `(1 x _S_)` with sample names as columns
#' names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  temp <- count_data[colnames(count_data)!='gene'] %>%
    summarise(across(everything(), ~ sum(.)))
    
  return(temp)
}


#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {
  gene <- count_data$gene
  library_size <- get_library_size(count_data[c(-1)]) %>%
    dplyr::slice(1) %>% unlist(., use.names=FALSE)
  cpm <- as_tibble(t(apply(count_data[-1],1,function(x) x/library_size*10^6)))
  return(cbind(gene,cpm))
}

#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts

#' @param meta_data tibble: sample-level information tibble corresponding to the
#' count matrix columns
#'
#' @return tibble: DESeq2 normalized count matrix
#' @export
#'
#' @examples
#' `deseq_normalize(count_data, meta_data)`
deseq_normalize <- function(count_data, meta_data) {
  count_mat <- as.matrix(count_data[-1])
  row.names(count_mat) <- count_data$gene
  
  dds <- DESeqDataSetFromMatrix(
    countData=count_mat,
    colData=tibble(sample_name=colnames(count_data[-1])),
    design=~1 # no formula is needed for normalization, ~1 produces a trivial design matrix
  )
  # compute normalization factors
  dds <- estimateSizeFactors(dds)
  
  # extract the normalized counts
  deseq_norm_counts <- as_tibble(counts(dds,normalized=TRUE)) %>%
    mutate(gene=count_data$gene) %>%
    relocate(gene)
  
  return(deseq_norm_counts)
}


#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
  pca <- prcomp(t(data))
  plot_data <- meta
  plot_data$PC1 <- pca$x[ , 1]
  plot_data$PC2 <- pca$x[ , 2]
  percent_var <- pca$sdev^2 / sum( pca$sdev^2 )
  pca_plot <- ggplot(plot_data, aes(x=PC1, y=PC2, col=timepoint)) +
    geom_point() +
    xlab(paste0("PC1: ",round(percent_var[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percent_var[2] * 100),"% variance")) +
    ggtitle(title)
  return(pca_plot)
}


#' Plot gene count distributions for each sample using boxplots.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
  temp <- tidyr::pivot_longer(data, cols=colnames(data), names_to='sample', values_to='counts') %>%
    mutate(sample=factor(sample,levels=colnames(data)))
  
  dist_plot <- ggplot(temp, aes(x=sample, y=counts, col=sample)) +
    geom_boxplot() +
    ggtitle(title)
  
  if (scale_y_axis) {
    dist_plot <- dist_plot + ggplot2::scale_y_log10()
  }
  return(dist_plot)
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis=FALSE, title="") {
  means <- apply(data, 1, mean)
  variances <- apply(data, 1, var)
  plot_data <- tibble(mean=means, variance=variances)
  plot_data$rank <- rank(plot_data$mean)
  mv_plot <- ggplot(plot_data, aes(x=rank, y=variance)) +
    geom_point(alpha=0.5) +
    geom_smooth(method='gam', formula = y ~ s(x, bs = "cs")) +
    xlab("Rank(Mean)") +
    ylab("Variance") +
    ggtitle(title)
  if (scale_y_axis) {
    mv_plot <- mv_plot + ggplot2::scale_y_log10()
  }
  return(mv_plot)
}
