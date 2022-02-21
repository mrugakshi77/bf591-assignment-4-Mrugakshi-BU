#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  read_tsv(filename)
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x 1+m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  # # non-tidyverse
  # #Method 1
  gene_vars <- apply(verse_counts[,-c(1)], 1, var)
  return(verse_counts[gene_vars!=0,])

  # #tidyverse
  # #Method 2
  # nonzero_gene_vars <- tidyr::pivot_longer(count_data, -c(gene),names_to="sample")  %>%
  #   dplyr::group_by(gene) %>%
  #   dplyr::summarize(var=var(value)) %>%
  #   dplyr::filter(var!=0) %>%
  #   return()
  # count_data %>%
  #   dplyr::filter(gene %in% nonzero_gene_vars$gene) %>%
  #   return()
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
  return(stringr::str_extract(x, "(?<=^v)(.*?)(?=_)"))
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
  stringr::str_split(x , "_")[[1]][2] %>%
    return()
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

  # #non-tidyverse
  # #Method 1
  # timepoints <- sapply(sample_names, timepoint_from_sample)
  # replicate <- sapply(sample_names, sample_replicate)
  # return(tibble::tibble(sample=sample_names, timepoint=timepoints, replicate=replicate))
  #
  #Method 2
  tibble(
    sample=sample_names,
    timepoint=vapply(sample_names, timepoint_from_sample, "x"),
    replicate=vapply(sample_names, sample_replicate, "x")
  ) %>%
    return()

  # # tidyverse
  #Method 3 -
  # tibble(
  #   sample=sample_names
  # ) %>%
  #   dplyr::mutate(temp=str_sub(sample,start=2)) %>%
  #   tidyr::separate(temp,into=c("timepoint","replicate"), sep="_") %>%
  #   return()
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble or named vector of read totals from each sample. Vectors must
#' be length `_S_ `, a tibble can be `(1 x _S_)` with sample names as columns
#' names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  # non-tidyverse
  # Method 1
  # count_data %>%
  #   apply(2,sum) %>%
  #   return()

  #Tidyverse
  #Method 2
  count_data[colnames(count_data)!='gene'] %>%
    dplyr::summarise(dplyr::across(everything(), ~ sum(.))) %>%
    return()

  # count_data %>%
  #   # pivot_longer(cols=where(is.numeric))%>%
  #   #dplyr::group_by(name) %>%
  #   dplyr::summarise(size = sum(value)) %>%
  #   return()

}


#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {

  gene <- count_data$gene
  size_factors <- get_library_size(count_data[c(-1)]) %>%
    dplyr::slice(1) %>% unlist(., use.names=FALSE)
  cpm <- as_tibble(t(apply(count_data[-1],1,function(x) x/size_factors*10^6)))
  return(cbind(gene,cpm))

  # note: the methods below have not been tested
  # #Method 2:
  # count_data %>%
  #   mutate(.[-c(1)]/(get_library_size(.[c(-1)])/10^6))) %>%
  #   return()

  # #Method 3:
  #size_factors <- get_library_size(count_data[c(-1)])/10^6
  #return(
  #  #cbind(count_data[1], count_data[-c(1)]/(size_factors[rep(1, nrow(count_data)),] ))
  #  cbind(count_data[1], count_data[-1]/(size_factors[rep(1, nrow(count_data)),] ))
  #)

  #Method 4: Tibble version of method 2
  #count_data %>%
  #  mutate((.[c(-1)])/(unlist(get_library_size(.[c(-1)])/10^6))) %>%
  #  return()

  # #Method 5: Tibble and a for loop
  # temp<- count_data
  # size_factors <- get_library_size(count_data[c(-1)])/10^6
  # for(column in colnames(size_factors)){
  #   temp[column] <- temp[column]/(size_factors[[column]][1])
  # }
  # return(temp)
}

#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
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
  genes <- count_data$gene
  count_data <- select(count_data, -c(gene))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=count_data,
    colData=meta_data,
    design=~1
  )
  dds <- estimateSizeFactors(dds)
  norm <- DESeq2::counts(dds, normalized=TRUE)
  norm <- tibble::as_tibble(norm) %>%
    mutate(gene=genes) %>%
    relocate(gene) %>% # move `gene` column to be first
    return()
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
  pca_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x=PC1, y=PC2, col=timepoint)) +
    ggplot2::geom_point() +
    ggplot2::xlab(paste0("PC1: ",round(percent_var[1] * 100),"% variance")) +
    ggplot2::ylab(paste0("PC2: ",round(percent_var[2] * 100),"% variance")) +
    ggplot2::ggtitle(title)
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

  long_counts <- tidyr::pivot_longer(data,
                                     cols=colnames(data),
                                     names_to='sample',
                                     values_to='counts') %>%
    mutate(sample=factor(sample,levels=colnames(data)))

  if (scale_y_axis) {
    # filter out counts == 0 so log10 below doesn't choke
    long_counts <- filter(long_counts,counts != 0)
  }

  dist_plot <- ggplot2::ggplot(long_counts, ggplot2::aes(x=sample, y=counts, col=sample)) +
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(title)

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
  plot_data <- tibble::tibble(mean=means, variance=variances)
  plot_data$rank <- rank(plot_data$mean)
  mv_plot <- ggplot2::ggplot(plot_data, aes(x=rank, y=variance)) +
    ggplot2::geom_point(alpha=0.5) +
    ggplot2::geom_smooth(method='gam', formula = y ~ s(x, bs = "cs")) +
    ggplot2::xlab("Rank(Mean)") +
    ggplot2::ylab("Variance") +
    ggplot2::ggtitle(title)
  if (scale_y_axis) {
    mv_plot <- mv_plot + ggplot2::scale_y_log10()
  }
  return(mv_plot)
}
