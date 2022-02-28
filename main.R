if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#Load Packages
library(DESeq2)
library(tidyverse)

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
  data <- read_tsv(filename, '\t')
  return(as_tibble(data))
}

verse_counts <- read_data('verse_counts.tsv')

#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x m) tibble of raw read counts
#'
#' @return tibble: a (n x m) tibble of raw reads with genes that have 
#' zero variance across samples removed
#' 
#' @note (g >= n)
#' 
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  filtered <- verse_counts[rowSums(verse_counts[,-1]) > 0,]
  return(filtered)
}

filt_v_counts <- filter_zero_var_genes(verse_counts)
#rownames(filt_v_counts)

#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point 
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(str) {
  tp <- substr(str,2,3)
  return(tp)
}

timepoint_from_sample('vAd_2')

#' Grab sample replicate number from sample name
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#' 
#' @note you may choose to return numeric values instead of strings here
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(str) {
  samp <- str_sub(str, -1)
  return(samp)
}

sample_replicate('vPO_1')

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
  tib <- as_tibble(sample_names)
  colnames(tib) <- 'sample'
  tib$timepoint <- timepoint_from_sample(tib$sample)
  tib$replicate <- sample_replicate(tib$sample)
  return(tib)
}

meta <- meta_info_from_labels(colnames(filt_v_counts)[colnames(filt_v_counts)!='gene'])
#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_
#' data tibble: a (n x m) tibble of raw read counts.
#'
#' @return named vector: numeric vector of read totals from each sample
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(data) {
  data[colnames(data)!='gene'] %>%
    summarise(across(everything(),~sum(.))) %>%
  return(data)
}

get_library_size(filt_v_counts)

#' Normalize raw count data to counts per million WITH pseudocounts using the 
#' following formula:
#'     (count + 1) / ((sample_library_size/10^6) + 1)
#'
#'
#' @param count_data tibble: a (n x m) matrix of raw read counts.
#'
#' @return tibble: a (n x m) matrix with read count normalized to counts
#' per million
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(data) { 
  gene <- data$gene
  size_factors <- get_library_size(data[c(-1)]) %>%
    slice(1) %>% unlist(., use.names=FALSE)
  cpm <- as_tibble(t(apply(data[-1],1,function(x) x/size_factors*10^6)))
  return(cbind(gene,cpm))
  
}

norm_cpm <- normalize_by_cpm(filt_v_counts)
norm_cpm

#' Normalize raw count data using DESeq2
#'
#'
#' @param count_data tibble: a (n x m) matrix of raw reads
#' @param meta_data tibble: sample-level information tibble containing time point,
#' sample, and replicate information.
#' @param design_formula formula: formula of comparision of interest
#'
#' @return tibble: a (n x m) tibble of DESeq2 normalized count data.
#'
#' @example ' `deseq_normalize(count_data, meta_data, ~ timepoint)`

deseq_normalize <- function(count_data, meta_data, design_formula) {
  genes <- count_data$gene
  count_data <- select(count_data, -c(gene))
  dds <- DESeqDataSetFromMatrix(
    countData=count_data,
    colData=meta_data,
    design=~1
  )
  dds <- estimateSizeFactors(dds)
  norm <- counts(dds, normalized=TRUE)
  norm <- as_tibble(norm) %>%
    mutate(gene=genes) %>%
    relocate(gene) %>% # move `gene` column to be first
    return()
}

norm_deseq <- deseq_normalize(filt_v_counts, meta)
norm_deseq

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

plot_pca <- function(data, meta, title) {
  pca <- prcomp(scale(t(data)), center=FALSE, scale=FALSE)
  plot_data <- meta
  plot_data$PC1 <- pca$x[ , 1]
  plot_data$PC2 <- pca$x[ , 2]
  percent_var <- pca$sdev^2 / sum( pca$sdev^2 )
  pca_plot <- ggplot(plot_data, ggplot2::aes(x=PC1, y=PC2, col=timepoint)) +
    geom_point() +
    xlab(paste0("PC1: ",round(percent_var[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percent_var[2] * 100),"% variance")) +
    ggtitle(title)
  return(pca_plot)
}

plot_pca(filt_v_counts, meta, "Raw Count PCA")

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
  long_counts <- pivot_longer(data,
                                  cols=colnames(data),
                                  names_to='sample',
                                  values_to='counts') %>%
    
    mutate(sample=factor(sample,levels=colnames(data)))
  
  if (scale_y_axis) {
    # filter out counts == 0 so log10 below doesn't choke
    long_counts <- filter(long_counts,counts != 0)
  }
  
  dist_plot <- ggplot(long_counts, aes(x=sample, y=counts, col=sample)) +
    geom_boxplot() +
    ggtitle(title)
  
  if (scale_y_axis) {
    dist_plot <- dist_plot + scale_y_log10()
  }
  return(dist_plot)
}

plot_sample_distributions(filt_v_counts, scale_y_axis = TRUE, title = 'Raw Count Distributions')

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
  return(NULL)
}
