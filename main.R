#Load Packages


#' Load data from a tsv located at a specific location
#'
#' @param filename (str):
#'
#' @return tibble: a (gene x sample) matrix with gene names as column names
#' 
#' @export
#'
#' @examples
#' `gene_data <- read_data('file/path/to/file.tsv')`
read_data <- function(filename){
  return(NULL)
}


#' Filter zero-variance genes from count matrix.
#'
#' @param count_data tibble: a (gene x sample) matrix of raw read counts
#'
#' @return tibble: a (gene x sample) matrix of raw reads with zero variance
#' genes removed.
#' @export
#'
#' @examples
#' `filtered_counts <- filter_zero_var_genes(count_data)`
filter_zero_var_genes <- function(count_data) {
  return(NULL)
}



#' Extract time point information from sample name
#'
#' @param x string: sample name from count matrix.
#'
#' @return string: string character representing sample time point 
#' @export
#'
#' @examples
#' `timepoint_from_sample("vAd_1")`
#' `"Ad"`
timepoint_from_sample <- function(x) {
  return(NULL)
}


#' Grab sample number from sample name
#'
#' @param x  string: sample name from count matrix.
#'
#' @return string: string character represent sample replicate number
#' @export
#'
#' @examples
#' `sample_number("vAd_1")`
#' `"1"`
sample_number <- function(x) {
  return(NULL)
}


#' Generate sample-level metadata from sample names
#'
#' @param sample_names vector: character vector of sample names from count matrix.
#'
#' @return tibble: tibble containing sample-level metadata including columns
#' "sample", "timepoint", and "replicate" that store sample names, sample
#' time points, and sample replicate number, respectively.
#' @export
#'
#' @examples
#' `meta_info_from_labels(colnames(count_matrix))`
meta_info_from_labels <- function(sample_names) {
  return(NULL)
}


#' Calculate total read counts for each sample in a count matrix.
#'
#' @param count_matrix tibble: a (gene x sample) matrix of raw read counts.
#'
#' @return vector: numeric vector of read totals from each sample
#' @export
#'
#' @examples
#' `get_library_size(count_matrix)`
get_library_size <- function(count_matrix) {
  return(NULL)
}


#' Normalize raw count matrix to counts per million.
#' 
#' Normalizes read counts using the following formula:
#' 
#' $C_{i, j} = \frac{X_{i, j}}{\frac{\sum \limits_{j}^N X_{i, j}}{10^6}}$
#'
#' @param count_matrix tibble: a (gene x sample) matrix of raw read counts.
#'
#' @return tibble: a (gene x sample) matrix with read count normalized to counts
#' per million 
#' @export
#'
#' @examples
#' `normalize_by_cpm(count_matrix)`
normalize_by_cpm <- function(count_matrix) { 
  return(NULL)
}


#' Normalize raw count matrix using DESeq2
#'
#' @param count_matrix tibble: a (gene x sample) matrix of raw reads
#' @param meta_data tibble: sample-level information tibble containing time point,
#' sample, and replicate information.
#' @param design_formula formula: formula of comparision of interest
#'
#' @return tibble: DESeq2 normalized count matrix.
#' @export
#'
#' @examples
#' `deseq_normalize(count_matrix, meta_data, ~ timepoint)`
deseq_normalize <- function(count_matrix, meta_data, design_formula) {
  return(NULL)
}


#' Perform and plot PCA using a provided count matrix.
#' 
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#' @param count_matrix tibble: count matrix to perform PCA over.
#' @param meta tibble: sample-level meta information
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs. 
#' @export
#'
#' @examples
#' `plot_pca(count_matrix, meta, "Raw Count PCA")`
plot_pca <- function(count_matrix, meta, title="") {
  return(NULL)
}


#' Plot gene count distributions for each sample using boxplots.
#' 
#'
#' @param count_matrix tibble: a (gene x sample) count matrix
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#' @export
#'
#' @examples
#' `plot_sample_distributions(count_matrix, scale_y_axis=TRUE, title='Raw Count Distributions')`
plot_sample_distributions <- function(count_matrix, scale_y_axis=FALSE, title="") {
  return(NULL)
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param count_matrix tibble: a (gene x sample) data matrix containing raw or 
#' normalized count values. 
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#' @export
#'
#' @examples
#' `plot_variance_vs_mean(count_matrix, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`
plot_variance_vs_mean <- function(count_matrix, scale_y_axis=FALSE, title="") {
  return(NULL)
}
