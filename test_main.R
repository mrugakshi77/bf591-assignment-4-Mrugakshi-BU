#!/usr/bin/Rscript
source("main.R")
library(testthat)

test_names <- c("a","b","c","d")
test_meta <- tibble(
  sample_name=test_names,
  timepoint=c("Ad","Ad","P0","P0"),
  replicate=c("1","2","1","2")
)
test_obj <- tibble(
  gene=paste0('g',seq(1,4)),
  a=c(1,2,1,0),
  b=c(2,3,1,0),
  c=c(3,4,1,0),
  d=c(4,5,1,0)
)
 
test_that("read_data loads correctly", {
  res <- read_data("verse_counts.tsv")
  
  expect_equal(dim(res), c(55416, 9))
  expect_true(is.data.frame(res))
  
})

test_that("filter_zero_var_genes works properly", {
 
  filtered <- filter_zero_var_genes(test_obj)
  
  expect_equal(nrow(filtered),2)
  expect_equal(ncol(filtered),5)
  expect_equal(filtered$gene,c("g1","g2"))
  
})

test_that("timepoint_from_sample", {
  expect_equal(timepoint_from_sample("vAd_1"),"Ad") 
})

test_that("sample_replicate", {
  expect_equal(sample_replicate("vAd_1"),"1") 
})

test_that("meta_info_from_labels", {
  test_names <- c("vAD_1","vAD_2","vP0_1","vP0_2")
  meta <- meta_info_from_labels(test_names)
  
  expect_equal(as.vector(meta$sample),test_names)
  expect_equal(as.vector(meta$timepoint),c("AD","AD","P0","P0"))
  expect_equal(as.vector(meta$replicate),c("1","2","1","2"))
})

test_that("get_library_size",{
  sizes <- get_library_size(test_obj) 
  expect_equal(colnames(sizes),test_names)
  expect_equal(dplyr::slice(sizes,1) %>% unlist(., use.names=FALSE),c(4,6,8,10))
})

test_that("normalize_by_cpm",{
  
  sizes <- get_library_size(test_obj) 
  
  cpm <- normalize_by_cpm(test_obj)
  
  expect_equal(nrow(cpm),4)
  expect_equal(ncol(cpm),5)
  expect_equal(cpm$gene,c("g1","g2","g3","g4"))
  
  expect_equal(cpm$a, test_obj$a/sizes$a*10^6)
  expect_equal(cpm$b, test_obj$b/sizes$b*10^6)
})

test_that("deseq_normalize",{
  res <- read_data("verse_counts.tsv")
  filtered <- filter_zero_var_genes(res)
  meta <- meta_info_from_labels(colnames(res[-1]))
  
  # this should really be done on a synthetic matrix where we know the answer by
  # construction
  norm_counts <- deseq_normalize(filtered,meta) 
  
  # there are 32613 filtered genes
  expect_equal(dim(norm_counts), c(32613, 9))
  expect_true(is.data.frame(norm_counts))
  
  # I looked this value up
  expect_equal(norm_counts$vP0_1[1],15.9359561)
  
})

test_that("plot_pca",{
  
  plot <- plot_pca(test_obj[-1],test_meta)
  
  # testing all geoms in the ggplot object
  geoms <- vapply(plot$layers,function(x) class(x$geom)[1],"a")
  expect_setequal(geoms,c("GeomPoint"))
})

test_that("plot_sample_distribution",{
  
  plot <- plot_sample_distributions(test_obj[-1])
  
  # testing all geoms in the ggplot object
  geoms <- vapply(plot$layers,function(x) class(x$geom)[1],"a")
  expect_setequal(geoms,c("GeomBoxplot"))
})

test_that("plot_variance_vs_mean",{
  
  plot <- plot_variance_vs_mean(test_obj[-1])
  
  # testing all geoms in the ggplot object
  geoms <- vapply(plot$layers,function(x) class(x$geom)[1],"a")
  expect_setequal(geoms,c("GeomPoint","GeomSmooth"))
})