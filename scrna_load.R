library(tidyverse)

setwd("/ix/djishnu/Aaron")

library(reticulate)

reticulate::use_condaenv(condaenv = "analysis_env", required = T)

library(anndata)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(SeuratWrappers)

library(SingleCellExperiment)
sc = reticulate::import("scanpy")


data_path = "/ix/djishnu/Aaron/scrnaseq/covid19atlas/"
pbmc_data = "meyer_nikolic_covid_pbmc.cellxgene.20210813.h5ad"

pbmc = anndata::read_h5ad(paste0(data_path, pbmc_data))

sce = SingleCellExperiment(
  assays = list(logcounts = t(pbmc$X)),
  colData = pbmc$obs,
  rowData = pbmc$var,
  reducedDims = pbmc$obsm
)

# import genes of interest

genes = read_tsv("/ix/djishnu/Aaron/data/gene_list.tsv")

gene_locations = which(pbmc$var_names %in% genes$unique.results.distinct.gene_name.)


# subset genes
pbmc_desired = pbmc[,gene_locations]

anndata::write_h5ad(pbmc_desired, filename = "/ix/djishnu/Aaron/data/genelist_counts.h5ad")


# # get log counts
# 
# lgc = logcounts(sce)
# 
# 
# # get rows with genes of interest
# 
# gene_expr = lgc[gene_locations,]
# 
# 
# 
# # conver to anndata
# 
# exprs = assay(sce, "logcounts")
# col_data = as.data.frame(colData(sce))
# row_data = as.data.frame(rowData(sce))
# embedding = as.data.frame(reducedDims(sce))
# 
# anndata_output = AnnData(X = exprs, obs = col_data, var = row_data,
#                             obsm = embedding)




