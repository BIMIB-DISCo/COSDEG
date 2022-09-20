### cosdeg.R

#webpath <- "https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/"
#dir.create("./data/raw", recursive = T)

#file_list <- c("Normal_PBMC_13.h5", "Normal_PBMC_14.h5", "Normal_PBMC_5.h5", "nCoV_PBMC_15.h5", "nCoV_PBMC_17.h5", "nCoV_PBMC_1.h5")
#for (i in file_list) {
#  download.file(url = paste0(webpath, i), destfile = paste0("./data/raw/", i))
#}

#install_github("https://github.com/chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
#install.packages('assertthat')


suppressMessages(require(assertthat))
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
suppressMessages(require(DoubletFinder))


#dataset_folders <- list.dirs(path = "data/", full.names = F, recursive = F)
#dataset_folders <- dataset_folders[-5] # remove raw folder
#data_dir <- list()
#
#for (d in dataset_folders) { 
#  data_dir[[d]] <- file.path("data", d, "filtered_feature_bc_matrix")
#}  
#
##load data
#s_data <- list()
#for (d in dataset_folders) { 
#  print(data_dir[[d]])
#  tmp_expression <- Seurat::ReadMtx(file.path(data_dir[[d]], "matrix.mtx.gz"), features = file.path(data_dir[[d]], "features.tsv.gz"), cells = file.path(data_dir[[d]], "barcodes.tsv.gz"))
#  s_data[[d]] <- CreateSeuratObject(counts = tmp_expression, project = d)
#}




filter_offgenes <- function (s_data, summary_qc, s_data_status) {
  s_data_dim <- list()
  data.filt <- list()
  selected_f0 <- list()
  
  print("Remove genes with zero value in all samples")
  dataset_folders <- names(s_data)
  
  for (d in dataset_folders) {
    selected_f0[[d]] <- rownames(s_data[[d]])[Matrix::rowSums(s_data[[d]]) > 0]
    data.filt[[d]] <- subset(s_data[[d]], features = selected_f0[[d]])
    
    s_data_dim[[d]] <- list(dim(data.filt[[d]]))
    summary_qc["offgenes",d] <- s_data_dim[d]
    
    print(dim(data.filt[[d]]))
    
    s_data_status[d, "offgenes"] <- 1
    
    #browser()
    #s_data_status[[d]][["offgenes"]] <- TRUE
  }
  return(list(data.filt = data.filt, selected_f0 = selected_f0, summary_qc = summary_qc, status=s_data_status))
} 


filter_emptylets <- function (s_data, summary_qc, perc_zeros, s_data_status) {
  s_data_dim <- list()
  data.filt <- list()
  selected_cE <- list()
  
  print("Remove cells with insufficient number of reads")
  dataset_folders <- names(s_data)
  
  for (d in dataset_folders) {
    #browser()
    selected_cE[[d]] <- colnames(s_data[[d]])[s_data[[d]]$nFeature_RNA/dim(s_data[[d]])[1]>=1-perc_zeros]
    data.filt[[d]] <- subset(s_data[[d]], cells = selected_cE[[d]])
    
    s_data_dim[[d]] <- list(dim(data.filt[[d]]))
    summary_qc["emptylets",d] <- s_data_dim[d]
    
    print(dim(data.filt[[d]]))
    
    #s_data_status[[d]][["emptylets"]] <- perc_zeros
    s_data_status[d, "emptylets"] <- perc_zeros
    data.filt[[d]][["selected_cE"]] <- selected_cE[[d]]
  }
  return(list(data.filt = data.filt, selected_cE = selected_cE, summary_qc = summary_qc, status=s_data_status))
}


filter_doublets <- function (s_data, summary_qc, multiplet_rate, s_data_status) {
  
  s_data_dim <- list()
  data.filt <- list()
  selected_cD <- list()
  
  sweep.res <- list()
  sweep.stats <- list()
  bcmvn <- list()
  pK_choose <- list()
  
  pca_npcs <- 20
  umap_npcs <- 10
  doublets_npcs <- 10
  doublets_pN <- 0.25
  
  plots <- list()
  
  print("Remove doublets")
  dataset_folders <- names(s_data)
  
  plots$AUC <- list()
  plots$BCmetric <- list()
  plots$parameter_estimation <- list()
  plots$doublets <- list()
  plots$violin <- list()
  
  for (d in dataset_folders) {
    print(d)
    data.filt[[d]] <- s_data[[d]]
    data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^MT-", col.name = "percent_mito")
    data.filt[[d]] <- NormalizeData(data.filt[[d]])
    data.filt[[d]] <- FindVariableFeatures(data.filt[[d]], verbose = F)
    data.filt[[d]] <- ScaleData(data.filt[[d]], vars.to.regress = c("nFeature_RNA", "percent_mito"),verbose = F)
    data.filt[[d]] <- RunPCA(data.filt[[d]], verbose = F, npcs = pca_npcs)
    data.filt[[d]] <- RunUMAP(data.filt[[d]], dims = 1:umap_npcs, verbose = F)
    
    # Can run parameter optimization with paramSweep, but skip for now.
    sweep.res[[d]] <- paramSweep_v3(data.filt[[d]], PCs = 1:doublets_npcs) 
    sweep.stats[[d]] <- summarizeSweep(sweep.res[[d]], GT = FALSE) 
    bcmvn[[d]] <- find.pK(sweep.stats[[d]])
    #browser()
    x <- plot(x = bcmvn[[d]]$ParamID, y = bcmvn[[d]]$MeanAUC, pch = 18, 
              col = "black", cex = 0.75, xlab = NA, ylab = NA)
    x <- lines(x = bcmvn[[d]]$ParamID, y = bcmvn[[d]]$MeanAUC, col = "black", lty = 2)
    plots$AUC[[d]] <- x 
    
    par(new = TRUE)
    x <- plot(x = bcmvn[[d]]$ParamID, y = bcmvn[[d]]$BCmetric, pch = 16, col = "#41b6c4", cex = 0.75)
    axis(side = 4)
    x <- lines(x = bcmvn[[d]]$ParamID, y = bcmvn[[d]]$BCmetric, col = "#41b6c4")
    plots$BCmetric[[d]] <- x
    
    plots$parameter_estimation[[d]] <- barplot(bcmvn[[d]]$BCmetric, names.arg = bcmvn[[d]]$pK, las=2)
  }
  #browser()
  for (d in dataset_folders) {
    pK <- as.numeric(as.character(bcmvn[[d]]$pK))
    BCmetric <- bcmvn[[d]]$BCmetric
    pK_choose[[d]] <- pK[which(BCmetric%in%max(BCmetric))]
  }
  #browser()
  for (d in dataset_folders) {
    # define the expected number of doublet cells.
    nExp <- round(ncol(data.filt[[d]]) * multiplet_rate)  # expect 4% doublets
    
    data.filt[[d]] <- doubletFinder_v3(data.filt[[d]], pN = doublets_pN, pK = pK_choose[[d]], nExp = nExp, PCs = 1:doublets_npcs)
    
    # name of the DF prediction can change, so extract the correct column name.
    DF.name = colnames(data.filt[[d]]@meta.data)[grepl("DF.classification", colnames(data.filt[[d]]@meta.data))]
    
    plots$doublets[[d]] <- cowplot::plot_grid(ncol = 2, DimPlot(data.filt[[d]], group.by = "orig.ident") + NoAxes(), DimPlot(data.filt[[d]], group.by = DF.name) + NoAxes())
    plots$violin[[d]] <-  VlnPlot(data.filt[[d]], features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  }
  
  for (d in dataset_folders) {
    DF.name = colnames(data.filt[[d]]@meta.data)[grepl("DF.classification", colnames(data.filt[[d]]@meta.data))]
    selected_cD[[d]] <- colnames(data.filt[[d]])[data.filt[[d]][[DF.name[1]]]=="Singlet"]
    data.filt[[d]] <- subset(data.filt[[d]], cells = selected_cD[[d]])
    
    s_data_dim[[d]] <- list(dim(data.filt[[d]]))
    summary_qc["doublets",d] <- s_data_dim[d]
    
    print(dim(data.filt[[d]]))
    
    #s_data_status[[d]][["doublets"]] <- multiplet_rate
    s_data_status[d, "doublets"] <- multiplet_rate
    data.filt[[d]][["selected_cD"]] <- selected_cD[[d]]
  }
  
  browser()
  return(list(data.filt = data.filt, selected_cD = selected_cD, 
              summary_qc = summary_qc, plots = plots, sweep.res = sweep.res, 
              sweep.stats = sweep.stats, bcmvn = bcmvn, pK_choose = pK_choose, status=s_data_status))
}



filter_lysed <- function (s_data, summary_qc, sigma_mito, s_data_status) {
  
  s_data_dim <- list()
  data.filt <- list()
  perc_mitoc <- list()
  selected_cM <- list()
  
  print("Remove lysed cells")
  dataset_folders <- names(s_data)
  
  for (d in dataset_folders) {
    data.filt[[d]] <- s_data[[d]]
    data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^MT-", col.name = "percent_mito")
    data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^RP[SL]", col.name = "percent_ribo")
    # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
    data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^HB[^(P)]", col.name = "percent_hb")
    data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "PECAM1|PF4", col.name = "percent_plat")
    
    perc_mitoc[[d]] <- mean(data.filt[[d]]$percent_mito/100.) + sigma_mito * sd(data.filt[[d]]$percent_mito/100.)
    
    selected_cM[[d]] <- colnames(data.filt[[d]])[data.filt[[d]]$percent_mito/100 <= perc_mitoc[[d]]]
    
    data.filt[[d]] <- subset(data.filt[[d]], cells = selected_cM[[d]])
    
    s_data_dim[[d]] <- list(dim(data.filt[[d]]))
    summary_qc["mitochodrial",d] <- s_data_dim[d]
    
    print(dim(data.filt[[d]]))
    
    #s_data_status[[d]][["mitochodrial"]] <- sigma_mito
    s_data_status[d, "mitochodrial"] <- sigma_mito
    data.filt[[d]][["selected_cM"]] <- selected_cM[[d]]
  }
  return(list(data.filt = data.filt, selected_cM = selected_cM, summary_qc = summary_qc, perc_mitoc = perc_mitoc, status=s_data_status))
}


pre_qc <- function (s_data, summary_qc) {
  
  s_data_dim <- list()
  
  print("Pre QC")
  dataset_folders <- names(s_data)
  
  for (d in dataset_folders) {  
    s_data_dim[[d]] <- list(dim(s_data[[d]]))
    summary_qc["pre_QC",d] <- s_data_dim[d]
    
    print(dim(s_data[[d]]))
  }
  
  return(list(summary_qc = summary_qc))
}

quality_check <- function(s_data, perc_zeros, sigma_mito, multiplet_rate) { 


  assert_that(is.numeric(perc_zeros))
  assert_that(is.numeric(sigma_mito))
  assert_that(is.numeric(multiplet_rate))
  assert_that(is.list(s_data))
  
  for (d in names(s_data)) {
    #s_data_dim[[d]] is a seurat object
    assert_that(class(s_data[[1]])[1]=="Seurat")
  }
  
  
  dataset_folders <- names(s_data)
  
  summary_qc <- data.frame()
  data.filt <- list()
  s_data_dim <- list()
  s_data_status <- data.frame()
  
  result <- pre_qc(s_data, summary_qc)
  summary_qc <- result$summary_qc
  
  
  # for (d in dataset_folders) {
  #   s_data_dim[[d]] <- list(dim(s_data[[d]]))
  #   summary_qc["pre_QC",d] <- s_data_dim[d]
  #   
  #   #s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "^MT-", col.name = "percent_mito")
  #   #s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "^RP[SL]", col.name = "percent_ribo")
  #   ## Percentage hemoglobin genes - includes all genes starting with HB except HBP.
  #   #s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "^HB[^(P)]", col.name = "percent_hb")
  #   #s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "PECAM1|PF4", col.name = "percent_plat")
  # }
  
  summary_qc
  
  result <- filter_offgenes(s_data, summary_qc, s_data_status)
  data.filt <- result$data.filt
  selected_f0 <- result$selected_f0
  summary_qc <- result$summary_qc
  s_data_status <- result$status
  
  # print("Remove genes with zero value in all samples")
  # selected_f0 <- list()
  # for (d in dataset_folders) {
  #   selected_f0[[d]] <- rownames(s_data[[d]])[Matrix::rowSums(s_data[[d]]) > 0]
  #   data.filt[[d]] <- subset(s_data[[d]], features = selected_f0[[d]])
  #   
  #   s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  #   summary_qc["zero_genes",d] <- s_data_dim[d]
  #   
  #   print(dim(data.filt[[d]]))
  # }
  
  summary_qc
  
  result <- filter_emptylets(s_data, summary_qc, perc_zeros, s_data_status)
  data.filt <- result$data.filt
  selected_cE <- result$selected_cE
  summary_qc <- result$summary_qc
  s_data_status <- result$status
  
  
  #print("Remove cells with insufficient number of reads")
  ##selecting cell with less than 90% of zeros
  ##alldata[1:26,1:20]@assays$RNA@counts
  ##colSums(alldata[1:26,1:20]@assays$RNA@counts==0)/26
  #selected_cE <- list()
  #for (d in dataset_folders) {
  #  selected_cE[[d]] <- colnames(data.filt[[d]])[data.filt[[d]]$nFeature_RNA/dim(data.filt[[d]])[1]>=1-perc_zeros]
  #  data.filt[[d]] <- subset(data.filt[[d]], cells = selected_cE[[d]])
  #  
  #  s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  #  summary_qc["empty_cells",d] <- s_data_dim[d]
  #  
  #  print(dim(data.filt[[d]]))
  #}
  
  summary_qc
  
  result <- filter_doublets(s_data, summary_qc, multiplet_rate, s_data_status)
  data.filt <- result$data.filt
  selected_cD <- result$selected_cD
  summary_qc <- result$summary_qc
  plots$doublets_filter <- result$plots
  pK_choose <- result$pK_choose
  s_data_status <- result$status
  
  # print("Remove doublets")
  # 
  # sweep.res <- list()
  # sweep.stats <- list()
  # bcmvn <- list()
  # 
  # for (d in dataset_folders) {
  #   print(d)
  #   data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^MT-", col.name = "percent_mito")
  #   data.filt[[d]] = NormalizeData(data.filt[[d]])
  #   data.filt[[d]] = FindVariableFeatures(data.filt[[d]], verbose = F)
  #   data.filt[[d]] = ScaleData(data.filt[[d]], vars.to.regress = c("nFeature_RNA", "percent_mito"),verbose = F)
  #   data.filt[[d]] = RunPCA(data.filt[[d]], verbose = F, npcs = 20)
  #   data.filt[[d]] = RunUMAP(data.filt[[d]], dims = 1:10, verbose = F)
  #   
  #   # Can run parameter optimization with paramSweep, but skip for now.
  #   sweep.res[[d]] <- paramSweep_v3(data.filt[[d]], PCs = 1:10) 
  #   sweep.stats[[d]] <- summarizeSweep(sweep.res[[d]], GT = FALSE) 
  #   bcmvn[[d]] <- find.pK(sweep.stats[[d]])
  #   #barplot(bcmvn[[d]]$BCmetric, names.arg = bcmvn[[d]]$pK, las=2)
  #   print(paste(d,"END"))
  # }
  # 
  # pK_choose <- list()
  # for (d in dataset_folders) {
  #   pK=as.numeric(as.character(bcmvn[[d]]$pK))
  #   BCmetric=bcmvn[[d]]$BCmetric
  #   pK_choose[[d]]=pK[which(BCmetric%in%max(BCmetric))]
  # }
  # 
  # for (d in dataset_folders) {
  #   # define the expected number of doublet cellscells.
  #   nExp <- round(ncol(data.filt[[d]]) * multiplet_rate)  # expect 4% doublets
  #   data.filt[[d]] <- doubletFinder_v3(data.filt[[d]], pN = 0.25, pK = pK_choose[[d]], nExp = nExp, PCs = 1:10)
  #   
  #   # name of the DF prediction can change, so extract the correct column name.
  #   DF.name = colnames(data.filt[[d]]@meta.data)[grepl("DF.classification", colnames(data.filt[[d]]@meta.data))]
  #   
  #   #cowplot::plot_grid(ncol = 2, DimPlot(data.filt[[d]], group.by = "orig.ident") + NoAxes(), DimPlot(data.filt[[d]], group.by = DF.name) + NoAxes())
  #   #VlnPlot(data.filt[[d]], features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  # }
  # 
  # selected_cD <- list()
  # for (d in dataset_folders) {
  #   DF.name = colnames(data.filt[[d]]@meta.data)[grepl("DF.classification", colnames(data.filt[[d]]@meta.data))]
  #   selected_cD[[d]] <- colnames(data.filt[[d]])[data.filt[[d]][[DF.name[1]]]=="Singlet"]
  #   data.filt[[d]] <- subset(data.filt[[d]], cells = selected_cD[[d]])
  #   
  #   s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  #   summary_qc["doublets",d] <- s_data_dim[d]
  # }
   
  
  result <- filter_lysed(s_data, summary_qc, sigma_mito, s_data_status)
  data.filt <- result$data.filt
  selected_cM <- result$selected_cM
  summary_qc <- result$summary_qc
  
  summary_qc
  
  # print("Remove lysed cells")
  # perc_mitoc <- list()
  # selected_cM <- list()
  # for (d in dataset_folders) {
  #   data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^MT-", col.name = "percent_mito")
  #   data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^RP[SL]", col.name = "percent_ribo")
  #   # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
  #   data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^HB[^(P)]", col.name = "percent_hb")
  #   data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "PECAM1|PF4", col.name = "percent_plat")
  #   
  #   perc_mitoc[[d]] <- mean(data.filt[[d]]$percent_mito/100.) + sigma_mito * sd(data.filt[[d]]$percent_mito/100.)
  #   
  #   selected_cM[[d]] <- colnames(data.filt[[d]])[data.filt[[d]]$percent_mito/100 <= perc_mitoc[[d]]]
  #   
  #   data.filt[[d]] <- subset(data.filt[[d]], cells = selected_cM[[d]])
  #   
  #   s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  #   summary_qc["mitochodrial",d] <- s_data_dim[d]
  #   print(dim(data.filt[[d]]))
  # }
  
  summary_qc
  
  selected <- list("selected_f0"=selected_f0, "selected_cE"=selected_cE, "selected_cD"=selected_cD, "selected_cM"=selected_cM)
  
  plots_handles <- list()
  
  return(list("data.filt"=data.filt, "bcmvn"=bcmvn, "summary_qc"=summary_qc, "selected"=selected, "plots_handles"=plots_handles))
}


perc_zeros <- 0.90
sigma_mito <- 1
multiplet_rate <- 0.04


#qc_results <- quality_check(s_data, perc_zeros, sigma_mito, multiplet_rate)
 

# ####################
# 
# # Merge datasets into one single seurat object
# alldata <- merge(s_data[[1]], s_data[-1], add.cell.ids = names(s_data))
# 
# feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
# VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
#   NoLegend()
# 
# # Merge datasets into one single seurat object
# alldata <- merge(data.filt[[1]], data.filt[-1], add.cell.ids = names(data.filt))
# 
# feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
# VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
#   NoLegend()
# 
# 
# for (d in dataset_folders) {
#   # Compute the relative expression of each gene per cell Use sparse matrix
#   # operations, if your dataset is large, doing matrix divisions the regular way
#   # will take a very long time.
#   par(mar = c(4, 8, 2, 1))
#   C <- data.filt[[d]]@assays$RNA@counts
#   C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
#   most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
#   boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
#           col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
# }
# 
# ####################

cluster_seurat_obj = function(alldata) {
  alldata = NormalizeData(alldata)
  alldata = FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(alldata)
  alldata <- ScaleData(alldata, features = all.genes)
  
  alldata <- RunPCA(alldata)
  
  alldata <- FindNeighbors(alldata, dims = 1:10)
  alldata <- FindClusters(alldata, resolution = 0.5)
  
  alldata <- RunUMAP(alldata, dims = 1:10)
  
  p = DimPlot(alldata, reduction = "umap")
  return(list(data.clustered = alldata, umap_plot = p))
}


############## Take alldata and find DEGs
in_out_fraction = function(seu_obj, sample_name, group, reference, interesting_genes, threshold=0, plus1=1) {
  "
    Input: 
    
      seu_obj is a anndata with the lognormalized counts
      reference is the string for the reference sample. 
      group is the string corresponding to the analysed subset of adata. 
      sample_name
      interesting_genes : is a list of var names corresponding to features in adata.
      plus1: quantity to add to numerator and denominator when calculating lof2fc 
      (should be set to 1 in case of low expressed genes)
      threshold: value below which a gene is considered not expressed
      
      return a pandas. For each gene it contains: log2 FC, 
                                                  fraction of cells in each group with expression > threshold,
                                                  mean for each group
    "
  
  print(group)
  
  fraction=data.frame(row.names = interesting_genes)
  
  # fraction.index.name='genes'
  
  if (reference != "rest") {
    gene_expression_ref = subset(seu_obj, idents = reference, features = interesting_genes)
  } else
    gene_expression_ref = subset(seu_obj, idents = setdiff(levels(seu_obj), group), features = interesting_genes)
  gene_expression_ref = gene_expression_ref@assays$RNA@data
  
  print("Reference:")
  print(dim(gene_expression_ref))
  ncells_ref = rep(ncol(gene_expression_ref), length(interesting_genes))
  #[gene_expression_ref.shape[0]] * len(interesting_genes)    
  
  gene_expression_group = subset(seu_obj, idents = group, features = interesting_genes)
  gene_expression_group = gene_expression_group@assays$RNA@data
  
  
  ncells_group = rep(ncol(gene_expression_group), length(interesting_genes))
  print("Group:")
  print(dim(gene_expression_group))
  fraction[,paste0('ratio_in_', group)] = rowSums(gene_expression_group > threshold) / ncells_group
  fraction[,paste0('ratio_in_', reference)] = rowSums(gene_expression_ref > threshold) / ncells_ref
  
  fraction[,paste0('mean_in_', group)] = rowMeans(expm1(gene_expression_group))
  fraction[,paste0('mean_in_', reference)] = rowMeans(expm1(gene_expression_ref))
  
  # Calculate lfc over the normalized gene expression 
  # (need to invert log with expm1)
  fraction[,'lfc'] = log2((rowMeans(expm1(gene_expression_group)) + plus1) / 
                            (rowMeans(expm1(gene_expression_ref)) + plus1)) 
  
  fraction['abs_lfc'] = abs(fraction['lfc'])
  
  return(fraction)
  
}


deg_analysis = function(seu_obj, key, group, reference) {
  # """
  #   Function that computes Differential Expression analysis and returns a DataFrame
  #   containing the result.
  # 
  #   Input:
  # 
  #   seu_obj: AnnData object containing log-normalized counts
  #   key: key in seu_obj.obs to consider 
  #   group: string containing the group on which we want to perform DEA
  #   reference: string containing the observation against which we want to compare group
  #   """
  if (reference != 'rest') {
    seu_obj_curr = subset(seu_obj, idents = c(group, reference))
  } else {
    seu_obj_curr = seu_obj
  }
  
  print(dim(seu_obj_curr))
  genesKeep = rownames(seu_obj_curr)[rowSums(seu_obj_curr@assays$RNA@counts) > 0]
  seu_obj_curr <- subset(seu_obj_curr, features = genesKeep)
  print(dim(seu_obj_curr))
  
  # sc.tl.rank_genes_groups(seu_obj_curr, 
  #                         groupby=key,
  #                         groups=[group], 
  #                         use_raw=False, 
  #                         rankby_abs=True, 
  #                         n_genes=seu_obj_curr.shape[1], 
  #                         method='wilcoxon',
  #                         corr_method='bonferroni',
  #                         pts = True,
  #                         tie_correct = True)
  # min.pct: nly test genes that are detected in a minimum fraction of min.pct cells 
  # in either of the two populations. Meant to speed up the function by not testing genes 
  # that are very infrequently expressed. Default is 0.1
  df = Seurat::FindMarkers(seu_obj_curr, 
                           ident.1 = group,
                           ident.2 = reference,
                           # cells.1 = WhichCells(seu_obj_curr, idents = group),
                           # cells.2 = WhichCells(seu_obj_curr, idents = reference),
                           logfc.threshold = 0, test.use = 'wilcox', min.pct = 0.2)
  
  #df = sc.get.rank_genes_groups_df(seu_obj_curr, group=group)
  colnames(df) = gsub(pattern = '1$', replacement = group, x = colnames(df))
  colnames(df) = gsub(pattern = '2$', replacement = reference, x = colnames(df))
  df_2 = in_out_fraction(seu_obj_curr, 
                         sample_name=key, 
                         group=group,
                         reference=reference, 
                         interesting_genes=rownames(seu_obj_curr),
                         plus1=1e-9)
  
  # df.set_index("names", inplace=True)
  ok_genes = intersect(rownames(df), rownames(df_2))
  df_final = cbind(df[ok_genes,], df_2[ok_genes,])
  df_final = df_final %>% arrange(desc("abs_lfc"))
  return(df_final)
}





meta_vars <- function(metadata_df, var_comparison_sel, var_strat_sel) {
  browser()
  var_comparison <- colnames(metadata_df)
  var_comparison
  
  #var_comparison_sel_values <- list()
  #if(!is.null(var_comparison_sel)) {
  #  if (var_comparison_sel=="")
  #    var_comparison_sel_values <- list()
  #  else  
      var_comparison_sel_values <- unique(metadata_df[,var_comparison_sel])
  #}
  
  all_var_strat <- lapply(as.list(metadata_df), function(x) x[!duplicated(x)])
  all_var_strat
  ##
  
  df0 <- metadata_df[,setdiff(var_comparison,var_comparison_sel), drop=FALSE]
  df0
  
  var_strat <- lapply(as.list(df0), function(x) x[!duplicated(x)]) #not needed
  var_strat
  
  var_non_strat <- colnames(df0)
  var_non_strat
  
  #if (is.null(var_strat_sel))
  #  df1 <- df0
  #else
  
  #var_strat_sel <- setdiff(var_strat_sel, var_strat)
  
  var_strat_remove <- var_strat_sel[var_strat_sel %in% var_comparison_sel_values ]
  var_strat_sel <- setdiff(var_strat_sel,var_strat_remove)  #if by chance the var_comparison_sel is changed so to include a stratified variable then it is removed from var_strat_sel
  
  
    df1 <- df0[rowSums(sapply(df0, `%in%`, var_strat_sel)) >= length(var_strat_sel), , drop=FALSE]
  df1 #df1 contains only rows with all the var_strat_sel
  
  
  # var_non_strat_null contains groups of not stratifiable variables because there is only one value associated to those groups
  var_non_strat_null <- colnames(df1[,sapply(df1, function(x) length(x[!duplicated(x)]))==1, drop=FALSE])
  var_non_strat_null
  

  # var_non_strat are the variables belonging to group names not yet stratified 
  var_non_strat <- setdiff(var_non_strat, var_non_strat_null)
  var_non_strat
  
  df2 <- df1[, var_non_strat, drop=FALSE]
  df2 #df2 contains rows with non stratified variables
  
  var_strat <- lapply(as.list(df2), function(x) x[!duplicated(x)])
  var_strat <- c(var_strat,var_strat_sel)
  
  external_res <- expand.grid(c(list(unique(metadata_df[,var_comparison_sel])), var_strat_sel))
  if (is.null(var_strat_sel))
    external_res <- data.frame(list(unique(metadata_df[,var_comparison_sel,drop=FALSE])))
  if (is.null(var_comparison_sel))  
    external_res <- data.frame(list(var_strat_sel))
  
    
  if (!is.null(var_comparison_sel))  
    colnames(external_res) <- c(var_comparison_sel, var_strat_sel)
  external_res
  
  if (length(var_strat_sel)+length(var_comparison_sel)>=2)
    com <- as.data.frame(combn(as.data.frame(t(as.matrix(external_res))),2))
  else
    com <- external_res
  
  
  xx <- c(metadata_df= list(df2), var_comparison = list(var_comparison), var_comparison_sel=var_comparison_sel, var_comparison_sel_values=list(var_comparison_sel_values), all_var_strat=list(all_var_strat), var_strat=list(var_strat), var_strat_sel=list(var_strat_sel), var_non_strat=list(var_non_strat))
  return(xx)
  
}

# # remove all objects that will not be used.
# rm(s_data,tmp_expression)
# 
# # run garbage collect to free up memory
# gc()


### end of file -- cosdeg.R
