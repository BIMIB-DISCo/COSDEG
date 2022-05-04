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
  douplets_pN <- 0.25
  
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
    
    x <- plot(x = bc.mvn$ParamID, y = bc.mvn$MeanAUC, pch = 18, 
              col = "black", cex = 0.75, xlab = NA, ylab = NA)
    x <- lines(x = bc.mvn$ParamID, y = bc.mvn$MeanAUC, col = "black", lty = 2)
    plots$AUC[[d]] <- x 
    
    par(new = TRUE)
    x <- plot(x = bc.mvn$ParamID, y = bc.mvn$BCmetric, pch = 16, col = "#41b6c4", cex = 0.75)
    axis(side = 4)
    x <- lines(x = bc.mvn$ParamID, y = bc.mvn$BCmetric, col = "#41b6c4")
    plots$BCmetric[[d]] <- x
    
    plots$parameter_estimation[[d]] <- barplot(bcmvn[[d]]$BCmetric, names.arg = bcmvn[[d]]$pK, las=2)
  }
  
  for (d in dataset_folders) {
    pK <- as.numeric(as.character(bcmvn[[d]]$pK))
    BCmetric <- bcmvn[[d]]$BCmetric
    pK_choose[[d]] <- pK[which(BCmetric%in%max(BCmetric))]
  }
  
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







# # remove all objects that will not be used.
# rm(s_data,tmp_expression)
# 
# # run garbage collect to free up memory
# gc()


### end of file -- cosdeg.R
