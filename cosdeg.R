webpath <- "https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/"
dir.create("./data/raw", recursive = T)

file_list <- c("Normal_PBMC_13.h5", "Normal_PBMC_14.h5", "Normal_PBMC_5.h5", "nCoV_PBMC_15.h5",
               "nCoV_PBMC_17.h5", "nCoV_PBMC_1.h5")
for (i in file_list) {
  download.file(url = paste0(webpath, i), destfile = paste0("./data/raw/", i))
}

install_github("https://github.com/chris-mcginnis-ucsf/DoubletFinder", upgrade = F)




suppressMessages(require(Seurat))
suppressMessages(require(Matrix))

suppressMessages(require(DoubletFinder))

dataset_folders <- list.dirs(path = "data/", full.names = F, recursive = F)
dataset_folders <- dataset_folders[-5] # remove raw folder

summary_qc <- data.frame()

s_data <- list()
data.filt <- list()
s_data_dim <- list()

for (d in dataset_folders) { 
  data_dir <- file.path("data", d, "filtered_feature_bc_matrix")
  print(data_dir)
  tmp_expression <- Seurat::ReadMtx(file.path(data_dir,"matrix.mtx.gz"), features = file.path(data_dir,"features.tsv.gz"), cells = file.path(data_dir,"barcodes.tsv.gz"))
  s_data[[d]] <- CreateSeuratObject(counts = tmp_expression, project = d)
  s_data_dim[[d]] <- list(dim(s_data[[d]]))
  summary_qc["preQC",d] <- s_data_dim[d]
  
  # Way1: Doing it using Seurat function
  s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "^MT-", col.name = "percent_mito")
  s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "^RP[SL]", col.name = "percent_ribo")
  # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
  s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "^HB[^(P)]", col.name = "percent_hb")
  s_data[[d]] <- PercentageFeatureSet(s_data[[d]], "PECAM1|PF4", col.name = "percent_plat")
  
}

summary_qc

# Merge datasets into one single seurat object
alldata <- merge(s_data[[1]], s_data[-1], add.cell.ids = names(s_data))

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()






selected_f <- list()
for (d in dataset_folders) {
  selected_f[[d]] <- rownames(s_data[[d]])[Matrix::rowSums(s_data[[d]]) > 0]
  data.filt[[d]] <- subset(s_data[[d]], features = selected_f[[d]])
  #s_data_dim[[d]] <- length(selected_f[[d]])
  s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  summary_qc["zero_genes",d] <- s_data_dim[d]
  
  print(dim(data.filt[[d]]))
}

summary_qc


#selecting cell with less than 90% of zeros
#alldata[1:26,1:20]@assays$RNA@counts
#colSums(alldata[1:26,1:20]@assays$RNA@counts==0)/26
selected_c <- list()
perc_zeros <- 0.90
for (d in dataset_folders) {
  selected_c[[d]] <- colnames(data.filt[[d]])[data.filt[[d]]$nFeature_RNA/dim(data.filt[[d]])[1]>=1-perc_zeros]
  data.filt[[d]] <- subset(data.filt[[d]], cells = selected_c[[d]])
  #s_data_dim[[d]] <- length(selected_c[[d]])
  s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  summary_qc["empty_cells",d] <- s_data_dim[d]
  
  print(dim(data.filt[[d]]))
}

summary_qc

perc_mitoc <- list()
for (d in dataset_folders) {
  # Way1: Doing it using Seurat function
  data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^MT-", col.name = "percent_mito")
  data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^RP[SL]", col.name = "percent_ribo")
  # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
  data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "^HB[^(P)]", col.name = "percent_hb")
  data.filt[[d]] <- PercentageFeatureSet(data.filt[[d]], "PECAM1|PF4", col.name = "percent_plat")
  
  perc_mitoc[[d]] <- mean(data.filt[[d]]$percent_mito/100.) + sd(data.filt[[d]]$percent_mito/100.)
  
  selected_c[[d]] <- colnames(data.filt[[d]])[data.filt[[d]]$percent_mito/100 <= perc_mitoc[[d]]]
  
  data.filt[[d]] <- subset(data.filt[[d]], cells = selected_c[[d]])
  #s_data_dim[[d]] <- length(selected_c[[d]])
  s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  summary_qc["mitochodrial",d] <- s_data_dim[d]
  print(dim(data.filt[[d]]))
}

summary_qc

for (d in dataset_folders) {
  # Compute the relative expression of each gene per cell Use sparse matrix
  # operations, if your dataset is large, doing matrix devisions the regular way
  # will take a very long time.
  par(mar = c(4, 8, 2, 1))
  C <- data.filt[[d]]@assays$RNA@counts
  C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
  boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
          col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
}

# Merge datasets into one single seurat object
alldata <- merge(data.filt[[1]], data.filt[-1], add.cell.ids = names(data.filt))

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()



for (d in dataset_folders) {
  print(d)
  # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  suppressMessages(require(DoubletFinder))
  
  data.filt[[d]] = NormalizeData(data.filt[[d]])
  
  data.filt[[d]] = FindVariableFeatures(data.filt[[d]], verbose = F)
  print(d)
  data.filt[[d]] = ScaleData(data.filt[[d]], vars.to.regress = c("nFeature_RNA", "percent_mito"),verbose = F)
  print(d)
  data.filt[[d]] = RunPCA(data.filt[[d]], verbose = F, npcs = 20)
  print(d)
  data.filt[[d]] = RunUMAP(data.filt[[d]], dims = 1:10, verbose = F)
  print(d)


  # Can run parameter optimization with paramSweep, but skip for now.
  sweep.res[[d]] <- paramSweep_v3(data.filt[[d]]) 
  sweep.stats[[d]] <- summarizeSweep(sweep.res[[d]], GT = FALSE) 
  bcmvn[[d]] <- find.pK(sweep.stats[[d]])
  barplot(bcmvn[[d]]$BCmetric, names.arg = bcmvn[[d]]$pK, las=2)
}

for (d in dataset_folders) {
  # define the expected number of doublet cellscells.
  nExp <- round(ncol(data.filt[[d]]) * 0.16)  # expect 4% doublets
  data.filt[[d]] <- doubletFinder_v3(data.filt[[d]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)


  # name of the DF prediction can change, so extract the correct column name.
  DF.name = colnames(data.filt[[d]]@meta.data)[grepl("DF.classification", colnames(data.filt[[d]]@meta.data))]



  cowplot::plot_grid(ncol = 2, DimPlot(data.filt[[d]], group.by = "orig.ident") + NoAxes(),
                     DimPlot(data.filt[[d]], group.by = DF.name) + NoAxes())
  
  VlnPlot(data.filt[[d]], features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  
  s_data_dim[[d]] <- list(dim(data.filt[[d]]))
  summary_qc["doublets",d] <- s_data_dim[d]
}


summary_qc

#########




# add metadata
#sdata.cov1$type = "Covid"
#sdata.cov15$type = "Covid"
#sdata.cov17$type = "Covid"
#sdata.ctrl5$type = "Ctrl"
#sdata.ctrl13$type = "Ctrl"
#sdata.ctrl14$type = "Ctrl"


# Merge datasets into one single seurat object
alldata <- merge(s_data[[1]], s_data[-1], add.cell.ids = names(s_data))

# remove all objects that will not be used.
rm(s_data,tmp_expression)

# run garbage collect to free up memory
gc()


# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")


# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")

alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")


feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()



FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)



#selecting expressed genes
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 0]


#selecting cell with less than 90% of zeros
#alldata[1:26,1:20]@assays$RNA@counts
#colSums(alldata[1:26,1:20]@assays$RNA@counts==0)/26
perc_zeros <- 0.95
#selected_c <- colnames(alldata)[Matrix::colSums(alldata@assays$RNA@counts==0)/dim(alldata)[1]<perc_zeros]
selected_c <- colnames(alldata)[alldata$nFeature_RNA/dim(alldata)[1]>=1-perc_zeros]

data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)


