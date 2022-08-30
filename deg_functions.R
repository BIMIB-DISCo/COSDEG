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
  #   reference: string containg the observation against which we want to compare group
  #   """
  if (reference != 'rest') {
    current_cells = WhichCells(alldata, idents = "H130001.out")
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








