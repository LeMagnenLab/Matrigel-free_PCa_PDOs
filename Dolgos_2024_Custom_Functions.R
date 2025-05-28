# Load needed files from previous experiment 
get_exp_file_path = function(root = "/scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/Projects", organ, project, samples_ID, prev_exp, pattern, last = T){
  
  sample_dir <- paste0(root, "/", organ, "/", project, "/exp/", samples_ID, "/")
  prev_exp_dir <- paste0(sample_dir, prev_exp)
  
  print(paste0("Searching file in ", prev_exp_dir))
  
  list_files = list.files(
    path = prev_exp_dir, 
    pattern = pattern, 
    full.names = T)
  
  if (last == T) {
    
    file = list_files[length(list_files)]
    
    return(file)
    
  }else{ return(list_files) }
  
  # Return the file path 
  # Then user choose the appropriate function to load the file depending on the type
  
}

## Create output path in exp folder
create_exp_folder <- function(  root = "/scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/Projects/", organ, project, samples_ID, exp) {
  
  sample_dir <- paste0(root, "/", organ, "/", project, "/exp/", samples_ID, "/")
  exp_dir <- paste0(sample_dir, exp, "/")
  
  if (dir.exists(sample_dir)) { # Test if sample folder is already existing
    if (dir.exists(exp_dir)) { # If yes test if experiment folder is already existing
      print(paste("The exp", exp, "already exists for the sample", samples_ID)) # If yes says so
    } else {
      dir.create(path = exp_dir) # If no creates it
    }
  } else {
    dir.create(path = sample_dir) # If no sample folder existing, creates it, then creates the experiment folder
    dir.create(path = exp_dir)  }
  
  return(exp_dir)
}

# Time stamp 
time_stamp = function(x){
  
  time = as.character(Sys.time())
  
  y = substr(time,3,4)
  m = substr(time,6,7)
  d = substr(time,9,10)
  h = substr(time,12,13)
  min = substr(time,15,16)
  
  return(paste0(y,m,d,"_",h,"h",min,"_"))
  
}

########################################
### Relative to scRNA-seq data wrangling
#######################################

## Harmonize features between sce from a list
common_features_intersection <- function(sce_list) {
  Reduce(intersect, lapply(sce_list, function(x) rownames(x)))
}

union_rowData_genes <- function(sce_list) {
  unique(unlist(lapply(sce_list, function(x) rowData(x)$SYMBOL)))
}

restrict_reorder_features_names <- function(sce_list) {
  common_genes <- common_features_intersection(sce_list)
  lapply(sce_list, function(x) {
    x[common_genes, ]
  })
}

extend_and_reorder_SCE <- function(sce_list) {
  
  # Get the union of all genes across SCE objects, maintaining the order
  union <- union_rowData_genes(sce_list)
  
  # Prepare loop objects
  list_sce_extended = list()
  sce_id = 1
  
  for (sce in sce_list) {
    
    print(paste("extending sce :", names(sce_list)[sce_id]))
    
    
    # Find which genes are missing in the current SCE object
    missing_genes <- setdiff(union, rowData(sce)$SYMBOL)
    
    # If there are missing genes, add them with 0 values
    if (length(missing_genes) > 0) {
      # Create a matrix of 0s for missing genes
      # Colnames are the names of the cells
      # rownames are the names of missing genes
      missing_data <- matrix(0, nrow = length(missing_genes), ncol = ncol(assay(sce)),
                             dimnames = list(missing_genes, colnames(assay(sce))))
      
      # Combine the original data with the missing data
      new_data <- rbind(counts(sce), missing_data)
    } else {
      new_data <- counts(sce) # If no missing genes then no chane to counts(sce)
    }
    
    # Reorder the combined data to match the order of all_genes
    # This ensures that the rows are in the same order as the union of all genes
    # As match returns a vector of the positions of (first) matches of its first argument in its second.
    new_data <- new_data[match(union, rownames(new_data)), ]
    
    # Combine the slots to create a new combined sce object
    # Updating it by just adding the new counts doesn't work as it doesn't match rowData length
    # No need to keep reduced_dims
    sce_extended <- SingleCellExperiment( 
      assays = list(counts = new_data),  
      rowData = data.frame(SYMBOL = rownames(new_data)), # RowData is similar among all elements 
      colData = colData(sce)
    ) 
    
    list_sce_extended[[sce_id]] = sce_extended
    
    print(paste("sce :", names(sce_list)[sce_id], "added to the list_sce_extended"))
    
    # Return the updated and reordered SCE object
    sce_id = sce_id + 1    
    
  }
  
  names(list_sce_extended) = names(sce_list)
  return(list_sce_extended)
  
}

########################################
### Relative to scRNA-seq signature scores
#######################################

# Extract all the gene names comprised in a signature presen in msigDB
######################################################################

# User give a proxy to look for a signature
# If many signature match witht the proxy, all matches will be displayed and user needs to refine

get_signature_data = function(pathway_proxy, MSigDB_category){
  
  library(msigdbr)
  library(clusterProfiler)
  
  h_t2g <- msigdbr(species = "Homo sapiens", category = MSigDB_category) %>% 
    dplyr::select(gs_name, entrez_gene)
  
  h_t2g = as_tibble(h_t2g)
  
  # Proxy for the targeted pathway
  h_t2g_hit = stringr::str_detect(h_t2g$gs_name, pathway_proxy)
  
  # Select only 
  h_t2g_hit = h_t2g %>%
    dplyr::filter(h_t2g_hit)
  
  # Number of SYMBOL ID is smaller than ENTREZ as there is replicates in h_t2g_
  # If some genes are not found there is a warning message
  
  h_t2g_hit_names = unique(h_t2g_hit$gs_name)
  
  if (length(h_t2g_hit_names) > 1) {
    
    warning(paste0("Multiple hit were detected with the pathway proxy: ", pathway_proxy))
    warning(paste0("Hit with proxy used :", h_t2g_hit_names, "\n \t "))
    
    stop(paste0("Refine your proxy to match only one the above mentionned hits"))
    
  }
  
  signature_data <- bitr(h_t2g_hit$entrez_gene, fromType = "ENTREZID",
                         toType = c("SYMBOL"),# UNIPROTID will be used later with KEGG
                         OrgDb = org.Hs.eg.db) %>%
    mutate(NAME =  h_t2g_hit_names)
  
  
  return(signature_data)  
  
}

##### Calculate signature module score and add assays to seurat object ##########
#################################################################################

## Reminder to access specific assay slot : seurat_object[["assay_name"]]@data

calculate_signature_score = function(seurat_object, signature_name, signature_genes, assay_name){
  
  DefaultAssay(seurat_object) = assay_name
  
  ## Add module 
  seurat_object <- AddModuleScore(seurat_object,
                                  assay = assay_name,
                                  features = list(signature_genes), 
                                  name = "signature_score" , 
                                  search = T)
  
  ## The calculated score is stored in the meta.data of the seurat object under signature_score1 column name
  # seurat_object@meta.data$signature_score1
  
  ## Creates a new assay where the only feature is the score of the above mentioned signature
  new_assay <- paste0(signature_name, "_score_only")
  
  seurat_object[[new_assay]] = CreateAssayObject(
    data = t(FetchData(
      object = seurat_object, 
      vars = 'signature_score1'))
  )
  
  
  # Creates a new assay where only the genes of the calculated signature are present
  new_assay <- paste0(signature_name, "_genes_only")
  
  seurat_object[[new_assay]] = CreateAssayObject(
    data = GetAssayData(seurat_object[signature_genes], assay_name))
  
  DefaultAssay(seurat_object) = assay_name
  
  return(seurat_object)
  
}

##### Add signature gene list (in rowData) and score (colData) to sce object  ###
#################################################################################

add_signature_score_to_sce = function(sce_object, seurat_object, signature_name, signature_genes, seurat_assay){
  
  ### Store the list of genes in the rowData matrix in the main_sce object
  column_name = paste0(signature_name, "_gene")
  rowSubset(sce_object, column_name) = signature_genes
  
  ### Check if genes have been added to the sce object
  check_genes = as_tibble(rowData(sce_object)[column_name])
  ind = which(check_genes == T)
  rowData(sce_object)$SYMBOL[ind]
  
  ### Store the log normed matrix as an alternative experiment in the sce object
  ## First extract the matrix from seurat object
  # Even though it's not a count matrix, it's necessary to add something in the count slot otherwise the conversion to seurat object will fail at the next step
  
  coldata_colname = paste0(signature_name, "_score")
  score = seurat_object[[seurat_assay]]@data
  sce_object$score = score # Need to create a temporary column with $ access otherwise it doesn't work
  colData(sce_object)[coldata_colname] = sce_object$score
  
  # Remove the temporary column names "score"
  ind_score_col = which(colnames(colData(sce_object)) == "score")
  colData(sce_object) = colData(sce_object)[,-ind_score_col]
  
  return(sce_object)
  
}

### Remove signature score outliers for color scale #####
#########################################################

get_signature_score_outliers = function(sce_object, sce_signature_column, signature_name, nb_outliers){
  
  # Histogram of signature score
  coldata_signature_ind = which(colnames(colData(sce_object)) == sce_signature_column)
  score = colData(sce_object)[coldata_signature_ind][,1] # Retrieves score for each cells (1st column)
  
  hist(score, breaks = 100, 
       main = paste0("Signature score distribution",
                     "\n",
                     unique(signature_name)))
  
  # Discard the outliers cells 
  five_min = head(order(score), n = nb_outliers)
  five_max = tail(order(score), n = nb_outliers)
  
  outliers = list(
    index = c(five_min = head(order(score), n = nb_outliers), five_max = tail(order(score), n = nb_outliers)),
    values = c(min_score = min(score[-five_min]), max_score = max(score[-five_max]))
  ) 
  
  return(outliers)
  
}


# Vilin plot with stat comparison on the mean
Violin_signature_score_stats <- function(sce_object, genes, signature_name, group1, group2, group1_color, group2_color, min_max_outliers, 
                                         out_path, file_suffix, 
                                         compute_stats = FALSE,
                                         stat_comparisons) {
  
  score_metadata_col <- paste0(signature_name, "_score")
  
  # Extracting signature scores from colData
  if (!is.null(group1) && !is.null(group2)) {
    
    col_nb <- c(which(colnames(colData(sce_object)) == score_metadata_col), 
                which(colnames(colData(sce_object)) == group1),
                which(colnames(colData(sce_object)) == group2))
    
    gene_exprs <- as_tibble(colData(sce_object)[col_nb])
    
    # Basic plot
    plot <- ggplot(gene_exprs[-min_max_outliers$index,], aes(x = .data[[group2]], y = .data[[score_metadata_col]], fill = .data[[group1]])) +
      geom_violin(alpha = 0.8, position = position_dodge(width = 0.8), scale = "width") +
      geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.5) +
      scale_fill_manual(values = group1_color, name = group1) +
      guides(fill = "none") +
      theme_minimal() + 
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),  # Customize axis title appearance
        axis.text = element_text(size = 10),    # Customize axis text appearance
        axis.line = element_line(size = 0.5),   # Customize axis line appearance
        axis.ticks = element_line(size = 0.5)) +
      labs(
        title = signature_name,
        x = NULL,
        y = "Signature score")
    
    # Add statistical comparisons if compute_stats is TRUE
    if (compute_stats) {
      plot <- plot +
        stat_compare_means(
          comparisons = stat_comparisons, # Define the pairwise comparison
          method = "wilcox.test",                    # Statistical test to use
          label = "p.signif",                   # Display significance stars
          size = 5        
        )
    }
    
  } else if (is.null(group2)) {
    
    col_nb <- c(which(colnames(colData(sce_object)) == score_metadata_col), 
                which(colnames(colData(sce_object)) == group1))
    
    gene_exprs <- as_tibble(colData(sce_object)[col_nb])
    
    # Basic plot
    plot <- ggplot(gene_exprs[-min_max_outliers$index,], aes(x = .data[[group1]], y = .data[[score_metadata_col]], fill = .data[[group1]])) +
      geom_violin(alpha = 0.8, position = position_dodge(width = 0.8), scale = "width") +
      geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.5) +
      scale_fill_manual(values = group1_color, name = group1) +
      guides(fill = "none") +
      theme_minimal() + 
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),  # Customize axis title appearance
        axis.text = element_text(size = 10),    # Customize axis text appearance
        axis.line = element_line(size = 0.5),   # Customize axis line appearance
        axis.ticks = element_line(size = 0.5)) +
      labs(
        title = signature_name,
        x = NULL,
        y = "Signature score")
    
    # Add statistical comparisons if compute_stats is TRUE
    if (compute_stats) {
      plot <- plot +
        stat_compare_means(
          comparisons = stat_comparisons, # Define the pairwise comparison
          method = "wilcox.test",                    # Statistical test to use
          label = "p.signif",                   # Display significance stars
          size = 5        
        )
    }
    
  }
  
  # Save the plot to a PDF file
  ggsave(
    plot = plot,
    filename = paste0(out_path, time_stamp(), "plots_Violin_Scores", file_suffix, ".pdf"),
    device = "pdf",
    width = 210,
    height = 150,
    units = "mm"
  )
}


###############################
### Check gene name aliases ###
###############################

check_aliase = function(sce, gene){
  
  library(org.Hs.eg.db)
  
  # All possible aliases existing
  all_aliases = AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = "ALIAS")
  
  # Testing first if the gene is in the dataset
  if (gene %in% rownames(sce)) {
    
    #print(paste("The gene", gene, "has been found in the dataset with this symbol."))
    
  } else{ # If not, say it and first test if there is any aliase existing 
    
    print(paste("The gene", gene, "HAS NOT BEEN FOUND in the dataset with this symbol."))
    
    ## Test if there are no aliases existing, if no are existing check spelling of the gene
    if (sum(gene %in% all_aliases$ALIAS) == 0){
      
      print(
        paste(
          "No aliase of", gene, "were found. Check the spelling of your gene"
        )
      )
    }
    
    else{# If not, store all the aliases existing in a vector
      
      symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = gene, keytype = "ALIAS", columns = "SYMBOL")
      symbol = symbol$SYMBOL[1] # Rarely, many symbols are associated to one aliase, juste take the first symbol
      
      print(paste("The gene", gene, "main name is", symbol, ". Searchnig now for all the aliases of this gene."))
      
      aliases <- AnnotationDbi::select(org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = c("ALIAS"))
      aliases = aliases$ALIAS
      
      print(paste("Aliase found for", gene, ":", aliases))
      
      # Test if any of those aliases is in our dataset, if not tell the user that this gene is not expressed
      if (sum(aliases %in% rownames(sce)) == 0) {
        
        print(paste("Unfortunately, none of the aliases were found in the dataset", gene, "is certainly not expressed."))
        
      }else{ # If yes, tell the user and change the name
        
        aliases_in =  aliases %in% rownames(sce)
        
        alias_in = aliases[which(aliases_in == TRUE)]
        
        if (length(alias_in) > 1) { 
          
          print(paste0("More than one aliase has been found in the dataset, please choose manually and change the name of the gene in the gene_set_table"))
          
        }else{
          
          print(
            paste(
              "The alias", alias_in, "has been found in the dataset.", 
              gene,"is now changed for", alias_in, "for the rest of the analysis"
            )
          ) 
          
          gene = alias_in
          
        }
        
        
      }
      
    }
    
    
  }
  
  # If no aliase found, then gene remanin the same, an extra step is applied in the main script to remove it from the analysis
  return(gene)
  
  
}


###########################################
#### violin split function
###########################################



GeomSplitViolin <- ggplot2::ggproto(`_class` = "GeomSplitViolin",`_inherit` =  GeomViolin,
                                    draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                      data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                      grp <- data[1, "group"]
                                      newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                      newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                      newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                                      
                                      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                                  1))
                                        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                        aesthetics$alpha <- rep(1, nrow(quantiles))
                                        both <- cbind(quantiles, aesthetics)
                                        quantile_grob <- GeomPath$draw_panel(both, ...)
                                        ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                      }
                                      else {
                                        ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                      }
                                    })



geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


################################
### Downsample a sce object ###
###############################

# Downsample function
downsample_sce <- function(sce, group_column, n_per_group) {
  
  # Check if group_column exists in colData
  if (!group_column %in% colnames(colData(sce))) {
    stop("The specified group_column does not exist in the metadata.")
  }
  
  # Extract metadata and create a downsampling index
  metadata_df <- as.data.frame(colData(sce))
  
  # Downsample the cells per group using dplyr
  sampled_metadata <- metadata_df %>%
    mutate(row_index = row_number()) %>%
    group_by(!!sym(group_column)) %>%
    sample_n(size = min(n_per_group, n()), replace = FALSE) %>%
    arrange(desc(Culture_Condition)) %>%
    ungroup()
  
  # Subset the SCE object based on the sampled cells
  downsampled_sce <- sce[, sampled_metadata$row_index]
  
  return(downsampled_sce)
}
