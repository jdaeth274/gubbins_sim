###############################################################################
## Exploring the point predictive value in gubbins ############################
###############################################################################

require(dplyr)
require(stringr)
require(ggplot2)
require(ape)
require(ggpubr)
require(data.table)
require(treespace)
require(viridis)
require(tictoc)

###############################################################################
## Functions ##################################################################
###############################################################################

recombination_gff_cleaner <- function(in_gff_pmen3){
  #browser()
  ## Takes delim file read in and cleans it up 
  pmen3_reccy_data <- data.frame(data = (matrix(data = NA, nrow = nrow(in_gff_pmen3),
                                                ncol = 9)))
  colnames(pmen3_reccy_data) <- c("isolate","number_isolates", "snp_count",
                                  "length","density", "start", "end", "start_node",
                                  "end_node")
  
  for( k in 1:nrow(in_gff_pmen3)){
    #browser()
    current_row <- in_gff_pmen3[k,]
    current_atty <- current_row$attribute
    current_atty_taxa_snp <- str_split_fixed(current_atty,"taxa=",2)[,2]
    taxa <- str_split_fixed(current_atty_taxa_snp,";snp_count=",2)[,1]
    nodes <- str_split_fixed(sub("node=","",str_split_fixed(current_atty, ";neg",2)[,1]), "->", 2)
    taxa_split <- str_split_fixed(taxa, "\\s",10)
    taxa_first <- taxa_split[1,which(taxa_split[1,] != "")[1]]
    #taxa_first <- sub("\\s*","",taxa_first)
    snp_count <- str_split_fixed(current_atty, ";snp_count=",2)[,2]
    snp_count <- as.numeric(sub(";","",snp_count))
    num_taxa <- str_count(taxa, pattern = "GCA") # This has to change wrt to dataset (can't think of a way round atm)
    length_of_event <- current_row$end - current_row$start + 1
    snp_density <- snp_count / length_of_event
    
    pmen3_reccy_data$isolate[k] <- taxa_first
    pmen3_reccy_data$number_isolates[k] <- num_taxa
    pmen3_reccy_data$snp_count[k] <- snp_count
    pmen3_reccy_data$length[k] <- length_of_event
    pmen3_reccy_data$density[k] <- snp_density
    pmen3_reccy_data$start[k] <- current_row$start
    pmen3_reccy_data$end[k] <- current_row$end 
    pmen3_reccy_data$start_node[k] <- nodes[1,1]
    pmen3_reccy_data$end_node[k] <- nodes[1,2]
    
  }
  return(pmen3_reccy_data)
}

delim_reader <- function(delim_file){
  #browser()
  ## read in the gff recombination events file from gubbins
  reccy_csv <- "No"
  try({
    reccy_csv <- read.delim(delim_file,header = FALSE, comment.char = "#")
    colnames(reccy_csv) <- c("type","prog","class","start","end","trent","alexander","arnold","attribute")
  }, silent = TRUE)
  return(reccy_csv)
}

typing_lapply_func <- function(row, reccy_preds){
  snp_type <- "S"
  potential_rec <- reccy_preds %>%
    filter((start_node == row[1]) &  (end_node == row[2]))
  if(nrow(potential_rec) > 0){
    
    if(any(data.table::between(as.integer(row[5]), potential_rec$start, potential_rec$end)))
      snp_type <- "r"
  }
  
  return(snp_type)
}

typing_gubbins_rec_apply <- function(reccy_preds, branch_base){
    snp_types <- apply(X = branch_base, MARGIN = 1, FUN = typing_lapply_func, 
                            reccy_preds = reccy_preds)
    branch_base$Type <- snp_types
    return(branch_base)
}



typing_gubbins_rec <- function(reccy_preds, branch_base){
  ## Function to take in branch base and classify each snp,
  ## either r for in a putative recombination event or S for a clonal frame mutation
  #browser()
  
    branch_base$Type <- "S"
    cat("\n")
    num_zeros <- nchar(nrow(branch_base))
    num_rows <- nrow(branch_base)
    for(k in 1:nrow(branch_base)){
      nchar_k <- nchar(k)
      nchar_0 <- num_zeros - nchar_k
      cat("\r", "Completed ", rep(0, nchar_0), k, " of ", num_rows, " SNPs", sep = "")
      current_snp <- branch_base[k,]
      potential_rec <- reccy_preds %>%
        filter((start_node == current_snp$start_node) &  (end_node == current_snp$end_node))
      if(nrow(potential_rec) > 0){
        
        if(any(data.table::between(current_snp$base_number, potential_rec$start, potential_rec$end)))
           branch_base$Type[k] <- "r"
      }
      
    }
  
  return(branch_base)
}

tip_nodes <- function(tree, node_index){
  ## Function to get all the tip nodes from a start node 
  #browser()
  tot_paths <- ape::nodepath(tree)
  tips <- NULL
  for(k in 1:length(tot_paths)){
    if(node_index %in% tot_paths[[k]]){
      tips <- append(tips, tot_paths[[k]][length(tot_paths[[k]])])
    }
  }
  return(tips)
}


append_ancestors <- function(gub_tree, branch_base){
  ## Function to find all the ancestors of an end Node and append these to a branch_base to match the 
  ## output from the simulator 
  
  branch_base$Taxa <- ""
  branch_base$end_node_index <- NA
  cat("\n")
  num_zeros <- nchar(nrow(branch_base))
  num_rows <- nrow(branch_base)
  tot_nodes <- c(gub_tree$tip.label, gub_tree$node.label)
  for(k in 1:nrow(branch_base)){
    nchar_k <- nchar(k)
    nchar_0 <- num_zeros - nchar_k
    cat("\r", "Completed ", rep(0, nchar_0), k, " of ", num_rows, " SNPs", sep = "")
    current_snp <- branch_base[k,]
    if(grepl("taxon", current_snp$end_node)){
      branch_base[k,"Taxa"] <- sub("taxon_","",current_snp$end_node)
    }else{

      start_node <- which(tot_nodes == current_snp$end_node)
      tips <- tip_nodes(gub_tree, start_node)
      tips_format <- paste(tips, collapse = ",")
      branch_base[k,"Taxa"] <- sub("taxon_","",tips_format)
    }
    node_index <- which(tot_nodes == current_snp$end_node)
    branch_base$end_node_index[k] <- node_index
  }
    return(branch_base)

}

taxa_ancestors <- function(tree, taxa_names){
  # Function to take in a list of taxa and find their MRCA
  
  root_index <- ape::nodepath(tree)[[1]][1]
  initial_path <- ape::nodepath(tree, root_index, taxa_names[1])
  for(k in 2:length(taxa_names)){
    next_path <- ape::nodepath(tree, root_index, taxa_names[k])
    overlaps <- which(initial_path %in% next_path)
    initial_path <- initial_path[overlaps]
  }
  
  mrca_index <- initial_path[length(initial_path)]
  
  return(mrca_index)
  
}

simul_node_finder <- function(simul_tree, simul_data){
  ## Function to turn the simul taxa with mutation into a node identifier for 
  ## comparison to Gubbins results 
  
  #browser()
  simul_data$end_node_index <- NA
  cat("\n")
  num_zeros <- nchar(nrow(simul_data))
  num_rows <- nrow(simul_data)
  tot_nodes <- c(simul_tree$tip.label,as.character(seq(Ntip(simul_tree) + 1, length.out = Nnode(simul_tree))))
  for(k in 1:nrow(simul_data)){
    nchar_k <- nchar(k)
    nchar_0 <- num_zeros - nchar_k
    cat("\r", "Completed ", rep(0, nchar_0), k, " of ", num_rows, " SNPs", sep = "")
    
    current_snp <- simul_data[k,]
    num_taxa <- str_count(current_snp$Taxa,",")
    if(num_taxa == 0){
      
      taxa_num <- current_snp$Taxa
      simul_data$end_node_index[k] <- which(tot_nodes == paste("taxon_",taxa_num, sep = ""))
      
    }else{
      taxa_list <- str_split_fixed(current_snp$Taxa, pattern = ",", n = num_taxa + 1)[1,]
      taxa_names <- paste("taxon_",taxa_list, sep = "")
      tip_indexes <- which(tot_nodes %in% taxa_names)
      mrca_index <- taxa_ancestors(simul_tree, tip_indexes)
      simul_data$end_node_index[k] <- mrca_index
      
    }
  }
  
  return(simul_data)
}

node_dictionary <- function(gubbins_tree, simul_tree){
  ## Function to create a mapping of the gubbins nodes over to the 
  ## simul nodes in the tree, trees here have to be identical in 
  ## topology of nodes 
  #browser()
  gubbins_nodes_df <- as.data.frame(matrix(nrow = Nnode(gubbins_tree) + Ntip(gubbins_tree), ncol = 3))
  colnames(gubbins_nodes_df) <- c("node_name_gubb","node_index_gubb","taxa")
  gubbins_nodes_df$node_name_gubb <- c(gubbins_tree$tip.label,gubbins_tree$node.label)
  gubbins_nodes_df$node_index_gubb <- seq(1, length.out = (Nnode(gubbins_tree) + Ntip(gubbins_tree)))
  gubbins_tot <- c(gubbins_tree$tip.label, gubbins_tree$node.label)
  for(k in 1:nrow(gubbins_nodes_df)){
    current_row <- gubbins_nodes_df[k,]
    if(k <= Ntip(gubbins_tree)){
      gubbins_nodes_df$taxa[k] <- gubbins_nodes_df$node_name_gubb[k]
    }else{
      tips <- tip_nodes(gubbins_tree, current_row$node_index)
      tips_taxon <- sort(gubbins_tot[tips])
      taxon_col <- paste(tips_taxon, collapse = "-")
      gubbins_nodes_df$taxa[k] <- taxon_col
    }
  }
  if(is.null(simul_tree$node.label))
    simul_tree <- makeNodeLabel(simul_tree)
  simul_nodes_df <- as.data.frame(matrix(nrow = Nnode(simul_tree) + Ntip(simul_tree), ncol = 3))
  colnames(simul_nodes_df) <- c("node_name_simul","node_index_simul","taxa")
  simul_nodes_df$node_name_simul <- c(simul_tree$tip.label, simul_tree$node.label)
  simul_nodes_df$node_index_simul <- seq( 1, length.out = Nnode(simul_tree) + Ntip(simul_tree))
  simul_tot <- c(simul_tree$tip.label, simul_tree$node.label)
  for(k in 1:nrow(simul_nodes_df)){
    current_row <- simul_nodes_df[k,]
    if(k <= Ntip(simul_tree)){
      simul_nodes_df$taxa[k] <- simul_nodes_df$node_name_simul[k]
    }else{
      tips <- tip_nodes(simul_tree, current_row$node_index)
      tips_taxon <- sort(simul_tot[tips])
      taxon_col <- paste(tips_taxon, collapse = "-")
      simul_nodes_df$taxa[k] <- taxon_col
    }
  }
  
  ## total df with cross links
  tot_df <- left_join(gubbins_nodes_df, simul_nodes_df, by = "taxa")
  
  return(tot_df)
}

data_cleaner <- function(gubbins_gff, gubbins_branch_base_csv, gubbins_tree_file,
                         simul_summary_file, simul_tree_file){
  ## Function to load up the neccessary data to check for the Gubbins
  ## SNPs predictions 
  browser()
  gubbins_reccy_gff <- delim_reader(gubbins_gff)
  gubbins_reccy_csv <- recombination_gff_cleaner(gubbins_reccy_gff)
  
  gubbins_branch_base <- read.csv(gubbins_branch_base_csv,
                                  stringsAsFactors = FALSE)
  gubbins_tree <- read.tree(file = gubbins_tree_file)
  
  branch_base_typed <- typing_gubbins_rec(gubbins_reccy_csv, gubbins_branch_base)
  
  branch_base_typed_taxa <- append_ancestors(gubbins_tree, branch_base_typed)
  
  ###############################################################################
  ## Now we have to read in the simul data and link it to start and end nodes ###
  ###############################################################################
  
  simul_summary <- read.table(simul_summary_file, sep = "\t",
                              comment.char = "", header = TRUE)
  
  simul_tree <- read.tree(simul_tree_file)
  
  simul_data_nodes <- simul_node_finder(simul_tree = simul_tree,
                                        simul_data = simul_summary)
  
  tree_distances <- treeDist(gubbins_tree, simul_tree, lambda = 0)
  
  if(tree_distances != 0){
    print("Warning the Tree topologies are not identical, SNP predictions will not 
          work properly")
  }
  ## Convert the gubbins node_indexes to match the simul data 
  
  node_diccers <- node_dictionary(gubbins_tree, simul_tree)
  
  branch_base_typed_taxa_simmed <- branch_base_typed_taxa %>%
    left_join(node_diccers %>% select(node_index_gubb, node_index_simul),
              by = c("end_node_index" = "node_index_gubb")) %>%
    select(-end_node_index) %>%
    rename(end_node_index = node_index_simul) %>%
    as.data.frame()
  
  
  simul_clonal_snps <- simul_data_nodes #%>% filter(Type == "S")
  
  return(list(gubbins_snps = branch_base_typed_taxa_simmed,
              simul_snps = simul_clonal_snps,
              gubbins_tree = gubbins_tree,
              simul_tree = simul_tree))
  
}

simul_apply_func <- function(row, gubbins_snps, taxa_step){
  #browser()
  tic(msg = "This long for one apply run")
  sim_gubbins <- "No"
  snp_diff <- -1
  gubbins_sim <- "No"
  gubbins_index <- 0
  if(taxa_step == "taxa"){
    gubbins_taxa_division <- gubbins_snps[gubbins_snps$Taxa == row[5],]
  }else{
    gubbins_taxa_division <- gubbins_snps[gubbins_snps$end_node_index == as.numeric(row[6]),]
  }
  if(nrow(gubbins_taxa_division) > 0){
    ## Now check for position of the SNP, allow for a SNP diff of +- 5bp
    gubbins_taxa_pos_division <- gubbins_taxa_division %>%
      filter((base_number >= (as.numeric(row[3]) - 5)) & (base_number <= (as.numeric(row[3]) + 5)))
    if(nrow(gubbins_taxa_pos_division) > 0){
      if(nrow(gubbins_taxa_pos_division) > 1){
        diff_vals <- abs(as.numeric(row[3]) - gubbins_taxa_pos_division$base_number)
        closest <- which.min(diff_vals)
        gubbins_recon <- gubbins_taxa_pos_division[closest,]
        sim_gubbins <- "Yes"
        snp_diff <- as.numeric(row[3]) - gubbins_recon$base_number
        gubbins_sim <-  "Yes"
        gubbins_index <- gubbins_recon$index
      }else{
        gubbins_recon <- gubbins_taxa_pos_division
        sim_gubbins <- "Yes"
        snp_diff <- as.numeric(row[3]) - gubbins_recon$base_number
        gubbins_sim <- "Yes" 
        gubbins_index <- gubbins_recon$index
      }
    }
    
  }
  toc()
  return(paste(sim_gubbins, as.character(snp_diff), gubbins_sim, as.character(gubbins_index), sep = "~"))
  
}



simul_looper <- function(gubbins_snps, simul_snps, taxa_step = "taxa",
                         snps = NULL, rec_rate, branch_rate){
  ## Function to loop through the simul snps testing if the taxa identified are the 
  ## same as the simul, if so see if the gubbins matches a snp here. 
#  browser()
  if(is.null(snps)){
    gubbins_snps <- gubbins_snps %>% 
      mutate(simul = "No") %>%
      mutate(index = row_number()) %>%
      mutate(snp_diff = NaN) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      mutate(gubbins = "No") %>%
      mutate(snp_diff = NaN) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    
  }else{
    gubbins_snps <- gubbins_snps %>% 
      filter(Type == snps) %>%
      mutate(simul = "No") %>%
      mutate(index = row_number()) %>%
      mutate(snp_diff = NaN) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      filter(Type == snps) %>%
      mutate(gubbins = "No") %>%
      mutate(snp_diff = NaN) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
  }
  cat("\n")
  num_zeros <- nchar(nrow(simul_snps))
  num_rows <- nrow(simul_snps)
  
  for(k in 1:nrow(simul_snps)){
    nchar_k <- nchar(k)
    tic(msg = "Took this long for the loop")
    nchar_0 <- num_zeros - nchar_k
    cat("\r", "Completed ", rep(0, nchar_0), k, " of ", num_rows, " SNPs", sep = "")
    current_simul <- simul_snps[k,]
    ## Check any gubbins rows with similar taxa
    if(taxa_step == "taxa"){
      gubbins_taxa_division <- gubbins_snps[gubbins_snps$Taxa == current_simul$Taxa,]
    }else{
      gubbins_taxa_division <- gubbins_snps[gubbins_snps$end_node_index == current_simul$end_node_index,]
    }
    if(nrow(gubbins_taxa_division) > 0){
      ## Now check for position of the SNP, allow for a SNP diff of +- 5bp
      gubbins_taxa_pos_division <- gubbins_taxa_division %>%
        filter((base_number >= (current_simul$Start - 5)) & (base_number <= (current_simul$Start + 5)))
      if(nrow(gubbins_taxa_pos_division) > 0){
        if(nrow(gubbins_taxa_pos_division) > 1){
          diff_vals <- abs(current_simul$Start - gubbins_taxa_pos_division$base_number)
          closest <- which.min(diff_vals)
          gubbins_recon <- gubbins_taxa_pos_division[closest,]
          simul_snps[k,"gubbins"] <- "Yes"
          simul_snps[k,"snp_diff"] <- current_simul$Start - gubbins_recon$base_number
          gubbins_snps[gubbins_snps$index == gubbins_recon$index, "simul"] <- "Yes"
          gubbins_snps[gubbins_snps$index == gubbins_recon$index, "snp_diff"] <- current_simul$Start - gubbins_recon$base_number
        }else{
          gubbins_recon <- gubbins_taxa_pos_division
          simul_snps[k,"gubbins"] <- "Yes"
          simul_snps[k,"snp_diff"] <- current_simul$Start - gubbins_recon$base_number
          gubbins_snps[gubbins_snps$index == gubbins_recon$index, "simul"] <- "Yes"
          gubbins_snps[gubbins_snps$index == gubbins_recon$index, "snp_diff"] <- current_simul$Start - gubbins_recon$base_number
        }
      }
      
    }
    toc()
  }
  
  
  
  ## Calculate PPV as true pos by gubbins over true pos by gubbins and false pos by gubbins
  true_pos <- nrow(gubbins_snps %>% filter(simul == "Yes"))
  false_pos <- nrow(gubbins_snps %>% filter(simul == "No"))
  gubbins_ppv <- round((true_pos / (true_pos + false_pos)) * 100,digits = 3)
  
  ## Calcultate Sensitivity from simul dataset, gubbins annotated over total snps
  simul_pos <- nrow(simul_snps %>% filter(gubbins == "Yes"))
  gubbins_sen <- round((simul_pos / nrow(simul_snps)) * 100, digits = 3)
  cat("\n")  
  cat(paste("Gubbins PPV = ", as.character(gubbins_ppv), "%", sep = ""), "\n")
  cat(paste("Gubbins sensitivity = ", as.character(gubbins_sen), "%", sep = ""), "\n")
  dist_plot <- ggplot(data = simul_snps %>%
                        mutate(x_col = "Simulated")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = gubbins)) +
    labs(x = "", y = "Number of descendant tips from Node", colour = "Reconstructed by Gubbins")
  ppv_plot <- ggplot(data = gubbins_snps %>%
                       mutate(x_col = "Gubbins")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = simul)) +
    labs(x = "", y = "Number of descendant tips from Node", colour = "Within simulated data")
  
  ## Create summary df
  summary_dataset <- as.data.frame(matrix(nrow = 1, ncol = 5))
  colnames(summary_dataset) <- c("ppv","sensitivity","rec_rate","branch_rate","rec-branch")
  summary_dataset[1,1] <- gubbins_ppv
  summary_dataset[1,2] <- gubbins_sen
  summary_dataset[1,3] <- rec_rate
  summary_dataset[1,4] <- branch_rate
  summary_dataset[1,5] <- paste(as.character(rec_rate), as.character(branch_rate), sep = "-")
  
  
  
  return(list(simul_df = simul_snps, gubbins_snps = gubbins_snps,
              sens_plot = dist_plot, ppv_plot = ppv_plot,
              summary = summary_dataset))
  
}

simul_looper_apply <- function(gubbins_snps, simul_snps, taxa_step = "taxa",
                         snps = NULL, rec_rate, branch_rate){
  ## Function to loop through the simul snps testing if the taxa identified are the 
  ## same as the simul, if so see if the gubbins matches a snp here. 

  if(is.null(snps)){
    gubbins_snps <- gubbins_snps %>% 
      mutate(simul = "No") %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    
  }else{
    gubbins_snps <- gubbins_snps %>% 
      filter(Type == snps) %>%
      mutate(index = row_number()) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      filter(Type == snps) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
  }
  cat("\n")
  tic(msg = "Time for apply runs")
  tot_res <- apply(X = simul_snps, MARGIN = 1, FUN = simul_apply_func, 
                   gubbins_snps = gubbins_snps,
                   taxa_step = taxa_step)
  toc()
  matty_res <- data.frame(matrix(unlist(strsplit(tot_res,split = "~")), ncol = 4, byrow = TRUE)) %>%
    rename(gubbins = 1, snp_diff = 2, simul = 3, index = 4) %>%
    mutate(snp_diff = as.numeric(snp_diff),
           index = as.numeric(index))
  simul_snps <- simul_snps %>%
    mutate(gubbins = matty_res$gubbins[row_number()],
           snp_diff = matty_res$snp_diff[row_number()])
  gubbins_snps <- gubbins_snps %>% 
    left_join(matty_res %>% select(snp_diff, simul, index), by = c("index" = "index")) %>%
    mutate(simul = ifelse(is.na(simul), "No", simul))
  
  
  ## Calculate PPV as true pos by gubbins over true pos by gubbins and false pos by gubbins
  true_pos <- nrow(gubbins_snps %>% filter(simul == "Yes"))
  false_pos <- nrow(gubbins_snps %>% filter(simul == "No"))
  gubbins_ppv <- round((true_pos / (true_pos + false_pos)) * 100,digits = 3)
  
  ## Calcultate Sensitivity from simul dataset, gubbins annotated over total snps
  simul_pos <- nrow(simul_snps %>% filter(gubbins == "Yes"))
  gubbins_sen <- round((simul_pos / nrow(simul_snps)) * 100, digits = 3)
  cat("\n")  
  cat(paste("Gubbins PPV = ", as.character(gubbins_ppv), "%", sep = ""), "\n")
  cat(paste("Gubbins sensitivity = ", as.character(gubbins_sen), "%", sep = ""), "\n")
  dist_plot <- ggplot(data = simul_snps %>%
                        mutate(x_col = "Simulated")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = gubbins)) +
    labs(x = "", y = "Number of descendant tips from Node", colour = "Reconstructed by Gubbins")
  ppv_plot <- ggplot(data = gubbins_snps %>%
                       mutate(x_col = "Gubbins")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = simul)) +
    labs(x = "", y = "Number of descendant tips from Node", colour = "Within simulated data")
  
  ## Create summary df
  summary_dataset <- as.data.frame(matrix(nrow = 1, ncol = 5))
  colnames(summary_dataset) <- c("ppv","sensitivity","rec_rate","branch_rate","rec-branch")
  summary_dataset[1,1] <- gubbins_ppv
  summary_dataset[1,2] <- gubbins_sen
  summary_dataset[1,3] <- rec_rate
  summary_dataset[1,4] <- branch_rate
  summary_dataset[1,5] <- paste(as.character(rec_rate), as.character(branch_rate), sep = "-")
  
  
  
  return(list(simul_df = simul_snps, gubbins_snps = gubbins_snps,
              sens_plot = dist_plot, ppv_plot = ppv_plot,
              summary = summary_dataset))
  
}


###############################################################################
## Load up the Gubbins recombinations data ####################################
###############################################################################

gubbins_reccy_gff <- delim_reader("fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.recombination_predictions.gff")
gubbins_reccy_csv <- recombination_gff_cleaner(gubbins_reccy_gff)
gubbins_branch_base <- read.csv("fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.embl.csv.csv",
                                stringsAsFactors = FALSE)
#gubbins_tree <- read.tree(file = gubbins_tree_file)


## Lets run this on the fasttree-iqtree-joint rec 0.1 branch 0.1 dataset initially

setwd("~/Dropbox/phd/gubbins_testing/gubbins_ppv_data/")

snp_data_jar <- data_cleaner(gubbins_gff = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.recombination_predictions.gff",
                         gubbins_branch_base_csv = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.embl.csv.csv",
                         gubbins_tree_file = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.node_labelled.final_tree.tre",
                         simul_summary_file = "sim-branch-0.1-rec-0.1.summary",
                         simul_tree_file = "sim-branch-0.1-rec-0.1.tree")
snp_data_mar <- data_cleaner(gubbins_gff = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.recombination_predictions.gff",
                             gubbins_branch_base_csv = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.embl.csv.csv",
                             gubbins_tree_file = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.node_labelled.final_tree.tre",
                             simul_summary_file = "sim-branch-0.1-rec-0.1.summary",
                             simul_tree_file = "sim-branch-0.1-rec-0.1.tree")

snp2_data <- data_cleaner(gubbins_gff = "./sim-0.1-snp2/ft-iq-jar-sim-branch-0.1-rec-0.1.recombination_predictions.gff",
                          gubbins_branch_base_csv = "./sim-0.1-snp2/ft-iq-jar-sim-branch-0.1-rec-0.1.embl_branch.csv.csv",
                          gubbins_tree_file = "./sim-0.1-snp2/ft-iq-jar-sim-branch-0.1-rec-0.1.node_labelled.final_tree.tre",
                          simul_summary_file = "sim-branch-0.1-rec-0.1.summary",
                          simul_tree_file = "sim-branch-0.1-rec-0.1.tree")

snp2_data_0.5 <- data_cleaner(gubbins_gff = "./sim-0.5-snp2/ft-iq-jar-sim-branch-0.5-rec-0.5.recombination_predictions.gff",
                              gubbins_branch_base_csv = "./sim-0.5-snp2/ft-iq-jar-sim-branch-0.5-rec-0.5.embl_branch.csv.csv",
                              gubbins_tree_file = "./sim-0.5-snp2/ft-iq-jar-sim-branch-0.5-rec-0.5.node_labelled.final_tree.tre",
                              simul_summary_file = "sim-rec-0.5_ft_iq_data/sim-branch-0.5-rec-0.5.summary",
                              simul_tree_file = "sim-rec-0.5_ft_iq_data/sim-branch-0.5-rec-0.5.tree")
snp_data_0.5 <- data_cleaner(gubbins_gff = "./sim-rec-0.5_ft_iq_data/fasttree-iqtree-joint.recombination_predictions.gff",
                             gubbins_branch_base_csv = "./sim-rec-0.5_ft_iq_data/fasttree-iqtree-joint.embl_branch_base.csv.csv",
                             gubbins_tree_file = "./sim-rec-0.5_ft_iq_data/fasttree-iqtree-joint.node_labelled.final_tree.tre",
                             simul_summary_file = "./sim-rec-0.5_ft_iq_data/sim-branch-0.5-rec-0.5.summary",
                             simul_tree_file = "./sim-rec-0.5_ft_iq_data/sim-branch-0.5-rec-0.5.tree")

plotTreeDiff(simul_tree, gubbins_tree, treesFacing = TRUE)

pdf(file = "./snp_0.5_tree_dists.pdf",width = 20, height = 15, paper = "a4r")
plotTreeDiff(snp_data_0.5$gubbins_tree,snp_data_0.5$simul_tree, type="cladogram", use.edge.length=FALSE, 
             treesFacing = TRUE, edge.width=2, colourMethod = "pallette",
             palette = funky, baseCol = "midnightblue", cex = 0.5)

dev.off()
## Lets loop through the simul_summary, see if there is a SNP reconstructed by GUBBINS
## with the same taxa in the tree 



system.time(test_out <- simul_looper(snp_data_jar$gubbins_snps, snp_data_jar$simul_snps, taxa_step = "index",
                         branch_rate = 0.1, rec_rate = 0.1, snps = "r"))
system.time(test_out_apply <- simul_looper_apply(snp_data_jar$gubbins_snps, snp_data_jar$simul_snps, taxa_step = "index",
                         branch_rate = 0.1, rec_rate = 0.1, snps = "r"))

test_out$sens_plot
test_out$ppv_plot

test_out_snp2 <- simul_looper(snp2_data$gubbins_snps, snp2_data$simul_snps, taxa_step = "index",
                              branch_rate = 0.1, rec_rate = 0.1)

test_out_0.5 <- simul_looper(snp_data_0.5$gubbins_snps, snp_data_0.5$simul_snps, taxa_step = "index",
                             branch_rate = 0.5, rec_rate = 0.5)
test_out_0.5$sens_plot
test_out_0.5$ppv_plot

test_out_0.5_snp2 <- simul_looper(snp2_data_0.5$gubbins_snps, snp2_data_0.5$simul_snps, taxa_step = "index",
                             branch_rate = 0.5, rec_rate = 0.5)
test_out_0.5_snp2$sens_plot
test_out_0.5_snp2$ppv_plot


## Lets plot out the success by taxa depth 
gubbins_res_depth <- gubbins_res %>%
  mutate(reconned = ifelse(simul == "Yes",1,0)) %>%
  mutate(row_val = 1) %>% group_by(taxa_depth) %>%
  summarise(num_reconstructed = sum(reconned), 
            num_total = sum(row_val)) %>%
  mutate(success = (num_reconstructed/ num_total) * 100) %>%
  as.data.frame()

simul_res_depth <- simul_res %>%
  mutate(reconned = ifelse(gubbins == "Yes",1,0)) %>%
  mutate(row_val = 1) %>% group_by(taxa_depth) %>%
  summarise(num_reconstructed = sum(reconned), 
            num_total = sum(row_val)) %>%
  mutate(success = (num_reconstructed/ num_total) * 100) %>%
  as.data.frame()

ggplot(data = simul_res_depth) + 
  geom_col(aes(x = taxa_depth, y = success))

example_depth <- simul_res %>% filter(taxa_depth > 5) %>% head()

## Lets look into an example of the taxa two isolates 





