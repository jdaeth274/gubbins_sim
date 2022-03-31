###############################################################################
## Thesis runs through on the PPV and gubbins sensitivity for the repped data #
###############################################################################

require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(stringr, quietly = TRUE, warn.conflicts = FALSE)
require(ape, quietly = TRUE, warn.conflicts = FALSE)
require(data.table, quietly = TRUE, warn.conflicts = FALSE)
require(devtools, quietly = TRUE, warn.conflicts = FALSE)
# if(!("treespace" %in% rownames(installed.packages())))
#   devtools::install_github("thibautjombart/treespace")
# require(treespace, quietly = TRUE, warn.conflicts = FALSE)
require(snow, quietly = TRUE, warn.conflicts = FALSE)
require(argparse, quietly = TRUE, warn.conflicts = FALSE)

###############################################################################
## Functions ##################################################################
###############################################################################

get_input <- function(){
  parser <- ArgumentParser(description='Calculate the PPV and sensitivity of gubbins runs')
  parser$add_argument('--embl-dir', type="character", required = TRUE,
                      help='Directory to classified snps csvs ',
                      dest = "embl_dir")
  parser$add_argument('--tree-dir', dest='tree_dir', required = TRUE,
                      help='Directory to trees, orig and gubbins should be here')
  parser$add_argument('--summary-dir', type="character", required = TRUE,
                      help='Directory of the summary files for the true data',
                      dest='summary_dir')
  parser$add_argument('--threads', type="integer", default = 1,
                      help='Number of threads to use for calculations',
                      dest='threads')
  parser$add_argument('--snps', type = "character", default = NULL,
                      help = "Whether to type snps by recombination type, default is null, put any word to reconstruct separately (should change this to boolean)",
                      dest = "snps")
  parser$add_argument('--out', dest = 'out', type = "character", required = TRUE,
                      help = "Location for out csv")
  
  return(parser$parse_args())
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
  tot_nodes <- c(gub_tree$tip.label, gub_tree$node.label)
  for(k in 1:nrow(branch_base)){
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
  tot_nodes <- c(simul_tree$tip.label,as.character(seq(Ntip(simul_tree) + 1, length.out = Nnode(simul_tree))))
  for(k in 1:nrow(simul_data)){
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

data_cleaner <- function(gubbins_snps, gubbins_tree_file,
                         simul_summary_file, simul_tree_file){
  ## Function to load up the neccessary data to check for the Gubbins
  ## SNPs predictions 
  #browser()
  gubbins_tree <- read.tree(file = gubbins_tree_file)
  
  branch_base_typed_taxa <- append_ancestors(gubbins_tree, gubbins_snps)
  
  ###############################################################################
  ## Now we have to read in the simul data and link it to start and end nodes ###
  ###############################################################################
  
  simul_summary <- read.table(simul_summary_file, sep = "\t",
                              comment.char = "", header = TRUE)
  
  simul_tree <- read.tree(simul_tree_file)
  
  simul_data_nodes <- simul_node_finder(simul_tree = simul_tree,
                                        simul_data = simul_summary)
  # treespace install very tricky on the cluster  
  # tree_distances <- treeDist(gubbins_tree, simul_tree, lambda = 0)
  # 
  # if(tree_distances != 0){
  #   print("Warning the Tree topologies are not identical, SNP predictions will not 
  #         work properly")
  # }
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
  
  return(paste(sim_gubbins, as.character(snp_diff), gubbins_sim, as.character(gubbins_index), sep = "~"))
  
}


simul_looper_dplyr <- function(gubbins_snps, simul_snps, taxa_step = "taxa",
                         snps = NULL, rec_rate, branch_rate, data_rep,
                         first_mod, main_mod, recon){
  ## Function to loop through the simul snps testing if the taxa identified are the 
  ## same as the simul, if so see if the gubbins matches a snp here. 

  if(is.null(snps)){
    gubbins_snps <- gubbins_snps %>% 
      mutate(simul = "No") %>%
      mutate(index = row_number()) %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      mutate(snp_index = row_number()) %>%
      mutate(gubbins = "No") %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    
  }else{
    
    gubbins_snps <- gubbins_snps %>% 
      filter(Type == snps) %>%
      mutate(simul = "No") %>%
      mutate(index = row_number()) %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1)) 
    simul_snps <- simul_snps %>%
      mutate(snp_index = row_number()) %>%
      filter(Type == snps) %>%
      mutate(gubbins = "No") %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1)) %>%
      mutate(Start = Start + 1)
  }
  #browser()
  if(taxa_step == "taxa"){
    sim_typed <- simul_snps %>%
      left_join(gubbins_snps %>% select(c(Taxa, base_number, simul, index)), by = c("Taxa" = "Taxa", "Start" = "base_number")) %>%
      mutate(gubbins = ifelse(is.na(simul), "No", "Yes"))
    
    gubb_typed <- gubbins_snps %>%
      mutate(gubbins = ifelse(index %in% sim_typed$index, "Yes","No"))
  }else{
    #gubbins_taxa_division <- gubbins_snps[gubbins_snps$end_node_index == current_simul$end_node_index,]
    test_full_join <- simul_snps %>%
      # Crossing
      full_join(gubbins_snps, by = c("end_node_index")) %>%
      # Filter
      filter(base_number >= Start - 5,
             base_number <= Start + 5) %>%
      rename(Type = 2, Taxa = 5) %>%
      mutate(testy = "yes") %>%
      mutate(match_difference = Start - base_number) %>% arrange(index, match_difference) %>%
      filter(duplicated(index) == FALSE) %>% 
      arrange(snp_index, match_difference) %>% filter(duplicated(snp_index) == FALSE)
    sim_typed <- simul_snps %>%
      mutate(gubbins = ifelse(snp_index %in% test_full_join$snp_index, "Yes","No"))
    gubb_typed <- gubbins_snps %>%
      mutate(simul = ifelse(index %in% test_full_join$index, "Yes","No"))
    
    
  }
  
  test_one <- "Yes"
  ## Calculate PPV as true pos by gubbins over true pos by gubbins and false pos by gubbins
  true_pos <- nrow(gubb_typed %>% filter(simul == "Yes"))
  false_pos <- nrow(gubb_typed %>% filter(simul == "No"))
  gubbins_ppv <- round((true_pos / (true_pos + false_pos)) * 100, digits = 3)
  
  ## Calcultate Sensitivity from simul dataset, gubbins annotated over total snps
  simul_pos <- nrow(sim_typed %>% filter(gubbins == "Yes"))
  gubbins_sen <- round((simul_pos / nrow(sim_typed)) * 100, digits = 3)
  # cat("\n")  
  # cat(paste("Gubbins PPV = ", as.character(gubbins_ppv), "%", sep = ""), "\n")
  # cat(paste("Gubbins sensitivity = ", as.character(gubbins_sen), "%", sep = ""), "\n")
  # dist_plot <- ggplot(data = simul_snps %>%
  #                       mutate(x_col = "Simulated")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = gubbins)) +
  #   labs(x = "", y = "Number of descendant tips from Node", colour = "Reconstructed by Gubbins")
  # ppv_plot <- ggplot(data = gubbins_snps %>%
  #                      mutate(x_col = "Gubbins")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = simul)) +
  #   labs(x = "", y = "Number of descendant tips from Node", colour = "Within simulated data")
  # 
  ## Create summary df
  #browser()
  summary_dataset <- as.data.frame(matrix(nrow = 1, ncol = 9))
  colnames(summary_dataset) <- c("ppv","sensitivity","rec_rate","branch_rate","data_rep","first_mod","main_mod", "recon", "snp")
  summary_dataset[1,1] <- gubbins_ppv
  summary_dataset[1,2] <- gubbins_sen
  summary_dataset[1,3] <- rec_rate
  summary_dataset[1,4] <- branch_rate
  summary_dataset[1,5] <- data_rep
  summary_dataset[1,6] <- first_mod
  summary_dataset[1,7] <- main_mod
  summary_dataset[1,8] <- recon
  if(is.null(snps)){
    summary_dataset[1,9] <- "all"
  }else{
    summary_dataset[1,9] <- snps
  }
  return(list(simul_df = simul_snps, gubbins_snps = gubbins_snps,
              summary = summary_dataset))
  
}


simul_looper <- function(gubbins_snps, simul_snps, taxa_step = "taxa",
                         snps = NULL, rec_rate, branch_rate, data_rep,
                         first_mod, main_mod, recon){
  ## Function to loop through the simul snps testing if the taxa identified are the 
  ## same as the simul, if so see if the gubbins matches a snp here. 
  #browser()
  if(is.null(snps)){
    gubbins_snps <- gubbins_snps %>% 
      mutate(simul = "No") %>%
      mutate(index = row_number()) %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      mutate(gubbins = "No") %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    
  }else{
    
    gubbins_snps <- gubbins_snps %>% 
      filter(Type == snps) %>%
      mutate(simul = "No") %>%
      mutate(index = row_number()) %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
    simul_snps <- simul_snps %>%
      filter(Type == snps) %>%
      mutate(gubbins = "No") %>%
      mutate(snp_diff = 0) %>%
      mutate(taxa_depth = (str_count(Taxa, pattern = ",") + 1))
  }
  #browser()
  simul_rows <- nrow(simul_snps)
  simmer <- 1
  while(simmer <= simul_rows){
    current_simul <- simul_snps[simmer,]
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
           simul_snps[simmer,"gubbins"] <- "Yes"
           simul_snps[simmer,"snp_diff"] <- current_simul$Start - gubbins_recon$base_number
           gubbins_snps <- gubbins_snps %>% 
             mutate(simul = ifelse(index == gubbins_recon$index, "Yes", simul),
                    snp_diff = ifelse(index == gubbins_recon$index,(current_simul$Start - gubbins_recon$base_number), snp_diff))
         }else{
           gubbins_recon <- gubbins_taxa_pos_division
           simul_snps[simmer,"gubbins"] <- "Yes"
           simul_snps[simmer,"snp_diff"] <- current_simul$Start - gubbins_recon$base_number
           gubbins_snps <- gubbins_snps %>% 
            mutate(simul = ifelse(index == gubbins_recon$index, "Yes", simul),
                   snp_diff = ifelse(index == gubbins_recon$index,(current_simul$Start - gubbins_recon$base_number), snp_diff))
           # gubbins_snps[gubbins_snps$index == gubbins_recon$index, "simul"] <- "Yes"
           # gubbins_snps[gubbins_snps$index == gubbins_recon$index, "snp_diff"] <- (current_simul$Start - gubbins_recon$base_number)
         }
       }
    }
    simmer <- simmer + 1
  }
  
  test_one <- "Yes"
  ## Calculate PPV as true pos by gubbins over true pos by gubbins and false pos by gubbins
  true_pos <- nrow(gubbins_snps %>% filter(simul == "Yes"))
  false_pos <- nrow(gubbins_snps %>% filter(simul == "No"))
  gubbins_ppv <- round((true_pos / (true_pos + false_pos)) * 100, digits = 3)
  
  ## Calcultate Sensitivity from simul dataset, gubbins annotated over total snps
  simul_pos <- nrow(simul_snps %>% filter(gubbins == "Yes"))
  gubbins_sen <- round((simul_pos / nrow(simul_snps)) * 100, digits = 3)
  # cat("\n")  
  # cat(paste("Gubbins PPV = ", as.character(gubbins_ppv), "%", sep = ""), "\n")
  # cat(paste("Gubbins sensitivity = ", as.character(gubbins_sen), "%", sep = ""), "\n")
  # dist_plot <- ggplot(data = simul_snps %>%
  #                       mutate(x_col = "Simulated")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = gubbins)) +
  #   labs(x = "", y = "Number of descendant tips from Node", colour = "Reconstructed by Gubbins")
  # ppv_plot <- ggplot(data = gubbins_snps %>%
  #                      mutate(x_col = "Gubbins")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = simul)) +
  #   labs(x = "", y = "Number of descendant tips from Node", colour = "Within simulated data")
  # 
  ## Create summary df
  #browser()
  summary_dataset <- as.data.frame(matrix(nrow = 1, ncol = 9))
  colnames(summary_dataset) <- c("ppv","sensitivity","rec_rate","branch_rate","data_rep","first_mod","main_mod", "recon", "snp")
  summary_dataset[1,1] <- gubbins_ppv
  summary_dataset[1,2] <- gubbins_sen
  summary_dataset[1,3] <- rec_rate
  summary_dataset[1,4] <- branch_rate
  summary_dataset[1,5] <- data_rep
  summary_dataset[1,6] <- first_mod
  summary_dataset[1,7] <- main_mod
  summary_dataset[1,8] <- recon
  summary_dataset[1,9] <- snps
  
  return(list(simul_df = simul_snps, gubbins_snps = gubbins_snps,
              summary = summary_dataset))
  
}

simul_looper_apply <- function(gubbins_snps, simul_snps, taxa_step = "taxa",
                               snps = NULL, rec_rate, branch_rate, data_rep,
                               first_mod, main_mod, recon){
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
  #tic(msg = "Time for apply runs")
  tot_res <- apply(X = simul_snps, MARGIN = 1, FUN = simul_apply_func, 
                   gubbins_snps = gubbins_snps,
                   taxa_step = taxa_step)
  #toc()
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
  # cat("\n")  
  # cat(paste("Gubbins PPV = ", as.character(gubbins_ppv), "%", sep = ""), "\n")
  # cat(paste("Gubbins sensitivity = ", as.character(gubbins_sen), "%", sep = ""), "\n")
  # dist_plot <- ggplot(data = simul_snps %>%
  #                       mutate(x_col = "Simulated")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = gubbins)) +
  #   labs(x = "", y = "Number of descendant tips from Node", colour = "Reconstructed by Gubbins")
  # ppv_plot <- ggplot(data = gubbins_snps %>%
  #                      mutate(x_col = "Gubbins")) + geom_jitter(aes(x = x_col, y = taxa_depth, colour = simul)) +
  #   labs(x = "", y = "Number of descendant tips from Node", colour = "Within simulated data")
  # 
  ## Create summary df
  summary_dataset <- as.data.frame(matrix(nrow = 1, ncol = 9))
  colnames(summary_dataset) <- c("ppv","sensitivity","rec_rate","branch_rate","data_rep","first_mod","main_mod", "recon", "snp")
  summary_dataset[1,1] <- gubbins_ppv
  summary_dataset[1,2] <- gubbins_sen
  summary_dataset[1,3] <- rec_rate
  summary_dataset[1,4] <- branch_rate
  summary_dataset[1,5] <- data_rep
  summary_dataset[1,6] <- first_mod
  summary_dataset[1,7] <- main_mod
  summary_dataset[1,8] <- recon
  summary_dataset[1,9] <- snps
  
  
  
  
  return(list(simul_df = simul_snps, gubbins_snps = gubbins_snps,
              summary = summary_dataset))
  
}

file_df_creator <- function(file_list, file_suffix){
  ## Function to swap the list of file names over to a df for easier indexing
  #browser()
  bassy_names <- basename(file_list)
  models_only <- sub(file_suffix, "JC",bassy_names, perl = TRUE)
  model_df <- as.data.frame(matrix(unlist(str_split(models_only, "-")), byrow = TRUE, ncol = 12)) %>%
    rename(start_mod = 1, main_mod = 2, recon = 3, fluff1 = 4,
           data_rep = 5, fluff2 = 6, fluff3 = 7, branch = 8, 
           fluff4 = 9, rec = 10, sim_mod = 11, gubb_mod = 12) %>%
    mutate(tot_file = file_list[row_number()]) %>%
    select(-c(fluff1, fluff2, fluff3, fluff4))
  
  return(model_df)
  
  
  
}

sim_ppv_getter <- function(sim_df, gubbins_df, tree_df, snp_type = NULL){
  #browser()
  out_df <- NULL
  for(k in 1:(nrow(sim_df) - 1)){
    current_sim <- sim_df[k,]
    current_summary <- current_sim$tot_file
    ## Now we have to loop through the tree models for each run
    gubbins_mods <- gubbins_df %>%
      filter(data_rep == current_sim$data_rep, 
             rec == current_sim$rec,
             branch == current_sim$branch)
    start_time <- Sys.time()
    for(gubb in 1:nrow(gubbins_mods)){
      gubb_start <- Sys.time()
      gubbins_row <- gubbins_mods[gubb,]
      gubb_tree <- tree_df %>% filter(data_rep == current_sim$data_rep, 
                                      rec == current_sim$rec,
                                      branch == current_sim$branch,
                                      tot_mod == gubbins_row$tot_mod) %>%
        dplyr::pull(tot_file)
      sim_tree <- tree_df %>% filter(data_rep == current_sim$data_rep, 
                                     rec == current_sim$rec,
                                     branch == current_sim$branch,
                                     main_mod == "orig") %>%
        dplyr::pull(tot_file)
      gubbins_snps <- read.csv(gubbins_row$tot_file, stringsAsFactors = FALSE)
      
      cleaned_dat <- data_cleaner(gubbins_snps = gubbins_snps,
                                  gubbins_tree_file = gubb_tree,
                                  simul_summary_file = current_summary,
                                  simul_tree_file = sim_tree)
      ## Now we've cleaned the data we can run the ppv checker 
    # gubbins_S_data <- simul_looper(gubbins_snps = cleaned_dat$gubbins_snps,
    #               simul_snps = cleaned_dat$simul_snps,
    #               taxa_step = "index",
    #               snps = "S",
    #               branch_rate = current_sim$branch,
    #               rec_rate = current_sim$rec,
    #               main_mod = gubbins_row$main_mod,
    #               first_mod = gubbins_row$start_mod,
    #               recon = gubbins_row$recon,
    #               data_rep = gubbins_row$data_rep)
      testy <- "pops"
      if(is.null(snp_type)){
        gubbins_S_data <- simul_looper_dplyr(gubbins_snps = cleaned_dat$gubbins_snps,
                                             simul_snps = cleaned_dat$simul_snps,
                                             taxa_step = "index",
                                             snps = NULL,
                                             branch_rate = current_sim$branch,
                                             rec_rate = current_sim$rec,
                                             main_mod = gubbins_row$main_mod,
                                             first_mod = gubbins_row$start_mod,
                                             recon = gubbins_row$recon,
                                             data_rep = gubbins_row$data_rep)
        sim_data <- gubbins_S_data$summary
      }else{
        gubbins_S_data <- simul_looper_dplyr(gubbins_snps = cleaned_dat$gubbins_snps,
                                       simul_snps = cleaned_dat$simul_snps,
                                       taxa_step = "index",
                                       snps = "S",
                                       branch_rate = current_sim$branch,
                                       rec_rate = current_sim$rec,
                                       main_mod = gubbins_row$main_mod,
                                       first_mod = gubbins_row$start_mod,
                                       recon = gubbins_row$recon,
                                       data_rep = gubbins_row$data_rep)
        gubbins_r_data <- simul_looper_dplyr(gubbins_snps = cleaned_dat$gubbins_snps,
                                       simul_snps = cleaned_dat$simul_snps,
                                       taxa_step = "index",
                                       snps = "r",
                                       branch_rate = current_sim$branch,
                                       rec_rate = current_sim$rec,
                                       main_mod = gubbins_row$main_mod,
                                       first_mod = gubbins_row$start_mod,
                                       recon = gubbins_row$recon,
                                       data_rep = gubbins_row$data_rep)
        sim_data <- bind_rows(gubbins_S_data$summary, gubbins_r_data$summary)
      }
      out_df <- bind_rows(out_df, sim_data)
      end_diff <- Sys.time() - gubb_start
      cat(sprintf("Finished on model %s in sim %s in %s time \n", gubbins_row$tot_mod, 
                  sub("orig-orig-orig-","",basename(current_sim$tot_file)), end_diff))
      
    }
    
    cat(sprintf("Finished on sim %s in %s \n", basename(current_sim$tot_file), (Sys.time() - start_time)))
  }
  return(out_df)
}

sim_apply_function <- function(row_index, sim_df, gubbins_df, tree_df, snp_type){
  ## Function to apply the ppv calculator in parallel
#  browser()
  sim_rows <- sim_df[row_index,]
  ppv_df <- sim_ppv_getter(sim_rows, gubbins_df, tree_df, snp_type)
  
  return(ppv_df)
  
}

main_ppv_function <- function(sim_summaries, gubbins_snps, tree_files, threads = 1,
                              snp_type = NULL){
  ## Loop through the csvs for the summaries, then get them into the right format,
  ## find the corresponding snp file and then run through them both. 
  #browser()
  sim_df <- file_df_creator(sim_summaries, "JC\\..*$")
  gubbins_df <- file_df_creator(gubbins_snps,"JC\\..*$") %>% 
    mutate(tot_mod = paste(start_mod, main_mod, recon, sep = "-"))
  tree_df <- file_df_creator(tree_files, "JC\\..*$") %>% 
    mutate(tot_mod = paste(start_mod, main_mod, recon, sep = "-"))
  
  ## Now we've got the file df lets loop through the sim_df and get the appropriate 
  ## gubbins df to run through 
  ## Set up clustersim_df
  # browser()
  # test_dataset <- sim_apply_function(c(5,6), sim_df = sim_df,
  #                  gubbins_df = gubbins_df, 
  #                    tree_df = tree_df,
  #                  snp_type = snp_type)
  newEnv <- new.env()
  
  # Assigning variables
  newEnv$sim_df <- sim_df 
  newEnv$gubbins_df <- gubbins_df
  newEnv$tree_df <- tree_df
  newEnv$snp_type <- snp_type
  #newEnv$output_file <- output_file
  
  summary_rows <- seq(1, nrow(sim_df))
  # Set up the SNOW cluster
  cluster_function <- makeCluster(spec = threads, type = "SOCK")
  # Get the row split for the cluster 
  parra_input <- clusterSplit(cluster_function, summary_rows)
  # Get the data copied over 
  cat(sprintf("Copying over the data \n"))
  clusterExport(cluster_function, "sim_summaries", envir = newEnv)
  clusterExport(cluster_function, "gubbins_df", envir = newEnv)
  clusterExport(cluster_function, "tree_df", envir = newEnv)
  clusterExport(cluster_function, "snp_type", envir = newEnv)
  #clusterExport(cluster_function, "output_file", envir = newEnv)
  # Get the function copied over 
  cat(sprintf("Copying over the functions \n"))
  clusterExport(cluster_function, "sim_apply_function")
  clusterExport(cluster_function, "sim_ppv_getter")
  clusterExport(cluster_function, "simul_looper")
  clusterExport(cluster_function, "data_cleaner")
  clusterExport(cluster_function, "node_dictionary")
  clusterExport(cluster_function, "simul_node_finder")
  clusterExport(cluster_function, "taxa_ancestors")
  clusterExport(cluster_function, "append_ancestors")
  clusterExport(cluster_function, "tip_nodes")
  clusterExport(cluster_function, "simul_looper_dplyr")
  clusterEvalQ(cluster_function, library(dplyr))
  clusterEvalQ(cluster_function, library(ape))
  clusterEvalQ(cluster_function, library(stringr))
  # Now lets run the results 
  cat(sprintf("Running on the parrallel cores: "))
  job_start <- Sys.time()
  snp_data_parra <- clusterApply(cluster_function, parra_input, 
                                 fun = sim_apply_function,
                                 sim_df = sim_df,
                                 gubbins_df = gubbins_df, 
                                 tree_df = tree_df,
                                 snp_type = snp_type)
  job_end <- Sys.time()
  cat(sprintf("Done in: %s \n", (job_end - job_start)))
  stopCluster(cluster_function)
  ppv_data <- base::as.data.frame(dplyr::bind_rows(snp_data_parra)) %>%
    mutate(index = row_number())
  
  return(ppv_data)
  
  }

main_func <- function(){
  start_time <- Sys.time()
  input_data <- get_input()
  print(input_data$summary_dir)
  sim_summaries <- list.files(input_data$summary_dir,
                              full.names = TRUE, pattern = "*\\.summary")
  
  gubb_snps <- list.files(input_data$embl_dir,
                          full.names = TRUE, pattern = "*\\.csv")
  
  tree_df <- list.files(input_data$tree_dir,
                        full.names = TRUE, pattern = "*\\.tre")
  snp_type <- input_data$snps
  threads <- input_data$threads
  test_full <- main_ppv_function(sim_summaries, gubb_snps, tree_df,
                                 threads = threads, snp_type = snp_type)
  write.csv(test_full, input_data$out,
            row.names = FALSE, quote = FALSE)
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  cat(sprintf("Finished complete run in: %s", tot_time))
  
}

###############################################################################
## Main run ###################################################################
###############################################################################

main_func()
# sim_summaries <- list.files("~/Dropbox/phd/gubbins_testing/theses_analyses/sim_summaries_jc/",
#                             full.names = TRUE, pattern = "*\\.summary")
# 
# gubb_snps <- list.files("~/Dropbox/phd/gubbins_testing/theses_analyses/jc_rep_10_snps/",
#                             full.names = TRUE, pattern = "*\\.csv")
# 
# tree_df <- list.files("~/Dropbox/phd/gubbins_testing/theses_analyses/jc_rep_10_trees/",
#                         full.names = TRUE, pattern = "*\\.tre")
# test_full <- main_ppv_function(sim_summaries, gubb_snps, tree_df, threads = 8, snp_type = NULL)
# 
# out_data <- main_ppv_function(sim_summaries, gubb_snps, tree_df, threads = 8)
# 
# head(out_data)
# write.csv(out_data, "~/Dropbox/phd/gubbins_testing/theses_analyses/ppv_out_data_dplyr.csv",
#           row.names = FALSE, quote = FALSE)
# 
# ## Now for da plot plot 
# 
# 
# tot_data <- constant_branch_all$total_data
# setDT(tot_data)[, r_end := r + 5]
# tot_data


