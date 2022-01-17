###############################################################################
## Exploring the point predictive value in gubbins ############################
###############################################################################

require(dplyr)
require(stringr)
require(ggplot2)
require(ape)
require(ggpubr)
require(data.table)

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
  tot_paths <- ape::nodepath(tree)
  tips <- NULL
  for(k in 1:length(tot_paths)){
    if(node_index %in% tot_paths[[k]]){
      tips <- append(tips, paths[[k]][length(paths[[k]])])
    }
  }
  return(tips)
}


append_ancestors <- function(gub_tree, branch_base){
  ## Function to find all the ancestors of an end Node and append these to a branch_base to match the 
  ## output from the simulator 
  
  branch_base$Taxa <- ""
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
  }
    return(branch_base)
}

###############################################################################
## Load up the Gubbins recombinations data ####################################
###############################################################################

## Lets run this on the fasttree-iqtree-joint rec 0.1 branch 0.1 dataset initially

setwd("~/Dropbox/phd/gubbins_testing/gubbins_ppv_data/")

gubbins_reccy_gff <- delim_reader("fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.recombination_predictions.gff")
gubbins_reccy_csv <- recombination_gff_cleaner(gubbins_reccy_gff)

gubbins_branch_base <- read.csv("fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.embl.csv.csv",
                                stringsAsFactors = FALSE)
gubbins_tree <- read.tree(file = "fasttree-iqtree-joint-sim-branch-0.1-rec-0.1.node_labelled.final_tree.tre")

branch_base_typed <- typing_gubbins_rec(gubbins_reccy_csv, gubbins_branch_base)

branch_base_typed_taxa <- append_ancestors(gubbins_tree, branch_base_typed)


###############################################################################
## Now we have to read in the simul data and link it to start and end nodes ###
###############################################################################

simul_summary <- read.table("sim-branch-0.1-rec-0.1.summary", sep = "\t",
                            comment.char = "", header = TRUE)

simul_clonal_snps <- simul_summary %>% filter(Type == "S")
head(simul_clonal_snps)





