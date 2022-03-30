###############################################################################
## Creating the classified SNPs object from command line ######################
###############################################################################
require(argparse, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE,warn.conflicts = FALSE)
require(stringr, quietly = TRUE,warn.conflicts = FALSE)
require(data.table, quietly = TRUE, warn.conflicts = FALSE)
require(snow, quietly = TRUE, warn.conflicts = FALSE)
###############################################################################
## FUNCTIONS ##################################################################
###############################################################################

get_input <- function(){
  parser <- ArgumentParser(description='Create classified SNP object')
  parser$add_argument('--reccy-gff', type="character", required = TRUE,
                      help='recombination predictions in gff format ')
  parser$add_argument('--branch-base', dest='branch_base', required = TRUE,
                      help='embl branch_base file in csv format')
  parser$add_argument('--threads', type="integer", default = 1,
                      help='Number of threads to use for the SNP classification',
                      dest='threads')
  parser$add_argument('--out', dest = 'out', type = "character", required = TRUE,
                      help = "Location for out csv")
  
  return(parser$parse_args())
}


recombination_gff_cleaner <- function(in_gff_pmen3){
  
  ## Takes delim file read in and cleans it up 
  pmen3_reccy_data <- tryCatch({
    data.frame(data = (matrix(data = NA, nrow = nrow(in_gff_pmen3),
                                                ncol = 9)))
  }, error=function(cond){
    return(NULL)
  })
  if(is.null(pmen3_reccy_data)){
    return(NULL)
  }else{
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
}

delim_reader <- function(delim_file){
  #browser()
  ## read in the gff recombination events file from gubbins
  reccy_csv <- NULL
  try({
    reccy_csv <- read.delim(delim_file,header = FALSE, comment.char = "#")
    colnames(reccy_csv) <- c("type","prog","class","start","end","trent","alexander","arnold","attribute")
  }, silent = TRUE)
  return(reccy_csv)
}

typing_lapply_func <- function(row, reccy_preds, num_zeros, num_rows){
  nchar_k <- nchar(row[length(row)])
  nchar_0 <- num_zeros - nchar_k
  cat("\r", "Completed ", rep(0, nchar_0), row[length(row)], " of ", num_rows, " SNPs", sep = "")
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
  cat("\n")
  num_zeros <- nchar(nrow(branch_base))
  num_rows <- nrow(branch_base)
  branch_base <- branch_base %>% mutate(index = row_number())
  snp_types <- apply(X = branch_base, MARGIN = 1, FUN = typing_lapply_func, 
                     reccy_preds = reccy_preds, num_zeros = num_zeros,
                     num_rows = num_rows)
  branch_base$Type <- snp_types
  cat("\n")
  return(branch_base)
}

classifier_parrallely <- function(unique_rows, mutty_recon, clade_csv){
  ## This function acts as way into the cohort_allocator function
  ## getting the hesids from the cluster input and subsetting the hes-data on 
  ## the rows 
  
  input_mutty_rows <- mutty_recon[unique_rows,]
  typed_snps <- typing_gubbins_rec_apply(clade_csv, input_mutty_rows)
  
  return(typed_snps)
}


parrallel_typing <- function(clade_csv, mutty_recon, threads){
  #browser
  
  newEnv <- new.env()
  
  # Assigning variables
  newEnv$clade_csv <- clade_csv
  newEnv$mutty_recon <- mutty_recon
  newEnv$threads <- threads
  
  ## Set up the data for parrallelisation
  mutty_rows <- seq(1, nrow(mutty_recon))
  # Set up the SNOW cluster
  cluster_function <- makeCluster(spec = threads, type = "SOCK")
  # Get the row split for the cluster 
  parra_input <- clusterSplit(cluster_function, mutty_rows)
  # Get the data copied over 
  cat(sprintf("Copying over the data \n"))
  clusterExport(cluster_function, "clade_csv", envir = newEnv)
  clusterExport(cluster_function, "mutty_recon", envir = newEnv)
  # Get the function copied over 
  cat(sprintf("Copying over the functions \n"))
  clusterExport(cluster_function, "typing_gubbins_rec_apply")
  clusterExport(cluster_function, "typing_lapply_func")
  clusterExport(cluster_function, "classifier_parrallely")
  clusterEvalQ(cluster_function, library(dplyr))
  # Now lets run the results 
  cat(sprintf("Running on the parrallel cores: "))
  job_start <- Sys.time()
  snp_data_parra <- clusterApply(cluster_function, parra_input, 
                                 fun = classifier_parrallely,
                                 mutty_recon = mutty_recon,
                                 clade_csv = clade_csv)
  job_end <- Sys.time()
  cat(sprintf("Done in: %s \n", (job_end - job_start)))
  stopCluster(cluster_function)
  classified_snps <- base::as.data.frame(dplyr::bind_rows(snp_data_parra)) %>%
    mutate(index = row_number())
  
  return(classified_snps)
}

main_func <- function(){
  
  input_args <- get_input()

  clade_gff <- delim_reader(input_args$reccy_gff)
  if(is.null(clade_gff)){
    ## We'll still recreate the snps, we'll just set everything to be S!
    mutty_recon <- read.csv(input_args$branch_base,
                            stringsAsFactors = FALSE)
    classified_snps <- mutty_recon %>% mutate(index = row_number()) %>%
      mutate(Type = "S")
    
    
  }else{
    clade_csv <- recombination_gff_cleaner(clade_gff)
  
    mutty_recon <- read.csv(input_args$branch_base,
                            stringsAsFactors = FALSE)
    threads <- input_args$threads
    if(threads > 1){
      classified_snps <- parrallel_typing(clade_csv, mutty_recon, threads)
    }else{
      system.time(classified_snps <- typing_gubbins_rec_apply(clade_csv, mutty_recon))
    }
  }
    write.csv(classified_snps,
              file = input_args$out,
              row.names = FALSE,
              quote = FALSE)  
  
}

###############################################################################
## Run the script #############################################################
###############################################################################

main_func()

cat("\n", "Done")
