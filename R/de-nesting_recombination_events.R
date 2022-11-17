###############################################################################
## Testing for overlap in recombination events from the summarry files ########
###############################################################################

require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(data.table, quietly = TRUE, warn.conflicts = FALSE)
require(stringr, quietly = TRUE, warn.conflicts = FALSE)
require(argparse, quietly = TRUE, warn.conflicts = FALSE)
require(snow, quietly = TRUE, warn.conflicts = FALSE)


###############################################################################
## Functions ##################################################################
###############################################################################

taxa_0_adder <- function(row){
  separate_taxa <- unlist(strsplit(row[5], ","))
  for(taxa in 1:length(separate_taxa)){
    if(nchar(separate_taxa[taxa]) == 1){
      separate_taxa[taxa] <- paste(0, separate_taxa[taxa], sep = "")
    }
  }
  taxa_paste <- paste(separate_taxa, collapse = ",")
  return(taxa_paste)
}

typing_lapply_func <- function(row, rec_val){
  cat(sprintf("\rOn this row number %s ", row[7]))
  
  overlaps <- ""
  start_overs <- rec_val[which(data.table::between(as.integer(row[3]), rec_val$Start, rec_val$End)) ,] %>%
    filter(Event. != as.integer(row[1])) %>%
    filter(ntaxa <= as.integer(row[9]))
  end_overs <- rec_val[which(data.table::between(as.integer(row[ 4]), rec_val$Start, rec_val$End)) ,] %>%
    filter(Event. != as.integer(row[1])) %>%
    filter(ntaxa <= as.integer(row[9]))
  tot_overs <- bind_rows(start_overs, end_overs) %>% 
    distinct() 
  if(nrow(tot_overs) > 0){
    
    overlaps <- paste(tot_overs$Event., collapse = "-")
  }
  
  return(overlaps)  
  
}

mutate_taxa_func <- function(taxa, ref_taxa){
  taxa_split <- gsub(",","|", ref_taxa)
  grepl(taxa_split, taxa)
}

masking_detection_func <- function(row, rec_val){
  cat(sprintf("\rOn this row number %s ", row[7]))
  
  if(row[8] == ""){
    return("")
  }else{
    to_strip <- ""
    events_to_investigate <- as.integer(unlist(strsplit(row[8],"-")))
    events_string <- rec_val %>% filter(Event. %in% events_to_investigate) %>%
      mutate(taxa_overlap = ifelse(mutate_taxa_func(Taxa,row[5]), "Yes","No")) %>%
      filter(taxa_overlap == "Yes") %>%
      pull(Event.)
    if(length(events_string) > 0){
      
      to_strip <- paste(events_string, collapse = "-")
    }else{
      to_strip <- ""
    }
    
    return(to_strip)
  }
  
  
  
}

remove_overlaps <- function(rec_val, tot_df){
  ## Function to meld the tot_df removing overwritten recombination events, snps and clonal SNPS
  ## Go through the rec_val line by line, if line has worries deal with, for S just remove 
  ## from the tot_df and for recs, remove rs by nesting, going from the deepest nodes to the 
  ## ones closest to the tips 
  cat(sprintf("\n"))
  tot_df$nests <- 0
  for(k in 1:nrow(rec_val)){
    cat(sprintf("\rWorking on row: %s    of %s", k, nrow(rec_val)))
    if(rec_val$worries[k] == ""){
      next
    }else if(rec_val$Type[k] == "S"){
      ## Just removing the overwritten SNPs
      tot_df <- tot_df[-which(tot_df$Event. == rec_val$Event.[k]),]
    }else{
      
      ## So now we should just be in the recombinations that have some 
      ## overwritten segments 
      indiv_overlaps <- unlist(strsplit(rec_val[k, "worries"], "-"))
      overlaps <- rec_val[rec_val$Event. %in% indiv_overlaps,] %>%
        arrange(desc(ntaxa)) %>% as.data.frame()
      current_row <- rec_val[k,]
      nested <- 0
      slimmed <- "No"
      for(over in 1:nrow(overlaps)){
        current_over <- overlaps[over,]
        
        ## Need what sort of overlap there is 
        if(between(current_row$Start, current_over$Start, current_over$End) & between(current_row$End, current_over$Start, current_over$End)){
          ## remove all the r snps and rec event associated with the event 
          tot_df <- tot_df[-which(tot_df$Event. == rec_val$Event.[k]),]
          break
        }else{
          ## Check for overlap with start value 
          if(current_row$Start >= current_over$Start & current_row$Start <= current_over$End & current_row$End >= current_over$End){
            ## Here we've got the overlap across the start end:
            ## row             \-----------------\
            ## start     \-----------\
            ## end                   \-----------\
            ## So we'll remove all snps from the over end to the current start and 
            ## reshape the start to be the over end 
            remove_rows <- which((tot_df$Type == "r") & (tot_df$Event. == current_row$Event.) & (tot_df$Start <= current_over$End))
            if(length(remove_rows) > 0){
              tot_df <- tot_df[-remove_rows,]
            }
            current_row$Start <- current_over$End
            slimmed <- "Yes"
          }else if(current_row$End >= current_over$Start & current_row$End <= current_over$End & current_row$Start <= current_over$Start){
            ## Here we've got the overlap across the End end:
            ## row             \-----------------\
            ## start                      \---------------\
            ## end             \----------\
            ## So we'll remove all snps from the over end to the current start and 
            ## reshape the start to be the over end 
            remove_rows <- which((tot_df$Type == "r") & (tot_df$Event. == current_row$Event.) & (tot_df$Start >= current_over$Start))
            if(length(remove_rows) > 0){
              tot_df <- tot_df[-remove_rows,]
            }
            current_row$End <- current_over$Start
            slimmed <- "Yes"
          }else{
            ## Now we must have a new recombination nested within the larger one:
            ## row          \-------------------------------\
            ## start               \-----------------\
            ## end           \-----                   ------\
            ## For this then I'll keep the one event, the row as the same, but just remove the SNPs that have been overlapped
            ## I'll make a note of the nested within ones 
            nested <- nested + 1
            remove_rows <- which((tot_df$Type == "r") & (tot_df$Event. == current_row$Event.) & (tot_df$Start >= current_over$Start) & 
                                   (tot_df$Start <= current_over$End))
            if(length(remove_rows) > 0){
              tot_df <- tot_df[-remove_rows,]
            }
            slimmed <- "Yes"
            
          }
          
        }
        
      }
      ## Now I'll alter the main rec event if it still exists in the tot df to update with the slimmed down ranges
      if(slimmed == "Yes"){
        tot_df <- tot_df %>%
          mutate(Start = ifelse(Type == "R" & Event. == current_row$Event., current_row$Start, Start),
                 End = ifelse(Type == "R" & Event. == current_row$Event., current_row$End, End),
                 nests = ifelse(Type == "R" & Event. == current_row$Event., nested, nests))
      }
    }
    
    
  }
  cat(sprintf("\n"))  
  return(tot_df)
}


complete_run_through <- function(summary_df, out_loc){
  start_time_over <- Sys.time()
  reccs_only <- summary_df %>% filter(Type != "r") %>%
    mutate(nchar_taxa = nchar(Taxa),
           index = row_number(),
           overlaps = "",
           ntaxa = str_count(Taxa, pattern = ",") + 1)
  
  ## Add on the 0s for grep matching 
  retyped_taxa <- apply(reccs_only, 1, taxa_0_adder)
  reccs_only$Taxa <- retyped_taxa
  ## Get the general overlaps to sort through 
  cat(sprintf("Getting the general overlaps: \n"))
  start_tim <- Sys.time()
  overlap_col <- apply(reccs_only, 1, typing_lapply_func, rec_val = reccs_only)
  reccs_only$overlaps <- overlap_col
  cat(sprintf("\n Took this long: %s \n", Sys.time() - start_tim))
  ## Get the overlaps we want to worry about 
  cat(sprintf("Filtering the worrysome overlaps: \n"))
  start_tim <- Sys.time()
  worries <- apply(reccs_only, 1, masking_detection_func,
                   rec_val = reccs_only)
  reccs_only$worries <- worries
  cat(sprintf("\n Took this long: %s \n", Sys.time() - start_tim))
  ## Remove the worrysome overlaps 
  cat(sprintf("Removing the overlaps: \n"))
  start_time <- Sys.time()
  
  tot_df <- remove_overlaps(reccs_only, summary_df)
  
  count_compare <- dplyr::count(tot_df, Type) %>%
    rename(n_removed = n) %>% 
    left_join(dplyr::count(summary_df, Type)) %>%
    rename(n_original = n)
  print(count_compare)
  cat(sprintf("Took this long: %s \n", Sys.time() - start_time))
  cat(sprintf("Took this long overall: %s \n", Sys.time() - start_time_over))
  
  write.table(tot_df, out_loc,col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = "\t")
  
}

nest_snow_func <- function(summary_rows, list_of_summaries, out_dir){
  
  run_through_rows <- list_of_summaries[summary_rows]
  for(k in 1:length(run_through_rows)){
    current_summary <- run_through_rows[k]
    current_output <- paste(out_dir, basename(current_summary), sep = "")
    summary_df <- read.table(current_summary, sep = "\t", comment.char = "", header = TRUE) 
    ## Updated the summary output with pareto runs so need to remove these additional
    ## columns as the rest is based on indexes!
    if("orig_start" %in% colnames(summary_df)){
      summary_df <- summary_df %>%
        select(-c(orig_start, orig_end))
    }
    complete_run_through(summary_df, current_output)    
  }
}

main_nest_func <- function(list_of_summaries, threads, out_dir){
  newEnv <- new.env()
  # Assigning variables
  newEnv$list_of_summaries <- list_of_summaries 
  newEnv$out_dir <- out_dir
  summary_rows <- seq(1, length(list_of_summaries))
  # Set up the SNOW cluster
  cluster_function <- makeCluster(spec = threads, type = "SOCK")
  # Get the row split for the cluster 
  parra_input <- clusterSplit(cluster_function, summary_rows)
  # Get the data copied over 
  cat(sprintf("Copying over the data \n"))
  clusterExport(cluster_function, "list_of_summaries", envir = newEnv)
  clusterExport(cluster_function, "out_dir", envir = newEnv)
  # Get the function copied over 
  cat(sprintf("Copying over the functions \n"))
  clusterExport(cluster_function, "complete_run_through")
  clusterExport(cluster_function, "remove_overlaps")
  clusterExport(cluster_function, "masking_detection_func")
  clusterExport(cluster_function, "mutate_taxa_func")
  clusterExport(cluster_function, "typing_lapply_func")
  clusterExport(cluster_function, "taxa_0_adder")
  clusterExport(cluster_function, "nest_snow_func")
  clusterEvalQ(cluster_function, library(dplyr))
  clusterEvalQ(cluster_function, library(stringr))
  # Now lets run the results 
  cat(sprintf("Running on the parrallel cores: "))
  job_start <- Sys.time()
  snp_data_parra <- clusterApply(cluster_function, parra_input, 
                                 fun = nest_snow_func,
                                 list_of_summaries = list_of_summaries,
                                 out_dir = out_dir)
  job_end <- Sys.time()
  cat(sprintf("Done in: %s \n", (job_end - job_start)))
  stopCluster(cluster_function)
  
}

main_func <- function(){
  ## Get the input 
  
  summary_input <- get_input()
  sim_summaries <- list.files(summary_input$summary_dir,
                              full.names = TRUE, pattern = "*\\.summary")
  print(sim_summaries)
  out_dir <- summary_input$out
  last_char <- substr(out_dir, nchar(out_dir) - 1 + 1, nchar(out_dir))
  if(last_char != "/")
    out_dir <- paste(out_dir, "/", sep = "")
  threads <- summary_input$threads
  main_nest_func(sim_summaries, threads = threads, out_dir = out_dir)
  
  
}

get_input <- function(){
  parser <- ArgumentParser(description='Unnest the summary files ')
  parser$add_argument('--summary-dir', type="character", required = TRUE,
                      help='Directory of the summary files for the true data',
                      dest='summary_dir')
  parser$add_argument('--threads', type="integer", default = 1,
                      help='Number of threads to use for calculations',
                      dest='threads')
  parser$add_argument('--out', dest = 'out', type = "character", required = TRUE,
                      help = "Folder to put the out summaries in ")
  
  return(parser$parse_args())
}


################################################################################
## Run through #################################################################
################################################################################

main_func()

# summary_list <- list.files("~/Dropbox/phd/gubbins_testing/theses_analyses/sim_summaries_jc/",
#                            pattern = "\\.summary", full.names = TRUE)
# out_dir <- as.data.frame("~/Dropbox/phd/gubbins_testing/theses_analyses/un-nested-summaries/") %>%
#   rename(dir = 1)
# main_nest_func(summary_list, out_dir, 6)
# 
# 
# example_summary <-  read.table("~/Dropbox/phd/gubbins_testing/theses_analyses/sim_summaries_jc/orig-orig-orig-rep-1-sim-branch-0.1-rec-5.0-JC-JC.summary", sep = "\t",
#                                comment.char = "", header = TRUE)
# example_summary_rec_0.1 <-  read.table("~/Dropbox/phd/gubbins_testing/theses_analyses/sim_summaries_jc/orig-orig-orig-rep-1-sim-branch-0.1-rec-0.1-JC-JC.summary", sep = "\t",
#                                comment.char = "", header = TRUE)
# example_summary_rec_2 <-  read.table("~/Dropbox/phd/gubbins_testing/theses_analyses/sim_summaries_jc/orig-orig-orig-rep-1-sim-branch-0.1-rec-2.0-JC-JC.summary", sep = "\t",
#                                        comment.char = "", header = TRUE)
# 
# write.table(tot_df_2$tot_df,
#             "~/Dropbox/phd/gubbins_testing/theses_analyses/test_summo_out.summary",
#             col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
# tot_df <- complete_run_through(example_summary)
# tot_df_0.1 <- complete_run_through(example_summary_rec_0.1)
# tot_df_1 <- complete_run_through(example_summary_rec_1)
# tot_df_2 <- complete_run_through(example_summary_rec_1)
# reccs_only <- example_summary %>% filter(Type != "r") %>%
#   mutate(nchar_taxa = nchar(Taxa),
#          index = row_number(),
#          overlaps = "",
#          ntaxa = str_count(Taxa, pattern = ",") + 1)
# 
# 
# retyped_taxa <- apply(reccs_only, 1, taxa_0_adder)
# 
# reccs_only$Taxa <- retyped_taxa
# 
# ## Write this out to a list function!
# 
# overlap_col <- apply(reccs_only, 1, typing_lapply_func, rec_val = reccs_only)
# reccs_only$overlaps <- overlap_col
# 
# ## Get the overlaps we want to worry about 
# worries <- apply(reccs_only, 1, masking_detection_func,
#                  rec_val = reccs_only)
# reccs_only$worries <- worries
# 
# tot_df <- remove_overlaps(reccs_only, example_summary)
# dplyr::count(tot_df, nests)
# tot_df[tot_df$nests == 4,]
# tot_df %>% filter(Event. == 48 & Type == "r") %>% nrow()
# example_summary %>% filter(Event. == 48 & Type == "r") %>% nrow()
# head(tot_df)## Now to 
# 
# ## lets go through it as a list 
# 
# reccs_list <- split(reccs_only,seq(nrow(reccs_only)))   
# taxa_split <- function(list_item){
#   taxa_split <- unlist(strsplit(list_item$Taxa, ","))
#   list_item <- as.list(list_item)
#   list_item$taxa_indiv <- taxa_split
#   return(list_item)
# }
# reccs_list_taxa <- lapply(reccs_list, taxa_split)
# 


