# selection_etl_script
#
# ETL script to load the selection data warehouse
#
# PLAN:
#
#    1. Populate dimension lookup tables: dimPopData, dimExperiment, dimPos, dimGene
#            - dimPopData and dimExperiment are straight inserts
#            - dimPos is initially built from dimGene.csv, 
#                 will require loading into staging and then munging to create dimPos and dimGene

library(data.table)
library(dplyr)  # for the %>% function :)
setwd("~/Git_repos/SelectionDW_ETL/")
options(scipen = 999)

# Transformation Functions
#
# Each function defines a data transformation specific to the statistic results. 
# Must return a data.table (pop, chr, start, end, nzise, variable, value)
# where: 
#    variable: the statistic e.g. tajimad, Eta_E, FuLi_D etc.
#    value: the value that goes with the statistic, as the name states.
tajima_ <- function () {
  
  function (results) {
    colnames(results) <- c("chrom", "BIN_START", "NumSites_TajimasD", "TajimaD", "pop", "chr")
    tmp_ <- results[, .(pop, chr, BIN_START, NumSites_TajimasD, TajimaD)]
    tmp_ <- melt.data.table(tmp_, id.vars = names(tmp_)[-5])
    
    tmp_[, chrom_start := BIN_START]
    tmp_[, nsize := diff(results[2:3, BIN_START])]
    tmp_[, chrom_end := BIN_START + nsize - 1]
    
    return (tmp_)
  }
}
fanwu_ <- function () {
  
  function (results) {
    tmp_ <- results[, c(4,5,7,9:17), with = FALSE]
    colnames(tmp_) <- c("chrom_start", "chrom_end", "nsize", "S", "Eta", "Eta_E", "Pi", "FuLi_D", "FuLi_F", "FayWu_H", "pop", "chr")
    
    tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "chr", "chrom_start", "chrom_end", "nsize"))
    
    return (tmp_)
  }
}

kaks_ <- function () {
  function (results) {
    colnames(results) <- c("GeneID", "GeneName", "ka", "ks", "kcomputed", "pop", "chr")
    tmp_ <- melt.data.table(results, id.vars = c("GeneID", "GeneName"))
    
    return (tmp_)
  }
  # code here to map genes to positions.
}
af_ <- function () {
  function (results) {
    colnames(results) <- c("POS", "Ref", "Alt", "Anc", "MAF", "DAF", "pop", "chr")
    
    # Discussed with Murray, dropping allele information at the moment, and will load this
    # separately in the future.
    
    tmp_ <- results[, .(pop, chr, chrom_start = POS, chrom_end = POS + 1, nsize = 1, MAF, DAF)]
    tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "chr", "chrom_start", "chrom_end", "nsize"))
    
    return (tmp_)
    
  }
}
nsl_ <- function () {
  function (results) {
    colnames(results) <- c("locus_id", "chrom_start", "nsl_freq_1", "sL1","sL0", "unstd_nsl", "norm_nsl", "significant_nsl", "pop", "chr")
    
    tmp_ <- results[, .(pop, chr, chrom_start, chrom_end = chrom_start + 1, nsize = 1, nsl_freq_1, sL1, sL0, unstd_nsl, norm_nsl, significant_nsl)]
    tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "chr", "chrom_start", "chrom_end", "nsize"))
    
    return (tmp_)
  }
}
ihs_ <- function () {
  function (results) {
    colnames(results) <- c("locus_id", "chrom_start", "ihs_freq_1", "ihh1","ihh0", "unstd_ihs", "norm_ihs", "significant_ihs", "pop", "chr")
    
    tmp_ <- results[, .(pop, chr, chrom_start, chrom_end = chrom_start + 1, nsize = 1, ihs_freq_1, ihh1, ihh0, unstd_ihs, norm_ihs, significant_ihs)]
    tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "chr", "chrom_start", "chrom_end", "nsize"))
    
    return (tmp_)
  }
}

xpehh_ <- function () {
  # xpehh is cross-populational
  function (results) {
    colnames(results) <- c('pos','chrom_start','gpos', "popA_freq_1", "ihhA",'popB_freq_1',"ihhB", "unstd_xpehh", "norm_xpehh", "significant_xpehh","pop","pop2", 'chr')
    
    tmp_ <- results[, .(pop, pop2, chr, chrom_start, chrom_end = chrom_start + 1, nsize = 1, popA_freq_1, ihhA, popB_freq_1, ihhB, unstd_xpehh, norm_xpehh, significant_xpehh)]
    tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "pop2", "chr", "chrom_start", "chrom_end", "nsize"))
    
    return (tmp_)
    
  }
}

fst_ <- function () {
  # fst is cross-populational
  function(results){
    colnames(results) <- c('chrom','chrom_start','chrom_end','nvariants_FST','weighted_FST','mean_FST','pop', 'pop2','chr')
    tmp_ <- results[, .(pop, pop2, chr, chrom_start, chrom_end, nsize = chrom_end-chrom_start, nvariants_FST, weighted_FST, mean_FST )]
    tmp_[weighted_FST < 0] <- 0
    tmp_[mean_FST < 0] <- 0
    tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "pop2", "chr", "chrom_start", "chrom_end", "nsize"))
  
    return (tmp_)
    
  }
}

utility_switch <- function (file_) {
  
  util_ <- if (grepl(".taj_d$", file_)) list(id = "tajima", util_ = tajima_(), type='intra')
  else if (grepl(".faw$", file_)) list(id = "fanwu", util_ = fanwu_(), type='intra')
  else if (grepl(".fst$", file_)) list(id = "fst", util_ = fst_(), type='inter')
  #else if (grepl(".kaks$", file_)) list(id = "kaks", util_ = kaks_(), type='inter')
  else if (grepl(".af$", file_)) list(id = "af", util_ = af_(), type = 'intra')
  else if (grepl(".nsl.", file_)) list(id = "nsl", util_ = nsl_(), type = 'intra')
  else if (grepl(".ihs.", file_)) list(id = "ihs", util_ = ihs_(), type = 'intra')
  else if (grepl(".xpehh.", file_)) list(id = "xpehh", util_ = xpehh_(), type = 'inter')
  else NULL
  
  return (util_)
}

#### helper functions

db_execute <- function (query, user = 'murraycadzow', db = "selectiondw_test") {
  execute <- sprintf('psql -U %s -d %s -c "%s"',
                     user, db, query)
  
  return (system(execute))
}

#bulk_load <- function (filename, db_table) {
#    query <- sprintf("
#                     LOAD DATA LOCAL INFILE '%s' INTO TABLE %s
#                     WITH DELIMITER AS ','
#                     LINES TERMINATED BY '\n'
#                     IGNORE 1 LINES;
#                     ", filename, db_table)
#    return (query)
#}


bulk_load <- function (filename, db_table) {
  query <- sprintf("
                   COPY %s FROM '%s'
                   CSV
                   NULL AS 'NA'
                   HEADER
                   ", db_table, filename)
  return (query)
}

init_ <- function () {
  
  # populate dimPopData
  db_execute(bulk_load(paste0(getwd(),"/data/populations.csv"), "staging_pop"))
  db_execute("select populate_pop();")
  
  # populate dimExperiment
  db_execute("select populate_experiment();")
  
  # bulk load gene data
  db_execute(bulk_load(paste0(getwd(),"/data/dimGene_psql.csv"), "staging_gene"))
  
  # init dimPos
  db_execute("select populate_pos();")
  
  # init dimGene
  db_execute("select populate_gene();")
  db_execute("truncate staging_gene;")
  
  # init dimStat
  db_execute("select populate_stat();")
}

# Extract data from raw results files
# We (Murray and I) have agreed on the following:
#
#    1. the ETL pipeline will operate on a single, user-defined directory 
#       which contains all of the results files for a SINGLE statistic.
#    2. the statistic type will be identified by the file extensions
#    3. The user must supply the experiment id as an arguement to main()

# parse_
# extracts population and chromosome from a results file
parse_ <- function (x) {
  prefix_ <- strsplit(x, "\\.")[[1]][1]
  chr_ <- as.integer(substring(prefix_, 
                               first = 4, 
                               last = ifelse((nchar(prefix_) - 3) > 1, 5, 4)))  # some mumbo jumbo to handle inconsistent naming conventions
  pop_ <- substring(prefix_, first = 1, last = 3)
  
  return (data.table(chr = chr_, pop = pop_))
}

# parseInter_
# extracts populations and chromosome from a inter population file
# format: pop1_pop2_chr.suffix
# pop is 3 characters
parseInter_ <- function(x){
  prefix_ <- strsplit(x, "\\.")[[1]][1]
  chr_ <- as.integer(strsplit(prefix_, '_')[[1]][3])
  pop_ <- strsplit(prefix_, '_')[[1]][1]
  pop2_ <- strsplit(prefix_, '_')[[1]][2]
  return (data.table(chr = chr_, pop = pop_, pop2 = pop2_))
}

# extract_
# Given a list of files, extract the raw data and combine
# NOTE: assumes that all files are the output of a single statistic
extract_ <- function (files = NULL, utility = NULL) {
  
  header_ <- function () {
    nrows <- if (utility[["id"]] == "fanwu") 5
    else if (utility[["id"]] %in% c("ihs", "nsl")) 0
    else 1
  }
  
  raw <- rbindlist(
    lapply(files, function (f_) {
      if(utility[['type']] == 'intra'){
        meta <- parse_(f_)
      } else {
        meta <- parseInter_(f_)
      }
      
      
      tmp <- data.table(read.table(f_, skip = header_(), header = FALSE))
      tmp[, pop := meta[["pop"]]]
      
      if(utility[['type']] == 'inter'){
        tmp[, pop2 := meta[["pop2"]]]
      }
      
      tmp[, chr := meta[["chr"]]]
      
      return (tmp)
    })
  )
  
  return (raw)
}
# end result intra: (chr, start, end, pop, stat, nsize, slide, type)
# end result inter: ("pop", "pop2", "chr", "chrom_start", "chrom_end", "nsize", "variable", "value")
transform_ <- function (results, utility = NULL) {
  
  print("        Wrangling results files ...")
  if(utility[['type']] == 'intra' ){
    schema <- c("pop", "chr", "chrom_start", "chrom_end", "nsize", "variable", "value")
  } else {
    schema <- c("pop", "pop2", "chr", "chrom_start", "chrom_end", "nsize", "variable", "value")
  }
  tmp <- utility[["util_"]](results)
  
  return (tmp[, schema, with = FALSE])
  
}





load_ <- function (results, experiment_id = 1) {
  
  print("        Loading data warehouse...")
  # write results file into tmp file
  result_file <- "/tmp/temporary_results.csv"
  write.csv(results, result_file, row.names = FALSE, quote = FALSE)
  
  # bulk load result file
  db_execute(bulk_load(result_file, "staging_results"))
  
  # check if there are new variants
  db_execute("select check_variants();")
  
  # unstage the data from the staging table into intraSel fact table
  db_execute(sprintf("select unstage(%s)", experiment_id))
  
  # truncate the staging table
  db_execute("truncate staging_results;")
  
}


loadInter_ <- function (results, experiment_id = 1) {
  
  print("        Loading data warehouse...")
  # write results file into tmp file
  result_file <- "/tmp/temporary_results.csv"
  write.csv(results, result_file, row.names = FALSE, quote = FALSE)
  
  # bulk load result file
  db_execute(bulk_load(result_file, "staging_inter_results"))
  
  # check if there are new variants
  db_execute("select check_inter_variants();")
  
  # unstage the data from the staging table into intraSel fact table
  db_execute(sprintf("select inter_unstage(%s)", experiment_id))
  
  # truncate the staging table
  db_execute("truncate staging_inter_results;")
  
}



pipeline <- function (experiment_id, files, utility) {
  if(utility[['type']] == 'intra'){
    extract_(files, utility) %>% 
      transform_(utility = utility) %>% 
      load_(experiment_id = experiment_id)
  } else {
    extract_(files, utility) %>% 
      transform_(utility = utility) %>% 
      loadInter_(experiment_id = experiment_id)
  }
}



main <- function (experiment_id = 1, results_directory = NULL) {
  
  # Find results files
  setwd(results_directory)
  files_ <- list.files(".", full.names = FALSE, recursive = FALSE, pattern = "\\.")
  
  # Link current files with a statistic and appropriate utility function
  util_ <- utility_switch(files_[1])
  
  if (is.null(util_)) {
    print("Unable to recognise type of results files.")
    stop()
  }
  
  # skip fst, xpehh and kaks at this stage.
  if (util_[["id"]] %in% c("kaks")) {
    print("Not able to load kaks at this stage.")
    return (-1)
  }
  
  # Run pipeline per population
  # Note: cannot run over ALL populations at once, as this produces
  #       millions of rows of data.
  if(util_[['type']] == 'intra'){
    populations <- rbindlist(lapply(files_, parse_))
  } else {
    populations <- rbindlist(lapply(files_, parseInter_))
  }
  
  for (pop in unique(populations[["pop"]])) {
    
    print(sprintf("----    %s    ----", pop))
    
    popfiles <- files_[grepl(paste0("^",pop), files_)]
      if(util_[['type']] == 'inter'){
        for(pop2 in unique(populations[['pop2']])){
          if(pop != pop2){
            print(sprintf("----  subpop  %s    ----", pop2))
            pop2files <- popfiles[grepl(pop2, popfiles)]
		if(length(pop2files) > 0){ pipeline(experiment_id, pop2files, util_)}
          }
        }
      }else{
        pipeline(experiment_id, popfiles, util_)
      }
  }
  return (0)
}

### For testing only
test_pipeline <- function (directory = "~/Git_repos/SelectionDW_ETL/") {
  
  setwd(directory)
  
  print("Initialising the data warehouse")
  init_()
  
  main(experiment_id = 1, "./data/SelectionData/FAWH")
}

####    test_pipeline()

