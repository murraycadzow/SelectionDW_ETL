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
library(plyr)  # for the %>% function :)
setwd("/mnt/DataDrive/Scratch/SelectionDW/")


# Transformation Functions
#
# Each function defines a data transformation specific to the statistic results. These
# functions are utility functions passed to transform_().
# Each function operates on a results file, which I have specifically not supplied as
# a parameter. It will be inherited in the calling scope within the transform_() fucntion
tajima_ <- function () {
    
    function (results) {
        colnames(results) <- c("chr", "BIN_START", "N_SNPS", "TajimaD", "pop", "chr")
        tmp_ <- results[, .(pop, chr, BIN_START, N_SNPS, TajimaD)]
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

get_utility_function <- function (base_directory) {
    util_ <- if (grepl("TD", base_directory)) tajima_()
             else if (grepl("FAWH", base_directory)) fanwu_()
             else tajima_()
    return (util_)
}

#### helper functions
db_execute <- function (query, username = "etl_user", host = "localhost", db = "selectionDW") {
    execute <- sprintf('mysql -u %s -h %s -D %s -e "%s"',
                       username, host, db, query)
    
    return (system(execute))
}

bulk_load <- function (filename, db_table) {
    query <- sprintf("
                     LOAD DATA LOCAL INFILE '%s' INTO TABLE %s
                     FIELDS TERMINATED BY ','
                     LINES TERMINATED BY '\n'
                     IGNORE 1 LINES;
                     ", filename, db_table)
    return (query)
}

init_ <- function () {
    
    # populate dimPopData
    db_execute(bulk_load("./data/populations.csv", "staging_pop"))
    db_execute("call populate_pop();")
    
    # populate dimExperiment
    db_execute("call populate_experiment();")
    
    # bulk load gene data
    db_execute(bulk_load("./data/dimGene.csv", "staging_gene"))
    
    # init dimPos
    db_execute("call populate_pos();")
    
    # init dimGene
    db_execute("call populate_gene();")
    db_execute("truncate staging_gene;")
    
    # init dimStat
    db_execute("call populate_stat();")
}

# Extract data from raw results files
# This requires a little bit of thought.
#   - currently I have two directories: TD and FAWH
#   - each directory contains sub-folders for each population
#   - the folder names map to dimPopData.pop
#   - each folder contains 1 file per chromosome, these are consistenly formatted
#   - the raw format of TD and FAWH results files are very different
#   - however, they can both be munged into a consistent format (pop, chr, start, end, statistic, n(size))

# for population, X, read all results files
extract_ <- function (population = NULL, base_directory = NULL) {
    
    print("        Extracting results files ...")
    root_ <- sprintf("%s/%s/results/", base_directory, population)
    chromosome_files <- list.files(root_, full.names = FALSE)
    
    parse_ <- function (filename) {
        prefix_ <- strsplit(filename, "\\.")[[1]][1]
        chr_ <- as.integer(substring(prefix_, first = 4))
        
        return (chr_)
    }
    raw <- rbindlist(lapply(chromosome_files, 
                                function (x) {
                                    tmp <- data.table(read.table(paste0(root_, x), 
                                                 skip = ifelse(grepl("FAWH", x), 5, 1),
                                                 header = FALSE))
                                    tmp[, pop := population]
                                    tmp[, chr := parse_(x)]
                                    
                                    return (tmp)
                                }))
    return (raw)
}
# end result: (chr, start, end, pop, stat, nsize, slide, type)
transform_ <- function (results, utility = NULL) {
    
    print("        Wrangling results files ...")
    schema <- c("pop", "chr", "chrom_start", "chrom_end", "nsize", "variable", "value")
    tmp <- utility(results)
    
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
    db_execute("call check_variants();")
    
    # unstage the data from the staging table into intraSel fact table
    db_execute(sprintf("call unstage(%s)", experiment_id))
    
    # truncate the staging table
    db_execute("truncate staging_results;")
    
}

pipeline <- function (pop, base, utility, experiment_id) {
    extract_(population = pop, base_directory = base) %>% 
        transform_(utility = utility) %>% 
        load_(experiment_id = experiment_id)   
}

main <- function (experiment_id = 1, init = FALSE, dir_ = "/mnt/DataDrive/Scratch/SelectionDW/data/SelectionResults") {
    
    # initialise the DW
    if (init) {
        print("Initialising the data warehouse...")
        init_()
    }
    
    # begin to load data files
    base_directories <- list.dirs(dir_, full.names = TRUE, recursive = FALSE)
    
    for (base_ in base_directories) {
        
        print(sprintf("----    %s    ----", base_))
        util_ <- get_utility_function(base_)
        
        for (pop_ in list.dirs(base_, full.names = FALSE, recursive = FALSE)) {
            print(sprintf("    %s -- ", pop_))
            try(pipeline(pop_, base_, util_, experiment_id))
        }
    }
}

main(experiment_id = 1, init = TRUE)

# an example of loading experiment 2
# main(experiment_id = 2)