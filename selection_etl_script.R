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
        tmp_[, slide := nsize]
        
        return (tmp_)
    }
}
fanwu_ <- function () {
    
    function (results) {
        tmp_ <- results[, c(4,5,7,9:17), with = FALSE]
        colnames(tmp_) <- c("chrom_start", "chrom_end", "nsize", "S", "Eta", "Eta_E", "Pi", "FuLi_D", "FuLi_F", "FayWu_H", "pop", "chr")
        
        tmp_ <- melt.data.table(tmp_, id.vars = c("pop", "chr", "chrom_start", "chrom_end", "nsize"))
        tmp_[, slide := nsize]
        
        return (tmp_)
    }
}

utility_switch <- function (file_) {
    
    util_ <- if (grepl(".taj_d$", file_)) list(id = "tajima", util_ = tajima_())
             else if (grepl(".faw$", file_)) list(id = "fanwu", util_ = fanwu_())
             else if (grepl(".fst$", file_)) list(id = "fst", util_ = fst_())
             else if (grepl(".kaks$", file_)) list(id = "kaks", util_ = kaks_())
             else if (grepl(".af$", file_)) list(id = "af", util_ = af_())
             else if (grepl(".nsl.", file_)) list(id = "nsl", util_ = nsl_())
             else if (grepl(".ihs.", file_)) list(id = "ihs", util_ = ihs_())
             else NULL
    
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
    chr_ <- as.integer(substring(prefix_, first = 4))
    pop_ <- substring(prefix_, first = 1, last = 3)
    
    return (data.table(chr = chr_, pop = pop_))
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
            meta <- parse_(f_)
            
            tmp <- data.table(read.table(f_, skip = header_(), header = FALSE))
            tmp[, pop := meta[["pop"]]]
            tmp[, chr := meta[["chr"]]]
            
            return (tmp)
        })
    )
                                    
    return (raw)
}
# end result: (chr, start, end, pop, stat, nsize, slide, type)
transform_ <- function (results, utility = NULL) {
    
    print("        Wrangling results files ...")
    schema <- c("pop", "chr", "chrom_start", "chrom_end", "nsize", "variable", "value")
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
    db_execute("call check_variants();")
    
    # unstage the data from the staging table into intraSel fact table
    db_execute(sprintf("call unstage(%s)", experiment_id))
    
    # truncate the staging table
    db_execute("truncate staging_results;")
    
}

pipeline <- function (experiment_id, files, utility) {
    extract_(files, utility) %>% 
        transform_(utility = utility) %>% 
        load_(experiment_id = experiment_id)   
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
    
    # Run pipeline per population
    # Note: cannot run over ALL populations at once, as this produces
    #       millions of rows of data.
    populations <- rbindlist(lapply(files_, parse_))
    for (pop in unique(populations[["pop"]])) {
        
        print(sprintf("----    %s    ----", pop))
        
        popfiles <- files_[grepl(pop, files_)]
        pipeline(experiment_id, popfiles, util_)
    }
}

### For testing only
test_pipeline <- function (directory = "/mnt/DataDrive/Scratch/SelectionDW/") {
    
    setwd(directory)
    
    print("Initialising the data warehouse")
    init_()
    
    main(experiment_id = 1, "./data/SelectionResults/FAWH")
}

####    test_pipeline()
