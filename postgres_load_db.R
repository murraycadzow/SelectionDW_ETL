setwd('~/Git_repos/SelectionDW_ETL/')
source('selection_etl_script_postgres.R')

data_dir <- '/media/xsan/scratch/merrimanlab/murray/working_dir/coreExome_selection/NZ_coreExome_1kgp/data/dbload/'
exid <- 5
for(chr in 1:1){
  print(chr)

  print('tajd')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'tajd/chr',chr,'/'))
  print('fawh')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'fawh/chr',chr,'/'))
  print('ihs')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'ihs/norm/chr',chr,'/'))
  print('nsl')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'nsl/norm/chr',chr,'/'))
  print('daf')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'daf/chr',chr,'/'))
  print('fst')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'fst/chr',chr,'/'))
  print('xpehh')
  main(experiment_id = exid, results_directory = paste0(data_dir, 'xpehh/norm/chr',chr,'/'))
  
  print('_________________________________________')
  print("#########################################")
  print("#########################################")
  print('_________________________________________')
  
}

