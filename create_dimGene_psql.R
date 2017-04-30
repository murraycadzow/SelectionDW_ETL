# create dimGene_psql.csv
library(dplyr)
ensGene <- read.delim('data/ensembl_dimGene.txt', skip=1, header=TRUE, sep ='\t', stringsAsFactors = FALSE)

ensGene2 <- ensGene %>%mutate(gene_id = 1) %>%  select(gene_id ,"ensembl_gene_id"=GeneID, "hgnc_symbol"=HGNCsymbol, "chromosome_name" = Chromosome_scaffold.name, "start_position" = GeneStart, "end_position"= GeneEnd) 
ensGene3 <- ensGene2%>% filter(!is.na(ensembl_gene_id) & chromosome_name %in% 1:22)
write.table(ensGene3, file = 'data/dimGene_psql.csv', row.names = FALSE, col.names=TRUE, quote=FALSE, sep =',')
