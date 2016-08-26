--- DDL script for selection database
---
--- Murray Cadzow
--- June 2016
drop database if exists selectiondw;
create database selectiondw;
\connect selectiondw;

---#### STAGING AREA ####
drop table if exists staging_pop;
create table staging_pop
(
	pop varchar(10) not null,
    description varchar(140),
    superPop varchar(3),
    genotypePlatform varchar(32)
);

drop table if exists staging_gene;
create table staging_gene
(
	gene_id int,
    EnsGeneID varchar(32),
    GeneName varchar(32),
    chrom smallint,
    start_position int,
    end_position int
);

drop table if exists staging_results;
create table staging_results
(
	pop varchar(3),
    chrom smallint,
    chrom_start int,
    chrom_end int,
    nsize int,
    variable varchar(20),
    statValue real
);
create index idx_stage_results on staging_results (chrom, chrom_start, chrom_end);

---#### Main DW ####
---# position table
---# Contains the Chromosome and Allele information for each SNP
---# notes:  bp is the chromosomal position in bases
drop table if exists dimPos;
create table dimPos (
    posID serial not null,
    chrom smallint,
    chrom_start int,
    chrom_end int,
	primary key (posID)
);
create index idx_pos_chromosome on dimPos (chrom, chrom_start, chrom_end);
alter table dimPos
  add constraint uniq_coords unique (chrom, chrom_start, chrom_end);


---# dimPopData
---# contains the information about the populations
drop table if exists dimPopData;
create table dimPopData (
	popID serial not null,
    pop varchar(10) not null,
    description varchar(140),
    superPop varchar(3),
    genotypePlatform varchar(32),
	primary key (popID)
);

---# dimStat
---# selection statistic names
drop table if exists dimStat;
create table dimStat (
    statID serial not null,
    statName varchar(20) not null,
    statDescription varchar(128),
    primary key (statID)
);

---# dimGene
---# positions of genes by ensembl ID or gene name
drop table if exists dimGene;
create table dimGene (
    posID int not null,
    GeneName varchar(32),
    EnsGeneID varchar(32) not NULL,
    primary key (GeneName,EnsGeneID),
    foreign key (posID) references dimPos(posID)
);

---# dimExperiment
drop table if exists dimExperiment;
create table dimExperiment (
	experimentID serial not null,
    description varchar(200),
    primary key (experimentID)
);



---# intraSel
---# store the inter population stats
drop table if exists intraSel;
create table intraSel (
    posID int not null,
    popID int not null,
    statValue real,
    statID int not null,
    experimentID int not null,
    primary key (popID, posID, statID, experimentID),
    foreign key (posID) references dimPos(posID),
    foreign key (popID) references dimPopData(popID),
    foreign key (statID) references dimStat(statID),
    foreign key (experimentID) references dimExperiment(experimentID)
);

---# intraSel
---# store the inter population stats
drop table if exists interSel;
create table interSel (
    posID int not null,
    popID1 int not null,
    popID2 int not null,
    statValue real,
    statID int not null,
    experimentID int not null,
    foreign key (posID) references dimPos(posID),
    foreign key (popID1) references dimPopData(popID),
    foreign key (popID2) references dimPopData(popID),
    foreign key (statID) references dimStat(statID),
    foreign key (experimentID) references dimExperiment(experimentID)
);
