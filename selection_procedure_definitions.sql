use selectionDW;

# initial population of dimExperiment
drop procedure if exists populate_pop;
delimiter //
create procedure populate_pop()
begin
	INSERT INTO dimPopData (pop, description, superPop, genotypePlatform)
        SELECT * from staging_pop;
end //
delimiter ;

# initial population of dimExperiment
drop procedure if exists populate_experiment;
delimiter //
create procedure populate_experiment()
begin
	INSERT INTO dimExperiment (description)
        VALUES ('AXIOM_unimputed'), 
               ('OMNI_unimputed'), 
               ('AXIOM and OMNI imputed info score threshold 0.3'), 
               ('AXIOM and OMNI imputed info score threshold 0.8'), 
               ('coreExome_v24_unimputed');
end //
delimiter ;

# initial population of dimPos
drop procedure if exists populate_pos;
delimiter //
create procedure populate_pos()
begin
	INSERT INTO dimPos (chrom, chrom_start, chrom_end)
        SELECT DISTINCT chrom, start_position, end_position
        FROM staging_gene;
end //
delimiter ;

# initial population of dimGene
drop procedure if exists populate_gene;
delimiter //
create procedure populate_gene()
begin
	INSERT INTO dimGene (posID, GeneName, EnsGeneID)
        SELECT DISTINCT pos.posID, stage.GeneName, stage.EnsGeneID
        FROM staging_gene as stage
        INNER JOIN dimPos as pos
          ON (
				pos.chrom = stage.chrom 
			AND pos.chrom_start = stage.start_position
            AND pos.chrom_end = stage.end_position
		  );
end //
delimiter ;

# initial population of dimStat
drop procedure if exists populate_stat;
delimiter //
create procedure populate_stat()
begin
	INSERT INTO dimStat (statName, statDescription)
        VALUES ('TajimaD', 'Tajimas D'), 
               ('NumSites_TajimasD', 'Number of sites in window for Tajimas D'), 
               ('NumSites_FayWuH', 'Number of sites in window for Fay and Wus H'), 
               ('S', 'Segregating sites'), 
               ('Eta', NULL),
               ('Eta_e', NULL),
               ('Pi', NULL),
               ('FuLi_D', 'Fu and Lis D'),
               ('FuLi_F', 'Fu and Lis F'),
               ('FayWu_H','Fay and Wus H');
end //
delimiter ;

# populate dimPos with new variants
drop procedure if exists check_variants;
delimiter //
create procedure check_variants()
begin
	INSERT IGNORE INTO dimPos (chrom, chrom_start, chrom_end)
        SELECT DISTINCT chrom, chrom_start, chrom_end
        FROM staging_results;
end //
delimiter ;

# take results from staging and populate intraSel
# NOTE: takes an INT as input parameter - this is the experiment ID
#       this must be determined programmatically during the loading script
drop procedure if exists unstage;
delimiter //
create procedure unstage (IN experiment_id int) 
begin
	INSERT INTO intraSel (posID, popID, statValue, statID, experimentID)
		SELECT d2.posID, d1.popID, st.statValue, d3.statID, experiment_id 
        FROM staging_results as st
          INNER JOIN dimPopData as d1 on d1.pop = st.pop
          INNER JOIN dimPos as d2 on (
				d2.chrom = st.chrom
			AND d2.chrom_start = st.chrom_start
            AND d2.chrom_end = st.chrom_end)
		  INNER JOIN dimStat as d3 on d3.statName = st.variable;
end //
delimiter ;
