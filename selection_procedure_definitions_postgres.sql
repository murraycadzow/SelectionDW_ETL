\connect selectiondw;

---# initial population of dimExperiment
CREATE OR REPLACE FUNCTION public.populate_pop()
    RETURNS void
    LANGUAGE 'plpgsql'
    VOLATILE NOT LEAKPROOF 
AS $function$
BEGIN
	INSERT INTO dimPopData (pop, description, superPop, genotypePlatform)
        SELECT * from staging_pop;
END;
$function$;

ALTER FUNCTION public.populate_pop()
    OWNER TO postgres;

GRANT ALL ON FUNCTION public.populate_pop() TO PUBLIC;


---# initial population of dimExperiment
CREATE OR REPLACE FUNCTION public.populate_experiment()
    RETURNS void
    LANGUAGE 'plpgsql'
    NOT LEAKPROOF 
AS $function$
BEGIN
	INSERT INTO dimExperiment (description)
        VALUES ('AXIOM_unimputed'), 
               ('OMNI_unimputed'), 
               ('AXIOM and OMNI imputed info score threshold 0.3'), 
               ('AXIOM and OMNI imputed info score threshold 0.8'), 
               ('coreExome_v24_unimputed');
END
$function$;

ALTER FUNCTION public.populate_experiment()
    OWNER TO postgres;

GRANT ALL ON FUNCTION public.populate_experiment() TO PUBLIC;

---# initial population of dimPos
CREATE OR REPLACE FUNCTION public.populate_pos()
    RETURNS void
    LANGUAGE 'plpgsql'
    NOT LEAKPROOF 
AS $function$
begin
	INSERT INTO dimPos (chrom, chrom_start, chrom_end)
        SELECT DISTINCT chrom, start_position, end_position
        FROM staging_gene;
end 
$function$;

ALTER FUNCTION public.populate_pos()
    OWNER TO postgres;

GRANT ALL ON FUNCTION public.populate_pos() TO PUBLIC;


---# initial population of dimGene
CREATE OR REPLACE FUNCTION public.populate_gene()
    RETURNS void
    LANGUAGE 'plpgsql'
    NOT LEAKPROOF 
AS $function$
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
end
$function$;

ALTER FUNCTION public.populate_gene()
    OWNER TO postgres;

GRANT ALL ON FUNCTION public.populate_gene() TO PUBLIC;

---# initial population of dimStat
CREATE OR REPLACE FUNCTION public.populate_stat()
    RETURNS void
    LANGUAGE 'plpgsql'
    NOT LEAKPROOF 
AS $function$
begin
	INSERT INTO dimStat (statName, statDescription)
        VALUES ('TajimaD', 'Tajimas D'), 
               ('NumSites_TajimasD', 'Number of sites in window for Tajimas D'), 
               ('NumSites_FayWuH', 'Number of sites in window for Fay and Wus H'), 
               ('S', 'Segregating sites'), 
               ('Eta', NULL),
               ('Eta_E', NULL),
               ('Pi', NULL),
               ('FuLi_D', 'Fu and Lis D'),
               ('FuLi_F', 'Fu and Lis F'),
               ('FayWu_H','Fay and Wus H'),
               ('ka', 'KAKS statistic'),
               ('ks', 'KAKS statistic'),
               ('kcomputed', 'KAKS ka / ks + 1'),
               ('MAF', 'AF statistic'),
               ('DAF', 'AF statistic'),
               ('nsl_freq_1', NULL),
               ('sL1', NULL),
               ('sL0', NULL),
               ('unstd_nsl', NULL),
               ('norm_nsl', NULL),
               ('significant_nsl', NULL),
               ('ihs_freq_1', NULL),
               ('ihh1', NULL),
               ('ihh0', NULL),
               ('unstd_ihs', NULL),
               ('norm_ihs', NULL),
               ('significant_ihs', NULL);
end
$function$;

ALTER FUNCTION public.populate_stat()
    OWNER TO postgres;

GRANT ALL ON FUNCTION public.populate_stat() TO PUBLIC;

---# populate dimPos with new variants
CREATE OR REPLACE FUNCTION public.check_variants()
    RETURNS void
    LANGUAGE 'plpgsql'
    VOLATILE

AS $function$begin
	INSERT INTO dimpos (chrom, chrom_start, chrom_end)
        SELECT DISTINCT chrom, chrom_start, chrom_end
        FROM staging_results
	ON CONFLICT  (chrom, chrom_start, chrom_end)
	DO NOTHING;
end
$function$;


---# take results from staging and populate intraSel
---# NOTE: takes an INT as input parameter - this is the experiment ID
---#       this must be determined programmatically during the loading script
CREATE OR REPLACE FUNCTION public.unstage(IN experiment_id integer)
    RETURNS void
    LANGUAGE 'plpgsql'
    NOT LEAKPROOF 
AS $function$
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
end
$function$;

ALTER FUNCTION public.unstage()
    OWNER TO postgres;

GRANT ALL ON FUNCTION public.unstage() TO PUBLIC;
