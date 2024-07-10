-- N7 - A
create or replace view average_score_per_mol as
select source_chembl_id, avg(tanimoto_similarity)
from fact_molecule_similarities
group by source_chembl_id
;

-- N7 - B
create or replace view average_deviation_of_alogp as
with molecule_properties as (
	select chembl_id, alogp
	from dim_molecules
),
similarity_data AS (
    SELECT
        f.source_chembl_id,
        f.target_chembl_id,
        mp_source.alogp AS source_alogp,
        mp_target.alogp AS target_alogp,
        ABS(mp_target.alogp - mp_source.alogp) AS alogp_deviation
    FROM
        fact_molecule_similarities f
    INNER JOIN
        molecule_properties mp_source ON f.source_chembl_id = mp_source.chembl_id
    INNER JOIN
        molecule_properties mp_target ON f.target_chembl_id = mp_target.chembl_id
)
SELECT
    source_chembl_id,
    AVG(alogp_deviation) AS average_alogp_deviation
FROM
    similarity_data
GROUP BY
    source_chembl_id;



-- N8
-- A
create or replace view pivot_similarity_score as
(
select target_chembl_id,
       max(case when source_chembl_id = 'CHEMBL2' then tanimoto_similarity end)      as "CHEMBL2",
       max(case when source_chembl_id = 'CHEMBL112998' then tanimoto_similarity end) as "CHEMBL112998",
       max(case when source_chembl_id = 'CHEMBL267915' then tanimoto_similarity end) as "CHEMBL267915",
       max(case when source_chembl_id = 'CHEMBL6200' then tanimoto_similarity end)   as "CHEMBL6200",
       max(case when source_chembl_id = 'CHEMBL414181' then tanimoto_similarity end) as "CHEMBL414181",
       max(case when source_chembl_id = 'CHEMBL6363' then tanimoto_similarity end)   as "CHEMBL6363",
       max(case when source_chembl_id = 'CHEMBL6212' then tanimoto_similarity end)   as "CHEMBL6212",
       max(case when source_chembl_id = 'CHEMBL6366' then tanimoto_similarity end)   as "CHEMBL6366",
       max(case when source_chembl_id = 'CHEMBL269132' then tanimoto_similarity end) as "CHEMBL269132",
       max(case when source_chembl_id = 'CHEMBL6334' then tanimoto_similarity end)   as "CHEMBL6334"
from fact_molecule_similarities
where fact_molecule_similarities.source_chembl_id in
      ('CHEMBL2',
       'CHEMBL112998',
       'CHEMBL267915',
       'CHEMBL6200',
       'CHEMBL414181',
       'CHEMBL6363',
       'CHEMBL6212',
       'CHEMBL6366',
       'CHEMBL269132',
       'CHEMBL6334')
group by target_chembl_id
order by target_chembl_id
    );

-- B
create or replace view next_most_similar_molexule as
select
    source_chembl_id,
    target_chembl_id,
    tanimoto_similarity,
    lead(tanimoto_similarity, 1) over
    	(partition by target_chembl_id order by tanimoto_similarity desc) AS next_most_similar_target,
    nth_value(target_chembl_id, 2)
    	over (partition by source_chembl_id order by tanimoto_similarity desc rows between unbounded preceding and current row)
    	as second_most_similar_target
from fact_molecule_similarities
where source_chembl_id != target_chembl_id;

-- C


CREATE OR REPLACE VIEW grouped_averages AS
WITH grouped_data AS (
    SELECT
        ms.source_chembl_id,
        mp.aromatic_rings,
        mp.heavy_atoms,
        AVG(ms.tanimoto_similarity) AS avg,
        GROUPING(ms.source_chembl_id) AS source,
        GROUPING(mp.aromatic_rings) AS aromatic,
        GROUPING(mp.heavy_atoms) AS heavy
    FROM
        fact_molecule_similarities ms
    JOIN
        staging_compound_properties mp ON mp.chembl_id = ms.source_chembl_id
    GROUP BY
        GROUPING SETS (
            (ms.source_chembl_id),
            (mp.aromatic_rings, mp.heavy_atoms),
            (mp.heavy_atoms),
            ()
        )
)
SELECT
    CASE
        WHEN source = 0 THEN source_chembl_id
        ELSE 'TOTAL'
    END AS source_molecule,
    CASE
        WHEN aromatic = 0 THEN aromatic_rings::TEXT
        ELSE 'TOTAL'
    END AS aromatic_rings,
    CASE
        WHEN heavy = 0 THEN heavy_atoms::TEXT
        ELSE 'TOTAL'
    END AS heavy_atoms,
    avg
FROM
    grouped_data