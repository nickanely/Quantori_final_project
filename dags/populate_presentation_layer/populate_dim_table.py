POPULATE_DIM_TABLE = """
    INSERT INTO dim_molecules (
        chembl_id,
        molecule_type,
        mw_freebase,
        alogp,
        psa,
        cx_logp,
        molecular_species,
        full_mwt,
        aromatic_rings,
        heavy_atoms
    )
    SELECT DISTINCT
        scp.chembl_id,
        smd.molecule_type,
        scp.mw_freebase,
        scp.alogp,
        scp.psa,
        scp.cx_logp,
        scp.molecular_species,
        scp.full_mwt,
        scp.aromatic_rings,
        scp.heavy_atoms
    FROM staging_compound_properties scp
    INNER JOIN staging_molecule_dictionary smd ON scp.chembl_id = smd.chembl_id
    WHERE scp.chembl_id IN (
        SELECT DISTINCT source_chembl_id AS chembl_id
        FROM fact_molecule_similarities
        UNION
        SELECT DISTINCT target_chembl_id AS chembl_id
        FROM fact_molecule_similarities
    );
"""
