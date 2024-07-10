-- Create the molecules dimension table
CREATE TABLE dim_molecules (
    chembl_id VARCHAR(50) PRIMARY KEY,
    molecule_type VARCHAR(50),
    mw_freebase FLOAT,
    alogp FLOAT,
    psa FLOAT,
    cx_logp FLOAT,
    molecular_species VARCHAR(50),
    full_mwt FLOAT,
    aromatic_rings INT,
    heavy_atoms INT
);

-- Create the molecule_similarities fact table
CREATE TABLE fact_molecule_similarities (
    target_chembl_id VARCHAR(50),
    source_chembl_id VARCHAR(50),
    tanimoto_similarity FLOAT,
    has_duplicates_of_last_largest_score BOOLEAN,
    PRIMARY KEY (source_chembl_id, target_chembl_id)
    FOREIGN KEY (source_chembl_id) REFERENCES dim_molecules(chembl_id),
    FOREIGN KEY (target_chembl_id) REFERENCES dim_molecules(chembl_id)
);