-- Staging tables
CREATE TABLE staging_chembl_id_lookup
(
    chembl_id VARCHAR(50),
    entity_type VARCHAR(50),
    last_active FLOAT,
    resource_url VARCHAR(50),
    status VARCHAR(50)
);

CREATE TABLE IF NOT EXISTS naneli.staging_compound_properties
(
    chembl_id text COLLATE pg_catalog."default" NOT NULL,
    alogp text COLLATE pg_catalog."default",
    aromatic_rings double precision,
    cx_logd text COLLATE pg_catalog."default",
    cx_logp text COLLATE pg_catalog."default",
    cx_most_apka text COLLATE pg_catalog."default",
    cx_most_bpka text COLLATE pg_catalog."default",
    full_molformula text COLLATE pg_catalog."default",
    full_mwt text COLLATE pg_catalog."default",
    hba double precision,
    hba_lipinski double precision,
    hbd double precision,
    hbd_lipinski double precision,
    heavy_atoms double precision,
    molecular_species text COLLATE pg_catalog."default",
    mw_freebase text COLLATE pg_catalog."default",
    mw_monoisotopic text COLLATE pg_catalog."default",
    np_likeness_score text COLLATE pg_catalog."default",
    num_lipinski_ro5_violations double precision,
    num_ro5_violations double precision,
    psa text COLLATE pg_catalog."default",
    qed_weighted text COLLATE pg_catalog."default",
    ro3_pass text COLLATE pg_catalog."default",
    rtb double precision,
    CONSTRAINT compound_properties_pk_1 PRIMARY KEY (chembl_id)
);



CREATE TABLE IF NOT EXISTS naneli.staging_molecule_dictionary
(
    chembl_id text COLLATE pg_catalog."default" NOT NULL,
    availability_type double precision,
    black_box_warning bigint,
    chebi_par_id double precision,
    chemical_probe bigint,
    chirality bigint,
    dosed_ingredient boolean,
    first_approval double precision,
    first_in_class bigint,
    helm_notation text COLLATE pg_catalog."default",
    indication_class text COLLATE pg_catalog."default",
    inorganic_flag bigint,
    max_phase text COLLATE pg_catalog."default",
    molecule_type text COLLATE pg_catalog."default",
    natural_product bigint,
    oral boolean,
    orphan bigint,
    parenteral boolean,
    polymer_flag bigint,
    pref_name text COLLATE pg_catalog."default",
    prodrug bigint,
    structure_type text COLLATE pg_catalog."default",
    therapeutic_flag boolean,
    topical boolean,
    usan_stem text COLLATE pg_catalog."default",
    usan_stem_definition text COLLATE pg_catalog."default",
    usan_substem text COLLATE pg_catalog."default",
    usan_year double precision,
    withdrawn_flag boolean,
    CONSTRAINT staging_molecule_dictionary_pkey PRIMARY KEY (chembl_id)
);




CREATE TABLE IF NOT EXISTS naneli.staging_compound_structures
(
    chembl_id text COLLATE pg_catalog."default" NOT NULL,
    canonical_smiles text COLLATE pg_catalog."default",
    molfile text COLLATE pg_catalog."default",
    standard_inchi text COLLATE pg_catalog."default",
    standard_inchi_key text COLLATE pg_catalog."default",
    CONSTRAINT compound_structures_pk_1 PRIMARY KEY (chembl_id)
);