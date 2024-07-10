from typing import Optional

from pydantic import BaseModel


class ChemblIdLookup(BaseModel):
    chembl_id: str
    entity_type: str
    last_active: Optional[float]
    resource_url: str
    status: str


class MoleculeProperties(BaseModel):
    alogp: Optional[str] = None
    aromatic_rings: Optional[float] = None
    cx_logd: Optional[str] = None
    cx_logp: Optional[str] = None
    cx_most_apka: Optional[str] = None
    cx_most_bpka: Optional[str] = None
    full_molformula: Optional[str] = None
    full_mwt: Optional[str] = None
    hba: Optional[float] = None
    hba_lipinski: Optional[float] = None
    hbd: Optional[float] = None
    hbd_lipinski: Optional[float] = None
    heavy_atoms: Optional[float] = None
    molecular_species: Optional[str] = None
    mw_freebase: Optional[str] = None
    mw_monoisotopic: Optional[str] = None
    np_likeness_score: Optional[str] = None
    num_lipinski_ro5_violations: Optional[float] = None
    num_ro5_violations: Optional[float] = None
    psa: Optional[str] = None
    qed_weighted: Optional[str] = None
    ro3_pass: Optional[str] = None
    rtb: Optional[float] = None


class MoleculeStructures(BaseModel):
    canonical_smiles: Optional[str] = None
    molfile: Optional[str] = None
    standard_inchi: Optional[str] = None
    standard_inchi_key: Optional[str] = None


class MoleculeDictionary(BaseModel):
    molecule_chembl_id: str
    availability_type: Optional[float] = None
    black_box_warning: Optional[int] = None
    chebi_par_id: Optional[float] = None
    chemical_probe: Optional[int] = None
    chirality: Optional[int] = None
    dosed_ingredient: Optional[bool] = None
    first_approval: Optional[float] = None
    first_in_class: Optional[int] = None
    helm_notation: Optional[str] = None
    indication_class: Optional[str] = None
    inorganic_flag: Optional[int] = None
    max_phase: Optional[str] = None
    molecule_type: Optional[str] = None
    natural_product: Optional[int] = None
    oral: Optional[bool] = None
    parenteral: Optional[bool] = None
    polymer_flag: Optional[int] = None
    pref_name: Optional[str] = None
    prodrug: Optional[int] = None
    structure_type: Optional[str] = None
    therapeutic_flag: Optional[bool] = None
    topical: Optional[bool] = None
    usan_stem: Optional[str] = None
    usan_stem_definition: Optional[str] = None
    usan_substem: Optional[str] = None
    usan_year: Optional[float] = None
    withdrawn_flag: Optional[bool] = None
    molecule_properties: Optional[MoleculeProperties] = None
    molecule_structures: Optional[MoleculeStructures] = None
