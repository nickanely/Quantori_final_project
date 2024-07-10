from typing import List, Dict, Tuple

import pandas as pd

from build_module.chembl_data_processing.models.pydantic_model import (
    ChemblIdLookup,
    MoleculeDictionary,
)


def process_chembl_id_lookup_chunk(records: List[Dict]) -> pd.DataFrame:
    validated_records = [ChemblIdLookup(**record).model_dump() for record in records]
    return pd.DataFrame(validated_records)


def process_molecule_chunk(records: List[Dict]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    molecule_dictionary_data = []
    compound_properties_data = []
    compound_structures_data = []

    for record in records:
        molecule = MoleculeDictionary(**record)
        molecule_dictionary_data.append({
            'chembl_id': record['molecule_chembl_id'],
            **molecule.model_dump(exclude={'molecule_chembl_id', 'molecule_properties', 'molecule_structures'})
        })

        if molecule.molecule_properties:
            compound_properties_data.append({
                'chembl_id': record['molecule_chembl_id'],
                **molecule.molecule_properties.model_dump()
            })

        if molecule.molecule_structures:
            compound_structures_data.append({
                'chembl_id': record['molecule_chembl_id'],
                **molecule.molecule_structures.model_dump()
            })

    return (
        pd.DataFrame(molecule_dictionary_data),
        pd.DataFrame(compound_properties_data),
        pd.DataFrame(compound_structures_data),
    )
