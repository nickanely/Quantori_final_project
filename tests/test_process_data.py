import pandas as pd

from script.build_dwh.process_data_DEPRECATED import (
    process_chembl_id_lookup_chunk,
    process_molecule_chunk
)


def test_process_chembl_id_lookup_chunk():
    records = [{'chembl_id': 'CHEMBL1', 'entity_type': 'Small Molecule', 'last_active': '2023-01-01',
                'resource_url': 'http://example.com', 'status': 'active'}]
    df = process_chembl_id_lookup_chunk(records)
    assert isinstance(df, pd.DataFrame)
    assert df.iloc[0]['chembl_id'] == 'CHEMBL1'


def test_process_molecule_chunk():
    records = [{'molecule_chembl_id': 'CHEMBL1', 'molecule_properties': {'alogp': 1.0},
                'molecule_structures': {'canonical_smiles': 'C'}}]
    molecule_df, _, _ = process_molecule_chunk(records)
    assert isinstance(molecule_df, pd.DataFrame)
    assert molecule_df.iloc[0]['molecule_chembl_id'] == 'CHEMBL1'
