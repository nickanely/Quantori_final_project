import pandas as pd
import pyarrow.dataset as ds
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

from Final_project.Airflow_scripts.agile_project.dags.script.processor.process_source_molecules import (
    ProcessSourceMolecules
)


class MockS3Saver:
    def __init__(self):
        self.s3fs_client = None

    def list_files(self) -> list[str]:
        return ['test_file.parquet']


class MockDataset:
    def __init__(self):
        self.data = [
            {'chembl_id': 'CHEMBL1', 'fingerprint': 'AAABBB'},
            {'chembl_id': 'CHEMBL2', 'fingerprint': 'CCC'},  # Invalid fingerprint
            {'chembl_id': 'CHEMBL1', 'fingerprint': 'AAABBB'}  # Duplicate
        ]

    def to_batches(self, **kwargs):
        df = pd.DataFrame(self.data)
        yield df


@pytest.fixture
def source_processor():
    s3_saver = MockS3Saver()
    processor = ProcessSourceMolecules(s3_saver=s3_saver, s3_fingerprints_folder="test_folder")
    return processor


def test_load_and_precompute_fingerprints(monkeypatch, source_processor):
    monkeypatch.setattr(ds, 'dataset', lambda *args, **kwargs: MockDataset())
    df = source_processor.load_and_precompute_fingerprints()

    # Assertions
    assert len(df) == 1
    assert df['chembl_id'].iloc[0] == 'CHEMBL1'
    assert isinstance(df['computed_fingerprint'].iloc[0], str)


def test_create_bit_vect():
    smiles = 'CCO'
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
    fingerprint_str = DataStructs.BitVectToText(fingerprint)

    bit_vect = ProcessSourceMolecules._create_bit_vect(fingerprint_str)

    assert isinstance(bit_vect, str)


def test_create_bit_vect_invalid():
    bit_vect = ProcessSourceMolecules._create_bit_vect("invalid_fingerprint")
    assert bit_vect is None