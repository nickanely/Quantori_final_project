import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from Final_project.Airflow_scripts.agile_project.dags.script.processor.process_similarity_score import (
    ProcessSimilarityScore
)


class MockS3Saver:
    def save_similarity_results_in_parallel(self, similarity_results):
        pass


@pytest.fixture
def similarity_processor():
    # Sample DataFrames
    target_data = {'chembl_id': ['CHEMBL1', 'CHEMBL2'],
                   'smiles': ['CCO', 'c1ccccc1']}
    source_data = {'chembl_id': ['CHEMBL3', 'CHEMBL4', 'CHEMBL5'],
                   'computed_fingerprint': ['FP1', 'FP2', 'FP3']}

    target_df = pd.DataFrame(target_data)
    source_df = pd.DataFrame(source_data)
    return ProcessSimilarityScore(s3_saver=MockS3Saver(),
                                  target_mols_df=target_df,
                                  source_mols_df=source_df)


def test_parse_smiles():
    smiles = "CCO"
    mol = ProcessSimilarityScore._parse_smiles(smiles)
    assert mol is not None
    assert Chem.MolToSmiles(mol) == smiles


def test_parse_smiles_invalid():
    with pytest.raises(ValueError):
        ProcessSimilarityScore._parse_smiles("invalid-smiles")


def test_calculate_similarity_score(similarity_processor):
    # Create a sample target fingerprint
    target_smiles = 'CCO'
    target_mol = Chem.MolFromSmiles(target_smiles)
    target_fps = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, 1024)

    # Call the method
    similarity_df = similarity_processor._calculate_similarity_score(target_fps)

    assert not similarity_df.empty


def test_compute_top_10():
    # Example DataFrame
    data = {'chembl_id': ['CHEMBL1', 'CHEMBL2', 'CHEMBL3', 'CHEMBL4',
                          'CHEMBL5', 'CHEMBL6', 'CHEMBL7', 'CHEMBL8',
                          'CHEMBL9', 'CHEMBL10', 'CHEMBL11'],
            'similarity_score': [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
                                 0.1, 0.05, 0.05]}
    df = pd.DataFrame(data)
    top_10_df = ProcessSimilarityScore._compute_top_10(df)

    assert len(top_10_df) == 10
    assert top_10_df['has_duplicates_of_last_largest_score'].iloc[-1] == True


def test_compute_similarity_scores(similarity_processor):
    # This is more of an integration tests
    similarity_results = similarity_processor.compute_similarity_scores()

    assert len(similarity_results) == 2  # Two target molecules in the fixture

    for chembl_id, top_10_df in similarity_results:
        assert isinstance(chembl_id, str)
        assert isinstance(top_10_df, pd.DataFrame)
        assert len(top_10_df) <= 10
