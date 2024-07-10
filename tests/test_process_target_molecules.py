import pytest
from Final_project.Airflow_scripts.agile_project.dags.script.processor.process_target_molecules import (
    ProcessTargetMolecules
)


class MockS3Saver:
    def read_csv_from_s3(self, s3_path: str):
        test_data = """chembl_id,smiles
        CHEMBL123,CCO
        CHEMBL456,c1ccccc1
        CHEMBL789,CC(=O)O
        """
        return test_data


@pytest.fixture
def target_processor():
    return ProcessTargetMolecules(s3_saver=MockS3Saver(), input_folder="test_input")


def test_load_and_clean_csv_files(target_processor):
    file_list = ["test_file_1.csv", "test_file_2.csv"]
    df = target_processor.load_and_clean_csv_files(file_list)

    # Assertions
    assert len(df) == 3
    assert df['chembl_id'].tolist() == ['CHEMBL123', 'CHEMBL456', 'CHEMBL789']
    assert df['smiles'].tolist() == ['CCO', 'c1ccccc1', 'CC(=O)O']


def test_clean_csv(target_processor):
    content = """chembl_id,smiles
    CHEMBL1,C
    CHEMBL2,CC, extra comma
    CHEMBL3,C#C 
    CHEMBL4,
    CHEMBL5,CC(C)C(C)CαC
    """
    cleaned_df = target_processor._clean_csv(content, "test_file.csv")

    assert len(cleaned_df) == 3
    assert cleaned_df['chembl_id'].tolist() == ['CHEMBL1', 'CHEMBL3', 'CHEMBL5']
    assert cleaned_df['smiles'].tolist() == ['C', 'C#C', 'CC(C)C(C)CαC']
