import pytest

from script.build_dwh.services.db_handler import (
    process_molecule_chunk,
    process_chembl_id_lookup_chunk
)

MOCK_RECORD = {
    "chembl_id": "CHEMBL1",
    "entity_type": "COMPOUND",
    "last_active": 34,
    "resource_url": "/chembl/api/data/molecule/CHEMBL1",
    "status": "ACTIVE"
}

MOCK_MOLECULE_RECORD = {
    'molecule_chembl_id': 'CHEMBL414196',
    'availability_type': -1,
    'black_box_warning': 0,
    'chebi_par_id': None,
    'chemical_probe': 0,
    'chirality': -1,
    'dosed_ingredient': False,
    'first_approval': None,
    'first_in_class': -1,
    'helm_notation': None,
    'indication_class': None,
    'inorganic_flag': -1,
    'max_phase': None,
    'molecule_type': 'Small molecule',
    'natural_product': 0,
    'oral': False,
    'orphan': -1,
    'parenteral': False,
    'polymer_flag': 0,
    'pref_name': None,
    'prodrug': -1,
    'structure_type': 'MOL',
    'therapeutic_flag': False,
    'topical': False,
    'usan_stem': None,
    'usan_stem_definition': None,
    'usan_substem': None,
    'usan_year': None,
    'withdrawn_flag': False,
    'molecule_properties': {
        'alogp': '2.68',
        'aromatic_rings': 3,
        'cx_logd': '2.71',
        'cx_logp': '2.94',
        'cx_most_apka': 'NULL',
        'cx_most_bpka': '7.24',
        'full_molformula': 'C24H27N5O3',
        'full_mwt': '433.51',
        'hba': 7,
        'hba_lipinski': 8,
        'hbd': 1,
        'hbd_lipinski': 2,
        'heavy_atoms': 32,
        'molecular_species': 'NEUTRAL',
        'mw_freebase': '433.51',
        'mw_monoisotopic': '433.2114',
        'np_likeness_score': '-0.89',
        'num_lipinski_ro5_violations': 0,
        'num_ro5_violations': 0,
        'psa': '93.81',
        'qed_weighted': '0.66',
        'ro3_pass': 'N',
        'rtb': 5
    },
    'molecule_structures': {
        'canonical_smiles': 'COc1cc2nc(N3CCN(C(=O)C4CC4c4ccccc4)CC3)nc(N)c2cc1OC',
        'molfile': "\n     RDKit          2D\n\n 32 36  0  0  0  0  0  0  0  0999 V2000\n    5.5375   -0.1792    "
                   "0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9667   -1.4250    0.0000 C   0  0  0  0  0  "
                   "0  0  0  0  0  0  0\n    1.9750   -2.2542    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    "
                   "6.3625   -0.1667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.8250    0.2333    0.0000 "
                   "C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5417   -2.2542    0.0000 C   0  0  0  0  0  0  0  0 "
                   " 0  0  0  0\n    5.9417   -0.8917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2625   "
                   "-2.6667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2542   -1.0167    0.0000 N   0  0 "
                   " 0  0  0  0  0  0  0  0  0  0\n    0.5417   -1.4292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0 "
                   " 0\n    2.6792   -1.0125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1042   -0.1792   "
                   " 0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1708   -2.6750    0.0000 C   0  0  0  0  0  "
                   "0  0  0  0  0  0  0\n   -0.1708   -1.0167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   "
                   "-0.8833   -2.2625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8833   -1.4375    "
                   "0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.8167    1.0583    0.0000 O   0  0  0  0  0  "
                   "0  0  0  0  0  0  0\n    7.0750    0.2458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    "
                   "3.4000   -1.4167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6792   -0.1875    0.0000 "
                   "C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.3917    0.2333    0.0000 C   0  0  0  0  0  0  0  0 "
                   " 0  0  0  0\n    4.1125   -1.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2542   "
                   "-3.4917    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6000   -2.6750    0.0000 O   0  0 "
                   " 0  0  0  0  0  0  0  0  0  0\n   -1.5958   -1.0250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0 "
                   " 0\n    7.7875   -0.1667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.0625    1.0708   "
                   " 0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.3125   -1.4375    0.0000 C   0  0  0  0  0  "
                   "0  0  0  0  0  0  0\n   -2.3125   -2.2625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    "
                   "7.7750    1.4833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.4917    0.2500    0.0000 "
                   "C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.4917    1.0750    0.0000 C   0  0  0  0  0  0  0  0 "
                   " 0  0  0  0\n  2 11  1  0\n  3  2  1  0\n  4  1  1  0\n  5  1  1  0\n  6 10  1  0\n  7  1  1  0\n "
                   " 8  3  2  0\n  9  2  2  0\n 10  9  1  0\n 11 19  1  0\n 12  5  1  0\n 13  6  2  0\n 14 10  2  0\n "
                   "15 16  2  0\n 16 14  1  0\n 17  5  2  0\n 18  4  1  0\n 19 22  1  0\n 20 21  1  0\n 21 12  1  0\n "
                   "22 12  1  0\n 23  8  1  0\n 24 15  1  0\n 25 16  1  0\n 26 18  2  0\n 27 18  1  0\n 28 25  1  0\n "
                   "29 24  1  0\n 30 27  2  0\n 31 26  1  0\n 32 30  1  0\n  7  4  1  0\n 20 11  1  0\n 31 32  2  0\n "
                   " 6  8  1  0\n 13 15  1  0\nM  END\n\u003E \u003Cchembl_id\u003E\nCHEMBL414196\n\n\u003E "
                   "\u003Cchembl_pref_name\u003E\nNone",
        'standard_inchi': "InChI=1S/C24H27N5O3",
        'standard_inchi_key': "QSRCXSDOJVDQBI-UHFFFAOYSA-N"
    }
}


@pytest.mark.parametrize("records, expected_length", [
    ([MOCK_RECORD], 1),
    ([], 0)
])
def test_process_chembl_id_lookup_chunk(records, expected_length):
    df = process_chembl_id_lookup_chunk(records)
    assert len(df) == expected_length
    if records:
        assert df.iloc[0]['chembl_id'] == 'CHEMBL1'


@pytest.mark.parametrize("records, expected_lengths", [
    ([MOCK_MOLECULE_RECORD], (1, 1, 1)),
    ([], (0, 0, 0)),
    ([MOCK_MOLECULE_RECORD, MOCK_MOLECULE_RECORD], (2, 2, 2)),
])
def test_process_molecule_chunk(records, expected_lengths):
    molecule_dict_df, compound_properties_df, compound_structures_df = process_molecule_chunk(records)

    assert len(molecule_dict_df) == expected_lengths[0]
    assert len(compound_properties_df) == expected_lengths[1]
    assert len(compound_structures_df) == expected_lengths[2]

    if records:
        assert molecule_dict_df.iloc[0]['chembl_id'] == 'CHEMBL414196'
        if 'molecule_properties' in records[0]:
            assert compound_properties_df.iloc[0]['alogp'] == '2.68'
        if 'molecule_structures' in records[0]:
            assert compound_structures_df.iloc[0]['standard_inchi_key'] == "QSRCXSDOJVDQBI-UHFFFAOYSA-N"


def test_process_large_molecule_chunk():
    records = [MOCK_MOLECULE_RECORD] * 1000
    molecule_dict_df, compound_properties_df, compound_structures_df = process_molecule_chunk(records)

    assert len(molecule_dict_df) == 1000
    assert len(compound_properties_df) == 1000
    assert len(compound_structures_df) == 1000


def test_dataframe_structure():
    records = [MOCK_MOLECULE_RECORD]
    molecule_dict_df, compound_properties_df, compound_structures_df = process_molecule_chunk(records)

    # Verify columns in molecule_dictionary DataFrame
    expected_columns = [
        'chembl_id', 'availability_type', 'black_box_warning', 'chebi_par_id', 'chemical_probe',
        'chirality', 'dosed_ingredient', 'first_approval', 'first_in_class', 'helm_notation', 'indication_class',
        'inorganic_flag', 'max_phase', 'molecule_type', 'natural_product', 'oral', 'orphan', 'parenteral',
        'polymer_flag', 'pref_name', 'prodrug', 'structure_type', 'therapeutic_flag', 'topical', 'usan_stem',
        'usan_stem_definition', 'usan_substem', 'usan_year', 'withdrawn_flag'
    ]
    assert list(molecule_dict_df.columns) == expected_columns

    # Verify columns in compound_properties DataFrame
    expected_columns = [
        'chembl_id', 'alogp', 'aromatic_rings', 'cx_logd', 'cx_logp', 'cx_most_apka', 'cx_most_bpka',
        'full_molformula', 'full_mwt', 'hba', 'hba_lipinski', 'hbd', 'hbd_lipinski', 'heavy_atoms', 'molecular_species',
        'mw_freebase', 'mw_monoisotopic', 'np_likeness_score', 'num_lipinski_ro5_violations', 'num_ro5_violations',
        'psa', 'qed_weighted', 'ro3_pass', 'rtb'
    ]
    assert list(compound_properties_df.columns) == expected_columns

    # Verify columns in compound_structures DataFrame
    expected_columns = ['chembl_id', 'canonical_smiles', 'molfile', 'standard_inchi', 'standard_inchi_key']
    assert list(compound_structures_df.columns) == expected_columns
