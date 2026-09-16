from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def test_synthetic_table_shape():
    assert pd.read_excel(ROOT/"data/synthetic/PDAC_In_Silico_Simulation_Data(101).xlsx").shape==(100,9)
def test_marker_table_has_20_rows():
    assert len(pd.read_csv(ROOT/"data/derived/cluster_7_differential_markers.csv"))==20
def test_mmp9_not_in_deposited_top20():
    df=pd.read_csv(ROOT/"data/derived/cluster_7_differential_markers.csv"); assert "MMP9" not in set(df["gene"].astype(str))
