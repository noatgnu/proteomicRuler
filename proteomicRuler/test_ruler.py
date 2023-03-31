from unittest import TestCase

import pandas as pd

from proteomicRuler.ruler import Ruler, add_mw


class TestRuler(TestCase):
    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName)
        self.da = pd.read_csv(
            r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For copy number test\combined\txt\proteinGroups.txt",
            sep="\t")

    def test_init(self):
        ruler = Ruler(self.da)
        print(ruler.factor)
        print(ruler.organism)
        # ruler.load_genes()
        # ruler.df.to_csv("test.csv")
        # ruler.plot()


class Test(TestCase):
    def test_add_mw(self):
        df = pd.read_csv(r"report.pg_matrix.tsv", sep="\t")
        # used as unique index and to directly fetch mw data from UniProt
        mw_col = "Mass"
        # molecular weight column name
        ploidy = 2
        # ploidy number
        total_cellular_protein_concentration = 200
        intensity_columns = df.columns[5:]
        accession_id_col = "Protein.Ids"
        df = add_mw(df, accession_id_col)
        df = df[pd.notnull(df[mw_col])]
        df[mw_col] = df[mw_col].astype(float)
        ruler = Ruler(df, intensity_columns, mw_col, accession_id_col, ploidy, total_cellular_protein_concentration) #
        ruler.df.to_csv("output.txt", sep="\t", index=False)