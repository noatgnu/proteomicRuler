from unittest import TestCase

import pandas as pd

from proteomicRuler.ruler import Ruler


class TestRuler(TestCase):
    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName)
        self.da = pd.read_csv(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For copy number test\combined\txt\proteinGroups.txt", sep="\t")


    def test_init(self):
        ruler = Ruler(self.da)
        print(ruler.factor)
        print(ruler.organism)
        #ruler.load_genes()
        #ruler.df.to_csv("test.csv")
        #ruler.plot()