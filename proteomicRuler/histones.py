import pandas as pd
import numpy as np

class HistoneDB:
    def __init__(self):
        self.df = pd.DataFrame()

    def get_histones(self):
        with open(r"C:\Users\toanp\PycharmProjects\proteomicRuler\proteomicRuler\organisms.json", "rb") as jsonfile:
            df = pd.read_json(jsonfile).T
            df = df.explode("histone_ids").reset_index(drop=True)
            df.set_index("histone_ids", inplace=True)
            self.df = df

    def check_histones(self, ids, organism=None):
        share = np.intersect1d(ids, self.df.index)
        result = self.df.loc[share]
        if organism:
            result = result[result["name"] == organism]
        return result

    def get_organism(self, ids):
        result = self.check_histones(ids)
        data = result.groupby(["name", "genome_size"]).size()
        return data.sort_values().index[-1]



