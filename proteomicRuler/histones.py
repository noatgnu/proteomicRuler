import pandas as pd
import numpy as np
import os

class HistoneDB:
    def __init__(self):
        self.df = pd.DataFrame()

    def get_histones(self):
        # get histone protein and organism ids from json file in the same directory as this file
        d, _ = os.path.split(__file__)
        with open(os.path.join(d, "organisms.json"), "rb") as jsonfile:
            df = pd.read_json(jsonfile).T
            df = df.explode("histone_ids").reset_index(drop=True)
            df.set_index("histone_ids", inplace=True)
            self.df = df

    def check_histones(self, ids, organism=None):
        # check if ids are histones
        share = np.intersect1d(ids, self.df.index)
        result = self.df.loc[share]
        if organism:
            result = result[result["name"] == organism]
        return result

    def get_organism(self, ids):
        # get organism name from histone ids
        result = self.check_histones(ids)
        data = result.groupby(["name", "genome_size"]).size()
        return data.sort_values().index[-1]



