from proteomicRuler.ruler import Ruler, add_mw
import pandas as pd

if __name__ == "__main__":
    need_mw = True
    # used to indicate whether or not to directly fetch mw data from UniProt
    accession_id_col = "Protein IDs"
    # used as unique index and to directly fetch mw data from UniProt
    filename = r"C:\Users\Toan Phung\Downloads\proteinGroups.txt"
    df = pd.read_csv(filename, sep="\t")

    # used for loading the input dataframe
    intensity_columns = df.columns[57:57+16]
    df = df[(df["Only identified by site"]!="+")&(df["Only identified by site"]!="Reverse")&(df["Only identified by site"]!="Potential contaminant")]
    df["cumulative_intensity"] = df[intensity_columns].sum(axis=1)

    for i, r in df.iterrows():
        if r["cumulative_intensity"] != 0:
            for c in intensity_columns:
                df.at[i, c] = r[c]/r["cumulative_intensity"]*r["Intensity"]

    # string of column names containing samples' intensity
    mw_col = "Mass"
    # molecular weight column name
    ploidy = 2
    # ploidy number
    total_cellular_protein_concentration = 200
    # cellular protein concentration used for calculation of total volume

    if need_mw:
        df = add_mw(df, accession_id_col)
    df = df[pd.notnull(df[mw_col])]
    df[mw_col] = df[mw_col].str.replace(",", "")
    df[mw_col] = df[mw_col].astype(float)
    ruler = Ruler(df, intensity_columns, mw_col, accession_id_col, ploidy, total_cellular_protein_concentration)
    ruler.df.to_csv(filename+"output.tsv", sep="\t", index=False)
    ruler.plot()
    print("Finished")
