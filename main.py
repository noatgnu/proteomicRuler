from proteomicRuler.ruler import Ruler
import pandas as pd
from uniprot.parser import UniprotParser, UniprotSequence
from io import StringIO


def add_mw(df, accession_id_col):
    df1 = df.copy()
    for i, r in df1.iterrows():
        seq = UniprotSequence(r[accession_id_col], True)
        if seq.accession:
            df1.at[i, "Accession"] = str(seq)
    accessions = df1["Accession"].unique()
    parser = UniprotParser(accessions, True)
    data = []
    for i in parser.parse("tab", method="post"):
        frame = pd.read_csv(StringIO(i), sep="\t")
        frame = frame.rename(columns={frame.columns[-1]: "query"})
        data.append(frame)
    data = pd.concat(data, ignore_index=True)
    unmatched = []
    for a in accessions:
        if a not in data["query"].values:
            unmatched.append(a)
    if unmatched:
        print("Non-Uniprot ID found:", unmatched)
    data["Gene names"] = data["Gene names"].fillna("")
    data["gene_name_list"] = data["Gene names"].str.split(" ")
    data.loc[:, "Gene names"] = data["gene_name_list"].map(lambda x: x[0])
    data1 = data[["query",  "Entry name", "Mass", "Gene names"]]
    data1["queries"] = data1["query"].str.split(",")
    data1 = data1.explode("queries")
    res = df1.merge(data1, how="left", left_on="Accession", right_on="queries")
    res = res.drop(columns=["Accession", "queries", "query"])
    res = res.groupby(accession_id_col).head(1)
    return res

if __name__ == "__main__":
    need_mw = True
    # used to indicate whether or not to directly fetch mw data from UniProt
    accession_id_col = "Protein ID"
    # used as unique index and to directly fetch mw data from UniProt
    filename = r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For MS-Fragger\C21orf2_protein.tsv"
    # used for loading the input dataframe
    intensity_columns = ["WT-01", "WT-02", "WT-03", "WT-04", "WT-05", "C21Orf2-01", "C21Orf2-02", "C21Orf2-03", "C21Orf2-04", "C21Orf2-05"]
    # string of column names containing samples' intensity
    mw_col = "Mass"
    # molecular weight column name
    ploidy = 2
    # ploidy number
    total_cellular_protein_concentration = 200
    # cellular protein concentration used for calculation of total volume

    df = pd.read_csv(filename, sep="\t")
    if need_mw:
        df = add_mw(df, accession_id_col)
    df = df[pd.notnull(df[mw_col])]
    df[mw_col] = df[mw_col].str.replace(",", "")
    df[mw_col] = df[mw_col].astype(float)
    ruler = Ruler(df, intensity_columns, mw_col, accession_id_col, ploidy, total_cellular_protein_concentration)
    ruler.df.to_csv(filename+"output.tsv", sep="\t", index=False)
    print("Finished")
