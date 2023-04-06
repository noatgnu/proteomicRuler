import click
from proteomicRuler.ruler import Ruler, add_mw
import pandas as pd

# CLI for proteomicRuler
@click.command()
@click.option("--input", "-i", help="Input file containing intensity of samples and uniprot accession ids", type=click.File("rb"))
@click.option("--output", "-o", help="Output file", type=click.File("wb"))
@click.option("--ploidy", "-p", help="Ploidy of the organism", type=int, default=2)
@click.option("--total-cellular", "-t", help="Total cellular protein concentration", type=float, default=200)
@click.option("--mw-column", "-m", help="Molecular weight column name", type=str, default="Mass")
@click.option("--accession-id-col", "-a", help="Accession id column name", type=str, default="Protein.Ids")
@click.option("--intensity-columns", "-c", help="Intensity columns list delimited by commas", type=str, default=None)
@click.option("--get-mw", "-g", help="Get molecular weight from uniprot", is_flag=True)
def main(input, output, ploidy, total_cellular, get_mw, mw_column="Mass", accession_id_col="Protein.Ids", intensity_columns=None):
    # Read the input file
    df = pd.read_csv(input, sep="\t")
    print(df[pd.isnull(df[accession_id_col])])
    df = df[pd.notnull(df[accession_id_col])]
    # used as unique index and to directly fetch mw data from UniProt
    # molecular weight column name
    # ploidy number
    # total_cellular_protein_concentration = 200
    # intensity_columns = df.columns[5:]
    # accession_id_col = "Protein.Ids"
    if get_mw:
        df = add_mw(df, accession_id_col)
    if intensity_columns is None:
        intensity_columns = []
    else:
        intensity_columns = intensity_columns.split(",")
    ruler = Ruler(df, intensity_columns, mw_column, accession_id_col, ploidy, total_cellular)
    ruler.df.to_csv(output, sep="\t", index=False)