Proteomic Ruler
--

An implementation of the same algorithm from Perseus `Wiśniewski, J. R., Hein, M. Y., Cox, J. and Mann, M. (2014) A “Proteomic Ruler” for Protein Copy Number and Concentration Estimation without Spike-in Standards. Mol Cell Proteomics 13, 3497–3506.` used for estimation of protein copy number from deep profile experiment.

Requirements
--
Python >= 3.9

Installation
--
```bash
pip install proteomicruler
```

Usage
--

In order to use the package, it is required that the input data is loaded into a `pandas.DataFrame` object. The following
basic parameters are also required:
- `accession_id_col` - column name that contains protein accession ids
- `mw_col` - column name that contains molecular weight of proteins
- `ploidy` - ploidy number
- `total_cellular_protein_concentration` - total cellular protein concentration used for calculation of total volume
- `intensity_columns` - list of column names that contain sample intensities

```python
import pandas as pd

accession_id_col = "Protein IDs"
# used as unique index and to directly fetch mw data from UniProt

mw_col = "Mass"
# molecular weight column name

ploidy = 2
# ploidy number

total_cellular_protein_concentration = 200
# cellular protein concentration used for calculation of total volume

filename = r"example_data\example_data.tsv" # example data from Perseus
df = pd.read_csv(filename, sep="\t")

# selecting intensity columns
intensity_columns = df.columns[57:57+16] # select 16 columns starting from column 57th that contain sample intensity



```

If the data does not contain molecular weight information, it is required to fetch it from UniProt.

```python
from proteomicRuler.ruler import add_mw

df = add_mw(df, accession_id_col)
df = df[pd.notnull(df[mw_col])]
df[mw_col] = df[mw_col].str.replace(",", "")
df[mw_col] = df[mw_col].astype(float)
```

The RuleR object can be created by passing the `DataFrame` object and the required parameters.

```python
from proteomicRuler.ruler import Ruler

ruler = Ruler(df, intensity_columns, mw_col, accession_id_col, ploidy, total_cellular_protein_concentration) #
ruler.df.to_csv("output.txt", sep="\t", index=False)
```
