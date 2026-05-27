# Examples

## `replogle2022_walkthrough.ipynb`

A complete pipeline walkthrough on the package's built-in synthetic
Perturb-seq fixture (2000 cells, 10 guides, 5 perturbations) designed to
mirror a Replogle 2022 essential-gene CRISPRi screen at small scale.

Covers:

1. Loading 10x H5 + guide call CSVs
2. Guide assignment with the three-rule procedure
3. Per-cell, per-guide, per-perturbation QC
4. Cell-state landscape (UMAP)
5. Mixscape KO vs NP classification
6. Pseudobulk DE per perturbation
7. Cell-state effect map
8. Pathway scoring (decoupler)
9. Rendering the HTML report

### Running it

```bash
pip install -e ".[all,dev]"

jupyter nbconvert --to notebook --execute \
  examples/replogle2022_walkthrough.ipynb \
  --output replogle2022_walkthrough_executed.ipynb
```

Or open it in Jupyter Lab and execute cell-by-cell.

### Adapting it to real Replogle 2022 data

The notebook's first cell points at the `build_synthetic` test fixture.
To swap in the actual public data:

```python
import perturbflow as pf

adata = pf.read_10x_h5('replogle_2022/filtered_feature_bc_matrix.h5')
guide_calls = pf.read_guide_calls('replogle_2022/guide_calls.csv')
guide_metadata = pf.read_guide_metadata('replogle_2022/guide_library.csv')
```

Public data sources:

- **GEO accession:** [GSE177150](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177150)
  (raw and processed)
- **figshare:** [Genome-scale Perturb-seq](https://plus.figshare.com/articles/dataset/Genome-scale_Perturb-seq/20029387)
  (preprocessed AnnDatas)

For the Dixit 2016 original Perturb-seq dataset:

- **GEO accession:** [GSE90063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90063)
