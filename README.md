# Jahan BioOmics

Jahan BioOmics is a Streamlit application for exploratory analysis of label-free DIA proteomics reports. It brings data import, quality filtering, missing-value handling, visualization, differential analysis, and pathway enrichment into one interactive workflow.

## Workflow

1. Import a Spectronaut-style CSV or TSV report.
2. Filter by run-wise q-value and optionally remove single-hit entries.
3. Pivot protein quantities into a protein-by-sample matrix.
4. Apply a log2 transform and choose KNN, minimum-value, or zero imputation.
5. inspect PCA and sample-to-sample correlations.
6. Compare two conditions with Welch's t-test and an interactive volcano plot.
7. Query GO Biological Process and KEGG enrichment through g:Profiler.

## Input data

The application expects these columns:

- `PG.ProteinAccessions`
- `PG.Genes`
- `PG.ProteinDescriptions`
- `PG.Quantity`
- `R.Condition`
- `R.Replicate`
- `PG.QValue (Run-Wise)` for q-value filtering
- `PG.IsSingleHit` for optional single-hit filtering

No research dataset is included in this repository.

## Installation

```bash
git clone https://github.com/zanax1990/Jahan-BioOmics.git
cd Jahan-BioOmics
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Usage

```bash
streamlit run app.py
```

Upload a report through the sidebar. For local development, an optional default dataset can be provided without editing the source:

```bash
export BIOOMICS_DEFAULT_DATA=/path/to/report.csv
streamlit run app.py
```

## Methods

- KNN imputation uses three nearest neighbors.
- Differential testing uses Welch's unequal-variance t-test.
- Fold change is calculated as the difference between group means after log2 transformation.
- Enrichment queries are separated into upregulated and downregulated gene lists.

## Limitations

This is an exploratory analysis tool, not a validated clinical pipeline. The current differential analysis does not apply multiple-testing correction, and PCA is run on the imputed matrix without feature standardization. Input schema validation and automated tests are not yet included. Pathway enrichment requires network access to g:Profiler.

## Repository structure

```text
.
├── app.py
├── requirements.txt
└── README.md
```
