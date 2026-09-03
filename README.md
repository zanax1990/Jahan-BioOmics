# Jahan BioOmics

Jahan BioOmics is a Streamlit application for exploratory analysis of label-free DIA proteomics reports. It combines schema validation, quality filtering, missing-value handling, visualization, differential testing, and pathway enrichment in one reproducible workflow.

## Workflow

1. Import a comma- or tab-separated Spectronaut-style report.
2. Validate the required columns and reject invalid negative quantities.
3. Filter by run-wise q-value and optionally remove single-hit rows.
4. Aggregate duplicate protein/sample rows by their mean quantity.
5. Apply a log2 transform and select KNN, global-minimum, or zero imputation.
6. Inspect PCA and sample correlations.
7. Compare two conditions with protein-wise Welch tests and Benjamini-Hochberg FDR correction.
8. Query GO Biological Process and KEGG enrichment through g:Profiler.

## Input data

Required columns:

- `PG.ProteinAccessions`
- `PG.Genes`
- `PG.ProteinDescriptions`
- `PG.Quantity`
- `R.Condition`
- `R.Replicate`

Optional quality-control columns:

- `PG.QValue (Run-Wise)`
- `PG.IsSingleHit`

Research data are not included. `examples/synthetic_proteomics.csv` is a small structural check and does not represent an experiment.

## Installation

```bash
git clone https://github.com/zanax1990/Jahan-BioOmics.git
cd Jahan-BioOmics
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt -r requirements-dev.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Usage

```bash
streamlit run app.py
```

Upload a report through the sidebar. For local development, a default input can be configured without editing the source:

```bash
export BIOOMICS_DEFAULT_DATA=examples/synthetic_proteomics.csv
streamlit run app.py
```

## Statistical details

- Quantities equal to zero are treated as missing before log2 transformation.
- KNN imputation operates on the protein-by-sample matrix with up to three neighboring proteins.
- Duplicate protein/sample observations are aggregated by the arithmetic mean and reported in the interface.
- Differential testing uses Welch's unequal-variance t-test and requires at least two replicate samples in each group.
- Log2 fold change is `mean(group B) - mean(group A)`.
- Benjamini-Hochberg correction is applied across the tested proteins; volcano-plot categories use both the selected FDR and fold-change thresholds.
- PCA is calculated from the imputed log2 matrix without feature standardization.

## Tests

```bash
pytest
```

Tests cover schema validation, negative-value rejection, deterministic matrix construction, transformation and imputation, a known FDR example, fold-change direction, and replicate requirements. GitHub Actions runs the suite on every push and pull request.

## Repository structure

```text
.
├── app.py                       # Streamlit interface
├── bioomics.py                  # validated processing and statistics
├── examples/                    # synthetic input for structural checks
├── tests/                       # unit tests
└── .github/workflows/tests.yml  # continuous integration
```

## Limitations

This is an exploratory research application, not a validated clinical pipeline. Imputation can affect downstream variance and statistical significance; results should be checked against an analysis plan and independent software. The application does not model paired or blocked designs, batch effects, covariates, or protein-level mixed effects. Pathway enrichment requires network access to g:Profiler and depends on the selected organism and identifier mapping.
