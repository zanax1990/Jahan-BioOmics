🧬 Jahan BioOmics
Python Version Streamlit License Status

Jahan BioOmics is a robust, interactive bioinformatics dashboard for end-to-end analysis of label-free DIA (Data-Independent Acquisition) proteomics data. Built with Python and Streamlit, it streamlines the workflow from raw report processing to biological interpretation.

Note: Optimized for mouse datasets (Mus musculus), with explicit handling for zero-variance edge cases.

🚀 Key Features

🧹 Automated Quality Control
Real-time Filtering: Adjust Q-Value thresholds dynamically.
Data Cleaning: Option to remove single-hit proteins for higher confidence.
Pivot Logic: Automatically transforms long-format reports into expression matrices.

🤖 Imputation
KNN Imputation: Uses k-Nearest Neighbors to fill missing values based on expression patterns.
Fallback Methods: Options for "Min Value" or "Zero Fill" for comparison.
Log2 Transformation: Includes infinity/NaN protection for stable downstream analysis.

📊 Advanced Visualization
Interactive PCA: 2D Principal Component Analysis with dynamic grouping.
Correlation Heatmaps: Visualizes sample-to-sample reproducibility.
Volcano Plots: Interactive scatter plots with adjustable Log2FC and P-value cutoffs.

🧮 Statistics
Welch’s T-Test: Fast, vectorized implementation for group comparisons.
Edge Case Handling: Supports presence/absence scenarios where variance is zero in one or both groups.

🧬 Biological Insight
Pathway Enrichment: Integration with g:Profiler.
Gene Ontology & KEGG: Fetches enriched terms for upregulated vs. downregulated sets.

🛠 Installation

git clone https://github.com/zanax1990/Jahan-BioOmics.git
cd Jahan-BioOmics
pip install -r requirements.txt


Usage

streamlit run app.py


Dashboard opens at: http://localhost:8501

📂 Input Data Format
Designed to parse Spectronaut DIA-style reports (CSV or TSV). Required columns:

PG.ProteinAccessions
PG.Genes
PG.Quantity
R.Condition
R.Replicate
PG.QValue (Run-Wise)
PG.IsSingleHit

Dependencies
streamlit
pandas
numpy
scipy
scikit-learn
plotly
gprofiler-official

👨‍🔬 Author
Jahanbakhsh Ghasemi
Ph.D. Candidate, University of Connecticut (2025)
