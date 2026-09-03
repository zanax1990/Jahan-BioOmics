import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from gprofiler import GProfiler
from sklearn.decomposition import PCA

from bioomics import (
    build_expression_matrix,
    differential_test,
    filter_quality,
    impute_expression,
    log2_transform,
    sample_columns_for_condition,
    validate_input_table,
)


st.set_page_config(
    page_title="Jahan BioOmics",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def read_report(source) -> pd.DataFrame:
    """Read a comma- or tab-separated report from a path or upload."""
    data = pd.read_csv(source)
    if data.shape[1] < 2:
        if hasattr(source, "seek"):
            source.seek(0)
        data = pd.read_csv(source, sep="\t")
    return data


@st.cache_data(ttl=3600)
def run_pathway_enrichment(gene_list: list[str], organism: str = "mmusculus"):
    """Query GO Biological Process and KEGG terms through g:Profiler."""
    if not gene_list:
        return None
    profiler = GProfiler(return_dataframe=True)
    results = profiler.profile(
        organism=organism,
        query=gene_list,
        sources=["GO:BP", "KEGG"],
    )
    if results.empty:
        return None
    columns = ["source", "native", "name", "p_value", "intersection_size"]
    return results.loc[:, columns].sort_values("p_value").head(20)


st.sidebar.title("Jahan BioOmics")
st.sidebar.caption("Exploratory DIA proteomics analysis")
st.sidebar.subheader("1. Data import")
uploaded_file = st.sidebar.file_uploader("Upload report", type=["csv", "txt", "tsv"])
default_path = os.getenv("BIOOMICS_DEFAULT_DATA")

st.title("Jahan BioOmics")
st.markdown("Interactive quality control, exploratory analysis, differential testing, and pathway enrichment")

data = None
try:
    if uploaded_file is not None:
        data = read_report(uploaded_file)
    elif default_path:
        data = read_report(Path(default_path))
        st.success("Loaded the dataset configured by BIOOMICS_DEFAULT_DATA.")
    else:
        st.info("Upload a report to begin, or set BIOOMICS_DEFAULT_DATA for local development.")
except (OSError, ValueError, pd.errors.ParserError) as error:
    st.error(f"Could not read the report: {error}")

if data is not None:
    try:
        data = validate_input_table(data)
    except ValueError as error:
        st.error(str(error))
        st.stop()

    st.sidebar.subheader("2. Quality control")
    q_threshold = None
    if "PG.QValue (Run-Wise)" in data.columns:
        q_threshold = st.sidebar.slider("Maximum run-wise q-value", 0.0, 0.05, 0.01, 0.001)
    remove_single_hits = False
    if "PG.IsSingleHit" in data.columns:
        remove_single_hits = st.sidebar.checkbox("Remove single hits", value=True)

    try:
        filtered = filter_quality(
            data,
            max_q_value=q_threshold,
            remove_single_hits=remove_single_hits,
        )
        expression, protein_metadata, sample_groups, duplicate_count = build_expression_matrix(filtered)
    except ValueError as error:
        st.error(str(error))
        st.stop()

    transformed = log2_transform(expression)
    imputation_label = st.sidebar.radio("Missing-value method", ["KNN", "Minimum", "Zero"])
    imputation_method = {"KNN": "knn", "Minimum": "minimum", "Zero": "zero"}[imputation_label]
    try:
        imputed = impute_expression(transformed, imputation_method)
    except ValueError as error:
        st.error(str(error))
        st.stop()

    overview_tab, imputation_tab, pca_tab, differential_tab, pathway_tab = st.tabs(
        ["QC & overview", "Imputation", "PCA & correlation", "Differential analysis", "Pathway enrichment"]
    )

    with overview_tab:
        column_a, column_b, column_c = st.columns(3)
        column_a.metric("Input rows", len(data))
        column_b.metric("Rows after QC", len(filtered))
        column_c.metric("Samples", expression.shape[1])
        if duplicate_count:
            st.warning(f"Mean aggregation combined {duplicate_count} duplicate protein/sample row(s).")
        st.dataframe(expression.head(), use_container_width=True)

    with imputation_tab:
        figure = go.Figure()
        figure.add_trace(
            go.Histogram(x=transformed.to_numpy().ravel(), name="Log2 input", opacity=0.7, marker_color="gray")
        )
        figure.add_trace(
            go.Histogram(x=imputed.to_numpy().ravel(), name="Imputed", opacity=0.5, marker_color="steelblue")
        )
        figure.update_layout(barmode="overlay", title="Distribution before and after imputation")
        st.plotly_chart(figure, use_container_width=True)

    with pca_tab:
        if imputed.shape[0] < 2 or imputed.shape[1] < 2:
            st.warning("PCA requires at least two proteins and two samples.")
        else:
            pca = PCA(n_components=2)
            components = pca.fit_transform(imputed.T)
            pca_data = pd.DataFrame(components, columns=["PC1", "PC2"], index=imputed.columns)
            pca_data["Sample"] = pca_data.index
            pca_data["Group"] = [sample_groups[sample] for sample in pca_data.index]
            pca_figure = px.scatter(
                pca_data,
                x="PC1",
                y="PC2",
                color="Group",
                text="Sample",
                title=(
                    f"PCA: PC1 {pca.explained_variance_ratio_[0] * 100:.1f}%, "
                    f"PC2 {pca.explained_variance_ratio_[1] * 100:.1f}%"
                ),
            )
            pca_figure.update_traces(textposition="top center")
            st.plotly_chart(pca_figure, use_container_width=True)
            correlation = px.imshow(
                imputed.corr(),
                text_auto=".2f",
                aspect="auto",
                color_continuous_scale="RdBu_r",
                zmin=-1,
                zmax=1,
                title="Sample correlation",
            )
            st.plotly_chart(correlation, use_container_width=True)

    with differential_tab:
        conditions = sorted(set(sample_groups.values()))
        if len(conditions) < 2:
            st.warning("Differential analysis requires at least two conditions.")
        else:
            selector_a, selector_b = st.columns(2)
            condition_a = selector_a.selectbox("Group A", conditions, index=0)
            condition_b = selector_b.selectbox("Group B", conditions, index=1)
            group_a = sample_columns_for_condition(sample_groups, condition_a)
            group_b = sample_columns_for_condition(sample_groups, condition_b)
            try:
                statistics = differential_test(imputed, group_a, group_b).reset_index(drop=True)
            except ValueError as error:
                st.warning(str(error))
            else:
                results = pd.concat([protein_metadata.reset_index(drop=True), statistics], axis=1)
                results["NegLogFDR"] = -np.log10(results["FDR"].clip(lower=np.finfo(float).tiny))
                fold_change_cutoff = st.slider("Absolute log2 fold-change cutoff", 0.0, 4.0, 1.0, 0.1)
                fdr_cutoff = st.slider("FDR cutoff", 0.0, 0.1, 0.05, 0.001)
                categories = [
                    (results["FDR"] < fdr_cutoff) & (results["Log2FC"] > fold_change_cutoff),
                    (results["FDR"] < fdr_cutoff) & (results["Log2FC"] < -fold_change_cutoff),
                ]
                results["Category"] = np.select(
                    categories,
                    ["Upregulated", "Downregulated"],
                    default="Not significant",
                )
                st.session_state["stats_result"] = results

                volcano = px.scatter(
                    results,
                    x="Log2FC",
                    y="NegLogFDR",
                    color="Category",
                    hover_data=["PG.Genes", "PG.ProteinDescriptions", "P_Value", "FDR"],
                    color_discrete_map={
                        "Upregulated": "#B23A48",
                        "Downregulated": "#3B6EA8",
                        "Not significant": "#B8B8B8",
                    },
                    title=f"Volcano plot: {condition_b} relative to {condition_a}",
                    opacity=0.8,
                )
                volcano.add_hline(y=-np.log10(fdr_cutoff), line_dash="dash")
                volcano.add_vline(x=fold_change_cutoff, line_dash="dash")
                volcano.add_vline(x=-fold_change_cutoff, line_dash="dash")
                st.plotly_chart(volcano, use_container_width=True)

                significant = results.loc[results["Category"] != "Not significant"].sort_values("FDR")
                st.write(f"{len(significant)} proteins meet both thresholds.")
                st.dataframe(
                    significant[["PG.Genes", "Log2FC", "P_Value", "FDR", "Category"]],
                    use_container_width=True,
                )
                st.download_button(
                    "Download differential results",
                    data=results.to_csv(index=False).encode("utf-8"),
                    file_name="differential_results.csv",
                    mime="text/csv",
                )

    with pathway_tab:
        st.caption("Enrichment uses proteins meeting the selected fold-change and FDR thresholds.")
        if "stats_result" not in st.session_state:
            st.info("Run differential analysis first.")
        else:
            results = st.session_state["stats_result"]
            for category in ["Upregulated", "Downregulated"]:
                st.subheader(category)
                genes = (
                    results.loc[results["Category"] == category, "PG.Genes"]
                    .dropna()
                    .astype(str)
                    .drop_duplicates()
                    .tolist()
                )
                if len(genes) <= 5:
                    st.info("At least six genes are required for this enrichment query.")
                    continue
                try:
                    pathways = run_pathway_enrichment(genes)
                except Exception as error:
                    st.warning(f"The g:Profiler query failed: {error}")
                    continue
                if pathways is None:
                    st.info("No enriched pathways were returned.")
                    continue
                pathways = pathways.assign(neg_log_p=-np.log10(pathways["p_value"]))
                chart = px.bar(
                    pathways,
                    x="neg_log_p",
                    y="name",
                    orientation="h",
                    color="source",
                    title="Top enrichment terms",
                    height=420,
                )
                chart.update_layout(yaxis={"categoryorder": "total ascending"})
                st.plotly_chart(chart, use_container_width=True)
