"""Core data-processing functions for Jahan BioOmics."""

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import KNNImputer


REQUIRED_COLUMNS = (
    "PG.ProteinAccessions",
    "PG.Genes",
    "PG.ProteinDescriptions",
    "PG.Quantity",
    "R.Condition",
    "R.Replicate",
)


def validate_input_table(data: pd.DataFrame) -> pd.DataFrame:
    """Validate the minimum input schema and return a defensive copy."""
    missing = [column for column in REQUIRED_COLUMNS if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    if data.empty:
        raise ValueError("The input table is empty")

    validated = data.copy()
    validated["PG.Quantity"] = pd.to_numeric(validated["PG.Quantity"], errors="coerce")
    if validated["PG.Quantity"].notna().sum() == 0:
        raise ValueError("PG.Quantity does not contain numeric values")
    if (validated["PG.Quantity"].dropna() < 0).any():
        raise ValueError("PG.Quantity contains negative values")
    if validated["R.Condition"].isna().any() or validated["R.Replicate"].isna().any():
        raise ValueError("Condition and replicate identifiers must not be missing")
    if validated["PG.ProteinAccessions"].isna().any():
        raise ValueError("Protein accessions must not be missing")
    return validated


def filter_quality(
    data: pd.DataFrame,
    *,
    max_q_value: float | None = 0.01,
    remove_single_hits: bool = True,
) -> pd.DataFrame:
    """Apply optional run-wise q-value and single-hit filters."""
    filtered = validate_input_table(data)
    if max_q_value is not None and "PG.QValue (Run-Wise)" in filtered.columns:
        q_values = pd.to_numeric(filtered["PG.QValue (Run-Wise)"], errors="coerce")
        filtered = filtered.loc[q_values <= max_q_value].copy()
    if remove_single_hits and "PG.IsSingleHit" in filtered.columns:
        values = filtered["PG.IsSingleHit"]
        if values.dtype == bool:
            keep = ~values
        else:
            keep = ~values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
        filtered = filtered.loc[keep].copy()
    if filtered.empty:
        raise ValueError("No rows remain after quality filtering")
    return filtered


def build_expression_matrix(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str], int]:
    """Create a protein-by-sample matrix using mean aggregation for duplicates."""
    validated = validate_input_table(data)
    validated["SampleName"] = (
        validated["R.Condition"].astype(str)
        + "__replicate_"
        + validated["R.Replicate"].astype(str)
    )
    sample_groups = (
        validated[["SampleName", "R.Condition"]]
        .drop_duplicates()
        .set_index("SampleName")["R.Condition"]
        .astype(str)
        .to_dict()
    )
    index_columns = ["PG.ProteinAccessions", "PG.Genes", "PG.ProteinDescriptions"]
    duplicate_count = int(validated.duplicated(index_columns + ["SampleName"]).sum())
    pivot = validated.pivot_table(
        index=index_columns,
        columns="SampleName",
        values="PG.Quantity",
        aggfunc="mean",
        dropna=False,
    )
    pivot = pivot.dropna(how="all")
    if pivot.empty or pivot.shape[1] < 2:
        raise ValueError("At least one protein and two samples are required")
    metadata = pivot.index.to_frame(index=False)
    return pivot, metadata, sample_groups, duplicate_count


def log2_transform(expression: pd.DataFrame) -> pd.DataFrame:
    """Log2-transform positive quantities and represent zeros as missing."""
    numeric = expression.apply(pd.to_numeric, errors="coerce")
    if (numeric < 0).any().any():
        raise ValueError("Expression matrix contains negative quantities")
    transformed = np.log2(numeric.where(numeric > 0))
    return transformed.replace([np.inf, -np.inf], np.nan)


def impute_expression(
    expression: pd.DataFrame,
    method: str,
    *,
    n_neighbors: int = 3,
) -> pd.DataFrame:
    """Impute a log2 expression matrix with an explicitly selected method."""
    if expression.empty:
        raise ValueError("Expression matrix is empty")
    if expression.isna().all(axis=0).any():
        missing_samples = expression.columns[expression.isna().all(axis=0)].tolist()
        raise ValueError(f"Samples contain no finite values: {', '.join(missing_samples)}")

    normalized_method = method.strip().lower()
    if normalized_method == "knn":
        neighbors = min(n_neighbors, max(1, len(expression) - 1))
        values = KNNImputer(n_neighbors=neighbors).fit_transform(expression)
        return pd.DataFrame(values, index=expression.index, columns=expression.columns)
    if normalized_method == "minimum":
        minimum = float(np.nanmin(expression.to_numpy(dtype=float)))
        return expression.fillna(minimum)
    if normalized_method == "zero":
        return expression.fillna(0.0)
    raise ValueError("method must be one of: knn, minimum, zero")


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    """Adjust finite p-values with the Benjamini-Hochberg procedure."""
    values = np.asarray(list(p_values), dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if finite_indices.size == 0:
        return adjusted
    finite = np.clip(values[finite_indices], 0.0, 1.0)
    order = np.argsort(finite)
    ranked = finite[order]
    scale = finite.size / np.arange(1, finite.size + 1)
    ranked_adjusted = np.minimum.accumulate((ranked * scale)[::-1])[::-1]
    ranked_adjusted = np.clip(ranked_adjusted, 0.0, 1.0)
    restored = np.empty_like(ranked_adjusted)
    restored[order] = ranked_adjusted
    adjusted[finite_indices] = restored
    return adjusted


def differential_test(
    expression: pd.DataFrame,
    group_a_columns: Sequence[str],
    group_b_columns: Sequence[str],
) -> pd.DataFrame:
    """Run protein-wise Welch tests and report raw p-values, FDR, and log2FC.

    Log2 fold change is mean(group B) minus mean(group A). At least two
    replicate columns are required in each group.
    """
    if len(group_a_columns) < 2 or len(group_b_columns) < 2:
        raise ValueError("Welch's t-test requires at least two replicates per group")
    missing = [
        column
        for column in [*group_a_columns, *group_b_columns]
        if column not in expression.columns
    ]
    if missing:
        raise ValueError(f"Unknown sample columns: {', '.join(missing)}")

    data_a = expression.loc[:, group_a_columns].to_numpy(dtype=float)
    data_b = expression.loc[:, group_b_columns].to_numpy(dtype=float)
    variance_a = np.nanvar(data_a, axis=1, ddof=1)
    variance_b = np.nanvar(data_b, axis=1, ddof=1)
    both_constant = (variance_a == 0) & (variance_b == 0)
    p_values = np.ones(data_a.shape[0], dtype=float)
    test_rows = ~both_constant
    with np.errstate(invalid="ignore", divide="ignore"):
        _, p_values[test_rows] = stats.ttest_ind(
            data_b[test_rows],
            data_a[test_rows],
            axis=1,
            equal_var=False,
            nan_policy="omit",
        )
    p_values = np.nan_to_num(p_values, nan=1.0, posinf=1.0, neginf=1.0)
    log2_fold_change = np.nanmean(data_b, axis=1) - np.nanmean(data_a, axis=1)
    return pd.DataFrame(
        {
            "Log2FC": log2_fold_change,
            "P_Value": p_values,
            "FDR": benjamini_hochberg(p_values),
        },
        index=expression.index,
    )


def sample_columns_for_condition(sample_groups: dict[str, str], condition: str) -> list[str]:
    """Return sample columns assigned to a condition."""
    return [sample for sample, group in sample_groups.items() if group == condition]
