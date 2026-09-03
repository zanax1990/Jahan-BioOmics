import numpy as np
import pandas as pd
import pytest

from bioomics import (
    benjamini_hochberg,
    build_expression_matrix,
    differential_test,
    impute_expression,
    log2_transform,
    validate_input_table,
)


def input_table() -> pd.DataFrame:
    rows = []
    for protein, gene, quantities in [
        ("P1", "GENE1", [4.0, 5.0, 16.0, 20.0]),
        ("P2", "GENE2", [8.0, 8.0, 8.0, 8.0]),
    ]:
        for condition, replicate, quantity in zip(["A", "A", "B", "B"], [1, 2, 1, 2], quantities):
            rows.append(
                {
                    "PG.ProteinAccessions": protein,
                    "PG.Genes": gene,
                    "PG.ProteinDescriptions": f"Synthetic {protein}",
                    "PG.Quantity": quantity,
                    "R.Condition": condition,
                    "R.Replicate": replicate,
                }
            )
    return pd.DataFrame(rows)


def test_validation_reports_missing_columns():
    with pytest.raises(ValueError, match="Missing required columns"):
        validate_input_table(pd.DataFrame({"PG.Quantity": [1.0]}))


def test_validation_rejects_negative_quantity():
    data = input_table()
    data.loc[0, "PG.Quantity"] = -1
    with pytest.raises(ValueError, match="negative"):
        validate_input_table(data)


def test_build_expression_matrix_is_deterministic():
    expression, metadata, groups, duplicates = build_expression_matrix(input_table())
    assert expression.shape == (2, 4)
    assert metadata.shape[0] == 2
    assert set(groups.values()) == {"A", "B"}
    assert duplicates == 0


def test_log2_transform_and_minimum_imputation_are_finite():
    expression = pd.DataFrame({"S1": [1.0, 0.0], "S2": [4.0, 16.0]})
    transformed = log2_transform(expression)
    imputed = impute_expression(transformed, "minimum")
    assert np.isfinite(imputed.to_numpy()).all()
    assert imputed.loc[1, "S1"] == pytest.approx(0.0)


def test_benjamini_hochberg_known_values():
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.002])
    np.testing.assert_allclose(adjusted, [0.02, 0.04, 0.04, 0.008])


def test_differential_test_reports_expected_direction_and_fdr():
    expression, _, _, _ = build_expression_matrix(input_table())
    transformed = log2_transform(expression)
    results = differential_test(
        transformed,
        ["A__replicate_1", "A__replicate_2"],
        ["B__replicate_1", "B__replicate_2"],
    )
    assert results.iloc[0]["Log2FC"] > 0
    assert results.iloc[1]["Log2FC"] == pytest.approx(0.0)
    assert ((results["FDR"] >= 0) & (results["FDR"] <= 1)).all()


def test_differential_test_requires_replicates():
    expression = pd.DataFrame({"A": [1.0], "B": [2.0]})
    with pytest.raises(ValueError, match="at least two replicates"):
        differential_test(expression, ["A"], ["B"])
