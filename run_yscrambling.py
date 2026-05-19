#!/usr/bin/env python3
"""
Phase 3 - Step 3.6: Y-Scrambling chance-correlation test.

This script tests whether the locked Champion Model (`Hierarchical_MC_MLR`)
depends on a real descriptor-response relationship rather than chance
correlation. It uses only the training set, repeatedly permutes `logKoc`,
refits the same MLR model, and records the scrambled R2 distribution.

Run from the project root in the `qsar_ml` conda environment:

    conda activate qsar_ml
    python run_yscrambling.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


PROJECT_ROOT = Path(__file__).resolve().parent
FEATURE_DIR = PROJECT_ROOT / "data" / "features"

TRAIN_PATH = FEATURE_DIR / "filtered_train_pyqsar3.csv"
RAW_OUTPUT_PATH = FEATURE_DIR / "yscrambling_100_iterations.csv"
SUMMARY_OUTPUT_PATH = FEATURE_DIR / "yscrambling_summary.json"

TARGET_COL = "logKoc"
CHAMPION_MODEL = "Hierarchical_MC_MLR"
CHAMPION_ALGORITHM = "MLR"
CHAMPION_FEATURES = [
    "ABC",
    "BCUTs-1h",
    "C1SP2",
    "ETA_shape_y",
    "FilterItLogS",
    "NdS",
    "SlogP_VSA1",
    "SlogP_VSA2",
]

N_SCRAMBLING_ITERATIONS = 100
BASE_RANDOM_STATE = 20260519


def load_training_data() -> tuple[pd.DataFrame, pd.Series]:
    """Load the filtered training set and return locked descriptors and target."""
    if not TRAIN_PATH.exists():
        raise FileNotFoundError(f"Required training file not found: {TRAIN_PATH}")

    df = pd.read_csv(TRAIN_PATH)
    missing_features = [feature for feature in CHAMPION_FEATURES if feature not in df.columns]
    if missing_features:
        raise KeyError(f"Missing champion descriptor columns: {missing_features}")
    if TARGET_COL not in df.columns:
        raise KeyError(f"Missing target column: {TARGET_COL}")

    X = df.loc[:, CHAMPION_FEATURES].copy()
    y = df[TARGET_COL].copy()
    return X, y


def fit_original_model(X: pd.DataFrame, y: pd.Series) -> float:
    """Fit the unscrambled MLR model and return original training R2."""
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    return float(r2_score(y, y_pred))


def run_one_scramble(
    iteration: int,
    X: pd.DataFrame,
    y: pd.Series,
) -> Dict[str, float | int]:
    """Permute y, refit MLR, and return scrambled training R2."""
    rng = np.random.default_rng(BASE_RANDOM_STATE + iteration)
    y_scrambled = rng.permutation(np.asarray(y, dtype=float))

    model = LinearRegression()
    model.fit(X, y_scrambled)
    scrambled_pred = model.predict(X)
    scrambled_r2 = float(r2_score(y_scrambled, scrambled_pred))

    return {
        "Iteration": iteration,
        "Random_State": BASE_RANDOM_STATE + iteration,
        "Scrambled_R2": scrambled_r2,
    }


def summarize_results(
    original_r2_train: float,
    results_df: pd.DataFrame,
) -> Dict[str, object]:
    """Create JSON-serializable summary metrics and conclusion."""
    avg_r2_ysc = float(results_df["Scrambled_R2"].mean())
    sd_r2_ysc = float(results_df["Scrambled_R2"].std(ddof=1))
    max_r2_ysc = float(results_df["Scrambled_R2"].max())
    min_r2_ysc = float(results_df["Scrambled_R2"].min())

    conclusion = (
        "Chance correlation rejected: the original model R2 is substantially "
        "higher than the average scrambled R2."
        if original_r2_train > avg_r2_ysc
        else "Warning: scrambled models approach or exceed the original model R2."
    )

    return {
        "champion_model": CHAMPION_MODEL,
        "algorithm": CHAMPION_ALGORITHM,
        "n_features": len(CHAMPION_FEATURES),
        "features": CHAMPION_FEATURES,
        "n_scrambling_iterations": int(len(results_df)),
        "base_random_state": BASE_RANDOM_STATE,
        "original_R2_train": original_r2_train,
        "average_R2_y_sc": avg_r2_ysc,
        "sd_R2_y_sc": sd_r2_ysc,
        "min_R2_y_sc": min_r2_ysc,
        "max_R2_y_sc": max_r2_ysc,
        "conclusion": conclusion,
    }


def print_summary(summary: Dict[str, object]) -> None:
    print("=" * 80)
    print("Y-Scrambling Chance-Correlation Test")
    print(f"Champion Model: {summary['champion_model']}")
    print(f"Algorithm: {summary['algorithm']}")
    print(f"Features: {summary['features']}")
    print(f"Original R2_train: {summary['original_R2_train']:.6f}")
    print(
        "Average R2_y-sc: "
        f"{summary['average_R2_y_sc']:.6f} +/- {summary['sd_R2_y_sc']:.6f}"
    )
    print(f"Min R2_y-sc: {summary['min_R2_y_sc']:.6f}")
    print(f"Max R2_y-sc: {summary['max_R2_y_sc']:.6f}")
    print(f"Conclusion: {summary['conclusion']}")
    print(f"Raw iterations saved to: {RAW_OUTPUT_PATH}")
    print(f"Summary saved to: {SUMMARY_OUTPUT_PATH}")
    print("=" * 80)


def main() -> None:
    X, y = load_training_data()
    original_r2_train = fit_original_model(X, y)

    rows: List[Dict[str, float | int]] = [
        run_one_scramble(iteration, X, y)
        for iteration in range(1, N_SCRAMBLING_ITERATIONS + 1)
    ]
    results_df = pd.DataFrame(rows)
    results_df.to_csv(RAW_OUTPUT_PATH, index=False)

    summary = summarize_results(original_r2_train, results_df)
    SUMMARY_OUTPUT_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print_summary(summary)


if __name__ == "__main__":
    main()
