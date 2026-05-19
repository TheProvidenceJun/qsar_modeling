#!/usr/bin/env python3
"""
Phase 3 - Step 3.5: Integrated Best Model Selection & MCCV.

This script re-inspects the 16 evaluated production models, selects the
Champion strictly by the highest Q2_cv, locks its 8 descriptors and estimator,
and performs 100 Monte Carlo Cross-Validation splits on the merged 642-compound
filtered dataset.

Run from the project root in the `qsar_ml` conda environment:

    conda activate qsar_ml
    python run_mccv.py
"""

from __future__ import annotations

import ast
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR


PROJECT_ROOT = Path(__file__).resolve().parent
FEATURE_DIR = PROJECT_ROOT / "data" / "features"

LINEAR_METRICS_PATH = FEATURE_DIR / "12_model_extended_metrics.csv"
NONLINEAR_METRICS_PATH = FEATURE_DIR / "nonlinear_model_metrics.csv"
LINEAR_FEATURES_PATH = FEATURE_DIR / "final_selected_features_8vars.json"
BEST_NONLINEAR_CONFIG_PATH = FEATURE_DIR / "best_nonlinear_config.json"
FILTERED_TRAIN_PATH = FEATURE_DIR / "filtered_train_pyqsar3.csv"
FILTERED_TEST_PATH = FEATURE_DIR / "filtered_test_pyqsar3.csv"

MCCV_RAW_PATH = FEATURE_DIR / "mccv_100_iterations.csv"
MCCV_SUMMARY_PATH = FEATURE_DIR / "mccv_summary.json"

ID_COL = "SMILES"
TARGET_COL = "logKoc"
N_MCCV_ITERATIONS = 100
TEST_SIZE = 0.20
BASE_RANDOM_STATE = 20260519


@dataclass(frozen=True)
class Champion:
    model_id: str
    algorithm: str
    source: str
    q2_cv: float
    features: List[str]
    hyperparameters: Dict[str, Any]


def rmse(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


def parse_feature_list(value: Any) -> List[str]:
    """Parse feature lists stored as JSON strings, Python literals, or lists."""
    if isinstance(value, list):
        return [str(item) for item in value]
    if not isinstance(value, str):
        raise TypeError(f"Cannot parse feature list from {type(value).__name__}")

    value = value.strip()
    for parser in (json.loads, ast.literal_eval):
        try:
            parsed = parser(value)
            if isinstance(parsed, list):
                return [str(item) for item in parsed]
        except Exception:
            continue
    raise ValueError(f"Could not parse feature list: {value[:120]}")


def parse_params(value: Any) -> Dict[str, Any]:
    """Parse hyperparameters stored as JSON strings, Python literals, or dicts."""
    if isinstance(value, dict):
        return dict(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    if not isinstance(value, str):
        return {}

    value = value.strip()
    if not value:
        return {}
    for parser in (json.loads, ast.literal_eval):
        try:
            parsed = parser(value)
            if isinstance(parsed, dict):
                return dict(parsed)
        except Exception:
            continue
    return {}


def normalize_param_keys(params: Mapping[str, Any]) -> Dict[str, Any]:
    """Remove Pipeline prefixes such as `model__` for direct estimators."""
    normalized: Dict[str, Any] = {}
    for key, value in params.items():
        if key.startswith("model__"):
            normalized[key.split("__", 1)[1]] = value
        else:
            normalized[key] = value
    return normalized


def load_linear_features() -> Dict[str, List[str]]:
    payload = json.loads(LINEAR_FEATURES_PATH.read_text(encoding="utf-8"))
    if "all_model_selected_features" in payload:
        return {
            model_id: [str(feature) for feature in features]
            for model_id, features in payload["all_model_selected_features"].items()
        }
    if "selected_features" in payload:
        return {
            model_id: [str(feature) for feature in features]
            for model_id, features in payload["selected_features"].items()
        }
    raise KeyError("Linear feature JSON lacks selected feature mappings.")


def collect_candidate_models() -> List[Dict[str, Any]]:
    linear_df = pd.read_csv(LINEAR_METRICS_PATH)
    nonlinear_df = pd.read_csv(NONLINEAR_METRICS_PATH)
    linear_features = load_linear_features()

    candidates: List[Dict[str, Any]] = []

    for _, row in linear_df.iterrows():
        model_id = str(row["Model_ID"])
        regressor = str(row["Regressor"])
        if model_id not in linear_features:
            raise KeyError(f"Missing selected features for linear model {model_id}")
        candidates.append(
            {
                "model_id": model_id,
                "algorithm": regressor,
                "source": "linear",
                "Q2_cv": float(row["Q2_cv"]),
                "features": linear_features[model_id],
                "hyperparameters": {},
            }
        )

    for _, row in nonlinear_df.iterrows():
        model_id = str(row["model"])
        candidates.append(
            {
                "model_id": model_id,
                "algorithm": str(row["model_type"]),
                "source": "nonlinear",
                "Q2_cv": float(row["Q2_cv"]),
                "features": parse_feature_list(row["selected_features"]),
                "hyperparameters": parse_params(row.get("best_params")),
            }
        )

    if len(candidates) != 16:
        raise RuntimeError(f"Expected 16 candidate models, found {len(candidates)}")
    return candidates


def identify_champion() -> Champion:
    """Select the Champion strictly by highest Q2_cv."""
    candidates = collect_candidate_models()
    winner = max(candidates, key=lambda item: item["Q2_cv"])

    # `best_nonlinear_config.json` is read to satisfy the audit trail and to
    # verify availability, but it is not used for Champion selection because it
    # was generated by external-score ranking during Step 3.4.
    if BEST_NONLINEAR_CONFIG_PATH.exists():
        json.loads(BEST_NONLINEAR_CONFIG_PATH.read_text(encoding="utf-8"))

    features = [str(feature) for feature in winner["features"]]
    if len(features) != 8:
        raise ValueError(f"Champion must have exactly 8 features; found {len(features)}")

    return Champion(
        model_id=str(winner["model_id"]),
        algorithm=str(winner["algorithm"]),
        source=str(winner["source"]),
        q2_cv=float(winner["Q2_cv"]),
        features=features,
        hyperparameters=dict(winner["hyperparameters"]),
    )


def make_model(champion: Champion) -> Any:
    """Instantiate the locked Champion estimator."""
    algorithm = champion.algorithm.upper()
    params = normalize_param_keys(champion.hyperparameters)

    if algorithm == "MLR":
        return LinearRegression(**params)
    if algorithm == "PLS":
        params.setdefault("n_components", min(3, len(champion.features)))
        return PLSRegression(**params)
    if algorithm == "SVR":
        params.setdefault("kernel", "rbf")
        return Pipeline([
            ("scale", StandardScaler()),
            ("model", SVR(**params)),
        ])
    if algorithm == "RF":
        params.setdefault("random_state", BASE_RANDOM_STATE)
        params.setdefault("n_jobs", 1)
        return Pipeline([
            ("scale", StandardScaler()),
            ("model", RandomForestRegressor(**params)),
        ])
    raise ValueError(f"Unsupported champion algorithm: {champion.algorithm}")


def load_full_dataset(champion: Champion) -> Tuple[pd.DataFrame, pd.Series]:
    train_df = pd.read_csv(FILTERED_TRAIN_PATH)
    test_df = pd.read_csv(FILTERED_TEST_PATH)
    full_df = pd.concat([train_df, test_df], axis=0, ignore_index=True)

    missing = [feature for feature in champion.features if feature not in full_df.columns]
    if missing:
        raise KeyError(f"Champion features absent from filtered dataset: {missing}")
    if TARGET_COL not in full_df.columns:
        raise KeyError(f"Missing target column: {TARGET_COL}")

    X = full_df.loc[:, champion.features].copy()
    y = full_df[TARGET_COL].copy()
    if len(X) != 642:
        raise ValueError(f"Expected merged dataset N=642, found N={len(X)}")
    return X, y


def run_one_iteration(
    iteration: int,
    X: pd.DataFrame,
    y: pd.Series,
    champion: Champion,
) -> Dict[str, Any]:
    random_state = BASE_RANDOM_STATE + iteration
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=TEST_SIZE,
        random_state=random_state,
        shuffle=True,
    )

    model = make_model(champion)
    model.fit(X_train, y_train)
    y_pred = np.asarray(model.predict(X_test)).reshape(-1)

    return {
        "iteration": iteration,
        "random_state": random_state,
        "n_train": int(len(X_train)),
        "n_test": int(len(X_test)),
        "R2_ext": float(r2_score(y_test, y_pred)),
        "RMSE_ext": rmse(y_test, y_pred),
        "MAE_ext": float(mean_absolute_error(y_test, y_pred)),
    }


def summarize_results(champion: Champion, results_df: pd.DataFrame) -> Dict[str, Any]:
    metric_cols = ["R2_ext", "RMSE_ext", "MAE_ext"]
    summary = {
        "champion_model": champion.model_id,
        "algorithm": champion.algorithm,
        "source": champion.source,
        "selection_basis": "highest Q2_cv across 12 linear and 4 nonlinear evaluated models",
        "Q2_cv": champion.q2_cv,
        "n_features": len(champion.features),
        "features": champion.features,
        "hyperparameters": champion.hyperparameters,
        "mccv": {
            "n_iterations": int(len(results_df)),
            "test_size": TEST_SIZE,
            "base_random_state": BASE_RANDOM_STATE,
            "metrics": {},
        },
    }

    for col in metric_cols:
        summary["mccv"]["metrics"][col] = {
            "mean": float(results_df[col].mean()),
            "sd": float(results_df[col].std(ddof=1)),
            "min": float(results_df[col].min()),
            "max": float(results_df[col].max()),
        }
    return summary


def main() -> None:
    n_workers = max(1, (os.cpu_count() or 1) - 2)
    champion = identify_champion()
    X, y = load_full_dataset(champion)

    print("Champion selected strictly by Q2_cv")
    print(f"  Model: {champion.model_id}")
    print(f"  Algorithm: {champion.algorithm}")
    print(f"  Q2_cv: {champion.q2_cv:.12f}")
    print(f"  Features: {champion.features}")
    print(f"Running {N_MCCV_ITERATIONS} MCCV iterations with {n_workers} workers...")

    rows = Parallel(n_jobs=n_workers, backend="loky")(
        delayed(run_one_iteration)(iteration, X, y, champion)
        for iteration in range(1, N_MCCV_ITERATIONS + 1)
    )

    results_df = pd.DataFrame(rows).sort_values("iteration").reset_index(drop=True)
    results_df.to_csv(MCCV_RAW_PATH, index=False)

    summary = summarize_results(champion, results_df)
    MCCV_SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Saved raw MCCV iterations: {MCCV_RAW_PATH}")
    print(f"Saved MCCV summary: {MCCV_SUMMARY_PATH}")
    print(json.dumps(summary["mccv"]["metrics"], indent=2))


if __name__ == "__main__":
    main()
