#!/usr/bin/env python3
"""
Phase 3 - Step 3.4: Server-side non-linear QSAR modeling.

Run from the project root in the `qsar_ml` conda environment:

    conda activate qsar_ml
    python run_nonlinear_models.py

The script evaluates four Hierarchical-track non-linear models:
SVR_GA, SVR_MC, RF_GA, and RF_MC. Each search is constrained to exactly
8 descriptors and every evaluated subset is scored by an inner 5-fold
GridSearchCV. Required outputs are written to `data/features/`.
"""

from __future__ import annotations

import json
import logging
import math
import os
import random
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
from deap import base, creator, tools
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR


# ------------------------------- Configuration -------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent
FEATURE_DIR = PROJECT_ROOT / "data" / "features"

TRAIN_PATH = FEATURE_DIR / "filtered_train_pyqsar3.csv"
TEST_PATH = FEATURE_DIR / "filtered_test_pyqsar3.csv"
CLUSTER_JSON_PATH = FEATURE_DIR / "feature_clusters_pyqsar3.json"

LOG_PATH = FEATURE_DIR / "nonlinear_search_log.txt"
BEST_CONFIG_PATH = FEATURE_DIR / "best_nonlinear_config.json"
METRICS_PATH = FEATURE_DIR / "nonlinear_model_metrics.csv"

ID_COL = "SMILES"
TARGET_COL = "logKoc"
CLUSTER_TRACK = "Hierarchical"
N_FEATURES = 8
RANDOM_STATE = 42
CV_SPLITS = 5

# Production defaults. Increase GA_N_GENERATIONS or MC_N_ITERATIONS on the
# server if wall time allows.
GA_POPULATION_SIZE = 48
GA_N_GENERATIONS = 30
GA_CXPB = 0.65
GA_MUTPB = 0.35
GA_TOURNAMENT_SIZE = 3
MC_N_ITERATIONS = 1000

SVR_PARAM_GRID = {
    "model__C": [1, 10, 100],
    "model__gamma": ["scale", 0.01, 0.1],
}
RF_PARAM_GRID = {
    "model__n_estimators": [100, 300],
    "model__max_depth": [None, 10, 20],
}


# ---------------------------------- Logging ----------------------------------

def configure_logger() -> logging.Logger:
    FEATURE_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("qsar_nonlinear")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(LOG_PATH, mode="a", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


LOGGER = configure_logger()


# -------------------------------- Data records --------------------------------

@dataclass(frozen=True)
class SearchResult:
    selector: str
    model_type: str
    features: Tuple[str, ...]
    best_cv_score: float
    best_params: Dict[str, Any]


@dataclass(frozen=True)
class MetricRow:
    model: str
    selector: str
    model_type: str
    n_features: int
    selected_features: str
    best_params: str
    R2_train: float
    CCC_tr: float
    RMSE_train: float
    Q2_cv: float
    CCC_cv: float
    RMSE_cv: float
    MAE_cv: float
    Q2_ext_F1: float
    Q2_ext_F2: float
    Q2_ext_F3: float
    CCC_ext: float
    RMSE_ext: float
    MAE_ext: float


# --------------------------------- Metrics -----------------------------------

def as_array(values: Sequence[float]) -> np.ndarray:
    return np.asarray(values, dtype=float).reshape(-1)


def rmse(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    return float(math.sqrt(mean_squared_error(y_true, y_pred)))


def ccc(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """Concordance Correlation Coefficient."""
    x = as_array(y_true)
    y = as_array(y_pred)
    if len(x) != len(y):
        raise ValueError("CCC requires equal-length arrays.")
    if len(x) < 2:
        return float("nan")

    sx = np.std(x, ddof=0)
    sy = np.std(y, ddof=0)
    if sx == 0 or sy == 0:
        return 0.0

    rho = np.corrcoef(x, y)[0, 1]
    numerator = 2.0 * rho * sx * sy
    denominator = sx**2 + sy**2 + (np.mean(x) - np.mean(y)) ** 2
    return float(numerator / denominator) if denominator != 0 else 0.0


def q2_ext_f1(
    y_train: Sequence[float],
    y_ext: Sequence[float],
    y_ext_pred: Sequence[float],
) -> float:
    y_train = as_array(y_train)
    y_ext = as_array(y_ext)
    y_ext_pred = as_array(y_ext_pred)
    numerator = np.sum((y_ext - y_ext_pred) ** 2)
    denominator = np.sum((y_ext - np.mean(y_train)) ** 2)
    return float(1.0 - numerator / denominator) if denominator != 0 else float("nan")


def q2_ext_f2(y_ext: Sequence[float], y_ext_pred: Sequence[float]) -> float:
    return float(r2_score(y_ext, y_ext_pred))


def q2_ext_f3(
    y_train: Sequence[float],
    y_ext: Sequence[float],
    y_ext_pred: Sequence[float],
) -> float:
    y_train = as_array(y_train)
    y_ext = as_array(y_ext)
    y_ext_pred = as_array(y_ext_pred)
    numerator = np.sum((y_ext - y_ext_pred) ** 2) / len(y_ext)
    denominator = np.sum((y_train - np.mean(y_train)) ** 2) / len(y_train)
    return float(1.0 - numerator / denominator) if denominator != 0 else float("nan")


# ------------------------------ Input handling -------------------------------

def normalize_hierarchical_mapping(payload: Mapping[str, Any]) -> Dict[str, List[str]]:
    """
    Return cluster_id -> descriptors for the Hierarchical feature clustering.

    The current project JSON stores mappings as descriptor -> cluster_id. This
    parser also supports regenerated JSON stored as cluster_id -> descriptor list.
    """
    try:
        raw_mapping = payload["mappings"][CLUSTER_TRACK]
    except KeyError as exc:
        raise KeyError(f"Missing mappings['{CLUSTER_TRACK}'] in cluster JSON.") from exc

    if not isinstance(raw_mapping, Mapping):
        raise TypeError("Hierarchical mapping must be a JSON object.")

    values = list(raw_mapping.values())
    cluster_to_features: Dict[str, List[str]] = {}
    if values and all(isinstance(value, list) for value in values):
        for cluster_id, features in raw_mapping.items():
            cluster_to_features[str(cluster_id)] = [str(feature) for feature in features]
    else:
        for feature, cluster_id in raw_mapping.items():
            cluster_to_features.setdefault(str(cluster_id), []).append(str(feature))

    cluster_to_features = {
        cluster_id: sorted(set(features))
        for cluster_id, features in cluster_to_features.items()
        if features
    }
    if len(cluster_to_features) < N_FEATURES:
        raise ValueError(
            f"Need at least {N_FEATURES} non-empty clusters; "
            f"found {len(cluster_to_features)}."
        )

    return dict(
        sorted(
            cluster_to_features.items(),
            key=lambda item: (0, int(item[0])) if item[0].isdigit() else (1, item[0]),
        )
    )


def load_data() -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, List[str]]]:
    for path in [TRAIN_PATH, TEST_PATH, CLUSTER_JSON_PATH]:
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    train_df = pd.read_csv(TRAIN_PATH)
    test_df = pd.read_csv(TEST_PATH)
    for name, df in [("train", train_df), ("test", test_df)]:
        missing = {ID_COL, TARGET_COL}.difference(df.columns)
        if missing:
            raise ValueError(f"{name} data missing required columns: {sorted(missing)}")

    with CLUSTER_JSON_PATH.open("r", encoding="utf-8") as handle:
        cluster_to_features = normalize_hierarchical_mapping(json.load(handle))

    train_features = set(train_df.columns).difference({ID_COL, TARGET_COL})
    mapped_features = {
        feature for features in cluster_to_features.values() for feature in features
    }
    missing_from_train = sorted(mapped_features.difference(train_features))
    missing_from_mapping = sorted(train_features.difference(mapped_features))
    missing_from_test = sorted(train_features.difference(test_df.columns))
    if missing_from_train:
        raise ValueError(f"Cluster features absent from train data: {missing_from_train[:20]}")
    if missing_from_mapping:
        raise ValueError(f"Train descriptors absent from cluster mapping: {missing_from_mapping[:20]}")
    if missing_from_test:
        raise ValueError(f"Train descriptors absent from test data: {missing_from_test[:20]}")

    return train_df, test_df, cluster_to_features


# ------------------------- Cluster-aware subset logic -------------------------

def all_feature_indices(
    cluster_to_features: Mapping[str, Sequence[str]],
) -> Tuple[List[str], Dict[int, str]]:
    features: List[str] = []
    index_to_cluster: Dict[int, str] = {}
    for cluster_id, cluster_features in cluster_to_features.items():
        for feature in cluster_features:
            index_to_cluster[len(features)] = str(cluster_id)
            features.append(feature)
    return features, index_to_cluster


def random_cluster_aware_subset(
    rng: random.Random,
    cluster_to_features: Mapping[str, Sequence[str]],
) -> Tuple[str, ...]:
    """Sample exactly 8 descriptors, preferentially one per cluster."""
    cluster_ids = list(cluster_to_features)
    rng.shuffle(cluster_ids)
    selected = [
        rng.choice(list(cluster_to_features[cluster_id]))
        for cluster_id in cluster_ids[:N_FEATURES]
    ]

    all_features = sorted({f for features in cluster_to_features.values() for f in features})
    while len(selected) < N_FEATURES:
        candidate = rng.choice(all_features)
        if candidate not in selected:
            selected.append(candidate)
    return tuple(sorted(selected[:N_FEATURES]))


def repair_indices(
    indices: Iterable[int],
    rng: random.Random,
    features: Sequence[str],
    index_to_cluster: Mapping[int, str],
) -> List[int]:
    """Repair a GA chromosome to exactly 8 unique feature indices."""
    seen = set()
    unique: List[int] = []
    for idx in indices:
        idx = int(idx) % len(features)
        if idx not in seen:
            seen.add(idx)
            unique.append(idx)

    if len(unique) > N_FEATURES:
        rng.shuffle(unique)
        retained: List[int] = []
        retained_clusters = set()
        deferred: List[int] = []
        for idx in unique:
            cluster_id = index_to_cluster[idx]
            if cluster_id not in retained_clusters and len(retained) < N_FEATURES:
                retained.append(idx)
                retained_clusters.add(cluster_id)
            else:
                deferred.append(idx)
        unique = (retained + deferred)[:N_FEATURES]

    while len(unique) < N_FEATURES:
        current_clusters = {index_to_cluster[idx] for idx in unique}
        pool = [
            idx for idx in range(len(features))
            if idx not in seen and index_to_cluster[idx] not in current_clusters
        ]
        if not pool:
            pool = [idx for idx in range(len(features)) if idx not in seen]
        candidate = rng.choice(pool)
        seen.add(candidate)
        unique.append(candidate)

    rng.shuffle(unique)
    return unique[:N_FEATURES]


# ------------------------------ Model evaluation -----------------------------

def make_estimator_and_grid(model_type: str) -> Tuple[Pipeline, Dict[str, List[Any]]]:
    if model_type == "SVR":
        return (
            Pipeline([("scale", StandardScaler()), ("model", SVR(kernel="rbf"))]),
            SVR_PARAM_GRID,
        )
    if model_type == "RF":
        return (
            Pipeline([
                ("scale", StandardScaler()),
                ("model", RandomForestRegressor(random_state=RANDOM_STATE, n_jobs=1)),
            ]),
            RF_PARAM_GRID,
        )
    raise ValueError(f"Unsupported model type: {model_type}")


def grid_score_subset(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    features: Sequence[str],
    model_type: str,
    n_jobs: int,
) -> Tuple[float, Dict[str, Any]]:
    estimator, param_grid = make_estimator_and_grid(model_type)
    cv = KFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    grid = GridSearchCV(
        estimator=estimator,
        param_grid=param_grid,
        scoring="r2",
        cv=cv,
        n_jobs=n_jobs,
        refit=True,
        error_score="raise",
    )
    grid.fit(X_train.loc[:, list(features)], y_train)
    return float(grid.best_score_), dict(grid.best_params_)


def evaluate_subset(
    selector: str,
    model_type: str,
    features: Sequence[str],
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_jobs: int,
) -> SearchResult:
    features = tuple(sorted(features))
    score, params = grid_score_subset(X_train, y_train, features, model_type, n_jobs)
    return SearchResult(selector, model_type, features, score, params)


def log_progress(result: SearchResult, iteration: str) -> None:
    LOGGER.info(
        "[Model: %s_%s] [Iter/Gen: %s] Evaluated Features: %s -> "
        "Best Q2_cv: %.4f -> Params: %s",
        result.model_type,
        result.selector,
        iteration,
        list(result.features),
        result.best_cv_score,
        result.best_params,
    )


# ------------------------------- GA selection --------------------------------

def run_ga(
    model_type: str,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    cluster_to_features: Mapping[str, Sequence[str]],
    n_workers: int,
) -> SearchResult:
    selector = "GA"
    rng = random.Random(RANDOM_STATE + (11 if model_type == "SVR" else 17))
    features, index_to_cluster = all_feature_indices(cluster_to_features)
    feature_to_idx = {feature: idx for idx, feature in enumerate(features)}

    fitness_name = f"FitnessMax_{model_type}_{selector}"
    individual_name = f"Individual_{model_type}_{selector}"
    if not hasattr(creator, fitness_name):
        creator.create(fitness_name, base.Fitness, weights=(1.0,))
    if not hasattr(creator, individual_name):
        creator.create(individual_name, list, fitness=getattr(creator, fitness_name))

    toolbox = base.Toolbox()

    def init_individual() -> Any:
        subset = random_cluster_aware_subset(rng, cluster_to_features)
        indices = [feature_to_idx[feature] for feature in subset]
        return getattr(creator, individual_name)(
            repair_indices(indices, rng, features, index_to_cluster)
        )

    def evaluate_individual(individual: List[int]) -> Tuple[float]:
        individual[:] = repair_indices(individual, rng, features, index_to_cluster)
        subset = tuple(sorted(features[idx] for idx in individual))
        result = evaluate_subset(selector, model_type, subset, X_train, y_train, n_workers)
        return (result.best_cv_score,)

    def mate(ind1: List[int], ind2: List[int]) -> Tuple[List[int], List[int]]:
        tools.cxTwoPoint(ind1, ind2)
        ind1[:] = repair_indices(ind1, rng, features, index_to_cluster)
        ind2[:] = repair_indices(ind2, rng, features, index_to_cluster)
        return ind1, ind2

    def mutate(individual: List[int]) -> Tuple[List[int]]:
        for _ in range(rng.randint(1, 2)):
            pos = rng.randrange(len(individual))
            current_clusters = {index_to_cluster[idx] for idx in individual}
            pool = [
                idx for idx in range(len(features))
                if idx not in individual and index_to_cluster[idx] not in current_clusters
            ]
            if not pool:
                pool = [idx for idx in range(len(features)) if idx not in individual]
            individual[pos] = rng.choice(pool)
        individual[:] = repair_indices(individual, rng, features, index_to_cluster)
        return (individual,)

    toolbox.register("individual", init_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", evaluate_individual)
    toolbox.register("mate", mate)
    toolbox.register("mutate", mutate)
    toolbox.register("select", tools.selTournament, tournsize=GA_TOURNAMENT_SIZE)

    LOGGER.info(
        "Starting %s_GA: population=%d, generations=%d, inner GridSearch jobs=%d",
        model_type,
        GA_POPULATION_SIZE,
        GA_N_GENERATIONS,
        n_workers,
    )

    population = toolbox.population(n=GA_POPULATION_SIZE)
    hall = tools.HallOfFame(1)

    for individual in population:
        individual.fitness.values = toolbox.evaluate(individual)
    hall.update(population)

    for generation in range(1, GA_N_GENERATIONS + 1):
        offspring = list(map(toolbox.clone, toolbox.select(population, len(population))))
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if rng.random() < GA_CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        for mutant in offspring:
            if rng.random() < GA_MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        invalid = [individual for individual in offspring if not individual.fitness.valid]
        for individual in invalid:
            individual.fitness.values = toolbox.evaluate(individual)

        population[:] = offspring
        hall.update(population)
        best_subset = tuple(sorted(features[idx] for idx in hall[0]))
        best_result = evaluate_subset(selector, model_type, best_subset, X_train, y_train, n_workers)
        log_progress(best_result, str(generation))

    best_subset = tuple(sorted(features[idx] for idx in hall[0]))
    return evaluate_subset(selector, model_type, best_subset, X_train, y_train, n_workers)


# ------------------------------- MC selection --------------------------------

def run_mc(
    model_type: str,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    cluster_to_features: Mapping[str, Sequence[str]],
    n_workers: int,
) -> SearchResult:
    selector = "MC"
    rng = random.Random(RANDOM_STATE + (101 if model_type == "SVR" else 107))
    subsets = [
        random_cluster_aware_subset(rng, cluster_to_features)
        for _ in range(MC_N_ITERATIONS)
    ]
    subsets = list(dict.fromkeys(tuple(sorted(subset)) for subset in subsets))

    LOGGER.info(
        "Starting %s_MC: subsets=%d, outer workers=%d, inner GridSearch jobs=1",
        model_type,
        len(subsets),
        n_workers,
    )
    results = Parallel(n_jobs=n_workers, backend="loky")(
        delayed(evaluate_subset)(selector, model_type, subset, X_train, y_train, 1)
        for subset in subsets
    )

    best = results[0]
    for i, result in enumerate(results, start=1):
        if result.best_cv_score > best.best_cv_score:
            best = result
        log_progress(best, f"{i}/{len(results)}")
    return best


# ----------------------------- Final evaluation ------------------------------

def final_metrics(
    result: SearchResult,
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
) -> MetricRow:
    model_name = f"{CLUSTER_TRACK}_{result.selector}_{result.model_type}"
    features = list(result.features)
    X_train = train_df.loc[:, features]
    y_train = train_df[TARGET_COL]
    X_test = test_df.loc[:, features]
    y_test = test_df[TARGET_COL]

    estimator, _ = make_estimator_and_grid(result.model_type)
    estimator.set_params(**result.best_params)
    estimator.fit(X_train, y_train)

    train_pred = estimator.predict(X_train).reshape(-1)
    test_pred = estimator.predict(X_test).reshape(-1)

    cv = KFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    cv_estimator, _ = make_estimator_and_grid(result.model_type)
    cv_estimator.set_params(**result.best_params)
    cv_pred = cross_val_predict(cv_estimator, X_train, y_train, cv=cv, n_jobs=1).reshape(-1)

    return MetricRow(
        model=model_name,
        selector=result.selector,
        model_type=result.model_type,
        n_features=len(features),
        selected_features=json.dumps(features),
        best_params=json.dumps(result.best_params, sort_keys=True),
        R2_train=float(r2_score(y_train, train_pred)),
        CCC_tr=ccc(y_train, train_pred),
        RMSE_train=rmse(y_train, train_pred),
        Q2_cv=float(r2_score(y_train, cv_pred)),
        CCC_cv=ccc(y_train, cv_pred),
        RMSE_cv=rmse(y_train, cv_pred),
        MAE_cv=float(mean_absolute_error(y_train, cv_pred)),
        Q2_ext_F1=q2_ext_f1(y_train, y_test, test_pred),
        Q2_ext_F2=q2_ext_f2(y_test, test_pred),
        Q2_ext_F3=q2_ext_f3(y_train, y_test, test_pred),
        CCC_ext=ccc(y_test, test_pred),
        RMSE_ext=rmse(y_test, test_pred),
        MAE_ext=float(mean_absolute_error(y_test, test_pred)),
    )


def save_outputs(rows: Sequence[MetricRow]) -> None:
    metrics_df = pd.DataFrame([asdict(row) for row in rows])
    metrics_df = metrics_df.sort_values("Q2_ext_F2", ascending=False).reset_index(drop=True)
    metrics_df.to_csv(METRICS_PATH, index=False)

    best = metrics_df.iloc[0].to_dict()
    best_config = {
        "selection_basis": "highest Q2_ext_F2 among Step 3.4 nonlinear candidates",
        "model": best["model"],
        "algorithm_type": best["model_type"],
        "selector": best["selector"],
        "n_features": int(best["n_features"]),
        "features": json.loads(best["selected_features"]),
        "hyperparameters": json.loads(best["best_params"]),
        "metrics": {
            key: float(best[key])
            for key in [
                "R2_train",
                "CCC_tr",
                "RMSE_train",
                "Q2_cv",
                "CCC_cv",
                "RMSE_cv",
                "MAE_cv",
                "Q2_ext_F1",
                "Q2_ext_F2",
                "Q2_ext_F3",
                "CCC_ext",
                "RMSE_ext",
                "MAE_ext",
            ]
        },
    }
    with BEST_CONFIG_PATH.open("w", encoding="utf-8") as handle:
        json.dump(best_config, handle, indent=2)

    LOGGER.info("Saved metrics: %s", METRICS_PATH)
    LOGGER.info("Saved best config: %s", BEST_CONFIG_PATH)


# ---------------------------------- Main -------------------------------------

def main() -> None:
    start = time.time()
    n_workers = max(1, (os.cpu_count() or 1) - 2)
    random.seed(RANDOM_STATE)
    np.random.seed(RANDOM_STATE)

    LOGGER.info("=" * 80)
    LOGGER.info("Starting Phase 3 Step 3.4 nonlinear modeling")
    LOGGER.info("Workers: %d", n_workers)
    LOGGER.info("Cluster track: %s", CLUSTER_TRACK)
    LOGGER.info("Feature constraint: exactly %d descriptors", N_FEATURES)

    train_df, test_df, cluster_to_features = load_data()
    feature_cols = [col for col in train_df.columns if col not in {ID_COL, TARGET_COL}]
    X_train = train_df.loc[:, feature_cols]
    y_train = train_df[TARGET_COL]

    LOGGER.info("Train shape: %s", train_df.shape)
    LOGGER.info("Test shape: %s", test_df.shape)
    LOGGER.info(
        "Hierarchical clusters: %d; mapped descriptors: %d",
        len(cluster_to_features),
        sum(len(features) for features in cluster_to_features.values()),
    )

    search_results = [
        run_ga("SVR", X_train, y_train, cluster_to_features, n_workers),
        run_mc("SVR", X_train, y_train, cluster_to_features, n_workers),
        run_ga("RF", X_train, y_train, cluster_to_features, n_workers),
        run_mc("RF", X_train, y_train, cluster_to_features, n_workers),
    ]

    metric_rows: List[MetricRow] = []
    for result in search_results:
        row = final_metrics(result, train_df, test_df)
        metric_rows.append(row)
        LOGGER.info(
            "[Final: %s] Q2_ext_F2=%.4f Q2_cv=%.4f RMSE_ext=%.4f Features=%s",
            row.model,
            row.Q2_ext_F2,
            row.Q2_cv,
            row.RMSE_ext,
            json.loads(row.selected_features),
        )

    save_outputs(metric_rows)
    LOGGER.info("Completed in %.2f minutes", (time.time() - start) / 60.0)
    LOGGER.info("=" * 80)


if __name__ == "__main__":
    main()
