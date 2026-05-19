# Project Handoff

## Current Status
**PHASE FIVE COMPLETED - MANUSCRIPT HARDENING & DOCUMENTATION ALIGNMENT PHASE INITIATED.**

The QSAR logKoc modeling pipeline, validation protocol, scientific manuscript, reference structure, and source-code appendix are complete. The current repository state has moved from manuscript compilation into peer-review hardening and front-facing documentation alignment.

## Active Project
- **Project:** QSAR logKoc Modeling Project
- **Dataset:** QDB.177, Gramatica et al. 2014
- **Target property:** logKoc
- **Champion model:** `Hierarchical_MC_MLR`
- **Champion algorithm:** Multiple Linear Regression
- **Champion descriptors:** `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, `SlogP_VSA2`
- **Original benchmark:** Gramatica 2014 Model 4 external `R²_ext = 0.794`
- **Final champion external performance:** `Q²_ext F2 = 0.814252`
- **Benchmark descriptor count:** Eight descriptors
- **Current documentation state:** Manuscript complete, audited, citation-locked, and aligned with the repository README.

## Current Phase
**Manuscript hardening and repository documentation alignment.**

## Next Action Items
- Execute the feature-clustering ablation study that compares the champion pipeline against a baseline workflow without feature clustering.
- Compile absolute error-reduction metrics and statistical significance testing for the performance gain over the Gramatica 2014 benchmark.
- Run a cross-model feature-importance consensus analysis comparing `Hierarchical_MC_MLR`, `Hierarchical_MC_PLS`, and adjacent high-ranking candidates.
- Preserve the historical-dataset defense by framing the project as controlled algorithm-topology benchmarking under the fixed QDB.177 and eight-descriptor constraints.

## Final Deliverables
- `docs/Paper_Draft.md`
- `docs/PROJECT.md`
- `docs/HANDOFF.md`
- `README.md`
- `notebooks/01_data_split.ipynb`
- `notebooks/02_mordred_extract.ipynb`
- `notebooks/03_pyqsar3_modeling.ipynb`
- `notebooks/04_result_metrics.ipynb`
- `run_nonlinear_models.py`
- `run_mccv.py`
- `run_yscrambling.py`
- `data/processed/train.csv`
- `data/processed/test.csv`
- `data/features/mordred_train.csv`
- `data/features/mordred_test.csv`
- `data/features/filtered_train_pyqsar3.csv`
- `data/features/filtered_test_pyqsar3.csv`
- `data/features/feature_clusters_pyqsar3.json`
- `data/features/feature_clusters_pyqsar3_hierarchical.cluster`
- `data/features/feature_clusters_pyqsar3_kmeans.cluster`
- `data/features/feature_clusters_pyqsar3_som.cluster`
- `data/features/12_model_extended_metrics.csv`
- `data/features/final_selected_features_8vars.json`
- `data/features/test_predictions_8vars.csv`
- `data/features/nonlinear_model_metrics.csv`
- `data/features/best_nonlinear_config.json`
- `data/features/nonlinear_search_log.txt`
- `data/features/mccv_100_iterations.csv`
- `data/features/mccv_summary.json`
- `data/features/yscrambling_100_iterations.csv`
- `data/features/yscrambling_summary.json`
- `data/features/figure_predicted_vs_experimental.png`
- `data/features/figure_williams_plot.png`
- `data/features/final_comprehensive_metrics.csv`
- `data/features/benchmark_oecd_comparison.csv`

## Completed Milestones
- Parsed the QDB.177 `M2.logKoc` endpoint and mapped compound identifiers to Daylight SMILES records.
- Curated molecular structures using RDKit validation, canonicalization, salt/disconnected-structure exclusion, and duplicate removal.
- Secured 642 unique valid compounds from 643 raw QDB endpoint entries.
- Generated the Y-ranking 80/20 data split with 514 training compounds and 128 external test compounds.
- Confirmed zero canonical SMILES overlap between the training and external test sets.
- Calculated raw two-dimensional Mordred descriptors independently for training and external test compounds.
- Preserved the full 1,613-descriptor raw feature space during Phase Two without missing-value, zero-variance, or correlation filtering.
- Completed native PyQSAR3-compatible preprocessing using training-set-only filtering decisions.
- Reduced the modeling descriptor matrix to 377 robust molecular descriptors plus `SMILES` and `logKoc` columns.
- Corrected the clustering architecture from sample-wise partitioning to descriptor-wise feature clustering.
- Generated hierarchical, K-Means, and SOM descriptor-cluster mappings using native PyQSAR3 feature-clustering modules.
- Saved feature-cluster JSON and native `.cluster` sidecar files for all three descriptor-clustering tracks.
- Executed the strict eight-descriptor 12-model linear baseline matrix across clustering, feature-selection, and regression tracks.
- Calculated comprehensive internal and external validation metrics for all 12 linear baseline models.
- Generated and executed the server-side non-linear modeling workflow for four hierarchical-track SVR and RF candidates.
- Integrated all 16 candidate models under the strict `Q²_cv` champion-selection policy.
- Selected `Hierarchical_MC_MLR` as the final champion model with `Q²_cv = 0.8119744761107354`.
- Locked the final eight-descriptor champion feature set: `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, and `SlogP_VSA2`.
- Completed 100-iteration MCCV for the locked champion model, yielding average `R²_ext = 0.806404` with SD = 0.032097.
- Completed 100-iteration Y-scrambling for the locked champion model, yielding average `R²_y-sc = 0.016196` with SD = 0.007508 and maximum scrambled `R² = 0.039415`.
- Generated the predicted-versus-experimental diagnostic plot for training and external test series.
- Generated the Williams Plot and verified applicability-domain compliance using standardized residual limits of `±3` and warning leverage threshold `h* = 0.0525`.
- Compiled the final comprehensive 16-model metric table with MCCV and Y-scrambling fields in `final_comprehensive_metrics.csv`.
- Generated the benchmark and OECD validation comparison file, confirming that `Hierarchical_MC_MLR` exceeded the Gramatica 2014 external benchmark under the identical eight-descriptor constraint.
- Completed the final manuscript draft with Introduction, Materials and Methods, Results and Discussion, Conclusion, References, and Appendix source-code integration.
- Completed the formal manuscript audit against the master metric tables, benchmark comparison file, MCCV summary, and Y-scrambling summary.
- Locked the 15-citation ACS-style reference section with software, dataset, OECD, and descriptor-origin citations.
- Added peer-review hardening recommendations for low Silhouette score defense, statistical significance of the performance delta, descriptor-consensus analysis, and historical dataset framing.
- Re-authored the root repository README to align the public-facing project description with the finalized manuscript.

## Peer-Review Hardening Tracker
- **Low Silhouette score defense:** Pending ablation study comparing the clustered pipeline against a non-clustered baseline.
- **Pipeline complexity justification:** Pending RMSE and MAE delta compilation with statistical significance testing.
- **Descriptor interpretation defense:** Pending cross-model feature-consensus analysis across top-ranked models.
- **Historical dataset framing:** Active repository and manuscript framing now identifies QDB.177 as a controlled benchmarking vehicle rather than a large-scale data-mining corpus.

## Final Scientific Position
The final `Hierarchical_MC_MLR` model outperformed the Gramatica et al. (2014) Model 4 historical baseline in external predictivity while preserving the identical eight-descriptor constraint. The project established a defined endpoint, an unambiguous algorithm, a formally evaluated applicability domain, rigorous goodness-of-fit and validation statistics, and chance-correlation rejection by Y-scrambling. The final workflow satisfies the five OECD principles for QSAR validation and now includes a documented hardening roadmap for peer-review defense.
