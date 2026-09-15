# Analysis code — oldest-old rehabilitation inpatients (non-inferiority + logistic regression)

R code for two analyses of the same convalescent rehabilitation cohort:

1. Non-inferiority of motor FIM in the oldest-old, tested with the
   Brunner–Munzel statistic and its relative effect, with Hodges–Lehmann shifts,
   inverted confidence intervals, Holm adjustment, and an ATS test for
   interaction with a marginal gate.
2. Multivariable logistic regression of good functional outcome, using
   restricted cubic splines (RCS) for non-linearity and
   substantive-model-compatible multiple imputation (`smcfcs`) for missing
   covariates.

The scripts implement the analysis plan.

> No patient data are included in this repository.
> The individual patient records underlying the analysis are held under
> institutional data-governance rules and are not redistributed. Every script
> here reads from a local `data/` directory that is deliberately excluded by
> `.gitignore`. Running the code therefore reproduces the *procedure*, not the
> published numbers.

---

## Repository layout

```
.
├── 01_labels.R … 28_table2_model_specification.R   analysis scripts, numbered
│                                                   in run order (see below)
├── LICENSE                MIT — applies to the code
├── LICENSE-docs.txt       CC BY 4.0 — applies to documentation, figures, tables
├── CITATION.cff
└── (created at run time, not tracked)
    ├── data/              inputs and intermediate .rds / .csv objects
    └── figures/           generated .png / .pdf / .tiff
```

Scripts locate the project root with the `here` package, so `01_labels.R`
must stay at the repository root. Open the folder as an RStudio project (or
create an empty `.here` file) before running anything.

## Requirements

The reported analysis was run on R 4.6.0; R ≥ 4.2 is expected to work.
Packages used:

`Hmisc`, `MASS`, `ResourceSelection`, `brglm2`, `car`, `data.table`,
`detectseparation`, `flextable`, `forcats`, `furrr`, `future`, `ggcorrplot`,
`ggplot2`, `gridExtra`, `here`, `logistf`, `mice`, `officer`, `pROC`,
`patchwork`, `ragg`, `rankFD`, `renv`, `smcfcs`, `tictoc`

```r
install.packages(c(
  "Hmisc", "MASS", "ResourceSelection", "brglm2", "car", "data.table",
  "detectseparation", "flextable", "forcats", "furrr", "future", "ggcorrplot",
  "ggplot2", "gridExtra", "here", "logistf", "mice", "officer", "pROC",
  "patchwork", "ragg", "rankFD", "renv", "smcfcs", "tictoc"
))
```

`06_session_info_license.R` writes the exact package versions used
(`sessionInfo.txt`, `package_versions.csv`) for the record.

---

## Scripts

The numbering is the run order. The three blocks are **common** (01–06),
**non-inferiority** (07–10) and **regression** (11–28). Every script is
self-contained: none of them sources another except for `01_labels.R`.

### 01–06 Common

Shared preparation and cohort description, used by both analyses.

| Script | Plan | Step |
|---|---|---|
| `01_labels.R` | — | English display labels for variables, levels and pipelines; sourced by every later script |
| `02_preprocess.R` | §6 | Entry-error correction, JCS dichotomisation, factor levels and references, outcome dichotomisation, subset construction, missingness typing |
| `03_table1_baseline.R` | §12.2 | Table 1 — baseline characteristics; also builds `BNB_smallplus` |
| `04_excluded_vs_included.R` | §10 | Excluded patients vs the analysed cohort — standardised mean differences and tipping-point analysis |
| `05_person_days.R` | §4.3 | Follow-up person-days and person-years (STROBE 14c) |
| `06_session_info_license.R` | §14 | `sessionInfo()`, package versions, licence files |

### 07–10 Non-inferiority

Brunner–Munzel non-inferiority of the motor FIM, and its interactions.

| Script | Plan | Step |
|---|---|---|
| `07_mfim_violin_scatter.R` | §12.2 | Violin plot of the admission motor FIM and the age × admission motor FIM scatter |
| `08_noninferiority_tests.R` | §8.1–8.5 | Brunner–Munzel non-inferiority tests (36 contrasts × 2 reference definitions), Hodges–Lehmann shifts with inverted confidence intervals, Holm adjustment, stratum QC |
| `09_noninferiority_forest.R` | §12.2 | Forest plots of the non-inferiority results, one per reference definition |
| `10_interaction_ats.R` | §8.6–8.7 | ATS interaction tests, family admission by marginal gate, bootstrap of Δ relative effect and Δ Hodges–Lehmann |

`09` and `10` both read the tables written by `08`; re-running `08` means
re-running both.

### 11–28 Regression

| Script | Plan | Step |
|---|---|---|
| `11_mnar_systems.R` | §7.2 | The two MNAR systems (worst-value, missing-category) |
| `12_centering_corr.R` | §9.2 | Median centring; Pearson / Spearman correlation pre-screening |
| `13_linearity_rcs.R` | §9.3 | Non-linearity tests and fitted curves with restricted cubic splines |
| `14_imputation.R` | §9.4 | Multiple imputation via `smcfcs` (four pipelines) |
| `15_imputation_diagnostics.R` | §9.4 | Convergence traces, observed vs imputed distributions, FMI and MCSE |
| `16_separation_estimation.R` | §9.5 | (Quasi-)complete separation detection; Firth vs ordinary estimation |
| `17_influence_diagnostics.R` | §9.6 | Pearson residuals, leverage, Cook's distance, dfbetas, averaged across imputations |
| `18_outlier_definition.R` | §9.6 | Outliers fixed by criterion 1 **and** criterion 2 |
| `19_variable_selection.R` | §9.7 | Backward elimination of main effects (pooled likelihood-ratio test, EPV ceiling) |
| `20_interaction_screening.R` | §9.7 | Forward screening of two-way interactions |
| `21_final_model.R` | §9.7–9.8 | Final model fit and Rubin pooling |
| `22_or_curves.R` | §9.8 | Adjusted odds-ratio curves for the continuous predictors |
| `23_interaction_plots.R` | §9.8 | Stratified display of the retained interactions |
| `24_model_performance.R` | §9.9 | Discrimination and calibration |
| `25_gvif_collinearity.R` | §9.9 | Generalised variance inflation factors |
| `26_bootstrap_optimism.R` | §9.9 | Bootstrap optimism correction with selection inside each resample; variable-selection frequencies |
| `27_reporting_supplements.R` | §11 | Riley minimum sample size, EPV, reporting supplement tables |
| `28_table2_model_specification.R` | §12.2 | Table 2 — full specification of the final model |

## Running order

```
01 → 02 → 03 → 04 → 05 → 06          common
          └→ 07 → 08 → 09, 10        non-inferiority
          └→ 11 → 12 → 13 → 14 → 15 → 16 → 17 → 18 → 19 → 20 → 21
                 → 22 … 28            regression (22–28 in any order)
```

`05_person_days.R` needs `BNB_smallplus`, which `03_table1_baseline.R` writes.
Within the regression block, `22`–`28` depend on `21` but not on one another.

Output file names follow each script's own numbering, which is NOT the
figure numbering of the manuscript.

## Note on comments

The source comments were removed before publication, so these files carry the
executed logic only. This README and the step tables above are the intended
documentation; the analysis plan carries the statistical rationale.

## Licence

- Code (`*.R`): MIT — see [`LICENSE`](LICENSE)
- Documentation, figures, tables, reporting checklists: CC BY 4.0 — see [`LICENSE-docs.txt`](LICENSE-docs.txt)
- Individual patient data: not covered by either licence and not redistributed

## Citation

See [`CITATION.cff`](CITATION.cff). When citing the code, use the version
DOI shown on the Zenodo record for the release you actually used.

---

## 日本語

回復期リハビリテーション入院患者の同一コホートを対象とした 2 つの解析のコードです。

1. 超高齢者における機能改善の非劣性：Brunner–Munzel 統計量と相対効果による
   非劣性検定、Hodges–Lehmann 推定と反転信頼区間、Holm 調整、および marginal gate
   を用いた ATS による交互作用検定。
2. 良好な機能転帰の多変量ロジスティック回帰：非線形性は制限付き 3 次スプライン
   （RCS）、欠測共変量は substantive-model-compatible な多重代入（`smcfcs`）。

解析計画書の §4〜§15 に対応します。

個票データは含まれていません。 各スクリプトはローカルの `data/` を読みますが、
このディレクトリは `.gitignore` で除外されています。コードは手続きを再現するもので
あり、公表値そのものを再現するものではありません。

スクリプトは番号順に実行してください。番号は 共通（01–06）→ 非劣性（07–10）→
回帰（11–28）の順に付けています。`here` パッケージでプロジェクトルートを解決する
ため、`01_labels.R` はリポジトリ直下に置いたままにしてください。

`05_person_days.R` は `03_table1_baseline.R` が書き出す `BNB_smallplus` を読みます。
`09` と `10` はいずれも `08` の出力表を入力とするため、`08` を再実行した場合は両方を
再実行してください。

コメントは公開前に全削除しています。上記のステップ表がドキュメントの役割を果たし、
統計的な根拠は解析計画書に記載しています。
