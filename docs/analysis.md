# Analysis guide

## Simulation designs and plans

`config/setting_metadata.json` contains the original 90 cohort/settings. The distinct expanded independent and central cohorts are recorded in `config/primary_cohort_targets.csv`. Together with the additional central and individual-binary confirmations, these comprise the 107 reported cohort/settings. Cohorts that share a generating architecture retain separate identities and seeds.

`config/simulation_evaluation_plan.csv` records the 44,300 canonical evaluation datasets. The original source-phase labels are retained. `config/precision_extension_plan.csv` records the separate added calibration and evaluation phases. `config/point_validation_plan.csv` separates 160 development datasets from 2,000 fresh validation datasets. `config/prospective_validation_plan.csv` describes the additional 4,000 dedicated single-exposure and joint-exposure datasets.

The public entry points use the same generating functions and numerical rules, with portable paths and simplified execution. They do not launch a background queue. A one-dataset run is the default. Requests for all datasets require `--limit=all`.

Some saved source objects include fitting results and metadata as well as generated data. Their original `input_sha256` values authenticate those objects, not a newly serialized RDS containing only regenerated summaries. Regeneration checks should compare the numerical inputs, generating seeds and declared design. A missing original hash remains missing.

## Calibration and inference

`analysis/freeze_calibration.R` accepts calibration rows only, uses the type-8 95th percentile of absolute usable-null score statistics and records hashes of the input and cutoff files. It rejects an existing cutoff file. Added evaluation inputs must not be generated until the intended calibration is completed and frozen.

The reported calibration and historical cutoff versions are retained in `data/reported_metrics/`. Historical single-exposure cutoffs are not used for the updated UVMR-SuSiE fits. Their application and simulation summaries use nominal score tests. No method setting is selected using evaluation truth.

The optional R support-union routine is a separate inference sensitivity. It assumes independent Gaussian exposure/outcome summary errors, homogeneous standard errors, supplied LD and a prespecified upper bound on invalid variants. It retains empty, disconnected and unbounded confidence sets. The routine profiles nuisance coefficients using residualized summaries to avoid subtracting large quadratic forms. Common-scale normalization and stable quadratic roots reduce overflow and cancellation. A nonpositive or numerically singular nuisance quadratic block returns an unbounded projection. Calculations outside the supported floating-point range are rejected or return a conservative unbounded projection. Its stated finite-sample guarantee does not extend to estimated LD or structural logistic outcomes.

## Estimands and metrics

Conditional coefficients, score coverage and conditional-effect MSE/RMSE are evaluated against their generating coefficients. Repeated UVMR-SuSiE and TWAS analyses of multivariable data are marginal. Their rejection summaries remain separate from conditional-effect power, FDR, coverage and RMSE.

Within a setting, estimation error and interval measures pool the eligible returned targets. Settings receive equal weights in cross-setting summaries. Exact paired comparisons use only targets available to both methods. Monte Carlo resampling uses whole datasets and preserves dependence among their coefficients and selected sets.

Detection and null rejection use all scheduled nonzero and null targets, respectively. An unavailable test contributes no rejection. Coefficient-testing FDR is the mean replicate false discovery proportion. Candidate-member FDR and credible-set coverage are separate localization endpoints. An empty candidate selection has false discovery proportion zero. Credible-set coverage is undefined when no role-specific set is returned. Mixed sets remain unresolved.

## Comparators and benchmark

The comparator adapter retains IVW, Egger, Lasso, median, robust, LD-aware IVW/Egger, MV-IWAS, IWAS-Egger projection, MVMR-cML-SuSiE and marginal TWAS summaries. The native MVMR-cML-SuSiE path uses Steps 1–4 with BIC and a top-alpha representative model. No data-perturbation implementation is included. The fit and native-comparator seeds are separate when both are recorded in a plan.

The True causal exposures benchmark supplies the generating exposure support. It retains the joint MR-SuSiE estimated invalid-variant support and reference LD. Nonzero effect sizes are estimated from the data. Its cross-LD nuisance construction retains directions associated with invalid variants outside a clumped panel. It does not receive the generating nonzero effect magnitudes.

Raw fitting objects retain warnings, failed identification, nonconvergence, unavailable estimates, unselected exposures and unresolved credible sets. The full original coefficient records accompany the manuscript as Data S1.

## Runtime comparison

`config/runtime_pairs.csv` identifies the ten recorded paired inputs and their prepared-input hashes. Both methods receive the same summaries and clumped panel. The native comparator uses BIC without data perturbation. Method order alternates across datasets. Each attempt is retained, including a failure or a return with no coefficients. Compilation and package initialization are excluded from the timed section. The supplied recorded results are unchanged. Running the script creates a separate output directory and new machine-specific timings.
