# Method interface

`make_mr_susie()` constructs the common model. Its `fit()` interface accepts an m-by-M exposure-association matrix `B`, the corresponding standard-error matrix `sX`, outcome associations `y`, outcome standard errors `sy`, signed LD `R`, and the M-by-M exposure-error correlation matrix `Sigma_exp`.

The selection-outcome associations, standard errors and sample-size parameter are supplied separately. The caller explicitly declares whether the selection and estimation outcomes reuse data. The default interface assumes zero exposure–outcome error covariance. Supplying an exposure-error correlation matrix does not model that overlap.

For continuous-outcome summaries, use `selector_family="gaussian"`. The application runner uses the binary score adapter and the supplied effective sample size. A missing sample-size value is passed as NULL, preserving the large-sample score convention. The implementation requires susieR 0.14.2 and uses the fixed localization and point-estimation controls in the analysis.

The output separates:

- `point`: covariance-aware posterior means.
- `score_inference`: corrected-score centers, standard errors, interval endpoints and p-values.
- `selection`: exposure, invalid-variant and unresolved credible-set memberships.
- `family_tests`: raw and multiplicity-adjusted p-values with unavailable tests retained.
- `warnings`, `status` and `raw_fit`: diagnostic and complete fitting records.

UVMR-SuSiE uses M=1 in the same engine. Repeated single-exposure fits to multivariable data estimate marginal working associations. Their estimates are not conditional effects adjusted for the remaining exposures.

The isotropic point prior depends on coefficient units. A separate rescaling of one exposure can alter the point map. Units should therefore match the analysis and be reported with the estimates.

Score intervals condition on localization and working covariance inputs. Empirical undercoverage is reported in the manuscript. Point shrinkage does not alter the score centers, intervals or tests. Unavailable coefficients are not filled with zeros. Zeros supplied by the True causal exposures benchmark are labeled Given and have no empirical test or interval.
