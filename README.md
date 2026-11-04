# MR-SuSiE analysis and reproducibility

MR-SuSiE jointly localizes candidate causal exposures and direct variant effects. UVMR-SuSiE is the single-exposure case. MVMR-SuSiE applies the same model jointly to multiple exposures.

This directory is a standalone analysis companion to the MR-SuSiE software. It contains the common implementation, synthetic worked examples, simulation generators and plans, application scripts, comparator adapters, metric summaries and the Table S1 workbook. It can be added as a separate analysis directory in the existing repository. It does not replace that repository's installed package namespace.

## Requirements

The analyses were checked with R 4.5.1. A C++ compiler and R development tools are needed for the numerical helpers.

### Third-party R packages

| Package | Checked version | Purpose | License |
|---|---|---|---|
| Rcpp | 1.1.2 | C++ numerical helpers | GPL (>= 2) |
| RcppArmadillo | 15.4.2-1 | Matrix computations in C++ | GPL (>= 2) |
| susieR | 0.14.2 | Single-effect updates | BSD_3_clause + file LICENSE |
| MASS | 7.3-66 | Linear algebra | GPL-2 \| GPL-3 |
| jsonlite | 2.0.0 | JSON inputs and outputs | MIT + file LICENSE |
| digest | 0.6.39 | Input and calibration integrity | GPL (>= 2) |
| data.table | 1.18.4 | Comparator dependencies | MPL-2.0 \| file LICENSE |
| readxl | 1.5.0 | Workbook export | MIT + file LICENSE |
| MendelianRandomization | 0.10.0 | Retained MVMR estimators | GPL-2 \| GPL-3 |
| MRcML | 0.0.0.9000 | Native cML screening | MIT + file LICENSE |
| MVMRcML | 0.1.0 | Native cML BIC estimation | MIT + file LICENSE |
| MVMRcMLSuSiE | 0.1.0 | Native cML-SuSiE localization | GPL (>= 3) |

The worked examples require Rcpp, RcppArmadillo, susieR, MASS, jsonlite and digest. The comparator and workbook packages are needed for their respective optional analyses.

To install the core versions:

```r
install.packages("remotes")
required <- c(
  Rcpp = "1.1.2",
  RcppArmadillo = "15.4.2-1",
  susieR = "0.14.2",
  MASS = "7.3-66",
  jsonlite = "2.0.0",
  digest = "0.6.39"
)
for (package in names(required)) {
  remotes::install_version(package, version = required[[package]], upgrade = "never")
}
```

For the optional comparators:

```r
remotes::install_version("MendelianRandomization", version = "0.10.0", upgrade = "never")
remotes::install_github("xue-hr/MRcML@58a3481976093175cbb82a5e6bf06c10e3c576f2", upgrade = "never")
remotes::install_github("ZhaotongL/MVMR-cML@20d09a54637e70a35f31a43d56c1e6276c99bb17", upgrade = "never")
remotes::install_github("lapsumchan/MVMR-cML-SuSiE@da71c5cd62bda8b0e89f8732ca6013585fb9584b", upgrade = "never")
```

Workbook export uses readxl 1.5.0. The optional support-union sensitivity uses the chi-square functions supplied with base R.

<details>
<summary>Complete recorded R package versions, including indirect dependencies</summary>

| Package | Version |
|---|---|
| arrangements | 1.1.10 |
| askpass | 1.2.1 |
| base64enc | 0.1-6 |
| bslib | 0.11.0 |
| cachem | 1.1.0 |
| cellranger | 1.1.0 |
| cli | 3.6.6 |
| codetools | 0.2-20 |
| cpp11 | 0.5.5 |
| crayon | 1.5.3 |
| crosstalk | 1.2.2 |
| curl | 7.1.0 |
| data.table | 1.18.4 |
| DEoptimR | 1.2-0 |
| digest | 0.6.39 |
| dplyr | 1.2.1 |
| evaluate | 1.0.5 |
| farver | 2.1.2 |
| fastmap | 1.2.0 |
| fontawesome | 0.5.3 |
| foreach | 1.5.2 |
| fs | 2.1.0 |
| generics | 0.1.4 |
| ggplot2 | 4.0.3 |
| glmnet | 5.0 |
| glue | 1.8.1 |
| gmp | 0.7-5.1 |
| gtable | 0.3.6 |
| highr | 0.12 |
| hms | 1.1.4 |
| htmltools | 0.5.9 |
| htmlwidgets | 1.6.4 |
| httr | 1.4.8 |
| irlba | 2.3.7 |
| isoband | 0.3.0 |
| iterators | 1.0.14 |
| iterpc | 0.4.2 |
| jquerylib | 0.1.4 |
| jsonlite | 2.0.0 |
| knitr | 1.51 |
| labeling | 0.4.3 |
| later | 1.4.8 |
| lattice | 0.22-9 |
| lazyeval | 0.2.3 |
| lifecycle | 1.0.5 |
| magrittr | 2.0.5 |
| MASS | 7.3-66 |
| Matrix | 1.7-6 |
| MatrixModels | 0.5-4 |
| matrixStats | 1.5.0 |
| memoise | 2.0.1 |
| MendelianRandomization | 0.10.0 |
| mime | 0.13 |
| mixsqp | 0.3-54 |
| MRcML | 0.0.0.9000 |
| MVMRcML | 0.1.0 |
| MVMRcMLSuSiE | 0.1.0 |
| numDeriv | 2016.8-1.1 |
| openssl | 2.4.2 |
| otel | 0.2.0 |
| pillar | 1.11.1 |
| pkgconfig | 2.0.3 |
| plotly | 4.12.1 |
| plyr | 1.8.9 |
| prettyunits | 1.2.0 |
| progress | 1.2.3 |
| promises | 1.5.0 |
| purrr | 1.2.2 |
| quantreg | 6.1 |
| R6 | 2.6.1 |
| rappdirs | 0.3.4 |
| RColorBrewer | 1.1-3 |
| Rcpp | 1.1.2 |
| RcppArmadillo | 15.4.2-1 |
| RcppEigen | 0.3.4.0.2 |
| readxl | 1.5.0 |
| rematch | 2.0.0 |
| reshape | 0.8.10 |
| rjson | 0.2.23 |
| rlang | 1.3.0 |
| rmarkdown | 2.31 |
| robustbase | 0.99-7 |
| S7 | 0.2.2 |
| sass | 0.4.10 |
| scales | 1.4.0 |
| shape | 1.4.6.1 |
| SparseM | 1.84-2 |
| stringi | 1.8.7 |
| stringr | 1.6.0 |
| survival | 3.8-9 |
| susieR | 0.14.2 |
| sys | 3.4.3 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| tidyselect | 1.2.1 |
| tinytex | 0.60 |
| utf8 | 1.2.6 |
| vctrs | 0.7.3 |
| viridisLite | 0.4.3 |
| withr | 3.0.3 |
| xfun | 0.60 |
| yaml | 2.3.12 |

</details>

## Worked examples

From this directory, run:

```sh
Rscript examples/run_examples.R
```

The script analyzes one exposure and three exposures and checks the posterior-mean points, corrected-score centers, standard errors, confidence limits and p-values against the supplied expected outputs. The first run compiles small numerical helpers. Generated files are kept under `runtime/` and `results/`.

For use from R:

```r
source("R/load.R")
model <- make_mr_susie()
x <- jsonlite::read_json("examples/joint_strong_0001.json", simplifyVector = TRUE)
fit <- model$fit(
  x$B, matrix(x$sx, nrow(x$B), ncol(x$B)),
  x$y, rep(x$sy, nrow(x$B)), x$R, x$Sigma,
  x$selection_y, rep(x$sy, nrow(x$B)), 10000,
  rep("Outcome SD per exposure SD", ncol(x$B))
)
fit$point
fit$score_inference
```

Point estimates are posterior means. Intervals and tests are centered on the separate corrected-score estimates. Exposure units must be specified consistently. See [the method interface](docs/method.md).

## Analyses

A single planned simulation can be run with:

```sh
Rscript analysis/run_simulations.R --plan=config/simulation_evaluation_plan.csv --limit=1 --out=results/simulation_example
Rscript analysis/evaluate.R --input=results/simulation_example/coefficient_results.csv --out=results/simulation_example/metrics.csv
```

`--limit=all` explicitly requests the complete supplied plan. Add `--comparators=true` to run the retained comparators. The native MVMR-cML-SuSiE adapter uses BIC without data perturbation. Full-LD and exposure-only-clumped panels are saved separately. Existing nonempty output directories are rejected to protect completed results.

The dedicated prospective validation has its own generator and plan:

```sh
Rscript analysis/run_prospective_validation.R --limit=1 --out=results/prospective_example
```

Application scripts accept harmonized inputs or the synthetic example:

```sh
Rscript analysis/run_applications.R --input=examples/joint_strong_0001.json --n_eff=10000 --me=true --out=results/application_example
```

This command illustrates the application workflow with synthetic data. The 76 study configurations are listed in `config/application_configurations.csv`. Instructions for obtaining and preparing the study inputs are in [data access](docs/data_access.md). The supplied source-data tables contain derived results, not raw GWAS or participant-level data.

## Support-union sensitivity

The optional Gaussian sensitivity is implemented in R/support_union.R. It allows a prespecified maximum number of invalid variants and retains empty, disconnected and unbounded confidence sets.

From R:

```r
source("R/support_union.R")
x <- jsonlite::read_json("examples/support_union_example.json", simplifyVector = TRUE)
result <- do.call(union_ar, x)
result$sets
result$joint_null_p
```

The same example can be run from the command line:

```sh
Rscript analysis/support_union.R --input=examples/support_union_example.json --out=results/support_union_example.rds
```

The saved RDS object preserves infinite interval endpoints. The Gaussian assumptions and treatment of singular nuisance blocks are described in docs/analysis.md.

The R implementation was checked against 5,000 saved sensitivity calculations and 1,064 numerical cases. The supplied tests cover interval geometry, chi-square calculations, extreme-scale calculations and input validation. Independent checks also compare projections with closed-form ellipsoids and verify invariance to rescaling and reordering.

## Tables and runtime

Table S1 is supplied in data/Table_S1.xlsx. Run Rscript analysis/export_table.R to export its sheets as CSVs. Additional derived numerical tables are under data/reported_metrics/.

The ten-pair timing inputs, roster and recorded results are supplied separately. `Rscript analysis/benchmark_runtime.R --run=true` requests exactly ten paired attempts with alternating method order. Package initialization and compilation precede timing. Failed and unavailable fits remain in the output. The recorded timings are specific to the original machine and software environment.

## Reproducibility

The simulation families, phases, seeds, calibration rules, estimands and treatment of missing fits are described in [the analysis guide](docs/analysis.md). Calibration and evaluation remain separate. Existing UVMR cutoffs are not applied to the updated single-exposure fits. Marginal UVMR and TWAS associations are not scored against conditional generating coefficients.

The checks can be run with:

```sh
Rscript tests/check_model.R
Rscript tests/check_analysis.R
Rscript tests/check_support_union.R
```

The first numerical test may compile the helpers.

## Citation and license

Please cite the accompanying manuscript, “MR-SuSiE: A sparse Bayesian method for identifying causal exposures and invalid instruments in univariable and multivariable Mendelian randomization,” by Siyi Chen. Citation metadata are provided in `CITATION.cff`.

The analysis code is distributed under the MIT license. External packages and source data retain their respective licenses and access conditions.

<details>
<summary>Third-party license and attribution notices</summary>

```text
Third-party dependencies

The analysis implementation is supplied under LICENSE. External packages remain under their own licenses. Their source code and binaries are not included in this directory. The following notices identify the versions used for verification.

Rcpp 1.1.2
License: GPL (>= 2)

RcppArmadillo 15.4.2-1
License: GPL (>= 2)

susieR 0.14.2
License: BSD_3_clause + file LICENSE
YEAR: 2017-2022
COPYRIGHT HOLDER: Gao Wang, Peter Carbonetto, Yuxin Zou, Kaiqian Zhang, Matthew Stephens
ORGANIZATION: Gao Wang, Peter Carbonetto, Yuxin Zou, Kaiqian Zhang, Matthew Stephens

MASS 7.3-66
License: GPL-2 | GPL-3

jsonlite 2.0.0
License: MIT + file LICENSE
YEAR: 2020
COPYRIGHT HOLDER: Jeroen Ooms

digest 0.6.39
License: GPL (>= 2)

data.table 1.18.4
License: MPL-2.0 | file LICENSE
Mozilla Public License Version 2.0
==================================

1. Definitions
--------------

1.1. "Contributor"
    means each individual or legal entity that creates, contributes to
    the creation of, or owns Covered Software.

1.2. "Contributor Version"
    means the combination of the Contributions of others (if any) used
    by a Contributor and that particular Contributor's Contribution.

1.3. "Contribution"
    means Covered Software of a particular Contributor.

1.4. "Covered Software"
    means Source Code Form to which the initial Contributor has attached
    the notice in Exhibit A, the Executable Form of such Source Code
    Form, and Modifications of such Source Code Form, in each case
    including portions thereof.

1.5. "Incompatible With Secondary Licenses"
    means

    (a) that the initial Contributor has attached the notice described
        in Exhibit B to the Covered Software; or

    (b) that the Covered Software was made available under the terms of
        version 1.1 or earlier of the License, but not also under the
        terms of a Secondary License.

1.6. "Executable Form"
    means any form of the work other than Source Code Form.

1.7. "Larger Work"
    means a work that combines Covered Software with other material, in 
    a separate file or files, that is not Covered Software.

1.8. "License"
    means this document.

1.9. "Licensable"
    means having the right to grant, to the maximum extent possible,
    whether at the time of the initial grant or subsequently, any and
    all of the rights conveyed by this License.

1.10. "Modifications"
    means any of the following:

    (a) any file in Source Code Form that results from an addition to,
        deletion from, or modification of the contents of Covered
        Software; or

    (b) any new file in Source Code Form that contains any Covered
        Software.

1.11. "Patent Claims" of a Contributor
    means any patent claim(s), including without limitation, method,
    process, and apparatus claims, in any patent Licensable by such
    Contributor that would be infringed, but for the grant of the
    License, by the making, using, selling, offering for sale, having
    made, import, or transfer of either its Contributions or its
    Contributor Version.

1.12. "Secondary License"
    means either the GNU General Public License, Version 2.0, the GNU
    Lesser General Public License, Version 2.1, the GNU Affero General
    Public License, Version 3.0, or any later versions of those
    licenses.

1.13. "Source Code Form"
    means the form of the work preferred for making modifications.

1.14. "You" (or "Your")
    means an individual or a legal entity exercising rights under this
    License. For legal entities, "You" includes any entity that
    controls, is controlled by, or is under common control with You. For
    purposes of this definition, "control" means (a) the power, direct
    or indirect, to cause the direction or management of such entity,
    whether by contract or otherwise, or (b) ownership of more than
    fifty percent (50%) of the outstanding shares or beneficial
    ownership of such entity.

2. License Grants and Conditions
--------------------------------

2.1. Grants

Each Contributor hereby grants You a world-wide, royalty-free,
non-exclusive license:

(a) under intellectual property rights (other than patent or trademark)
    Licensable by such Contributor to use, reproduce, make available,
    modify, display, perform, distribute, and otherwise exploit its
    Contributions, either on an unmodified basis, with Modifications, or
    as part of a Larger Work; and

(b) under Patent Claims of such Contributor to make, use, sell, offer
    for sale, have made, import, and otherwise transfer either its
    Contributions or its Contributor Version.

2.2. Effective Date

The licenses granted in Section 2.1 with respect to any Contribution
become effective for each Contribution on the date the Contributor first
distributes such Contribution.

2.3. Limitations on Grant Scope

The licenses granted in this Section 2 are the only rights granted under
this License. No additional rights or licenses will be implied from the
distribution or licensing of Covered Software under this License.
Notwithstanding Section 2.1(b) above, no patent license is granted by a
Contributor:

(a) for any code that a Contributor has removed from Covered Software;
    or

(b) for infringements caused by: (i) Your and any other third party's
    modifications of Covered Software, or (ii) the combination of its
    Contributions with other software (except as part of its Contributor
    Version); or

(c) under Patent Claims infringed by Covered Software in the absence of
    its Contributions.

This License does not grant any rights in the trademarks, service marks,
or logos of any Contributor (except as may be necessary to comply with
the notice requirements in Section 3.4).

2.4. Subsequent Licenses

No Contributor makes additional grants as a result of Your choice to
distribute the Covered Software under a subsequent version of this
License (see Section 10.2) or under the terms of a Secondary License (if
permitted under the terms of Section 3.3).

2.5. Representation

Each Contributor represents that the Contributor believes its
Contributions are its original creation(s) or it has sufficient rights
to grant the rights to its Contributions conveyed by this License.

2.6. Fair Use

This License is not intended to limit any rights You have under
applicable copyright doctrines of fair use, fair dealing, or other
equivalents.

2.7. Conditions

Sections 3.1, 3.2, 3.3, and 3.4 are conditions of the licenses granted
in Section 2.1.

3. Responsibilities
-------------------

3.1. Distribution of Source Form

All distribution of Covered Software in Source Code Form, including any
Modifications that You create or to which You contribute, must be under
the terms of this License. You must inform recipients that the Source
Code Form of the Covered Software is governed by the terms of this
License, and how they can obtain a copy of this License. You may not
attempt to alter or restrict the recipients' rights in the Source Code
Form.

3.2. Distribution of Executable Form

If You distribute Covered Software in Executable Form then:

(a) such Covered Software must also be made available in Source Code
    Form, as described in Section 3.1, and You must inform recipients of
    the Executable Form how they can obtain a copy of such Source Code
    Form by reasonable means in a timely manner, at a charge no more
    than the cost of distribution to the recipient; and

(b) You may distribute such Executable Form under the terms of this
    License, or sublicense it under different terms, provided that the
    license for the Executable Form does not attempt to limit or alter
    the recipients' rights in the Source Code Form under this License.

3.3. Distribution of a Larger Work

You may create and distribute a Larger Work under terms of Your choice,
provided that You also comply with the requirements of this License for
the Covered Software. If the Larger Work is a combination of Covered
Software with a work governed by one or more Secondary Licenses, and the
Covered Software is not Incompatible With Secondary Licenses, this
License permits You to additionally distribute such Covered Software
under the terms of such Secondary License(s), so that the recipient of
the Larger Work may, at their option, further distribute the Covered
Software under the terms of either this License or such Secondary
License(s).

3.4. Notices

You may not remove or alter the substance of any license notices
(including copyright notices, patent notices, disclaimers of warranty,
or limitations of liability) contained within the Source Code Form of
the Covered Software, except that You may alter any license notices to
the extent required to remedy known factual inaccuracies.

3.5. Application of Additional Terms

You may choose to offer, and to charge a fee for, warranty, support,
indemnity or liability obligations to one or more recipients of Covered
Software. However, You may do so only on Your own behalf, and not on
behalf of any Contributor. You must make it absolutely clear that any
such warranty, support, indemnity, or liability obligation is offered by
You alone, and You hereby agree to indemnify every Contributor for any
liability incurred by such Contributor as a result of warranty, support,
indemnity or liability terms You offer. You may include additional
disclaimers of warranty and limitations of liability specific to any
jurisdiction.

4. Inability to Comply Due to Statute or Regulation
---------------------------------------------------

If it is impossible for You to comply with any of the terms of this
License with respect to some or all of the Covered Software due to
statute, judicial order, or regulation then You must: (a) comply with
the terms of this License to the maximum extent possible; and (b)
describe the limitations and the code they affect. Such description must
be placed in a text file included with all distributions of the Covered
Software under this License. Except to the extent prohibited by statute
or regulation, such description must be sufficiently detailed for a
recipient of ordinary skill to be able to understand it.

5. Termination
--------------

5.1. The rights granted under this License will terminate automatically
if You fail to comply with any of its terms. However, if You become
compliant, then the rights granted under this License from a particular
Contributor are reinstated (a) provisionally, unless and until such
Contributor explicitly and finally terminates Your grants, and (b) on an
ongoing basis, if such Contributor fails to notify You of the
non-compliance by some reasonable means prior to 60 days after You have
come back into compliance. Moreover, Your grants from a particular
Contributor are reinstated on an ongoing basis if such Contributor
notifies You of the non-compliance by some reasonable means, this is the
first time You have received notice of non-compliance with this License
from such Contributor, and You become compliant prior to 30 days after
Your receipt of the notice.

5.2. If You initiate litigation against any entity by asserting a patent
infringement claim (excluding declaratory judgment actions,
counter-claims, and cross-claims) alleging that a Contributor Version
directly or indirectly infringes any patent, then the rights granted to
You by any and all Contributors for the Covered Software under Section
2.1 of this License shall terminate.

5.3. In the event of termination under Sections 5.1 or 5.2 above, all
end user license agreements (excluding distributors and resellers) which
have been validly granted by You or Your distributors under this License
prior to termination shall survive termination.

************************************************************************
*                                                                      *
*  6. Disclaimer of Warranty                                           *
*  -------------------------                                           *
*                                                                      *
*  Covered Software is provided under this License on an "as is"       *
*  basis, without warranty of any kind, either expressed, implied, or  *
*  statutory, including, without limitation, warranties that the       *
*  Covered Software is free of defects, merchantable, fit for a        *
*  particular purpose or non-infringing. The entire risk as to the     *
*  quality and performance of the Covered Software is with You.        *
*  Should any Covered Software prove defective in any respect, You     *
*  (not any Contributor) assume the cost of any necessary servicing,   *
*  repair, or correction. This disclaimer of warranty constitutes an   *
*  essential part of this License. No use of any Covered Software is   *
*  authorized under this License except under this disclaimer.         *
*                                                                      *
************************************************************************

************************************************************************
*                                                                      *
*  7. Limitation of Liability                                          *
*  --------------------------                                          *
*                                                                      *
*  Under no circumstances and under no legal theory, whether tort      *
*  (including negligence), contract, or otherwise, shall any           *
*  Contributor, or anyone who distributes Covered Software as          *
*  permitted above, be liable to You for any direct, indirect,         *
*  special, incidental, or consequential damages of any character      *
*  including, without limitation, damages for lost profits, loss of    *
*  goodwill, work stoppage, computer failure or malfunction, or any    *
*  and all other commercial damages or losses, even if such party      *
*  shall have been informed of the possibility of such damages. This   *
*  limitation of liability shall not apply to liability for death or   *
*  personal injury resulting from such party's negligence to the       *
*  extent applicable law prohibits such limitation. Some               *
*  jurisdictions do not allow the exclusion or limitation of           *
*  incidental or consequential damages, so this exclusion and          *
*  limitation may not apply to You.                                    *
*                                                                      *
************************************************************************

8. Litigation
-------------

Any litigation relating to this License may be brought only in the
courts of a jurisdiction where the defendant maintains its principal
place of business and such litigation shall be governed by laws of that
jurisdiction, without reference to its conflict-of-law provisions.
Nothing in this Section shall prevent a party's ability to bring
cross-claims or counter-claims.

9. Miscellaneous
----------------

This License represents the complete agreement concerning the subject
matter hereof. If any provision of this License is held to be
unenforceable, such provision shall be reformed only to the extent
necessary to make it enforceable. Any law or regulation which provides
that the language of a contract shall be construed against the drafter
shall not be used to construe this License against a Contributor.

10. Versions of the License
---------------------------

10.1. New Versions

Mozilla Foundation is the license steward. Except as provided in Section
10.3, no one other than the license steward has the right to modify or
publish new versions of this License. Each version will be given a
distinguishing version number.

10.2. Effect of New Versions

You may distribute the Covered Software under the terms of the version
of the License under which You originally received the Covered Software,
or under the terms of any subsequent version published by the license
steward.

10.3. Modified Versions

If you create software not governed by this License, and you want to
create a new license for such software, you may create and use a
modified version of this License if you rename the license and remove
any references to the name of the license steward (except to note that
such modified license differs from this License).

10.4. Distributing Source Code Form that is Incompatible With Secondary
Licenses

If You choose to distribute Source Code Form that is Incompatible With
Secondary Licenses under the terms of this version of the License, the
notice described in Exhibit B of this License must be attached.

Exhibit A - Source Code Form License Notice
-------------------------------------------

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

If it is not possible or desirable to put the notice in a particular
file, then You may include the notice in a location (such as a LICENSE
file in a relevant directory) where a recipient would be likely to look
for such a notice.

You may add additional accurate notices of copyright ownership.

Exhibit B - "Incompatible With Secondary Licenses" Notice
---------------------------------------------------------

  This Source Code Form is "Incompatible With Secondary Licenses", as
  defined by the Mozilla Public License, v. 2.0.

readxl 1.5.0
License: MIT + file LICENSE
YEAR: 2026
COPYRIGHT HOLDER: readxl authors

MendelianRandomization 0.10.0
License: GPL-2 | GPL-3

MRcML 0.0.0.9000
License: MIT + file LICENSE
YEAR: 2021
COPYRIGHT HOLDER: Haoran Xue; Xiaotong Shen; Wei Pan

MVMRcML 0.1.0
License: MIT + file LICENSE
YEAR: 2023
COPYRIGHT HOLDER: Zhaotong Lin

MVMRcMLSuSiE 0.1.0
License: GPL (>= 3)
```

</details>
