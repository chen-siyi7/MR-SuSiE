# Study data and input preparation

## Source resources

| Resource | Access |
|---|---|
| Global Lipids Genetics Consortium 2013 | https://csg.sph.umich.edu/willer/public/lipids2013/ |
| CARDIoGRAMplusC4D | https://cardiogramplusc4d.org/data-downloads/ |
| NHLBI GRASP full-results mirror | https://www.nhlbi.nih.gov/grasp/full-results |
| Oxford BIG imaging release | https://open.oxcin.ox.ac.uk/ukbiobank/big40/BIGv2/ |
| Bellenguez et al. Alzheimer disease stage I | https://www.ebi.ac.uk/gwas/studies/GCST90027158 |
| SCALLOP CVD1 | https://doi.org/10.5281/zenodo.2615265 |
| 1000 Genomes Phase 3 reference | https://www.internationalgenome.org/data-portal/data-collection/phase3/ |

Obtain these resources under their source access and reuse conditions. GRASP restricts reposting its catalog or significant subsets and commercial use without permission. Its terms are at https://www.nhlbi.nih.gov/grasp/terms. IGSR reuse information is at https://www.internationalgenome.org/IGSR_disclaimer/.

This upload directory contains synthetic examples and derived numerical results. Raw source GWAS, reference genotypes and participant-level data are not included.

## Harmonization

Match genome builds and variants before forming the analysis table. Express exposure and outcome effects on the same allele orientation. Exclude incompatible and palindromic variants. Selection and clumping use exposure associations only. Align signed LD to the reference allele, then multiply per-allele effects and standard errors by the reference-genotype SD.

`analysis/prepare_application_input.R` accepts an aligned wide summary CSV and a signed reference-LD CSV. Required columns are `variant`, `effect_allele`, `other_allele`, `reference_effect_allele`, `reference_other_allele`, `reference_sd`, `beta_outcome`, `se_outcome` and `beta_TRAIT`, `se_TRAIT` for each supplied trait. LD row and column names must be variant identifiers. The script excludes incompatible, palindromic and incomplete rows, aligns effect directions and applies reference-SD scaling.

For example:

```sh
Rscript analysis/prepare_application_input.R --summary=aligned_summary.csv --ld=reference_ld.csv --traits=LDL,HDL,TG --n1=100000 --n2=160000 --out=lipids.rds
Rscript analysis/run_applications.R --input=lipids.rds --n_eff=160000 --me=true --out=results/lipids
```

The numbers in this example illustrate the interface. Use the study-specific sample sizes and configurations in `config/application_configurations.csv` for the manuscript analyses.

## Study configurations

The lipid panel contains 149 variants. The imaging panels contain 96 and 330 variants at reference r² thresholds 0.01 and 0.05. The primary imaging/protein outcome is the complete Bellenguez stage-I release on the log-odds scale. Its effective sample size uses the working cohort-information approximation described in the manuscript, approximately 127,661. The listed sample-size sensitivities retain the original summary coefficients and standard errors.

The protein panels use r² thresholds 0.01, 0.05, 0.10, 0.20, 0.30 and 0.40. Eligible exposures have at least three own significant instruments after complete-row filtering. The primary r²=0.01 panel has 53 proteins. Larger panels have their own eligible-exposure counts and family thresholds.

Place the harmonized study RDS files under `data/application_inputs/` using the names in `config/application_configurations.csv`, or supply a manifest with `job_id`, `input_file`, `configuration` and `n_eff`. Then run:

```sh
Rscript analysis/run_applications.R --manifest=config/application_configurations.csv --limit=all --me=true --out=results/applications
```

Primary application estimation uses NOME. The ME option retains the same selected nuisance support and applies the supplied exposure-error covariance. Cross-GWAS exposure/outcome error covariance remains zero. Full-LD and exposure-only-clumped panels are distinct. Findings use p<0.05/M within the eligible exposure panel. Source overlap, selection and repeated sensitivity analyses remain outside that within-panel adjustment.
