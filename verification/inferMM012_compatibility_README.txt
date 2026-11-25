inferMM 0.1.2 reproduction compatibility (20 November 2026)

The archived full simulations were executed with inferMM 0.1.1. Their original
seeds, session information, fingerprints, outputs and CPU records are retained.
They have not been relabeled as 0.1.2 executions.

The current archive contains inferMM 0.1.2 with 192 SDL demonstration rows.
The initial candidate and preceding documentation-only builds retained all
12 R source files from 0.1.1; the historical checks below refer to those builds.
The current default-screen build changes fit_mm.R and methods.R. Ten other
R source files and the explicitly specified variance estimating functions
remain unchanged; see the separate current update record. SDL28/SDL56
simulation drivers explicitly retain the original seven-concentration subset.
Restoring the eighth observed concentration therefore does not change those
designs or the reported simulation reference values.

Checks
- Single50/single100: 2 replications per generating condition under each
  version, 48 rows per sample size; paired numerical equality.
- Alves 3/6-cluster modules: 2 replications per latent/error condition under
  each version, 40 rows per template; paired numerical equality.
- Additional benchmark: 10 replications per condition, 1,120 method rows;
  equality between versions and with the archived prefix.
- Matched comparison: 10 replications per condition, 460 method rows;
  equality between versions and with the archived prefix.
- Higher-moment comparison: 20 replications per condition, 240 method rows;
  equality with the archived prefix.
- Bootstrap smoke: documented n=20 and n=30 datasets, 399 resamples each;
  798 bootstrap refits reproduce their archived intervals.
- Projected-score illustration: profile/projected estimates, SEs and 95% limits
  reproduce the archived illustration.
- Corrected real data: all 5 tables match the corrected 0.1.1 production results
  exactly, including all 24 SDL model rows and unchanged Alves values.

Timing columns are excluded. Logical TRUE/FALSE and numeric 1/0 indicators
are compared by their common value. The paired empirical SDL template differs
only at floating-point rounding level (maximum relative difference 8.42e-16).
These bounded checks assess compatibility. They do not replace the archived
full simulation study or establish independence in the real experiment.

Files
inferMM012_estimation_R_sha256.csv: 12 pairs of estimation-source hashes.
inferMM012_numeric_checks.csv: check sizes and maximum differences.
inferMM012_demo_check.csv: final installed package and 192-row data assertions.
inferMM012_compatibility.json: final archive pin and detailed provenance.
inferMM012_package_check.json/txt: final package build/check report.

The numerical checks used the first 0.1.2 candidate. After documentation and
packaging cleanup, the final archive was reinstalled. Its 12 estimation R
files and demonstration CSV/RDA are byte-identical to the checked candidate.
Final installed demo assertions passed. packages/package_build_current.csv
and inferMM012_compatibility.json record the final archive SHA256.

Current execution examples (after package setup in README.txt)
  Rscript run.R real
  Rscript verify.R real
  Rscript run.R bootstrap_smoke
  Rscript examples/run_ma_genton_example.R
Use a separate extracted working copy for bounded simulation runs so that
partial-run outputs remain distinct from full-run outputs.

Later data-source reference update (20 November 2026)
The current_package pin now identifies the documentation-only reference build.
The original numeric-validation candidate, comparison records and saved
analysis values above are preserved. All 12 estimation R files, NAMESPACE
and R/CSV data are byte-identical to the preceding final 0.1.2 build.
See inferMM012_reference_update.json for the separate archive comparison and
installed-data assertions. The new package help PDF is in docs/.

METHOD-REFERENCE MANUAL UPDATE (20 November 2026)
The current 0.1.2 source package and PDF help manual now include references
for conditional-moment estimation, fixed-weight regression, analytic height
elimination and the wild-bootstrap multipliers. The package overview also
identifies related error-distribution background and the projected-score option.
screen_mm() and its AIC/BIC help are unchanged; no AIC/BIC references were added.
The exact current build is pinned in packages/package_build_current.csv.
The preceding documentation_reference_update and candidate-comparison fields
are historical records for their named earlier archives; the current method
reference update is recorded separately in inferMM012_method_reference_update.json.
All estimation R files, NAMESPACE, demo data, saved simulation/coverage/timing
results and prior numeric-validation records are preserved. No new simulation
or empirical refit was required for this documentation-only update.
The new ZIP is reproducibility_9078092_20261120_methodrefs.zip; earlier dated
archives remain available. See CHANGES_20261120_methodrefs_ko.txt.

CURRENT DEFAULT-SCREEN BUILD (20 November 2026)
The current package pin and check report identify the revised fitting
dispatcher and display methods, with default four-model screening.
Historical numeric comparisons, the candidate-vs-final packaging comparison,
bootstrap/coverage/timing records and prior reference-update fields retain
their original scope. The current update is recorded separately in
inferMM012_default_screen_update.json, including ten unchanged R files,
unchanged numerical-function hashes, installation checks and default/fixed
fit equivalence. Original code-identity CSV values were not overwritten.

CURRENT CRAN-PREFLIGHT METADATA BUILD (20 November 2026)
The current_package pin now identifies citation/DOI documentation corrections;
the estimation R sources, data and tests are identical to the checked default-
screen build. All prior numeric-check and build-specific records above are
historical and have not been overwritten. Current source/PDF preflight records
are in cran_preflight_20261120/, and cran_preflight_update gives their scope.
The one current check NOTE is a local HTML Tidy/V8 availability issue.
