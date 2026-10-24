# Current-source R numerical references

These fixtures were generated from clean source checkouts of `did`
`74b88eb07f1faa644df1271055a0f77848f6550c` (2.5.1.902) and `DRDID`
`eaac3b49cd225d0d7706c633858555a5f204e21e` (1.4.0). The commits were the
upstream heads verified on 2026-10-23. The word “current” identifies this
recorded verification date; the references do not update themselves.

`cases.json` defines 23 analytic configurations and two bootstrap
configurations. Python evaluates each in both execution modes. The analytic
cases cover weighted balanced panels, repeated cross sections, and unbalanced
panels; DR, regression, and IPW estimators; never-treated/varying-base and
not-yet-treated/universal-base settings; and the default, first-period,
base-period, and varying weight options. The two extra varying-weight cases
exercise the repeated-cross-section calculation on panel observations. The
unbalanced inputs deliberately omit some units at weight anchor periods.
Three separately declared analytic-cluster cases cover panel DR,
repeated-cross-section regression, and unbalanced-panel IPW with never-treated
controls and varying bases. These extend the original 22-case design without
changing its inputs or cases.

All cases include the group-time estimates, standard errors, critical values,
interval endpoints, Wald statistic and p-value, and full unit influence
functions. Simple, group, dynamic, and calendar aggregations include every
reported estimate, standard error and critical value, plus every aggregate
influence-function component. Normalization cells and the clustered Wald
suppression exercise missing-value masks. This is a representative regression
suite, not exhaustive coverage of every combination of arguments, singular
designs, or numerical conditioning.

## Regeneration

Use R with `pkgload`, `jsonlite`, `digest`, and the dependencies of both source
packages installed, then run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Rscript csdid/test_csdid/r_ref/current/generate.R \
  --did /path/to/did --drdid /path/to/DRDID
```

The script verifies the pinned commits and clean tracked source trees, loads
DRDID before did, and verifies that both namespaces come from those checkouts.
It reuses the committed input CSVs. `--regenerate-inputs` recreates them from
the fixed data-generation specification in the script; it is not needed to
refresh numerical outputs. `--output /path/to/directory` selects an alternate
output directory. Updating the reference commits is a deliberate code change.
Normal regeneration fails if any frozen input or multiplier CSV is missing;
it never recreates inputs implicitly. For a fresh output directory, copy the
five committed input/multiplier CSVs there first, or explicitly request
`--regenerate-inputs`.
Do not replace inputs, remove cases, or broaden tolerances to make a mismatch
pass.

R runs every case in both modes. `references.json.gz` contains the fast-mode
reference; `r_mode_checks.json` records each field's support check and maximum
absolute fast/standard difference. `manifest.json` records commits, SHA256
hashes of tracked R source files and fixture artifacts, loaded package
versions, environment information, and the tolerance policy declared before
comparison. `sessionInfo.txt` gives the complete R session details.

## Encoding and row identity

Each numeric payload has an explicit `shape`, row-major `values`, and separate
`nan`, `posinf`, and `neginf` arrays of zero-based flat positions. Nonfinite
positions have a zero placeholder in `values`. An absent R field remains JSON
`null`. Influence functions are always matrices, including one-column
components. Support and nonfinite masks must agree exactly; a numerical
tolerance cannot excuse a missing cell or mismatched missing value.

`unit_ids` lists the original shared-input identifiers in increasing order.
R's authoritative influence-function row names determine the mapping. For
true repeated cross sections, R replaces the id with a one-based input row
number; the exporter maps that number back to the original input `id`.
Python's corresponding synthetic row number is zero-based and needs the same
mapping. Aggregate influence functions use the identical row permutation.
Group-time columns must be compared by their `(group, t)` keys.

The predeclared componentwise tolerances are `atol=1e-7, rtol=1e-6` for public
outputs and influence functions, and `atol=1e-5, rtol=1e-5` for Wald outputs.
Existing broader differential audits found small cross-solver differences in
individual influence-function components. Any such mismatch in these fixed
cases must be retained and reported explicitly, rather than accommodated by
silently increasing a tolerance.

## Shared-draw bootstrap

The unclustered panel DR case and clustered panel regression case use 99
explicit frozen normal multiplier draws. The CSV column headers identify
units or clusters, and every bootstrap invocation reuses the same complete
matrix. R temporarily replaces its namespace function
`run_multiplier_bootstrap(inf.func, biters, pl, cores)` with the documented
matrix projection. The Python test replaces the corresponding helper.

The projection first aligns influence-function rows to the draw column IDs,
then computes `draws @ influence / number_of_rows`. The ordinary package
bootstrap code still performs its unit/cluster scaling, IQR standard-error
calculation, and confidence-band calculation. These fixtures test identical
draws and downstream inference; they do not claim that R and NumPy RNG seeds
produce identical multipliers or that 99 draws suffice for an empirical
analysis.
