# What package alignment establishes

The target is the six shared traditional DR, normalized IPW and regression
estimators, each for panel and repeated-cross-section data. It is not the entire
R DRDID API. The current PyPI releases remain `csdid==0.4.3` and `drdid==1.1.6`;
the corrections in CSDID #88 and DRDIDpy #21 are release candidates.

## Why CSDID could pass while standalone DRDID differed

| CSDID route | Python implementation actually used | R numerical counterpart |
|---|---|---|
| `dr`, panel/RC | Local `csdid.attgt_fnc.drdid_trim` | Traditional `DRDID::drdid_panel` / `drdid_rc` |
| `ipw`, panel/RC | Local `csdid.attgt_fnc.drdid_trim` | Normalized `DRDID::std_ipw_did_panel` / `std_ipw_did_rc` |
| `reg`, panel/RC | External `drdid.reg_did` | `DRDID::reg_did_panel` / `reg_did_rc` |

CSDID's local DR/IPW kernels already contained corrections missing from the
standalone distribution. Its public preprocessing also supplies weights and
constructs binary treatment indicators, bypassing some old standalone input
bugs. Therefore the earlier public CSDID-versus-R comparisons did not certify
the complete standalone Python package. Calling both packages fully compatible
would have exceeded that evidence.

The external DRDID package is still required at installation and import time:
CSDID loads its regression module when loading the estimator dispatch. A DR or
IPW fit uses the local kernels above, but removing the installed dependency
would still prevent CSDID from loading normally.

The new dtype audit also found genuine gaps in CSDID: integer subtraction could
overflow in the faster path, and direct float32 kernels did not always match
the identical represented values reduced in double precision. Built-in paths now
convert numeric arrays before arithmetic, matching R; raw user columns and
custom-callback inputs retain their prior contract. Closed-form and dtype
regressions cover these corrections.

## One reviewed installed pair

The proposed pair is `csdid==0.4.4` with `drdid==1.1.7`. CSDID now requires
`drdid==1.1.7`, so installs cannot silently retain the known-broken 1.1.6 or
select a future unreviewed backend. A backend upgrade must update the exact
dependency and peer lock together and pass the coordinated checks.
The exact reviewed peer commit/version is recorded in
[joint_sync/python_drdid.json](joint_sync/python_drdid.json).

Both mandatory Python CI jobs build the two exact Git revisions into wheels,
install them together in a fresh virtual environment, run both numerical suites
and the 72-case direct kernel bridge, and check module origins and every runtime
file's source/wheel/installation hash. The bridge enforces both raw-output and
predeclared original-unit tolerances. The full local run includes the two slow
tests; CI explicitly deselects them. JUnit, source IDs, wheel hashes, dependency
metadata and import-origin receipts are retained as artifacts.

The completed local pair (`e8354bc` / `eb89664`) passes **2,935 tests and 19
subtests**, with no failures, errors or skips, plus all 72 direct comparisons.
It also preserves exact numeric bytes, warnings and RNG state across 240
healthy workflows and 960 aggregates versus the prior CSDID baseline; see the
[installed-pair identity receipt](audit20261027/joint_installed_healthy_identity.json).

The separate required R job still generates current official R references.
Any unreviewed tracked R-tree change now requires review, including source-only
changes whose existing numerical cases still pass. Test changes and other source
changes remain distinct evidence. The runner still executes numerical checks
and preserves numerical-failure precedence, pins, fixtures and tolerances.

Publishing waits for the successful tested main revision. DRDID must be released
first. Before CSDID uploads its own release, its workflow requires the backend
wheel to exist on PyPI and match every reviewed Python source byte and the PyPI
digest. The existing publishing credentials are retained. No merge or upload is
performed by this audit.

For later backend upgrades, first validate both candidate commits. DRDID's
integration job also checks CSDID's default branch, so merge CSDID's reviewed
lock/pin update first; its publisher will block while the backend is absent.
Then rerun DRDID's checks against that CSDID branch, merge/release DRDID, and
retry CSDID publication on the tested main revision. This is a coordinated
upgrade, not an automatic adoption of a new backend.

## Boundaries retained explicitly

Equivalent results require the same retained sample, design/intercept, weights,
trimming, estimator and inference settings. R's high-level `drdid()` defaults to
improved/IPT; the shared functions here use traditional DR. Improved/IPT, TWFE,
alternative RC estimators, nonnormalized IPW, R formula wrappers and weighted
nuisance-refit bootstrap APIs are not added in this effort. Equal R and NumPy
seeds do not mean equal random draws; cross-language bootstrap checks use the
same frozen multipliers where tested.

There is also a retained numerical limitation. The adversarial outcome-scale
audit still has **36 raw aggregate-IF comparisons failing at outcome scale
`1e10`**, while ATT, SE, confidence limits, critical values and Wald outputs pass.
R's own fast/slow paths miss 42 raw IF fields in that audit. The largest recorded
Python/R violating difference is about `0.1288` in a matrix with entries on the
order of `1e11`, or `1.29e-11` in the original outcome units. The small
implementation discrepancy is already present in original units and is
amplified by the outcome scale; this is not solely an error introduced by
rescaling. The replay keeps exit 1; the 72
original-unit diagnostics passing do not erase the raw failures. See
[UPSTREAM_SYNC_20261026.md](UPSTREAM_SYNC_20261026.md). This is why neither a
passing fixed grid nor these releases establish universal numerical equivalence
within those fixed tolerances.

R-only pushes do not trigger the existing Python push workflow. Updating a
package still requires a coordinated reviewed port and release; CI detects
covered differences and unreviewed source changes rather than translating code.
