# csdid — Pre-existing Python Bug Fixes (unrelated to the R `did` v2.5.1 sync)

This document lists defects **in the Python `csdid` port itself** that were fixed
during this work but are **not** part of syncing with R `did` v2.5.1. They were
wrong (crashes or incorrect output) regardless of any R feature — i.e., they
would be bugs in csdid even with no R update. The R-feature work (input
validation, `notyettreated`, clustered inference, `compute_inffunc`,
`fix_weights`, `faster_mode`, etc.) is documented separately in
`CHANGELOG_v0.3.0.md`.

| # | Bug | Severity | Symptom | Fix |
|--:|-----|----------|---------|-----|
| 1 | `plot_attgt` single group | Medium | `TypeError` | `squeeze=False` |
| 2 | Panel post-treatment flag | High | wrong `post` column | pass `pst=post_treat` |
| 3 | Panel global-`n` overwrite | High (latent) | corrupt influence functions | use a cell-local `n_cell` |
| 4 | `allow_unbalanced_panel=False` branch | High | immediate `TypeError` | fix 3 typos |
| 5 | `raise "<string>"` in error paths | Medium | `TypeError` masks the real message | `raise ValueError(...)` |
| 6 | Results key `'post '` had a trailing space | Low | `results['post']` → `KeyError` | canonical `'post'` + back-compat `'post '` alias |

---

## 1. `plot_attgt()` crashes for a single group

**File:** `csdid/att_gt.py` · **Commit:** `72a20f4`

`plot_attgt(group=[g])` raised `TypeError: 'Axes' object is not subscriptable`.
`plt.subplots(nrows=len(group), ncols=1)` returns a **bare `Axes`** when
`len(group) == 1` (not an array), so the `axes[i]` access in the loop failed.

```python
# before
fig, axes = plt.subplots(nrows=len(group), ncols=1, figsize=(10, 5))
...
ax = axes[i]                      # TypeError when len(group) == 1

# after
fig, axes = plt.subplots(nrows=len(group), ncols=1, figsize=(10, 5), squeeze=False)
...
ax = axes[i, 0]                   # axes is always a 2-D array
```

Regression test: `test_ggdid.py::test_plot_attgt_with_single_group`.

---

## 2. Panel path: post-treatment indicator always 0

**File:** `csdid/attgt_fnc/compute_att_gt.py` · **Commit:** `5cf446a`

In the **panel** branch, the result was appended without the post-treatment
flag, so the `'post '` column returned for panel data was **always 0** (the
default). Aggregations that rely on the post indicator were therefore fed wrong
values for panel datasets.

```python
# before
add_att_data(att_gt, inf_f=inf_zeros)            # pst defaults to 0

# after
add_att_data(att_gt, pst=post_treat, inf_f=inf_zeros)
```

(The repeated-cross-section branch already passed `pst=post_treat`; only the
panel branch was affected.)

---

## 3. Panel path overwrote the global unit count `n`

**File:** `csdid/attgt_fnc/compute_att_gt.py` · **Commit:** `de8e0b5`

Inside the (g,t) loop the global unit count `n = dp['n']` was reassigned to the
per-cell `len(disdat)`. The influence-function vector is sized and scaled by
`n`, so a per-cell `n` produced cell-specific lengths/scaling. This was *latent*
for balanced panels (`len(disdat) == n`) but corrupted the influence functions —
and could crash `np.vstack(inf_func)` — for any panel where a cell had a
different unit count (e.g. unbalanced two-period slices).

```python
# before
n = len(disdat)                     # clobbers the global unit count
...
att_inf = n / n1 * att_inf_func

# after
n_cell = len(disdat)                # cell-local; global n is preserved
...
att_inf = n_cell / n1 * att_inf_func
inf_zeros = np.zeros(n)             # still the global n
```

---

## 4. `allow_unbalanced_panel=False` branch crashed immediately

**File:** `csdid/attgt_fnc/preprocess_did.py` · **Commit:** this change

The "force a balanced panel" branch (used when `allow_unbalanced_panel=False`)
contained three Python typos that made it crash before it could run. Because the
default is `allow_unbalanced_panel=True`, the branch was rarely exercised and the
bugs survived.

```python
# before
n = len(data[idname].unique)                 # len() of a bound method -> TypeError
print(n)                                     # stray debug print
n_keep = len(data.iloc[keepers, idname].unique())   # .iloc with a string column -> error
if len(data.loc[keepers] < len(data)):       # parens misplaced -> always truthy

# after
n = len(data[idname].unique())
n_keep = len(data.loc[keepers, idname].unique())
if len(data.loc[keepers]) < len(data):
```

Regression tests:
`test_review_fixes.py::TestUnbalancedPanel::test_allow_unbalanced_false_runs`
and `test_allow_unbalanced_false_balances` (the path now runs and, on
already-balanced data, matches the default to 1e-10).

---

## 5. `raise "<string>"` raised `TypeError`, masking the intended message

**File:** `csdid/attgt_fnc/preprocess_did.py` · **Commit:** this change

Three error paths used `raise "<message>"`. In Python 3 you cannot raise a
`str` — this produces `TypeError: exceptions must derive from BaseException`,
hiding the helpful message the author intended.

```python
# before
raise "All observations dropped due to missing data problems."
raise f"No valid groups. The variable in '{gname}' should be expressed as ..."
raise "Never-treated group is too small, try setting control_group='notyettreated'."

# after
raise ValueError("All observations dropped due to missing data problems.")
raise ValueError(f"No valid groups. The variable in '{gname}' should be ...")
raise ValueError("Never-treated group is too small, try setting control_group='notyettreated'.")
```

---

## 6. Results key `'post '` had a trailing space

**File:** `csdid/attgt_fnc/compute_att_gt.py`, `compute_att_gt2.py`, `csdid/att_gt.py`

The per-(g,t) results dict stored the post-treatment indicator under the key
`'post '` — **with a trailing space**. Anyone indexing the natural `'post'`
(no space) on `ATTgt(...).fit().results` got a `KeyError`. The trailing space
was an unintentional typo, not an API choice.

```python
# before  (producers)
output = {'group': group, 'year': year, "att": att_est, 'post ': post_array}
res.results['post']      # KeyError: 'post'

# after
output = {'group': group, 'year': year, "att": att_est, 'post': post_array}
# att_gt.py also keeps a backward-compatible alias so old code keeps working:
result['post '] = result['post']
res.results['post']      # ok
res.results['post ']     # still ok (deprecated alias)
```

`summ_attgt()` was hardened to select its display columns by name (instead of a
positional `drop('c')` + rename), so the extra alias key can never shift the
summary table. Regression tests:
`test_review_fixes.py::TestPostKeyAlias`.

---

_All fixes are covered by the suite (`csdid/test_csdid`, 561 passed)._
