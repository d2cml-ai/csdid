import warnings

import numpy as np
from scipy.stats import mstats, norm
from joblib import Parallel, delayed
import pandas as pd

from csdid.utils.bmisc import multiplier_bootstrap


def _cluster_sums_by_label(inf_func, cluster_labels):
    """Group IF rows once, retaining the original within-cluster sum order."""
    cluster_ids = np.unique(cluster_labels)
    cluster_if = np.zeros((len(cluster_ids), inf_func.shape[1]))
    # With few clusters, the original masks cost less than building sorted
    # groups. Object arrays also retain the old route: sorting object NaNs with
    # and without inverse codes can produce different unique-label sequences.
    if len(cluster_ids) <= 16 or np.asarray(cluster_labels).dtype.kind == "O":
        for i, label in enumerate(cluster_ids):
            cluster_if[i] = inf_func[cluster_labels == label].sum(axis=0)
        return cluster_if

    _, codes = np.unique(cluster_labels, return_inverse=True)
    order = np.argsort(codes, kind="stable")
    ends = np.cumsum(np.bincount(codes, minlength=len(cluster_ids)))
    start = 0
    for i, (label, end) in enumerate(zip(cluster_ids, ends)):
        rows = order[start:end]
        # Match the former equality mask, including its empty NaN/NaT groups.
        # Checking only this group's rows avoids rescanning all observations.
        rows = rows[cluster_labels[rows] == label]
        cluster_if[i] = inf_func[rows].sum(axis=0)
        start = end
    return cluster_if


def mboot(inf_func, DIDparams, pl=False, cores=1):
    # Setup needed variables
    data            = DIDparams['data'] 
    idname          = DIDparams['idname']
    clustervars     = DIDparams['clustervars']
    biters          = DIDparams['biters']
    tname           = DIDparams['tname']
    try:
        tlist           = np.sort(data[tname].unique())
    except Exception:
        tlist           = np.sort(data[tname].unique().to_numpy())
    alp             = DIDparams['alp']
    panel           = DIDparams['panel']
    true_repeated_cross_sections = DIDparams['true_repeated_cross_sections']

    # Get n observations (for clustering below)
    if panel:
        dta = data[ data[tname] == tlist[0] ]
    else:
        dta = data.copy()

    # Convert inf_func to matrix
    inf_func = np.asarray(inf_func)
    if inf_func.ndim == 1:
        inf_func = inf_func[:, np.newaxis]

    # Set correct number of units
    n = inf_func.shape[0]

    # Normalize clustervars: can be None, a string, or a list with one string
    if clustervars is not None:
        if isinstance(clustervars, str):
            clustervars = clustervars  # already a string
        elif isinstance(clustervars, (list, tuple)):
            if len(clustervars) == 0:
                clustervars = None
            elif len(clustervars) > 1:
                raise ValueError("Can't handle more than one cluster variable beyond idname.")
            else:
                clustervars = clustervars[0]
        # Drop if same as idname (clustering at unit level = no extra clustering)
        if clustervars == idname:
            clustervars = None

    # Validate cluster variable
    if clustervars is not None:
        if clustervars not in data.columns:
            # v7-DIV1: match R `did` (graceful degradation) -- warn and proceed with
            # non-clustered SEs rather than raising. R emits "reporting standard
            # errors that do NOT account for clustering" and falls back to the
            # unclustered bootstrap.
            warnings.warn(
                f"Cluster variable '{clustervars}' not found in data; reporting "
                f"standard errors that do NOT account for clustering."
            )
            clustervars = None
    if clustervars is not None:
        # Check that cluster variable does not vary over time within unit
        # Must check full data (not just first period) to catch time-varying clusters
        clust_tv = data.groupby(idname)[clustervars].nunique()
        if (clust_tv > 1).any():
            raise ValueError(
                f"Cluster variable '{clustervars}' varies over time within units. "
                f"Cluster variables must be time-invariant."
            )

    # Multiplier bootstrap
    n_clusters = n
    if clustervars is None:
        # Standard (non-clustered) bootstrap
        bres = np.sqrt(n) * run_multiplier_bootstrap(inf_func, biters, pl, cores)
    else:
        # Clustered bootstrap: draw one multiplier per cluster,
        # apply to cluster sums of the influence function (Callaway & Sant'Anna 2021, Remark 10)
        cluster_vec = dta[[idname, clustervars]].drop_duplicates().set_index(idname)[clustervars]
        # Align cluster assignments with IF rows (which are in unit order)
        unit_ids = dta[idname].unique()
        cluster_labels = cluster_vec.reindex(unit_ids).values
        
        # Aggregate IF to cluster sums
        cluster_if = _cluster_sums_by_label(inf_func, cluster_labels)
        n_clusters = len(cluster_if)
        
        bres = np.sqrt(n_clusters) * run_multiplier_bootstrap(cluster_if, biters, pl, cores)

    # Handle vector and matrix case differently to get nxk matrix
    if isinstance(bres, np.ndarray) and bres.ndim == 1:
        bres = bres.reshape(int(biters), -1)
    elif isinstance(bres, np.ndarray) and bres.ndim > 2:
        bres = bres.transpose()

    # Non-degenerate dimensions
    # A positive bootstrap dimension remains usable in any outcome units.
    # Avoid squaring here, which can underflow for small valid dimensions.
    ndg_dim = np.all(np.isfinite(bres), axis=0) & np.any(bres != 0, axis=0)
    bres = bres[:, ndg_dim]

    # All dimensions degenerate: no finite bootstrap statistics (matches R's ncol==0 guard).
    if bres.shape[1] == 0:
        se = np.full(ndg_dim.shape, np.nan)
        return {'bres': bres, 'V': np.nan, 'se': se, 'crit_val': np.nan}

    # Bootstrap variance matrix (this matrix can be defective because of degenerate cases)
    if len(bres) > 1:
        V = np.cov(bres, rowvar=False)
    else:
        V = np.full((bres.shape[1], bres.shape[1]), np.nan).squeeze()

    # Bootstrap standard error
    quantile_75 = np.quantile(bres, 0.75, axis=0, method="inverted_cdf")
    quantile_25 = np.quantile(bres, 0.25, axis=0, method="inverted_cdf")
    qnorm_75 = norm.ppf(0.75)
    qnorm_25 = norm.ppf(0.25)
    bSigma = (quantile_75 - quantile_25) / (qnorm_75 - qnorm_25)
    bSigma[~np.isfinite(bSigma) | (bSigma <= 0)] = np.nan

    # Critical value for uniform confidence band
    # O7-F1: biters=1 -> 25th and 75th percentiles coincide -> bSigma=0 -> bres/bSigma
    # all non-finite -> bT empty after isfinite filter -> np.quantile on empty array
    # crashes (IndexError). Suppress the divide-by-zero and guard the empty bT case so
    # the function degrades gracefully (crit_val=NaN), matching R's mboot and the
    # existing all-degenerate branch's contract, instead of crashing.
    with np.errstate(divide='ignore', invalid='ignore'):
        scaled = np.abs(bres / bSigma)
    scaled[np.isnan(scaled)] = -np.inf
    bT = np.max(scaled, axis=1)
    bT = bT[np.isfinite(bT)]
    if bT.size == 0:
        se = np.full(ndg_dim.shape, np.nan)
        se[ndg_dim] = bSigma * np.sqrt(n_clusters) / n
        return {'bres': bres, 'V': V, 'se': se, 'crit_val': np.nan}
    crit_val = np.quantile(bT, 1 - alp, method="inverted_cdf")
    
    # Standard error: R uses bSigma * sqrt(n_clusters) / n
    # For non-clustered (n_clusters == n), this equals bSigma / sqrt(n).
    # For clustered, the sqrt(n_clusters)/n scaling correctly accounts for
    # cluster-sum influence functions.
    se = np.full(ndg_dim.shape, np.nan)
    se[ndg_dim] = bSigma * np.sqrt(n_clusters) / n

    return {'bres': bres, 'V': V, 'se': se, 'crit_val': crit_val}

def run_multiplier_bootstrap(inf_func, biters, pl=False, cores=1):
    biters = int(biters)
    cores = max(1, int(cores))
    # Split biters into per-core chunks that always sum to biters exactly,
    # matching R's diff(round(seq(0, biters, length.out = cores + 1)))
    boundaries = np.round(np.linspace(0, biters, cores + 1)).astype(int)
    chunks = np.diff(boundaries)
    chunks = chunks[chunks > 0]  # drop empty chunks

    n = inf_func.shape[0]

    if n > 2500 and pl and cores > 1:
        # v10-F1: make the parallel path REPRODUCIBLE under np.random.seed while
        # keeping workers INDEPENDENT. loky workers do not inherit the parent's
        # global RNG, so we derive a deterministic per-chunk seed from the parent
        # global state: draw a uint32 entropy value off np.random (pinned by
        # np.random.seed), build a master SeedSequence, and spawn one independent
        # child SeedSequence per chunk. Distinct spawned sequences => independent
        # draw streams (no variance collapse); deterministic derivation => same
        # seed => same draws across runs.
        entropy = int(np.random.randint(0, 2 ** 31 - 1))
        live_chunks = [b for b in chunks if b > 0]
        children = np.random.SeedSequence(entropy).spawn(len(live_chunks))

        def parallel_function(b, ss):
            return multiplier_bootstrap(inf_func, b, seed=ss)

        results = Parallel(n_jobs=cores)(
            delayed(parallel_function)(b, ss) for b, ss in zip(live_chunks, children)
        )
        results = np.vstack(results)
    else:
        results = multiplier_bootstrap(inf_func, biters)

    return results



