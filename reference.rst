Reference
=========

This page documents every public function and method exposed by the
**csdid** package along with the meaning of every argument. Use it to
look up where an option goes (constructor, ``fit``, ``aggte``, etc.) and
what the default behavior is.

The package is organized around the ``ATTgt`` class. A typical workflow
is:

1. Build an estimator object with :class:`csdid.att_gt.ATTgt` —
   pass in the data and identifier columns.
2. Call :meth:`csdid.att_gt.ATTgt.fit` — choose the 2x2 DiD estimator
   and the base period.
3. Inspect the disaggregated group-time effects with
   :meth:`csdid.att_gt.ATTgt.summ_attgt` or plot them with
   :meth:`csdid.att_gt.ATTgt.plot_attgt`.
4. Aggregate to a smaller number of parameters with
   :meth:`csdid.att_gt.ATTgt.aggte` and plot with
   :meth:`csdid.att_gt.ATTgt.plot_aggte`.

.. contents::
   :local:
   :depth: 2


``csdid.att_gt.ATTgt``
----------------------

.. py:class:: ATTgt(yname, tname, idname, gname, data, control_group=['nevertreated', 'notyettreated'], xformla=None, panel=True, allow_unbalanced_panel=True, clustervar=None, weights_name=None, anticipation=0, cband=False, biters=1000, alp=0.05)

   Group-time average treatment effects estimator. Sets up the design
   from the long-format panel (or repeated cross-section), runs the
   pre-processing step, and stores everything needed for ``fit``.

   :param str yname: Name of the **outcome** column in ``data``.
   :param str tname: Name of the **time** column in ``data``.
   :param str idname: Name of the **unit identifier** column. Required
      for panel data; ignored when ``panel=False``.
   :param str gname: Name of the **group** column. This is the period in
      which a unit is first treated. Never-treated units should have
      ``gname = 0``.
   :param pandas.DataFrame data: The long-format dataset.
   :param control_group: Which units to use as controls. Either
      ``"nevertreated"`` (default) or ``"notyettreated"``. The
      not-yet-treated group is at least as large as the never-treated
      group and changes across time. A list ``['nevertreated',
      'notyettreated']`` resolves to the first entry.
   :type control_group: str or list of str
   :param str xformla: Patsy-style formula for the covariates, e.g.
      ``"lemp ~ lpop"``. The left-hand side is ignored; only the
      right-hand side is used to build the design matrix. ``None`` (or
      ``"y~1"``) requests **unconditional** parallel trends.
   :param bool panel: ``True`` (default) treats the data as a panel.
      ``False`` requests repeated cross-section mode (``idname`` is then
      ignored).
   :param bool allow_unbalanced_panel: When ``True`` (default), the
      package does **not** drop units with missing observations in some
      periods. Set to ``False`` to force a balanced panel.
   :param str clustervar: Name of an extra clustering variable. Up to
      two clustering variables are supported and one of them must be
      ``idname``.
   :param str weights_name: Name of a column with sampling weights. If
      ``None``, equal weights are used.
   :param int anticipation: Number of pre-treatment periods in which
      units are allowed to anticipate the treatment. Default ``0``.
   :param bool cband: When ``True``, store enough information for
      uniform confidence bands. Default ``False``.
   :param int biters: Number of multiplier-bootstrap iterations used for
      standard errors. Default ``1000``.
   :param float alp: Significance level for confidence bands. Default
      ``0.05``.

   .. py:method:: fit(est_method='dr', base_period='varying', bstrap=True)

      Estimate group-time average treatment effects.

      :param str est_method: The 2x2 DiD estimator used for each
         ``(g, t)`` cell. One of:

         * ``"dr"`` (default) — Sant'Anna and Zhao doubly-robust
           estimator.
         * ``"ipw"`` — inverse-probability-weighted estimator.
         * ``"reg"`` — outcome-regression estimator.

      :param str base_period: Reference period used to construct each
         :math:`ATT(g,t)`.

         * ``"varying"`` (default) — the reference period changes with
           ``t``. For pre-treatment periods it is ``t - 1``; for
           post-treatment periods it is the last period before
           treatment for group ``g``.
         * ``"universal"`` — a single, fixed reference period is used
           for each group (the last pre-treatment period, after
           anticipation). The entry at the base period itself is
           normalized to 0. With this option the result table contains
           one extra row per group.

      :param bool bstrap: When ``True`` (default), standard errors and
         critical values come from the multiplier bootstrap. When
         ``False``, asymptotic standard errors and a normal critical
         value (1.96) are used.

      :returns: ``self`` (the call is chainable). After fitting,
         ``self.results``, ``self.MP``, and ``self.did_object`` are
         populated.

   .. py:method:: summ_attgt(n=4)

      Build a tidy DataFrame of the ``ATT(g,t)`` table.

      :param int n: Number of decimals to round to. Default ``4``.
      :returns: ``self`` with the table available as ``self.summary2``.

   .. py:method:: aggte(typec='group', balance_e=None, min_e=-inf, max_e=inf, na_rm=False, bstrap=None, biters=None, cband=None, alp=None, clustervars=None)

      Aggregate the group-time effects into a smaller number of
      parameters. Returns ``self``; the aggregated object is stored on
      ``self.atte``.

      :param str typec: Type of aggregation:

         * ``"simple"`` — weighted average of all post-treatment
           ``ATT(g,t)`` (weights proportional to group size).
         * ``"dynamic"`` — average effects by length of exposure
           (event-study).
         * ``"group"`` (default) — average effect per treatment cohort.
         * ``"calendar"`` — average effect per calendar period.
      :param int balance_e: Only used when ``typec='dynamic'``. If set,
         drops cohorts that are not observed for at least ``balance_e +
         1`` post-treatment periods, so the composition is constant in
         event time. Default ``None`` (no balancing).
      :param float min_e: Smallest event time to include in the dynamic
         aggregation. Default ``-inf``.
      :param float max_e: Largest event time to include. Default
         ``inf``.
      :param bool na_rm: When ``True``, drop missing aggregated
         estimates. Default ``False``.
      :param bool bstrap: Override the ``bstrap`` setting stored on the
         fitted object. ``None`` (default) inherits from ``fit``.
      :param int biters: Override the number of bootstrap iterations.
      :param bool cband: Compute uniform confidence bands (requires
         ``bstrap=True``). ``None`` inherits.
      :param float alp: Override the significance level.
      :param list clustervars: Override the clustering variables.

   .. py:method:: plot_attgt(ylim=None, xlab=None, ylab=None, title='Group', xgap=1, ncol=1, legend=True, group=None, ref_line=0, theming=True, grtitle='Group')

      Faceted plot of the group-time treatment effects, one subplot per
      cohort.

      :param tuple ylim: ``(ymin, ymax)`` for the y-axis.
      :param str xlab: Label for the x-axis.
      :param str ylab: Label for the y-axis.
      :param str title: Subplot title prefix. Default ``"Group"``.
      :param int xgap: Spacing between x-axis ticks. Default ``1``.
      :param int ncol: Number of columns in the facet grid. Default
         ``1``.
      :param bool legend: Whether to draw the legend. Default ``True``.
      :param list group: Subset of groups to plot. ``None`` plots all.
      :param float ref_line: y-value of the reference line. Default
         ``0``.
      :param bool theming: Apply the package's default styling. Default
         ``True``.
      :param str grtitle: Per-subplot title prefix (used together with
         the group value). Default ``"Group"``.

   .. py:method:: plot_aggte(ylim=None, xlab=None, ylab=None, title='', xgap=1, legend=True, ref_line=0, theming=True, **kwargs)

      Plot the aggregated object stored on ``self.atte`` (produced by a
      prior call to :meth:`aggte`). For ``typec='group'`` it draws a
      ``splot``; for the other types it draws a ``gplot``.

      :param tuple ylim: y-axis limits.
      :param str xlab: x-axis label.
      :param str ylab: y-axis label.
      :param str title: Plot title. Empty defaults to ``"Average Effect
         by Group"`` (group aggregation) or ``"Average Effect by Length
         of Exposure"`` (dynamic / calendar).
      :param int xgap: Tick spacing on the x-axis.
      :param bool legend: Show legend.
      :param float ref_line: y-value of the reference line.
      :param bool theming: Apply default styling.



Notes on argument inheritance
-----------------------------

Some arguments can be specified in **more than one place** (for example
``bstrap``, ``biters``, ``cband``, ``alp``). The resolution order is:

1. Value passed to :meth:`ATTgt.aggte` overrides the value set on
   :meth:`ATTgt.fit`, which overrides the value set on
   :class:`ATTgt`.
2. Passing ``None`` to :meth:`ATTgt.aggte` means "inherit from the
   ``MP`` object built during ``fit``".

This matches the R ``did`` package convention.
