# permaBrioche
<table>
  <tr>
    <td style="vertical-align: top; padding-right: 20px;">

<!-- badges: start -->
<!-- badges: end -->

**permaBrioche** is a statistical framework for **confounding‑aware PERMANOVA**
and **interpretable distance‑based effect‑size estimation** in repeated‑measures
and other blocked study designs. It is particularly motivated by microbiome and
other high‑dimensional biological data, where subject‑level clustering and
longitudinal sampling are common.

permaBrioche addresses two well‑known limitations of standard PERMANOVA:

1. **Invalid permutation schemes** under subject‑level confounding  
   (e.g., longitudinal designs)
2. **Upward bias and poor interpretability of the PERMANOVA $R^2$ effect size**

The package implements:

- **Design‑aware permutation schemes** for invariant and variant covariates
- A **null‑centered $R^2$** for bias‑corrected variance explanation
- A **Hájek‑based distance effect size** with direct geometric interpretation
- Optional **location–dispersion decomposition** in the Euclidean case

    </td>
    <td style="vertical-align: top;">
      <img src="man/permaBriocheLogo.png" alt="permaBrioche logo" width="700">
    </td>
  </tr>
</table>

---

## Citation

This manuscript is currently under review.  
Software is available via the BioBakery ecosystem as **permaBrioche**.

---

## Inputs

**permaBrioche** requires two aligned objects:

### (a) Sample‑level metadata

A `data.frame` with **one row per sample**, containing at least:

*   `subject` – repeated‑measures (blocking) ID
*   `exposure` – the covariate of interest (e.g. binary treatment status or multi-level categorical factor)
*   **Sample Identifiers** provided in one of two ways:
    *   row names that uniquely identify samples and are used when constructing the feature matrix, or
    *   a column containing unique sample IDs (specified via the `sample_id` argument, e.g. `sample_id = "sample"`)

Each subject may appear in multiple rows.

### (b) Feature matrix (`bugs`)

A matrix of features (e.g., taxa abundances):

*   rows = samples
    *   (row names must use the same sample identifiers as the metadata)
*   columns = features

***

With these two objects, we can create the distance matrix for **permaBrioche**:

### (c) Distance matrix

Use vegdist to get a distance matrix with a chosen dissimilarity method:

```{r}
library(vegan)
dist_bc <- vegdist(bugs, method = "bray")
```
---

## Outputs

Depending on the method used, permaBrioche reports:

- **PERMANOVA results**
  - pseudo‑ $F$ statistic
  - empirical permutation *p*‑value
  - optional null‑centered $R^2$

- **Hájek distance‑based effect size**
  - effect on the distance scale
  - permutation‑based uncertainty
  - optional decomposition into:
    - **location** (mean / centroid shift)
    - **dispersion** (change in variability)

All outputs are **valid under repeated measures** and designed to be
**interpretable across studies**.

---

## Installation

permaBrioche is currently distributed via GitHub as part of the BioBakery
ecosystem. To install the package (including vignettes):

```r
library(devtools)

devtools::install_github(
  "biobakery/permabrioche",
  build_vignettes = TRUE
)
````

***

## Tutorial

**Read the vignette:**  
`browseVignettes("permabrioche")`

The vignette contains:

*   fully runnable examples
*   detailed explanations of permutation schemes
*   interpretation guidance for all effect sizes

The examples below are illustrative only.

***

## Example: Repeated‑measures PERMANOVA

```r
res_perm <- PERMANOVA_repeat_measures(
  formula           = dist_bc ~ exposure,
  data              = meta,
  sample_id         = "sample",
  blocking_variable = "subject",
  permutations      = 999
)
```

This function:

*   permutes **variant covariates within subjects**
*   permutes **invariant covariates across subjects**
*   returns valid inference under repeated measures

***

## Example: Null‑centered $R^2$

```r
res_perm_centered <- PERMANOVA_repeat_measures(
  formula           = dist_bc ~ exposure,
  data              = meta,
  sample_id         = "sample",
  blocking_variable = "subject",
  permutations      = 999,
  center_R2         = TRUE
)
```

**Interpretation**

> **Null‑centered $R^2$** = variation explained *beyond chance under the study
> design*

A value near zero indicates no meaningful effect.

***

## Example: Hájek distance‑based effect size

```r
res_hajek <- hajek_repeat_measures(
  formula           = dist_bc ~ exposure,
  data              = meta,
  bugs              = bugs,
  sample_id         = "sample",
  blocking_variable = "subject",
  covariate_name    = "exposure",
  permutations      = 999,
  method            = "bray"
)
```

Output:

*   `observed`: distance‑scale effect size
*   `pval`: permutation *p*‑value

This effect size is **directly interpretable**, unlike variance‑based summaries.

***

## Location–dispersion decomposition

In the **Euclidean case**, the Hájek effect decomposes as:

$$
\tau = \tau_{\text{location}} + \tau_{\text{dispersion}}
$$

*   **location**: mean (centroid) shift
*   **dispersion**: change in within‑group variability

This decomposition distinguishes systematic shifts from increased heterogeneity.

***

## When to use which method (permaBrioche function)

| Goal                      | Recommended method                        |
| ------------------------- | ----------------------------------------- |
| Hypothesis testing        | `PERMANOVA_repeat_measures()`             |
| Variance summary          | PERMANOVA + null‑centered $R^2$           |
| Interpretable effect size | `hajek_repeat_measures()`                 |
| Mean vs variability       | Hájek + location–dispersion decomposition |

***

## Support

*   GitHub issues: <https://github.com/biobakery/permabrioche>
*   BioBakery tools: <https://github.com/biobakery>

***
```
