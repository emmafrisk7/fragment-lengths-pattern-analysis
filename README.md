

# Fragment Length Analysis with PERMANOVA

An R script that compares fragment length distributions between two sample groups, feature by feature, using PERMANOVA (Bray-Curtis distance).

For each feature (e.g. a gene or region), it tests whether the distribution of fragment lengths differs between two groups (default: `OLD_NO` vs `OLD_YES`).

1. Reads all `*_counts.txt` files and merges them.
2. Joins group information from `sample_metadata.csv`.
3. Computes proportions:
   - `prop_sample`: `count` / total reads in the sample.
   - `prop_length`: `count` / total reads for the feature in the sample (proportion per fragment length).
4. Filters features: both groups must have at least `min_samples` samples with more than `min_reads` total reads for the feature.
5. For each feature (`run_permanova`):
   - builds a matrix with samples as rows and fragment lengths as columns (`prop_length`, missing values = 0),
   - computes Bray-Curtis distances,
   - runs PERMANOVA (`vegan::adonis2`, `dist ~ group`),
   - tests for differences in spread between groups (`vegan::betadisper` + `permutest`).
6. Adjusts p-values for multiple testing with Benjamini-Hochberg (BH), separately for PERMANOVA and dispersion.
7. Writes one CSV per comparison.

## Input

### Count files: `*_counts.txt`

Place these in the working directory. Whitespace-separated, **no header row**, four columns:

| Column | Content |
|---|---|
| 1 | `sample` (file name stem, optionally ending in `_lengths.txt`, which is stripped) |
| 2 | `feature` |
| 3 | `length` (fragment length) |
| 4 | `count` (number of reads) |

### Metadata: `sample_metadata.csv`

Must contain at least:

| Column | Content |
|---|---|
| `sample` | Same names as in the count files (after `_lengths.txt` is stripped), e.g. `AAG-5_S1` |
| `age` | e.g. `YOUNG` / `OLD` |
| `condition` | e.g. `YES` / `NO` |

The group is built as `paste(age, condition, sep = "_")`, e.g. `OLD_YES`.

## Dependencies

R packages: `dplyr`, `tidyr`, `stringr`, `purrr`, `vegan`, `tibble`, `readr`.

```r
install.packages(c("dplyr", "tidyr", "stringr", "purrr", "vegan", "tibble", "readr"))
```

## Settings

Set at the top of the script:

| Variable | Default | Meaning |
|---|---|---|
| `metadata_file` | `"sample_metadata.csv"` | Path to the metadata file |
| `output_dir` | `"results"` | Intended output folder (see *Known limitations*) |
| `min_reads` | `20` | Minimum reads (strictly greater than) for a sample to count |
| `min_samples` | `3` | Minimum number of samples per group that must pass `min_reads` |
| `selected_groups` | `c("OLD_NO", "OLD_YES")` | Groups to compare (all pairs are tested) |

Arguments to `run_permanova()`:

| Argument | Default | Meaning |
|---|---|---|
| `n_perm` | `9999` | Permutations for PERMANOVA |
| `n_disp_perm` | `999` | Permutations for the dispersion test |

## Running the analysis

1. Put the `*_counts.txt` files and `sample_metadata.csv` in the working directory.
2. Check the settings at the top of the script.
3. Make sure `run_permanova()` is defined **before** the loop that calls it.
4. Run the whole script.

The random seed (`set.seed(123)`) is set before the loop, so results are reproducible given the same input, R version and package versions.

## Output

One file per comparison: `permanova_<group1>_vs_<group2>.csv`, sorted by `p_adj`.

| Column | Meaning |
|---|---|
| `feature` | Feature name |
| `F` | Pseudo-F from PERMANOVA |
| `R2` | Proportion of variance explained by group |
| `p_value` | Permutation p-value, PERMANOVA |
| `dispersion_F` | F from the dispersion test |
| `dispersion_p` | p-value, dispersion test |
| `n_samples` | Number of samples included |
| `comparison` | Name of the comparison |
| `p_adj` | BH-adjusted `p_value` |
| `dispersion_p_adj` | BH-adjusted `dispersion_p` |

## Interpretation

- A significant PERMANOVA result means the distributions differ between groups, either in location (centroid) or in spread.
- If `dispersion_p_adj` is also significant, the PERMANOVA result may reflect differences in spread rather than location, so interpret it with caution.
- `R2` indicates effect size. A small p-value with a small `R2` is a weak but detectable difference.

## Excluded features

`run_permanova()` returns `NULL` (the feature is skipped) if any of the following hold:

- fewer than 4 samples,
- fewer than 2 fragment lengths,
- not exactly 2 groups,
- either group has fewer than 2 samples,
- any sample has only zeros.

## Known limitations

- **Missing combinations are not filled with 0.** A sample that is entirely absent for a feature is left out of that feature's matrix rather than counted as 0 reads.
- **Sample names must match exactly** between the count files and the metadata. Otherwise `sample_group` becomes `NA` and the sample is filtered out. Check with `sum(is.na(merged_counts$sample_group))`.
- **`output_dir` is currently unused.** CSV files are written to the working directory. Create the folder and use `file.path(output_dir, ...)` in `write_csv()` to save into `results/`.
- **Empty result table:** if no features pass validation, `p.adjust` can fail. Add `if (nrow(results) == 0) next` after `map_dfr`.
- **Duplicate rows** (the same sample, feature and length more than once) produce list-columns in `pivot_wider`. Add `values_fn = sum` if this can happen.
- Permutation tests are random: changing `set.seed()` or the number of permutations will shift p-values slightly.

## Suggested methods text

> Fragment length distributions per feature were compared between groups using PERMANOVA (`vegan::adonis2`, 9999 permutations) on Bray-Curtis distances computed from per-length proportions. Homogeneity of multivariate dispersion was tested with `betadisper` and `permutest` (999 permutations). Features were included if at least three samples per group had more than 20 reads. P-values were adjusted with Benjamini-Hochberg. Random seed: 123.
