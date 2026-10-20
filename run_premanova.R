run_permanova <- function(feat, dat, n_perm = 9999, n_disp_perm = 999) {
  
  wide <- dat %>%
    dplyr::filter(feature == feat) %>%
    tidyr::pivot_wider(
      id_cols     = c(sample, sample_group),
      names_from  = length,
      values_from = prop_length,
      values_fill = 0
    ) %>%
    dplyr::arrange(sample)
  
  group <- droplevels(factor(wide$sample_group))
  mat   <- as.matrix(dplyr::select(wide, -sample, -sample_group))
  
  # Skip features that can't support a two-group test
  if (nrow(mat) < 4 || ncol(mat) < 2 ||
      nlevels(group) != 2 || any(table(group) < 2) ||
      any(rowSums(mat) == 0)) {
    return(NULL)
  }
  
  d     <- vegan::vegdist(mat, method = "bray")
  perm  <- vegan::adonis2(d ~ group, data = data.frame(group),
                          permutations = n_perm)
  disp  <- vegan::permutest(vegan::betadisper(d, group),
                            permutations = n_disp_perm)
  
  tibble::tibble(
    feature      = feat,
    F            = perm$F[1],
    R2           = perm$R2[1],
    p_value      = perm$`Pr(>F)`[1],
    dispersion_F = disp$tab$F[1],
    dispersion_p = disp$tab$`Pr(>F)`[1],
    n_samples    = nrow(mat)
  )
}