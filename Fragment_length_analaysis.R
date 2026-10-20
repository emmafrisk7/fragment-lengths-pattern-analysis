#Package to
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(vegan)
library(tibble)
library(readr)



metadata_file <- "sample_metadata.csv"
output_dir <- "results"
min_reads <- 20
min_samples <- 3
selected_groups <- c("OLD_NO", "OLD_YES")


metadata <- read_csv(metadata_file,show_col_types = FALSE) %>%mutate(sample = as.character(sample),sample_group = paste(age, condition, sep = "_"))

files <- list.files(
  pattern = "_counts\\.txt$"
)

merged_counts <- map_dfr(files, function(f) {read.table(
    f,
    header = FALSE,
    col.names = c("sample", "feature", "length", "count")
  )}) %>%
  mutate(
    sample = str_remove(sample, "_lengths\\.txt$")
  )

merged_counts <- merged_counts %>%
  left_join(
    metadata %>% select(sample, sample_group),
    by = "sample"
  )

merged_counts <- merged_counts %>% group_by(sample) %>% mutate(
    total_reads_sample = sum(count),
    prop_sample = if_else(
      total_reads_sample > 0,
      count / total_reads_sample,
      0)) %>% ungroup() %>% group_by(sample, feature) %>% mutate(
    total_reads_feature = sum(count),
    prop_length = if_else(
      total_reads_feature > 0,
      count / total_reads_feature, 0) ) %>% ungroup()


comparisons <- combn(selected_groups,2,simplify = FALSE)

all_results <- list()

set.seed(123)

for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  comparison_name <- paste0(group1, "_vs_", group2)
  message("Running: ", comparison_name)
  df <- merged_counts %>% filter(sample_group %in% comparison)
  
  coverage <- df %>%
    distinct(feature, sample, sample_group, total_reads_feature) %>%
    group_by(feature, sample_group) %>%
    summarise(n_pass = sum(total_reads_feature > min_reads), .groups = "drop")
  
  # Behåll features där BÅDA grupperna har minst min_samples samples som klarar kravet
  valid_features <- coverage %>%
    group_by(feature) %>%
    filter(n_distinct(sample_group) == 2, all(n_pass >= min_samples)) %>%
    ungroup() %>%
    pull(feature) %>%
    unique()
  
  df <- df %>% filter(feature %in% valid_features)
  
  if (nrow(df) == 0) {
    message("No features passed filtering.")
    next
  }
  
  # Run PERMANOVA
  results <- map_dfr(unique(df$feature),run_permanova,dat = df)

  
  # Multiple testing correction
  results <- results %>%mutate(
      comparison = comparison_name,
      p_adj = p.adjust(
      p_value,
      method = "BH"),
      dispersion_p_adj = p.adjust(
      dispersion_p,
      method = "BH"
      )
    ) %>%arrange(p_adj)
  
# Export individual comparison
write_csv(results,paste0("permanova_",comparison_name,".csv"))
all_results[[comparison_name]] <- results}
message("Analysis complete!")

