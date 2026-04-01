library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(tidyr)
library(ggplot2)

insert_size_data <- list.files(pattern = 'stats$') |>
  tibble(file = _) |>
  mutate(sample_id = str_remove(file, '\\.stats$')) |>
  rowwise(sample_id) |>
  reframe(line = read_lines(file)) |>
  filter(str_detect(line, '^IS')) |>
  separate(line, into = c('IS', 'insert_size', 'pairs_total',
                          'inward_pairs', 'outward_pairs',
                          'other_pairs'),
           sep = '\t', convert = TRUE) |>
  select(-IS) |>
  filter(insert_size < 1000, pairs_total > 0) |>
  mutate(sample_id = fct_reorder(sample_id, insert_size,
                                 .fun = function(x, w) weighted.mean(x, w),
                                 w = pairs_total))

# Extract number of samples
n_samples <- n_distinct(insert_size_data$sample_id)

# Define height scaling function
calc_height <- function(n, base = 1.5, per_sample = 0.35, max_height = 20) {
  min(base + n * per_sample, max_height)
}

insert_size_violin <- ggplot(insert_size_data,
                             aes(y = sample_id,
                                 x = insert_size,
                                 weight = pairs_total)) +
  geom_violin() +
  scale_x_continuous(labels = scales::comma_format()) +
  labs(x = 'Insert Size (< 1,000)',
       y = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))

ggsave(
  filename = 'insert_size_violin.png',
  plot = insert_size_violin,
  height = calc_height(n_samples),
  width = 5
)