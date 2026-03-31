insert_size_violin <- list.files(pattern = 'stats$') |>
  tibble::tibble(file = _) |>
  dplyr::mutate(sample_id = stringr::str_remove(file, '.stats$')) |>
  dplyr::rowwise(sample_id) |>
  dplyr::reframe(line = readr::read_lines(file)) |>
  dplyr::filter(stringr::str_detect(line, '^IS')) |>
  tidyr::separate(line, into = c('IS', 'insert_size', 'pairs_total',
                                 'inward_pairs', 'outward_pairs',
                                 'other_pairs'),
                  sep = '\t', convert = TRUE) |>
  dplyr::select(-IS) |>
  dplyr::filter(insert_size < 1000) |>
  tidyr::uncount(pairs_total) |>
  dplyr::mutate(sample_id = forcats::fct_reorder(sample_id, insert_size)) |>
  ggplot2::ggplot(ggplot2::aes(y = sample_id, 
                               x = insert_size)) +
  ggplot2::geom_violin() +
  ggplot2::scale_x_continuous(labels = scales::comma_format()) +
  ggplot2::labs(x = 'Insert Size (< 1,000)',
                y = NULL) +
  ggplot2::theme_classic() +
  ggplot2::theme(panel.background = ggplot2::element_rect(colour = 'black'))

# Extract number of samples
n_samples <- insert_size_violin$data |>
  dplyr::pull(sample_id) |>
  dplyr::n_distinct()

# Define height scaling function
calc_height <- function(n, base = 1.5, per_sample = 0.35, max_height = 20) {
  height <- base + n * per_sample
  return(min(height, max_height))
}


ggplot2::ggsave(
  filename = 'insert_size_violin.png',
  plot = insert_size_violin,
  height = calc_height(n_samples),
  width = 5
)
