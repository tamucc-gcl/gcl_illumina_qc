library(tidyverse)

insert_data <- list.files(path = '~/../Downloads/insert_hist',
           full.names = TRUE) %>%
  read_delim(col_names = c('insert_size', 'number'),
             id = 'sample_id',
             show_col_types = FALSE) %>%
  mutate(sample_id = basename(sample_id) %>%
           str_remove('.insert_hist.tsv')) %>%
  filter(number > 0) %>%
  filter(insert_size < 1000)

mean_sizes <- dplyr::summarize(insert_data, 
                               total_n = sum(number),
                               mean_size = sum(insert_size * number) / sum(number),
                               sd_size = sqrt(sum(number * (insert_size - mean_size)^2) / (total_n - 1)),
                               .by = sample_id)

insert_data %>%
  ggplot(aes(x = insert_size)) +
  geom_histogram(aes(weight = number)) +
  geom_vline(data = mean_sizes,
             aes(xintercept = mean_size + -1 * sd_size),
             linetype = 'dashed') +
  geom_vline(data = mean_sizes,
             aes(xintercept = mean_size + 1 * sd_size),
             linetype = 'dashed') +
  geom_vline(data = mean_sizes,
             aes(xintercept = mean_size)) +
  geom_label(data = mean_sizes,
            aes(label = str_c("Mean Insert Size\n", scales::comma(mean_size)),
                x = mean_size + 1 * sd_size), 
            y = Inf, vjust = 1.5,
            hjust = 0, border.colour = 'white') +
  theme(x = 'Insert Size',
        y = 'Number') +
  facet_wrap(~sample_id,
             scales = 'free_y') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())



