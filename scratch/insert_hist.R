library(tidyverse)

insert_data <- list.files(path = '~/../Downloads/insert_hist',
           full.names = TRUE) %>%
  read_delim(col_names = c('insert_size', 'number'),
             id = 'sample_id',
             show_col_types = FALSE) %>%
  mutate(sample_id = basename(sample_id) %>%
           str_remove('.insert_hist.tsv')) %>%
  filter(number > 0,
         insert_size < 1000)

insert_data %>%
  ggplot(aes(x = insert_size, 
             weight = number)) +
  geom_histogram() +
  theme(x = 'Insert Size',
        y = 'Number')


dplyr::summarize(insert_data, mean_size = sum(insert_size * number) / sum(number))

read_delim('~/../Downloads/mean_insert_sizes.tsv',
           show_col_types = FALSE) %>%
  ggplot(aes(x = mean_insert_size)) +
  geom_histogram()
