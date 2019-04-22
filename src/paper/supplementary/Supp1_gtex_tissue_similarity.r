rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)

# utils functions
calc_jaccard_index = function(tissue_1, tissue_2, ...) {
  nrow(intersect(tissue_1, tissue_2)) /  nrow(union(tissue_1, tissue_2))
}


# load gtex data
gtex_data = list.files("data/TF_target_sources/inferredGTEx/aracne", full.names = T, 
                       recursive = T, pattern="network.txt") %>%
  map_df(function(f) {
    tissue = f %>%
      str_split(pattern="/") %>%
      pluck(1,5)
    read_delim(f, delim="\t") %>%
      mutate(tissue = tissue)
  })

# calculation of similarity across all tissue combinations
plan(multiprocess)
tissue_combination = gtex_data %>% 
  distinct(tissue) %>%
  mutate(tissue_2 = tissue) %>%
  expand(tissue, tissue_2) %>%
  rename(tissue_1 = tissue) %>%
  unite(comparison, tissue_1, tissue_2, sep = ":", remove = F) %>%
  gather(combination, tissue, -comparison) %>%
  arrange(comparison) %>%
  inner_join(nest(gtex_data, -tissue, .key="regulon"), by="tissue") %>%
  arrange(comparison) %>%
  select(-tissue) %>%
  spread(combination, regulon) %>%
  mutate(jaccard = future_pmap(., .f=calc_jaccard_index,.progress = T)) %>%
  select(comparison, jaccard) %>%
  gather(metric, value, jaccard) %>%
  unnest(value) %>%
  separate(comparison, into=c("tissue_1", "tissue_2"), sep = ":") %>%
  nest(-metric) %>%
  mutate(table = data %>% map(~spread(., tissue_2, value)))

# plot similarity heatmap
hmap = tissue_combination %>% 
  pull(table) %>% 
  pluck(1) %>%
  data.frame(row.names=1, check.names = F) %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n=9, name ="BuGn")))(100),
           cellwidth = 12, cellheight = 12, cluster_rows = T, cluster_cols = T, 
           fontsize = 12, treeheight_row = 15, treeheight_col = 15)



### number of samples per gtex tissue
s = list.files("data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix", 
               full.names = T) %>%
  map_df(function(i) {
    tissue = i %>%
      str_split(pattern="/") %>%
      pluck(1,5) %>%
      str_split(pattern = "[.]") %>%
      pluck(1,1)
    num_samples = ncol(read_delim(i, delim ="\t")) - 1
    tibble(tissue = tissue, num_samples = num_samples)
  })

num_samples = s %>% 
  ggplot(aes(x=fct_reorder(tissue, -num_samples), y=num_samples)) +
  geom_col() +
  coord_flip() +
  labs(x="Tissue", y="#Samples") +
  ylim(c(0,1500)) +
  geom_text(aes(label=num_samples), hjust=-0.2) +
  theme(aspect.ratio = c(1), plot.margin=grid::unit(c(0,0,0,0), "mm"))

### save plot
figs1 = plot_grid(hmap$gtable, num_samples, 
                  labels = c("A", "B"), 
                  ncol = 1)

ggsave(filename = "paper/figures/supplementary/FigureS1.png", 
       plot = figs1, width = 10, height = 12)

