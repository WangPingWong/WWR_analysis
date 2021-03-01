# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
##load packages
pacman::p_load(tidyverse, poppr, RColorBrewer, paletteer, scales)
pacman::p_load_gh("Mikata-Project/ggthemr")
##read pop data from 1000 genome
pop_data <- read_tsv("data/igsr_samples.tsv")

## sum by population
sum_pop <- pop_data %>% count(`Population name`) %>%
  arrange(desc(n))

##Select population of interset (EU and AU)
selected_pops <- pop_data %>%
  filter(`Superpopulation code`=="EUR" | `Population name`=="Australian")
selected_pops %>% count(`Population name`)


# Verify that the selected populations exist in each of the files
link_table <- read_tsv("data/igsr_30x_GRCh38.tsv")
count_sample <- function(x, sample_list){
  # x = link_table$Sample[1]; sample_list=selected_pops$`Sample name`
  samples <- unlist(str_split(x, ","))
  paste(sample_list[sample_list %in% samples], collapse = ",")
}
link_update_table <- link_table %>% select(url, Sample, Population) %>%
  mutate(our_samples=map_chr(Sample,
     ~count_sample(x=.x , sample_list =selected_pops$`Sample name` )))
samples_in_vcf <- str_split(unique(link_update_table$our_samples), ",") %>% unlist()
selected_pops <- selected_pops %>%
  mutate(In_vcf=`Sample name` %in% samples_in_vcf)
selected_pops %>% filter(grepl("aus", `Population name`, ignore.case = TRUE)) %>% 
  select(`Sample name`, `Population name`, `Data collections`,  In_vcf)

write_tsv(selected_pops, "EUR_AUS_1000G.tsv")

# PLINK PCA #####
# read in the eigenvectors, produced in PLINK
eigenvec <- read.table('data/plink.eigenvec', header = FALSE, skip=0, sep = ' ',
                       stringsAsFactors = FALSE) 
PCA_data <- eigenvec %>% .[2:ncol(.)] %>% as_tibble() %>% 
  setNames(c("Sample name", paste0("PC", 1:20))) %>% 
  right_join(pop_data, .) %>% mutate(Superpopulation=sub(",.+\\(SGDP\\)", "", sub(".+\\(SGDP\\),", "", `Superpopulation name`)),
    across(!starts_with("PC"), ~replace_na(.x, "Unknown")),
  across(contains("population"), 
         ~ifelse(grepl("UWW1", `Sample name`), "UWW1", .)))
# check
PCA_data %>% filter(grepl("UWW1", `Sample name`))

# levels of populations
unique(PCA_data$Superpopulation) #%>% sub(",.+\\(SGDP\\)", "", .) %>% sub(".+\\(SGDP\\),", "", .)


pca_var <- eigenvec[3:ncol(eigenvec)] %>% map_dbl(~var(.x))
# pca_var <-  PCA_data %>%  summarise(across(starts_with("PC"), ~var(.x)))
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
palettes_d_names %>% filter(type=="qualitative", length>=10, package!="palettetown")

superpop_pal <- unique(PCA_data$Superpopulation) %>% 
  setNames(c(as.character(paletteer_d("rcartocolor::Safe", n = length(.)-1, direction = -1)), "#794924FF"), .)
#ggsci::category10_d3
shapes <- c(rep(19, length(superpop_pal)-1), 17)
sizes <- c(rep(3, length(superpop_pal)-1), 4)
unknown <- PCA_data %>% filter(grepl("UWW1", `Sample name`))

### PCA plots super population ####
# Create the plot for C1 and C2
plot_comps <- 1:2
plot_axes <- paste0("PC", plot_comps)
pale_theme <- ggthemr("pale", text_size = 20, set_theme = FALSE)
ggplot(PCA_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), 
                               colour=Superpopulation, 
                               group=Superpopulation,
                               size=Superpopulation,
                               shape=Superpopulation)) + 
  geom_point(alpha = 0.65) + # size=3
  scale_size_manual(values = sizes) +
  scale_colour_manual(values = superpop_pal) + 
  scale_x_continuous(breaks = breaks_pretty(n=6)) +
  scale_y_continuous(breaks = breaks_pretty(n=6)) +
  scale_shape_manual(values = shapes) +

  labs( x=glue::glue("PC{plot_comps[1]}"),
        y=glue::glue("PC{plot_comps[2]}")) +
  pale_theme$theme
# Save plot to a pdf file
ggsave(filedate(glue::glue("UWW1_PCA_superpop_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here::here("output")), width = 10, height=8)

# Create the plot for C3 and C4
plot_comps <- 3:4
plot_axes <- paste0("PC", plot_comps)

ggplot(PCA_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), 
                               colour=Superpopulation, 
                               group=Superpopulation,
                               size=Superpopulation,
                               shape=Superpopulation)) + 
  geom_point(alpha = 0.65) + # size=3
  scale_size_manual(values = sizes) +
  scale_colour_manual(values = superpop_pal) + 
  scale_x_continuous(breaks = breaks_pretty(n=6)) +
  scale_y_continuous(breaks = breaks_pretty(n=6)) +
  scale_shape_manual(values = shapes) +
  
  labs( x=glue::glue("PC{plot_comps[1]}"),
        y=glue::glue("PC{plot_comps[2]}")) +
  pale_theme$theme
# Save plot to a pdf file
ggsave(filedate(glue::glue("UWW1_PCA_superpop_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here::here("output")), width = 10, height=8)



### PCA plots EUR/AUS population ####
EUR_PCA_data <- PCA_data %>% mutate(Population=sub(",.+$", "", `Population name`)) %>% 
  filter(`Superpopulation code`=="EUR" | `Population name`=="Australian" | `Population name`=="UWW1")
# Define colour palette
palettes_d_names %>% filter(type=="qualitative", length>=6, package!="palettetown")
pop_pal <- unique(EUR_PCA_data$Population) %>% 
  setNames(as.character(paletteer_d("awtools::mpalette", n = length(.))), .)
shapes <- c(rep(19, length(pop_pal)-1), 17)
sizes <- c(rep(3, length(pop_pal)-1), 4)
# Create the plot for C1 and C2
plot_comps <- 1:2
plot_axes <- paste0("PC", plot_comps)
pale_theme <- ggthemr("pale", text_size = 20, set_theme = FALSE)
ggplot(EUR_PCA_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), 
                               colour=Population, 
                               group=Population,
                               size=Population,
                               shape=Population)) + 
  geom_point(alpha = 0.65) + # size=3
  scale_size_manual(values = sizes) +
  scale_colour_manual(values = pop_pal) + 
  scale_x_continuous(breaks = breaks_pretty(n=6)) +
  scale_y_continuous(breaks = breaks_pretty(n=6)) +
  scale_shape_manual(values = shapes) +
  
  labs( x=glue::glue("PC{plot_comps[1]}"),
        y=glue::glue("PC{plot_comps[2]}")) +
  pale_theme$theme
# Save plot to a pdf file
ggsave(filedate(glue::glue("UWW1_PCA_EUR_pop_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here::here("output")), width = 8, height=6)
# Create the plot for PC3 and PC4
plot_comps <- 3:4
plot_axes <- paste0("PC", plot_comps)
pale_theme <- ggthemr("pale", text_size = 20, set_theme = FALSE)
ggplot(EUR_PCA_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), 
                                   colour=Population, 
                                   group=Population,
                                   size=Population,
                                   shape=Population)) + 
  geom_point(alpha = 0.65) + # size=3
  scale_size_manual(values = sizes) +
  scale_colour_manual(values = pop_pal) + 
  scale_x_continuous(breaks = breaks_pretty(n=6)) +
  scale_y_continuous(breaks = breaks_pretty(n=6)) +
  scale_shape_manual(values = shapes) +
  
  labs( x=glue::glue("PC{plot_comps[1]}"),
        y=glue::glue("PC{plot_comps[2]}")) +
  pale_theme$theme
# Save plot to a pdf file
ggsave(filedate(glue::glue("UWW1_PCA_EUR_pop_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here::here("output")), width = 8, height=6)
