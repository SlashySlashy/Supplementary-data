library(readxl)
library(lvplot)
library(showtext)
library(DescTools)
library(tidyverse)

#source("get-metadata.R")
source("dmg.R")
source("get_calculate_plot_grid.R")
#source("get_penalized_weighted_median_reads.R")
source("perk.R")
#source("damage_est_analysis.R")
source("perk_wrapper.R")
source("perk_wrapper_function.R")
source("get_dmg_decay_fit.R")


# Load metaDMG data


holi_data <- read_tsv("metadmg_august_irene.tsv") %>% filter(grepl("August", label))
colnames(holi_data)
head(holi_data)

# holi_data2 <- holi_data1 |> 
#   separate(sample_id, into = c("label", "prefix"), sep = "-", extra = "merge", fill = "right") |> 
#   mutate(label = ifelse(is.na(label), prefix, label)) # Handle cases without a dash
#holi_data <- holi_data |>
 # rename(label = library_id)

unique(holi_data$rank)

# set number of reads minimum
z <- 1000
# set taxonomic level 
y <- 'genus'

# Get genus and with at least 100 reads
holi_data_gen_euk_1000 <- holi_data |>
  filter(rank == y) |>
  filter(grepl("Eukaryota", taxa_path)) |>
  filter(nreads >= z) |>
  rename(tax_name = taxid, n_reads = nreads)


# Let's get the damage fits using CCC https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Lins_Concordance_Correlation_Coefficient.pdf 
samples <- holi_data$label |> unique()
dat1000 <- dmg_fwd_CCC(holi_data_gen_euk_1000, samples, ci = "asymptotic", nperm = 100, nproc = 10)

# Get good fits for the 100 reads
dat1000_filt <- dat1000 |> 
  ungroup() |>
  inner_join(holi_data_gen_euk_1000) |>
  mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
  mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))


# 100 reads, good fit
tax <- dat1000_filt |>
  ungroup() |>
  filter(fit == "good") |>
  group_by(label) |>
  slice_sample(n = 1000) |>
  ungroup()

tax

samples <- tax$label |> unique()
samples


plots100_good <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 15, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  data <- dat |> filter(label == x)
  grid_size <- calculate_plot_grid(length(data$tax_name))
  l <- lapply(data$name, function(X) {
    df1 <- data |>
      filter(name == X)
    p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
    p <- p + ggtitle(X)
    return(p)
  })
  plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
  tit <- paste0(x, " -- ", data$label_fig |> unique())
  ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                        color = "black", face = "bold", size = 12
  ))
}, .progress = TRUE)

names(plots100_good) <- samples

pdf("plots1000ggood.pdf", width = 20, height = 20)
print(plots100_good)
dev.off()


# 1000 reads, bad fit
tax <- dat1000_filt |>
  ungroup() |>
  filter(fit == "bad") |>
  group_by(label) |>
  slice_sample(n = 1000) |>
  ungroup()

samples <- tax$label |> unique()

plots100_bad <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 15, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  data <- dat |> filter(label == x)
  grid_size <- calculate_plot_grid(length(data$tax_name))
  l <- lapply(data$name, function(X) {
    df1 <- data |>
      filter(name == X)
    p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
    p <- p + ggtitle(X)
    return(p)
  })
  plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
  tit <- paste0(x, " -- ", data$label_fig |> unique())
  ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                        color = "black", face = "bold", size = 12
  ))
}, .progress = TRUE)

names(plots100_bad) <- samples

pdf("plots1000gbad.pdf", width = 20, height = 20)
print(plots100_bad)
dev.off()
