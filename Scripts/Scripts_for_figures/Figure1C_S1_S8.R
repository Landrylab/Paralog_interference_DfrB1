#############################################################
####               Figure 1E,S1B,S8                      ####
#############################################################

# Load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(viridis)
library(agricolae)
library(ggrepel)
library(bio3d)
library(readxl)

theme_set(theme_cowplot() +
            theme(panel.background = element_rect(fill = 'white'),
                  plot.background = element_rect(fill = 'white'),
                  axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12, face = 'bold'),
                  strip.text = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  axis.line = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = 'black', linewidth = 1) 
            )
)

panel_label_size = 14

# A function to draw heatmaps
draw_CHeatmap <- function(in_heatmap){
  return(grid.grabExpr(draw(in_heatmap)))
}

## Set the path to the main Github directory
# setwd('/path/to/Paralog_interference_DfrB1')
setwd('/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Angel/PhD_projects/Yacine/Manuscript/Github/Paralog_interference_DfrB1/')

## Define a function to calculate t_Int / t_NoInt
tint_tnoint <- function(ratio_mu, np_fix){
  return((1 + ratio_mu) / (np_fix + ratio_mu))
}


# Define values of ratio_mu and np_fix(F)
ratios_mu <- c(2, 3, 5, 10, 20)
values_npfix <- seq(from = 0, to = 1, by = 0.01)

## Define a dataframe and then populate it
df_final <- data.frame()

for(ratio_mu in ratios_mu){
  
  new_section <- data.frame(np_fix = values_npfix, 
                            ratio_mu = ratio_mu,
                            ratio_t = tint_tnoint(ratio_mu, values_npfix))
  df_final <- bind_rows(df_final, new_section)
  
}

## Draw the figure
p_fig1E <- df_final %>%
  mutate(ratio_mu = as.factor(ratio_mu)) %>%
  ggplot(aes(x = np_fix, y = ratio_t, colour = ratio_mu)) +
  geom_line(aes(group= ratio_mu)) +
  labs(
    x = 'Probability of fixation of interfering\nLOF relative to a neutral allele',
       y = expression(
         paste(
           bold('Redundancy time ratio ('),
           bold(T[Int]), bold('/'), bold(T[NoInt]), bold(')'),
           sep = '')
         ),
       colour = expression(paste(bold('(\u03bc'[A]), bold('+\u03bc'[E]), bold(')/'), 
                                 bold('\u03bc'[F]),
                                 sep = ''))
       ) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = 'inside',
        legend.position.inside = c(0.55, 0.75),
        legend.text = element_text(size = 9)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))
p_fig1E
ggsave(p_fig1E, device = cairo_pdf, width = 12, height = 8, dpi = 300, units = 'cm', 
       filename = 'Figures/Main_figures/Fig1E.pdf')
## The rest of the figure was assembled with Inkscape

#### Fig. S1B: Case when the loss of one copy is deleterious ####

# Set the mu_ratio to 2
mu_f <- 1
mu_e <- 1
mu_a <- 1
mu_ratio <- (mu_a + mu_e) / mu_f

# Use an array of values for NP_Fix(Loss)
list_p_fix_ae <- c(0.2, 0.4, 0.6, 0.8)
list_p_fix_f <- seq(from = 0.01, to = 1, by = 0.01)

## Define a function to calculate t_Int / t_NoInt when the loss of one copy is not neutral
tint_tnoint_not_neutral <- function(mu_f, mu_e, mu_a, p_fix_f, p_fix_ae){
  num_exp <- p_fix_ae*(mu_f + mu_e + mu_e)
  denom_exp <- p_fix_f*mu_f + p_fix_ae*(mu_a + mu_e)
  return(num_exp / denom_exp)
}

## Solve the equation for each of these parameters along a range of values for NP_Fix(F)
df_final_not_neutral <- data.frame()
for(p_fix_f in list_p_fix_f){
  
  new_section <- data.frame(p_fix_f = p_fix_f, 
                            p_fix_ae = list_p_fix_ae,
                            mu_f = mu_f, mu_e = mu_e, mu_a = mu_a, 
                            ratio_t = tint_tnoint_not_neutral(mu_f, mu_e, mu_a, p_fix_f, list_p_fix_ae))
  
  df_final_not_neutral <- bind_rows(df_final_not_neutral, new_section)
  
}

df_final_not_neutral %<>% arrange(desc(p_fix_f), p_fix_ae)

## Draw the curves
p_figS1B <- df_final_not_neutral %>%
  mutate(bool_t = ratio_t < 1) %>%
  mutate(p_fix_ae = factor(p_fix_ae, 
                           levels = c(0.2, 0.4, 0.6, 0.8)
                           )
         ) %>%
  ggplot(aes(x = p_fix_f, y = ratio_t, colour = p_fix_ae)) +
  geom_line(aes(group= p_fix_ae)) +
  labs(
    x = 'Probability of fixation of interfering\nLOF relative to a neutral allele',
    y = expression(
      paste(
        bold('Redundancy time ratio ('),
        bold(T[Int]), bold('/'), bold(T[NoInt]), bold(')'),
        sep = '')
    ),
       colour = expression(bold(NP['Fix(Loss)']))
  ) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = 'inside',
        legend.position.inside = c(0.55, 0.75),
        legend.text = element_text(size = 9)) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))
p_figS1B
ggsave(p_figS1B, device = cairo_pdf, width = 12, height = 8, dpi = 300, units = 'cm', 
       filename = 'Figures/Supp_figures/FigS1B.pdf')
## The rest of the figure was assembled with Inkscape

#### Figure S8 ####

#### Load the tables of interfering candidates and all DMS selection coefficients ####
dms_data <- read_delim('Data/DMS_bulk_competition_experiments/dms_selection_coefficients_Singleton_Duplicated.tsv', delim = '\t')

## Identify the WT codons at each position
dms_data_check_TAG <- dms_data %>% 
  separate(col = Genetic_code, into = c('WT_Codon', 'Position', 'Mutant_Codon'), sep = c(3, 5)) %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1), remove = F)

# Retrieve the genetic code
genetic_code <- dms_data_check_TAG %>% select(Mutant_Residue, Mutant_Codon) %>%
  arrange(Mutant_Residue, Mutant_Codon) %>% unique()

## Load the interfering mutations
interfering_muts <- read_delim('Supp_tables/TableS4_DominantNegativeCandidates.tsv', delim = '\t')

## Check the total number of codon mutations that would result in interference
interfering_muts_total_tmp <- left_join(x = interfering_muts, 
                                    y = genetic_code, 
                                    by = c('Mutant_Residue' = 'Mutant_Residue'), 
                                    relationship = 'many-to-many')

wt_codons <- dms_data_check_TAG %>% select(Position, WT_Codon) %>%
  mutate(Position = as.numeric(Position)) %>%
  unique()
# 69 codons

## Each codon has three positions, and each positions has three available mismatches
total_accessible_mismatch <- nrow(wt_codons)*3*3
# 621 mutations

## Get the total number of mutations sampled
dms_data_check_TAG_mut_count <- dms_data_check_TAG %>% 
  filter(Background == 'Duplicated', Mutant_Codon != WT_Codon, Mutant_Residue != '*') %>%
  select(Background, Position, Mutant_Codon) %>%
  unique()
# 1310 mismatch mutations sampled (residue level)
# 4116 mismatch mutations sampled (codon level)

interfering_muts_total <- left_join(x = interfering_muts_total_tmp %>% 
                                        mutate(Position = as.numeric(Position)), 
                                      y = wt_codons, 
                                      by = c('Position' = 'Position')
                                      )

count_mismatches <- function(codon1, codon2){
  mismatch_count <- 0
  for(i in 1:nchar(codon1)){
    if(str_sub(codon1, start = i, end =  i) != str_sub(codon2, start = i, end =  i)){
      mismatch_count <- mismatch_count + 1
    }
  }
  return(mismatch_count)
}

## Filter to count how many are one point mutation away from WT
interfering_muts_total %<>% rowwise() %>%
  mutate(num_mismatches = count_mismatches(WT_Codon, Mutant_Codon))
# 180 total interfering mutations (codon level)

interfering_muts_accessible <- interfering_muts_total %>% filter(num_mismatches == 1)
# 11 interfering mutations (accessible through one point mutation)

## Load biophysical effects to identify loss of affinity
biophys_effects <- read_delim('Supp_tables/TableS2_DMS_s_biophysical_effects.tsv', delim = '\t')

## Identify loss of affinity as mutations that destabilize the dimer interface
ddG_threshold <- 2
LOA_dimer <- biophys_effects %>% filter(Foldx_ddG_dimerization >= ddG_threshold)

# Add WT and mutant codons
LOA_dimer_tmp <- left_join(x = LOA_dimer %>% 
                             select(WT_Residue, Position, Mutant_Residue, Selection_coefficient_Duplicated,
                                    Foldx_ddG_dimerization
                                    ), 
                           y = genetic_code, 
                           by = c('Mutant_Residue' = 'Mutant_Residue'), 
                           relationship = 'many-to-many'
                           )


LOA_dimer_total <- left_join(x = LOA_dimer_tmp %>% 
                               mutate(Position = as.numeric(Position)), 
                             y = wt_codons, 
                             by = c('Position' = 'Position')
)
# 361 LOA codons

LOA_dimer_total %<>% rowwise() %>%
  mutate(num_mismatches = count_mismatches(WT_Codon, Mutant_Codon))

LOA_dimer_accessible <- LOA_dimer_total %>% filter(num_mismatches == 1)
# 34 accessible LOA mutations

## Get total and accessible stop codons (LOA) ##
stop_codons_count_sampled <- dms_data_check_TAG %>% 
  filter(Background == 'Duplicated', Mutant_Residue == '*') %>%
  select(Background, Position, WT_Codon, Mutant_Codon) %>%
  unique()
# 138 stop codons sampled

## Check count of accessible stop codons
stop_codons_count_sampled %<>% rowwise() %>%
  mutate(num_mismatches = count_mismatches(WT_Codon, Mutant_Codon))

accessible_stop_codons <- stop_codons_count_sampled %>% filter(num_mismatches == 1)
# 13 accessible stop codons

#### Use the estimates from the DMS data to calculate the increase           ####
#### in the residence time of redundant copies due to interference           ####

## Assume that LOA and LOE are neutral, sampled mutations ##

s_LOF <- -0.64 ## Based on the median effect of interfering mutations in the DMS data

total_mutations <- 4116

n_LOF <- 180
mu_LOF <- n_LOF / total_mutations
## Could also be this plus the proportion of mutations with strong destabilizing effects on tet. interface

n_LOA <- 361 + 138
mu_LOA <- n_LOA / total_mutations

## Calculate probability of fixation from selection coefficients
get_pfix <- function(sel_coeff, c, Neff){
  ## A function to get the probability of fixation of a mutation from its selection coefficient
  ## For haploid individuals:
  # c = 1 (one genome copy per individual)
  # p = 1 (number of individuals carrying the mutant)
  
  if(sel_coeff == 0){
    return(1/Neff)
  }else{
    p = 1/(c*Neff) # Assuming the mutation only appears in one individual
    numerator <- 1 - exp(-2*c*Neff*sel_coeff*p)
    denominator <- 1 - exp(-2*c*Neff*sel_coeff)
    
    return(numerator/denominator)
  }
}

Neff <- 1e6
get_pfix(-0.0001, 1, Neff)
Pfix_LOF <- get_pfix(s_LOF, 1, Neff) # For interfering mutations

residence_time1 <- function(Neff, Pfix_LOF, Pfix_LOA_LOE, mu_LOF, mu_LOA, mu_LOE){
  ## A function that compares the residence times with and without interference
  
  T_NoInt <- 1 / (mu_LOF*Neff*Pfix_LOA_LOE + mu_LOA*Neff*Pfix_LOA_LOE + mu_LOE*Neff*Pfix_LOA_LOE)
  T_Int <- 1 / (mu_LOF*Neff*Pfix_LOF + mu_LOA*Neff*Pfix_LOA_LOE + mu_LOE*Neff*Pfix_LOA_LOE)
  
  return(T_Int / T_NoInt)
}

Pfix_LOA_LOE <- get_pfix(0, 1, Neff) # Assuming the loss of one copy is neutral

## An example of the calculation
residence_time1(Neff, Pfix_LOF, Pfix_LOA_LOE, mu_LOF, mu_LOA, 0.1)
mu_LOE_values <- seq(from = 0.001, to = 0.7, by = 0.001)

Neff1 <- 1e6
data_model1 <- residence_time1(Neff1, Pfix_LOF, Pfix_LOA_LOE, mu_LOF, mu_LOA, mu_LOE_values)

## Load data from the Urtecho et al. 2023
data <-  read.table('Data/Ecoli_promoter_data/endo_scramble_expression_formatted_std.txt', 
                    header = T)
ttest <- read.table('Data/Ecoli_promoter_data/endo_scramble_ttests.txt', header = T)
data <- left_join(data, select(ttest, -tss_name), by = 'name')

## For each promoter type, count the total tested instances and the ones that decreased expression by 50%
data_all_entries <- data %>% ungroup() %>%
  group_by(tss_name) %>%
  summarise(total_tested = n())

data_loe <- data %>% ungroup() %>%
  filter(relative_exp <= 0.5, significant) %>%
  group_by(tss_name) %>%
  summarise(total_loe = n())

data_frac_loe <- left_join(x = data_all_entries, y = data_loe, 
                           by = c('tss_name' = 'tss_name')) %>%
  complete(fill = list(total_loe = 0)) %>%
  mutate(frac_loe = round(total_loe / total_tested, 4)) %>%
  filter(total_tested >= 15)

median_mu_LOE <-median(data_frac_loe$frac_loe)

data_fig_model1 <- data.frame(cbind(mu_LOE_values, data_model1))
colnames(data_fig_model1) <- c('mu_LOE', 'Ratio_residence_times')

p_analytical_sampled <- data_fig_model1 %>% ggplot(aes(x = mu_LOE, y = Ratio_residence_times)) +
  geom_line() +
  geom_vline(xintercept = median_mu_LOE, linetype = 'dashed') +
  labs(
    x = expression(bold(µ[E])),
    y = expression(
      paste(
        bold('Redundancy time ratio ('),
        bold(T[Int]), bold('/'), bold(T[NoInt]), bold(')'),
        sep = '')
    )
  ) +
  geom_rug(data = data_frac_loe, inherit.aes = F, aes(x = frac_loe))
p_analytical_sampled

mu_ratio <- (median_mu_LOE + mu_LOA) / mu_LOF
# 3.56

## Repeat for the accessible mutations ##

## Assume that LOA and LOE are neutral, accessible mutations ##
s_LOF <- -0.64 ## Based on the median effect of interfering mutations in the DMS data

total_mutations <- 621

n_LOF <- 11
mu_LOF <- n_LOF / total_mutations

n_LOA <- 34 + 13
mu_LOA <- n_LOA / total_mutations

Neff <- 1e6
get_pfix(-0.0001, 1, Neff)
Pfix_LOF <- get_pfix(s_LOF, 1, Neff) # For interfering mutations

Pfix_LOA_LOE <- get_pfix(0, 1, Neff) # Assuming the loss of one copy is neutral

## An example of the calculation
residence_time1(Neff, Pfix_LOF, Pfix_LOA_LOE, mu_LOF, mu_LOA, 0.1)
mu_LOE_values <- seq(from = 0.001, to = 0.7, by = 0.001)

Neff1 <- 1e6
data_model1 <- residence_time1(Neff1, Pfix_LOF, Pfix_LOA_LOE, mu_LOF, mu_LOA, mu_LOE_values)

data_fig_model1 <- data.frame(cbind(mu_LOE_values, data_model1))
colnames(data_fig_model1) <- c('mu_LOE', 'Ratio_residence_times')

p_analytical_accessible <- data_fig_model1 %>% ggplot(aes(x = mu_LOE, y = Ratio_residence_times)) +
  geom_line() +
  geom_vline(xintercept = median_mu_LOE, linetype = 'dashed') +
  labs(
    x = expression(bold(µ[E])),
    y = expression(
      paste(
        bold('Redundancy time ratio ('),
        bold(T[Int]), bold('/'), bold(T[NoInt]), bold(')'),
        sep = '')
    )
  ) +
  geom_rug(data = data_frac_loe, inherit.aes = F, aes(x = frac_loe))
p_analytical_accessible

mu_ratio <- (median_mu_LOE + mu_LOA) / mu_LOF
# 6.22 

## Add the panels for the Moran model ##

# Sampled
data_moran_sampled <- read_delim('Data/Moran_model/Sampled/MoranModel_0.2pct_number_events_LOA_LOE_Neutral_Sampled.tsv', delim ='\t')

data_moran_wide_sampled <- data_moran_sampled %>% 
  mutate(population = row_number()) %>%
  pivot_longer(cols = c('Num_events_int', 'Num_events_noint'), 
               names_to = 'Sim_type', values_to = 'Sim_value') %>%
  mutate(Sim_type = case_when(
    Sim_type == 'Num_events_int' ~ 'Interference', 
    Sim_type == 'Num_events_noint' ~ 'No interference'
  )
  )

p_moran_sampled <- data_moran_wide_sampled %>%
  ggplot(aes(x = Sim_value, fill = Sim_type)) +
  geom_histogram(alpha = 0.4, position = 'identity') +
  theme(legend.position = 'top', 
        legend.justification  = 'center') +
  xlab('Number of birth/death events\nbefore reaching 0.2N') + ylab('Number of populations') +
  labs(fill = '')
p_moran_sampled

## Calculate the mean increase in residence time
mean_noint <- mean(data_moran_sampled$Num_events_noint)
mean_int <- mean(data_moran_sampled$Num_events_int)
mean_increase <- mean_int / mean_noint

# Accessible
data_moran_accessible <- read_delim('Data/Moran_model/Accessible/MoranModel_0.2pct_number_events_LOA_LOE_Neutral_Accessible.tsv', delim ='\t')

data_moran_wide_accessible <- data_moran_accessible %>% 
  mutate(population = row_number()) %>%
  pivot_longer(cols = c('Num_events_int', 'Num_events_noint'), 
               names_to = 'Sim_type', values_to = 'Sim_value') %>%
  mutate(Sim_type = case_when(
    Sim_type == 'Num_events_int' ~ 'Interference', 
    Sim_type == 'Num_events_noint' ~ 'No interference'
  )
  )

p_moran_accessible <- data_moran_wide_accessible %>%
  ggplot(aes(x = Sim_value, fill = Sim_type)) +
  geom_histogram(alpha = 0.4, position = 'identity') +
  theme(legend.position = 'top', 
        legend.justification  = 'center') +
  xlab('Number of birth/death events\nbefore reaching 0.2N') + ylab('Number of populations') +
  labs(fill = '')
p_moran_accessible

## Calculate the mean increase in residence time
mean_noint <- mean(data_moran_accessible$Num_events_noint)
mean_int <- mean(data_moran_accessible$Num_events_int)
mean_increase <- mean_int / mean_noint

## Put the four panels together
p_figS8 <- plot_grid(p_analytical_sampled, p_moran_sampled, NULL, 
                             p_analytical_accessible, p_moran_accessible, NULL,
                             nrow = 2, labels = c('A','B', '','C','D', ''), label_size = panel_label_size,
                             label_fontface = 'bold', rel_widths = c(1, 1, 0.05))

ggsave(p_figS8, device = cairo_pdf, width = 20, height = 20, dpi = 300, units = 'cm', 
       filename = 'Figures/Supp_figures/FigS8.pdf')
