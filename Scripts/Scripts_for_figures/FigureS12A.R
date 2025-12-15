######################################################
####        Figure S12A: Tm curves                ####
######################################################

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

## Load the Prometheus data
data_tm_exp1 <- read_excel('Data/Protein_complex_stability/ProteinStability_Purif1.xlsx',
                      sheet = 'Overview', skip = 0)
data_tm_exp2 <- read_excel('Data/Protein_complex_stability/ProteinStability_Purif2.xlsx',
                           sheet = 'Overview', skip = 0)
data_tm_exp3 <- read_excel('Data/Protein_complex_stability/ProteinStability_Purif3.xlsx',
                           sheet = 'Overview', skip = 0)
data_tm_exp4 <- read_excel('Data/Protein_complex_stability/ProteinStability_Purif4.xlsx',
                           sheet = 'Overview', skip = 0)
data_tm_exp5 <- read_excel('Data/Protein_complex_stability/ProteinStability_Purif5.xlsx',
                           sheet = 'Overview', skip = 0)

#### Organize the above data to be able to link each capillary with a variant and a replicate ####

## Work with experiment 1 ##

# Standardize columns
data_tm_exp1 %<>%
  mutate(Replicate = rep(c(1, 2, 3), nrow(data_tm_exp1)/3)) %>%
  mutate(Sample_standard = `Sample ID`) %>%
  select(Capillary, Sample_standard, Replicate, `Onset #1 for Ratio`, 
         "Inflection Point #1 for Ratio", "Inflection Point #2 for Ratio")

colnames(data_tm_exp1) <- c('Capillary', 'Sample_standard', 'Replicate',
                            'Onset_1', 'Inflection_point_1', 'Inflection_point_2')


## Separate the replicate from the sample type
data_tm_exp1 %<>% separate(Sample_standard, into = c('Sample', 'Replicate'), sep = '_')

data_tm_exp1 %<>% filter(!(is.na(Sample))) %>%
  ## Standardize sample names 
  mutate(Sample_standard = str_replace(string = Sample, pattern = '[\\-]', replacement = '/')
  ) %>%
  mutate(
    Sample_standard = factor(Sample_standard, levels = c('WT', 'Q67C', 'S59Y', 'Q67C+S59Y',
                                                         'WT/WT',
                                                         'WT/Q67C', 'Q67C/WT', 'Q67C/Q67C',
                                                         'WT/S59Y', 'S59Y/WT', 'S59Y/S59Y',
                                                         'Q67C/S59Y', 'S59Y/Q67C'
    )
    )
  )

## Show only the data for the single copy variants
data_single_copy <- data_tm_exp1 %>% 
  filter(!(str_detect(string = Sample_standard, pattern = '[\\/]')))

## Work with experiment 2 ##

data_tm_exp2 %<>%
  separate(col = `Sample ID`, into = c('Sample_standard', 'Replicate'), sep = '_') %>%
  mutate(Sample_standard = str_replace(string = Sample_standard, pattern = '-', replacement = '+')) %>%
  select(Capillary, Sample_standard, Replicate, `Onset #1 for Ratio`, 
         "Inflection Point #1 for Ratio", "Inflection Point #2 for Ratio")

colnames(data_tm_exp2) <- c('Capillary', 'Sample_standard', 'Replicate',
                            'Onset_1', 'Inflection_point_1', 'Inflection_point_2')


## Work with experiments 3-5 ##

# Standardize columns
data_tm_exp3 %<>%
  mutate(Replicate = rep(c(1, 2, 3), nrow(data_tm_exp3)/3)) %>%
  mutate(Sample_standard = `Sample ID`) %>%
  select(Capillary, Sample_standard, Replicate, `Onset #1 for Ratio`, 
         "Inflection Point #1 for Ratio", "Inflection Point #2 for Ratio")

colnames(data_tm_exp3) <- c('Capillary', 'Sample_standard', 'Replicate',
                            'Onset_1', 'Inflection_point_1', 'Inflection_point_2')

data_tm_exp4 %<>%
  mutate(`Sample ID` = ifelse(or(`Sample ID` == 'WT1', `Sample ID` == 'WT2'),
                              'WT', `Sample ID`)) %>%
  mutate(Replicate = c(rep(c(1, 2, 3), (nrow(data_tm_exp4)-3)/3), 
                       c(4,5,6))
  ) %>%
  mutate(Sample_standard = `Sample ID`) %>%
  select(Capillary, Sample_standard, Replicate, `Onset #1 for Ratio`, 
         "Inflection Point #1 for Ratio")

colnames(data_tm_exp4) <- c('Capillary', 'Sample_standard', 'Replicate',
                            'Onset_1', 'Inflection_point_1')

data_tm_exp5 %<>%
  mutate(`Sample ID` = ifelse(or(`Sample ID` == 'WT_1', `Sample ID` == 'WT_2'),
                              'WT', `Sample ID`)) %>%
  mutate(Replicate = c(rep(c(1, 2, 3), (nrow(data_tm_exp5)-3)/3), 
                       c(4,5,6))
  ) %>%
  mutate(Sample_standard = `Sample ID`) %>%
  select(Capillary, Sample_standard, Replicate, `Onset #1 for Ratio`, 
         "Inflection Point #1 for Ratio", "Inflection Point #2 for Ratio")

colnames(data_tm_exp5) <- c('Capillary', 'Sample_standard', 'Replicate',
                            'Onset_1', 'Inflection_point_1', 'Inflection_point_2')

## Put all the data in the same table
data_tm_exp_all <- bind_rows(data_single_copy %>% 
                               select(Capillary, Sample_standard, Replicate, Onset_1, 
                                      Inflection_point_1,) %>%
                               mutate(Purification = '1', 
                                      Replicate = as.numeric(Replicate)), 
                             data_tm_exp2 %>%
                               select(Capillary, Sample_standard, Replicate, Onset_1, 
                                      Inflection_point_1) %>%
                               mutate(Purification = '2', 
                                      Replicate = as.numeric(Replicate)), 
                             data_tm_exp3 %>%
                               select(Capillary, Sample_standard, Replicate, Onset_1, 
                                      Inflection_point_1) %>%
                               mutate(Purification = '3'),
                             data_tm_exp4 %>%
                               select(Capillary, Sample_standard, Replicate, Onset_1, 
                                      Inflection_point_1) %>%
                               mutate(Purification = '4'),
                             data_tm_exp5 %>%
                               select(Capillary, Sample_standard, Replicate, Onset_1, 
                                      Inflection_point_1) %>%
                               mutate(Purification = '5'),
                             
)

#### Load the Tm curve data ####

## A function to load the signal ##
load_signal <- function(infile){
  data_signal <- read_excel(infile, sheet = 'Ratio', skip = 0)
  
  colnames(data_signal)[1:2] <- c('Time', 'Temperature')
  
  data_signal <- data_signal[3:(nrow(data_signal)), 1:(ncol(data_signal))]
  
  data_signal_long <- data_signal %>% ungroup() %>% 
    select(-Time) %>%
    group_by(Temperature) %>%
    pivot_longer(cols = colnames(data_signal)[3:ncol(data_signal)],
                 names_to = 'Capillary', values_to = 'Fluorescence') %>%
    mutate(Capillary = as.numeric(Capillary), 
           Fluorescence = as.numeric(Fluorescence), 
           Temperature = as.numeric(Temperature)) %>%
    arrange(Temperature, Capillary)
  
  return(data_signal_long)
}

data_tm_curve1_long <- load_signal('Data/Protein_complex_stability/ProteinStability_Purif1.xlsx')
data_tm_curve2_long <- load_signal('Data/Protein_complex_stability/ProteinStability_Purif2.xlsx')
data_tm_curve3_long <- load_signal('Data/Protein_complex_stability/ProteinStability_Purif3.xlsx')
data_tm_curve4_long <- load_signal('Data/Protein_complex_stability/ProteinStability_Purif4.xlsx')
data_tm_curve5_long <- load_signal('Data/Protein_complex_stability/ProteinStability_Purif5.xlsx')

## Merge with the annotations from the previous analysis ##
data_tm_curve_all <- bind_rows(
  data_tm_curve1_long %>% mutate(Purification = '1'), 
  data_tm_curve2_long %>% mutate(Purification = '2'), 
  data_tm_curve3_long %>% mutate(Purification = '3'), 
  data_tm_curve4_long %>% mutate(Purification = '4'), 
  data_tm_curve5_long %>% mutate(Purification = '5')
)

data_tm_curve_all_new <- inner_join(x = data_tm_curve_all %>% ungroup() %>%
                                     mutate(Capillary = as.numeric(Capillary), 
                                            Purification = as.numeric(Purification)), 
                                   y = data_tm_exp_all %>% ungroup() %>%
                                     select(-Onset_1, -Inflection_point_1) %>%
                                     mutate(Capillary = as.numeric(Capillary), 
                                            Purification = as.numeric(Purification)), 
                                   by = c('Capillary' = 'Capillary', 'Purification' = 'Purification')
)

samples_show <- c('WT', 'C47Y', 'W38L', 'Q39P', 'S59Y', 'H62N', 'Q67C', 'Y69D', 'Q67C+S59Y')

## Example of the Tm curves (need to put all of them together, indicate replicates, facet by purification)
p_signal <- data_tm_curve_all_new %>% rowwise() %>%
  filter(Sample_standard %in% samples_show) %>%
  mutate(Purification_new = str_c('Purification ', toString(Purification), sep = ''), 
         Sample_standard = factor(Sample_standard, levels = samples_show)) %>%
  ggplot(aes(x = Temperature, y = Fluorescence)) +
  geom_line(aes(group = interaction(Sample_standard, Capillary),
                colour = as.factor(Sample_standard))) +
  labs(colour = 'Sample') +
  facet_wrap(~Purification_new, nrow = 1) +
  theme(legend.position = 'top', 
        legend.justification = 'center') +
  ylab('Fluorescence\n') +
  scale_y_continuous(breaks = c(0.60, 0.80, 1.00, 1.20), labels = c('0.60', '0.80', '1.00', '1.20'))
p_signal

#### Load the first derivative data ####

# A function to load the first derivative
load_d1 <- function(infile){
  data_signal <- read_excel(infile, sheet = 'Ratio (1st deriv.)', skip = 0)
  
  colnames(data_signal)[1:2] <- c('Time', 'Temperature')
  
  data_signal <- data_signal[3:(nrow(data_signal)), 1:(ncol(data_signal))]
  
  data_signal_long <- data_signal %>% ungroup() %>% 
    select(-Time) %>%
    group_by(Temperature) %>%
    pivot_longer(cols = colnames(data_signal)[3:ncol(data_signal)],
                 names_to = 'Capillary', values_to = 'fluo_d1') %>%
    mutate(Capillary = as.numeric(Capillary), 
           fluo_d1 = as.numeric(fluo_d1), 
           Temperature = as.numeric(Temperature)) %>%
    arrange(Temperature, Capillary)
  
  return(data_signal_long)
}

data_tm_d1_purif1_long <- load_d1('Data/Protein_complex_stability/ProteinStability_Purif1.xlsx')
data_tm_d1_purif2_long <- load_d1('Data/Protein_complex_stability/ProteinStability_Purif2.xlsx')
data_tm_d1_purif3_long <- load_d1('Data/Protein_complex_stability/ProteinStability_Purif3.xlsx')
data_tm_d1_purif4_long <- load_d1('Data/Protein_complex_stability/ProteinStability_Purif4.xlsx')
data_tm_d1_purif5_long <- load_d1('Data/Protein_complex_stability/ProteinStability_Purif5.xlsx')

## Put all the data together and merge with the sample names
data_tm_d1_all <- bind_rows(
  data_tm_d1_purif1_long %>% mutate(Purification = '1'), 
  data_tm_d1_purif2_long %>% mutate(Purification = '2'), 
  data_tm_d1_purif3_long %>% mutate(Purification = '3'), 
  data_tm_d1_purif4_long %>% mutate(Purification = '4'), 
  data_tm_d1_purif5_long %>% mutate(Purification = '5')
)

data_tm_d1_all_new <- inner_join(x = data_tm_d1_all %>% ungroup() %>%
                                      mutate(Capillary = as.numeric(Capillary), 
                                             Purification = as.numeric(Purification)), 
                                    y = data_tm_exp_all %>% ungroup() %>%
                                      select(-Onset_1, -Inflection_point_1) %>%
                                      mutate(Capillary = as.numeric(Capillary), 
                                             Purification = as.numeric(Purification)), 
                                    by = c('Capillary' = 'Capillary', 'Purification' = 'Purification')
)

## Example of the Tm curves (need to put all of them together, indicate replicates, facet by purification)
p_d1 <- data_tm_d1_all_new %>% rowwise() %>%
  filter(Sample_standard %in% samples_show) %>%
  mutate(Purification_new = str_c('Purification ', toString(Purification), sep = ''), 
         Sample_standard = factor(Sample_standard, levels = samples_show)) %>%
  ggplot(aes(x = Temperature, y = fluo_d1)) +
  geom_line(aes(group = interaction(Sample_standard, Capillary),
                colour = as.factor(Sample_standard))) +
  labs(colour = 'Sample') +
  facet_wrap(~Purification_new, nrow = 1) +
  theme(legend.position = 'top', 
        legend.justification = 'center') +
  ylab('First derivative')
p_d1

## Extract legend from p_signal
p_legend <- get_legend(p_signal)

p_figS12A <- plot_grid(
  p_legend,
  p_signal + theme(legend.position = 'none'),
  p_d1 + theme(legend.position = 'none'), 
  nrow = 3, labels = c('', 'A', 'B'), 
  label_size = panel_label_size, label_fontface = 'bold', rel_heights = c(0.2, 1, 1)
  )

ggsave(p_figS12A, device = cairo_pdf, width = 25, height = 14, units = 'cm', dpi = 300, 
       filename = 'Figures/Supp_figures/FigS12A.pdf')
