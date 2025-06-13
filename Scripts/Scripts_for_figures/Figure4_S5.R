#####################################################################
####                       Figures_4_S5                          ####
#####################################################################

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

library(growthcurver)

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
setwd('/path/to/Paralog_interference_DfrB1')

#### Figures 4A-4B ####

## Load the DMS data
dms_data <- read_delim('Data/DMS_bulk_competition_experiments/dms_selection_coefficients_Singleton_Duplicated.tsv', delim = '\t')

## Confirm that TAG is not present
dms_data_check_TAG <- dms_data %>% 
  separate(col = Genetic_code, into = c('WT_Codon', 'Position', 'Mutant_Codon'), sep = c(3, 5))
all(dms_data_check_TAG$Mutant_Codon != 'TAG')
all(dms_data_check_TAG$Mutant_Codon != 'UAG')

## Save table S1 (selection coefficients for codons)
dms_data_codon <- dms_data %>% ungroup() %>%
  group_by(Background, Run, Mutation, Genetic_code) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient)) %>% ungroup() %>%
  group_by(Background, Mutation, Genetic_code) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

dms_data_codon_wide <- dms_data_codon %>% ungroup() %>%
  pivot_wider(names_from = Background, names_prefix = 'Selection_coefficient_',
              values_from = Selection_coefficient)
write.table(dms_data_codon_wide, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Supp_tables/TableS1_DMS_selectionCoefficients_codons.tsv')

## Initialize variables
dom_eff_abundance_all <- c()

data_annotation_dfrb1 <- read_delim('Data/DfrB1_annotation/data_annotation_DfrB1.txt', delim = '\t')

# Pivot to a longer format
data_annotation_new <- data_annotation_dfrb1 %>% rowwise() %>% mutate(
  `Only NADPH` = ifelse(and(NADPH == 1, 
                            sum(`A,C`, `A,D`, DHF, NADPH, Cat_residues, Buried, Disordered_region) == 1), 1, 0), 
  Unannotated = ifelse(sum(`A,C`, `A,D`, DHF, NADPH, Cat_residues, Buried, Disordered_region) == 0, 1, 0)
) %>%
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH', 'Cat_residues', 
                        'Only NADPH', 'Disordered_region', 'Buried', 'Unannotated'), 
               names_to = 'Site', values_to = 'Site_check') %>%
  filter(Site_check == 1) %>% 
  mutate(Site = case_when(
    Site == 'A,C' ~ 'Dimerization interface' , 
    Site == 'A,D' ~ 'Tetramerization interface', 
    Site == 'Cat_residues' ~ 'Catalytic residues', 
    Site =='Disordered_region'~ 'Disordered region', 
    Site == 'Buried' ~ 'Buried residues',
    Site == 'DHF' ~ 'DHF binding', 
    Site == 'NADPH' ~ 'NADPH binding',
    Site == 'Only NADPH' ~ 'Binding only NADPH',
    TRUE ~ Site
  )
  )

## Merge the categories for catalytic residues and DHF binding
data_annotation_new %<>% 
  filter(Site != 'Catalytic residues') %>%
  mutate(Site = ifelse(Site == 'DHF binding', 'Cat. residues/DHF binding', Site))

backgrounds <- unique(dms_data$Background)

## Define cutoff for standard deviations of stop codons used to identify interfering mutations
cutoff_sd <- 2.5

#### Panels of the negative interference effects in the single-copy and duplication backgrounds ####

## Single copy ##

## Average by amino acid but not by run or replicate
dms_data_cds <- dms_data %>% ungroup() %>%
  filter(Background == 'Singleton', !(is.na(Selection_coefficient))) %>%
  group_by(Run, Replicate, Mutation) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

mutation_list <- unique(dms_data_cds$Mutation)

dms_data_cds %<>% filter(!(is.na(Selection_coefficient)))

## Get a distribution of synonymous codons
synonymous_codons <- dms_data %>% ungroup() %>%
  filter(Background == 'Singleton', !(is.na(Selection_coefficient))) %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  separate(col = Genetic_code, into = c('WT_Codon', 'Position2', 'Mutant_Codon'), sep = c(3, -3)) %>%
  filter(and(WT_Residue == Mutant_Residue, WT_Codon != Mutant_Codon)) %>%
  group_by(WT_Residue, Position, Mutant_Residue) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

## Save the data for synonymous codons for EMPTY
dms_data_EMPTY <- dms_data_cds
synonymous_codons_EMPTY <- synonymous_codons
synonymous_codons_EMPTY_high <- quantile(synonymous_codons$Selection_coefficient, probs = 0.975)
synonymous_codons_EMPTY_low <- quantile(synonymous_codons$Selection_coefficient, probs = 0.025)

out_mut_list <- c()
p_val_list <- c()
estimate_list <- c()

for(mutation in mutation_list){
  
  ## Select only values for this mutation
  dms_mut_vals <- dms_data_cds %>% ungroup() %>%
    filter(Mutation == mutation)
  
  ## For t-test compared against the distribution of synonymous codons
  ttest_result <- t.test(dms_mut_vals$Selection_coefficient, synonymous_codons$Selection_coefficient)
  p_val <- ttest_result$p.value
  mean_val <- ttest_result$estimate[[1]][1]
  
  if(!(is.null(p_val))){
    out_mut_list <- c(out_mut_list, mutation)
    p_val_list <- c(p_val_list, p_val)
    estimate_list <- c(estimate_list, mean_val)
  }
  
}

df_pvals <- cbind(out_mut_list, estimate_list, p_val_list)
df_pvals <- as.data.frame(df_pvals)
colnames(df_pvals) <- c('Mutation', 'Mean_ttest', 'pval')

df_volcano_plot <- df_pvals %>% 
  mutate(Selection_coefficient = as.numeric(Mean_ttest), 
         pval = as.numeric(pval))

# Get the mean effect of stop codons
stop_codons_wt <- dms_data_cds %>% ungroup() %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  filter(Mutant_Residue == '*', !(is.na(Selection_coefficient)), 
         !(Position %in% c('15', '78'))
  ) %>%
  ## Average over replicates
  group_by(Run, WT_Residue, Position, Mutant_Residue) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient)) %>%
  ## Average over runs
  ungroup() %>% group_by(WT_Residue, Position, Mutant_Residue) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

mean_stop <- mean(stop_codons_wt$Selection_coefficient)
sem_stop <-sd(stop_codons_wt$Selection_coefficient) / sqrt(nrow(stop_codons_wt))
sd_stop <- sd(stop_codons_wt$Selection_coefficient)

## Save the list of stop codons for the single copy system
stop_codons_single_copy <- stop_codons_wt

## Save the data for stop codons in EMPTY
mean_stop_EMPTY <- mean_stop
sd_stop_EMPTY <- sd_stop

stop_EMPTY_high <- mean_stop + cutoff_sd*sd_stop_EMPTY
stop_EMPTY_low <- mean_stop -  cutoff_sd*sd_stop_EMPTY

## Average over run and replicate for all codons
dms_data_final <- dms_data_cds %>% ungroup() %>% 
  filter(!(is.na(Selection_coefficient))) %>% 
  ## Average over replicates
  group_by(Run, Mutation) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient)) %>%
  ## Average over runs
  ungroup() %>% group_by(Mutation) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

## Count how many mutations are classified as deleterious
dms_data_final_deleterious <- dms_data_final %>% 
  mutate(del_check = Selection_coefficient <= synonymous_codons_EMPTY_low, 
         stop_like = and(Selection_coefficient <= stop_EMPTY_high, 
                         Selection_coefficient >= stop_EMPTY_low),
         neg_interference = Selection_coefficient < stop_EMPTY_low)  %>%
  ## Remove stops
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  filter(WT_Residue != Mutant_Residue, Mutant_Residue != '*')

sum(dms_data_final_deleterious$del_check)
# 1055 deleterious mutations including stops (Singleton)
# 986 deleterious missense mutations (Singleton)

sum(dms_data_final_deleterious$stop_like) 
# 833 stop-like mutations including stops (Singleton)
# 766 stop-like missense mutations (Singleton)

sum(dms_data_final_deleterious$neg_interference) # 2 "interfering" (Singleton)

## Do Benjamini-Hochberg correction first
df_volcano_plot %<>% mutate(p_adj = p.adjust(pval, method = 'BH'))

df_volcano_plot %<>% ungroup() %>% rowwise() %>%
  mutate(log_pval = -log10(pval), 
         log_pval_adj = -log10(p_adj))

## Get the list of mutations from the barplot to label them here
mutations_label <- c('Q67C', 'S59A', 'S59C')

data_candidates_label <- df_volcano_plot %>%
  filter(Mutation %in% mutations_label)

df_volcano_plot_new <- df_volcano_plot %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  mutate(stop_check = ifelse(Mutant_Residue == '*', 'Stop', 'Missense'))

df_stop <- df_volcano_plot_new %>% filter(stop_check == 'Stop')

## Annotate the residues that are in contact with the substrates ##
df_volcano_plot_regions <- left_join(x = df_volcano_plot_new %>% mutate(Position = as.numeric(Position)), 
                                     y = data_annotation_new, 
                                     by = c('Position' = 'Position'), relationship = 'many-to-many')

substrate_residues <- df_volcano_plot_regions %>% filter(Site %in% c('DHF binding', 'Binding only NADPH'))

df_volcano_plot_regions %<>% 
  mutate(color_legend = case_when(
    stop_check == 'Stop' ~ 'Stop',
    Site == 'DHF binding' ~ 'DHF binding', 
    Site == 'Binding only NADPH' ~ 'Binding only NADPH',
    TRUE ~ 'Missense'
  )) %>%
  mutate(color_legend = factor(color_legend, levels = c('Missense', 'Stop', 'DHF binding', 'Only NADPH')))


## For some mutants I get very high -log(p_val), reduce them to at most a cutoff
# pval_cutoff <- 15
pval_cutoff <- 10

df_volcano_plot_regions %<>%
  mutate(shape_bool = ifelse(log_pval_adj > pval_cutoff, TRUE, FALSE)) %>%
  mutate(log_pval_adj = ifelse(log_pval_adj > pval_cutoff, pval_cutoff, log_pval_adj))

df_stop %<>%
  mutate(shape_bool = ifelse(log_pval_adj > pval_cutoff, TRUE, FALSE)) %>%
  mutate(log_pval_adj = ifelse(log_pval_adj > pval_cutoff, pval_cutoff, log_pval_adj))

data_candidates_label %<>%
  mutate(shape_bool = ifelse(log_pval_adj > pval_cutoff, TRUE, FALSE)) %>%
  mutate(log_pval_adj = ifelse(log_pval_adj > pval_cutoff, pval_cutoff, log_pval_adj))

dms_data_final_EMPTY <- dms_data_final

## Duplicated  (Duplicated) ##

## Average by amino acid but not by run or replicate
dms_data_cds <- dms_data %>% ungroup() %>%
  filter(Background == 'Duplicated', !(is.na(Selection_coefficient))) %>%
  group_by(Run, Replicate, Mutation) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

mutation_list <- unique(dms_data_cds$Mutation)

dms_data_cds %<>% filter(!(is.na(Selection_coefficient)))

out_mut_list <- c()
p_val_list <- c()
estimate_list <- c()

## Get a distribution of synonymous codons
synonymous_codons <- dms_data %>% ungroup() %>%
  filter(Background == 'Duplicated', !(is.na(Selection_coefficient))) %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  separate(col = Genetic_code, into = c('WT_Codon', 'Position2', 'Mutant_Codon'), sep = c(3, -3)) %>%
  filter(and(WT_Residue == Mutant_Residue, WT_Codon != Mutant_Codon)) %>%
  group_by(WT_Residue, Position, Mutant_Residue) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

## Save the data for synonymous codons for CDS
dms_data_CDS <- dms_data_cds
synonymous_codons_CDS <- synonymous_codons
synonymous_codons_CDS_high <- quantile(synonymous_codons$Selection_coefficient, probs = 0.975)
synonymous_codons_CDS_low <- quantile(synonymous_codons$Selection_coefficient, probs = 0.025)


for(mutation in mutation_list){
  
  ## Select only values for this mutation
  dms_mut_vals <- dms_data_cds %>% ungroup() %>%
    filter(Mutation == mutation)
  
  ## For t-test compared against the distribution of synonymous codons
  ttest_result <- t.test(dms_mut_vals$Selection_coefficient, synonymous_codons$Selection_coefficient)
  p_val <- ttest_result$p.value
  mean_val <- ttest_result$estimate[[1]][1]
  
  if(!(is.null(p_val))){
    out_mut_list <- c(out_mut_list, mutation)
    p_val_list <- c(p_val_list, p_val)
    estimate_list <- c(estimate_list, mean_val)
  }
  
}

df_pvals <- cbind(out_mut_list, estimate_list, p_val_list)
df_pvals <- as.data.frame(df_pvals)
colnames(df_pvals) <- c('Mutation', 'Mean_ttest', 'pval')

df_volcano_plot <- df_pvals %>% 
  mutate(Selection_coefficient = as.numeric(Mean_ttest), 
         pval = as.numeric(pval))

# Get the mean effect of stop codons
stop_codons_wt <- dms_data_cds %>% ungroup() %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  filter(Mutant_Residue == '*', !(is.na(Selection_coefficient)), 
         !(Position %in% c('15', '78'))
  ) %>%
  ## Average over replicates
  group_by(Run, WT_Residue, Position, Mutant_Residue) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient)) %>%
  ## Average over runs
  ungroup() %>% group_by(WT_Residue, Position, Mutant_Residue) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))


mean_stop <- mean(stop_codons_wt$Selection_coefficient)
sem_stop <-sd(stop_codons_wt$Selection_coefficient) / sqrt(nrow(stop_codons_wt))
sd_stop <- sd(stop_codons_wt$Selection_coefficient)

## Save the selection coefficients for stop codons in the expeirment with the duplication
stop_codons_dup <- stop_codons_wt

## Save data for stop codons
mean_stop_CDS <- mean_stop
sd_stop_CDS <- sd_stop

stop_CDS_high <- mean_stop_CDS + cutoff_sd*sd_stop_CDS
stop_CDS_low <- mean_stop_CDS - cutoff_sd*sd_stop_CDS

## Average over run and replicate for all codons
dms_data_final <- dms_data_cds %>% ungroup() %>% 
  filter(!(is.na(Selection_coefficient))) %>% 
  ## Average over replicates
  group_by(Run, Mutation) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient)) %>%
  ## Average over runs
  ungroup() %>% group_by(Mutation) %>%
  summarise(Selection_coefficient = mean(Selection_coefficient))

dms_data_final_deleterious_CDS <- dms_data_final %>% 
  mutate(del_check = Selection_coefficient <= synonymous_codons_CDS_low, 
         stop_like = and(Selection_coefficient <= stop_CDS_high, 
                         Selection_coefficient >= stop_CDS_low),
         neg_interference = Selection_coefficient < stop_CDS_low) %>%
  ## Remove stops
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  filter(WT_Residue != Mutant_Residue, Mutant_Residue != '*')

sum(dms_data_final_deleterious_CDS$del_check) 
# 1036 deleterious mutations with stops (Duplicated)
# 967 deleterious missense mutations (Duplicated)

sum(dms_data_final_deleterious_CDS$stop_like) 
# 619 stop-like mutations with stops (Duplicated)
# 554 stop-like missense mutations (Duplicated)

sum(dms_data_final_deleterious_CDS$neg_interference) # 63 "interfering" (Duplicated)

## Do Benjamini-Hochberg correction first
df_volcano_plot %<>% mutate(p_adj = p.adjust(pval, method = 'BH'))

df_volcano_plot %<>% ungroup() %>% rowwise() %>%
  mutate(log_pval = -log10(pval), 
         log_pval_adj = -log10(p_adj))


data_candidates_label <- df_volcano_plot %>%
  filter(Mutation %in% mutations_label)

df_volcano_plot_new <- df_volcano_plot %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  mutate(stop_check = ifelse(Mutant_Residue == '*', 'Stop', 'Missense'))

df_stop <- df_volcano_plot_new %>% filter(stop_check == 'Stop')

## Annotate the residues that are in contact with the substrates ##

df_volcano_plot_regions <- left_join(x = df_volcano_plot_new %>% mutate(Position = as.numeric(Position)), 
                                     y = data_annotation_new, 
                                     by = c('Position' = 'Position'), relationship = 'many-to-many')

substrate_residues <- df_volcano_plot_regions %>% filter(Site %in% c('DHF binding', 'Binding only NADPH'))

df_volcano_plot_regions %<>% 
  mutate(color_legend = case_when(
    stop_check == 'Stop' ~ 'Stop',
    Site == 'DHF binding' ~ 'DHF binding', 
    Site == 'Binding only NADPH' ~ 'Binding only NADPH',
    TRUE ~ 'Missense'
  )) %>%
  mutate(color_legend = factor(color_legend, levels = c('Missense', 'Stop', 'DHF binding', 'Binding only NADPH')))

dms_data_final_CDS <- dms_data_final

dms_data_final_both <- inner_join(x = dms_data_final_CDS %>% 
                                    mutate(Selection_coefficient_Duplicated = Selection_coefficient) %>%
                                    select(-Selection_coefficient), 
                                  y = dms_data_final_EMPTY %>% 
                                    mutate(Selection_coefficient_Singleton = Selection_coefficient) %>%
                                    select(-Selection_coefficient), 
                                  by = c('Mutation'))

dms_data_final_both_long <- dms_data_final_both %>% ungroup() %>% 
  group_by(Mutation) %>%
  pivot_longer(cols = c('Selection_coefficient_Duplicated', 'Selection_coefficient_Singleton'), 
               names_to = 'Condition', names_prefix = 'Selection_coefficient_',
               values_to = 'Selection_coefficient') %>%
  mutate(Condition = factor(Condition, levels = c('Singleton', 'Duplicated')))


#### Figure 4A: Boxplots for the replicates of stop codons and interfering mutations ####

## Gather the DMS data for Q67C and stop codons (versus WT+WT)
dms_data_Q67C <- dms_data %>% ungroup() %>%
  filter(Background == 'Duplicated',
         # Mutation %in% c('Q67C', 'S59A', 'S59C', 'S59Y'),
         Mutation %in% c('Q67C', 'S59Y'),
         !(is.na(Selection_coefficient)), 
  )

# Get the stop codons
stop_codons_cds_noavg <- dms_data %>% 
  filter(Background == 'Duplicated') %>%
  ungroup() %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  filter(Mutant_Residue == '*', !(is.na(Selection_coefficient)), 
         !(Position %in% c('15', '78'))
  )

## Concatenate the two tables and draw figure
data_fig4A <- bind_rows(
  dms_data_Q67C %>% 
    select(-Ranked_sel_coeff) %>% 
    separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1), remove = F) %>%
    mutate(Mutation = str_c('WT+',Mutation, sep = '')),
  stop_codons_cds_noavg %>% mutate(Mutation = 'WT+Stop')
)

## Rename the samples to match our naming convention
data_fig4A %<>% separate(col = Mutation, into = c('tmp1', 'tmp2'), sep = '\\+') %>%
  rowwise() %>%
  mutate(Mutation = str_c(tmp2, tmp1, sep = '/'))

## Perform statistical test manually to assign p-values
comps <- compare_means(Selection_coefficient~Mutation, data = data_fig4A %>% ungroup(), 
                       paired = FALSE, 
                       method = 't.test'
                       ) %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p >= 0.01, str_c('p = ', round(p, 2), sep = ''), 
                                  sprintf("p = %2.1e", as.numeric(p)))
  )
  )


comps_stop <- comps %>% 
  arrange(factor(group1, levels = c('Q67C/WT', 'S59Y/WT')), 
          factor(group2, levels = c('S59Y/WT', 'Stop/WT'))) %>%
  mutate(y_pos = c(0.25, 0.4, 0.1))

p_fig4A <- data_fig4A %>% 
  ggplot(aes(x = Mutation, y = Selection_coefficient)) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  geom_boxplot(alpha = 0.3, width = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_signif(data = as.data.frame(comps_stop), inherit.aes = FALSE,
              aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 4) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab('') + 
  ylab('Selection coefficient') +
  ggtitle('Bulk competition (Genotype vs WT/WT)')+
  ylim(-1, 0.6)
p_fig4A

#### Figure 4B: qPCR validations of interfering candidates ####

candidate_validations_qpcr_wtwt_fig4B <- read_delim('Data/qPCR_direct_competition/qPCR_data_fig4B.tsv', delim = '\t')

## Rename according to our convention
candidate_validations_qpcr_wtwt_fig4B %<>% 
  separate(col = Sample, into = c('tmp1', 'tmp2'), sep = '\\+') %>%
  rowwise() %>%
  mutate(Sample = str_c(tmp2, tmp1, sep = '/'))

## Perform statistical test manually to assign p-values
comps <- compare_means(sel_coeff~Sample, 
                       data = candidate_validations_qpcr_wtwt_fig4B %>% ungroup(), 
                       paired = FALSE, 
                       method = 't.test'
                       ) %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p >= 0.01, str_c('p = ', round(p, 2), sep = ''), 
                                  sprintf("p = %2.1e", as.numeric(p)))
  )
  )

comps %<>% mutate(y_pos = c(0.1, 0.25, 0.4))

p_fig4B <- candidate_validations_qpcr_wtwt_fig4B %>%
  mutate(Sample = factor(Sample, levels = c('Q67C/WT', 'S59Y/WT', 'M1*/WT'))) %>%
  ggplot(aes(x = Sample, y = sel_coeff)) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_signif(data = as.data.frame(comps), inherit.aes = FALSE,
              aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 4) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab('') + 
  ylab('Selection coefficient') +
  ggtitle('Pairwise competition (Genotype vs WT/WT)') +
  ylim(-1, 0.6)
p_fig4B

#### Figures 4C and 4D were added manually with Inkscape ####

#### Figure 4E: Tm analysis ####

## Load the Tm data
data_tm <- read_excel('Data/Protein_complex_stability/ProteinStability_Calorimetry_Purif1.xlsx',
                      sheet = 'Overview', skip = 0)

colnames(data_tm) <- c('Capillary', 'Sample_ID', 'Onset_1', 'Inflection_point_1', 'Inflection_point_2')

## Separate the replicate from the sample type
data_tm %<>% separate(Sample_ID, into = c('Sample', 'Replicate'), sep = '_')

data_tm %<>% filter(!(is.na(Sample))) %>%
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
data_single_copy <- data_tm %>% 
  filter(!(str_detect(string = Sample_standard, pattern = '[\\/]')))

## Load the data from the second experiment
data_tm_exp2 <- read_excel('Data/Protein_complex_stability/ProteinStability_Calorimetry_Purif2.xlsx',
                           sheet = 'Overview', skip = 0)

data_tm_exp2 %<>%
  separate(col = `Sample ID`, into = c('Sample_standard', 'Replicate'), sep = '_') %>%
  mutate(Sample_standard = str_replace(string = Sample_standard, pattern = '-', replacement = '+')) %>%
  select(Capillary, Sample_standard, Replicate, `Onset #1 for Ratio`, 
         "Inflection Point #1 for Ratio", "Inflection Point #2 for Ratio")

colnames(data_tm_exp2) <- c('Capillary', 'Sample_standard', 'Replicate',
                            'Onset_1', 'Inflection_point_1', 'Inflection_point_2')

## Put the two experiments together
data_tm_exp_both <- rbind(data_single_copy %>% 
                            select(Capillary, Sample_standard, Replicate, Onset_1, 
                                   Inflection_point_1, Inflection_point_2) %>%
                            mutate(Purification = '1'), 
                          data_tm_exp2 %>%
                            select(Capillary, Sample_standard, Replicate, Onset_1, 
                                   Inflection_point_1, Inflection_point_2) %>%
                            mutate(Purification = '2')
)

## Compare using an ANOVA ##
m <- aov(Inflection_point_1~Purification + Sample_standard, 
         # data=data_single_copy
         data=data_tm_exp_both
         )
anova_test <- anova(m)
write.table(anova_test, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Supp_tables/TableS6_ANOVA_Tm_both_exp.tsv')

tukey_test <- HSD.test(m, trt = 'Sample_standard')
tukey_test2 <- TukeyHSD(m, 'Sample_standard')
write.table(tukey_test2$Sample_standard, 
            file = 'Supp_tables/TableS6_Tukey_test_Tm_both_exp.tsv',
            row.names = T, col.names = T, sep = '\t', append = F)

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(62, 62, 62, 62))
tukey_annotation$Sample_standard <- rownames(tukey_annotation)

p_fig4E <- data_tm_exp_both %>% 
  ggplot(aes(x = Sample_standard, y = Inflection_point_1, 
             shape = as.factor(Purification))) +
  geom_jitter(width = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = 'top', 
        legend.justification = 'center'
        ) +
  xlab('') +
  labs(y = expression(
    paste(bold(T[m]), 
    bold(' (°C)')
    )
    ), 
    shape ='Purification'
  ) +
  scale_shape_manual(values = c(16, 17)) +
  geom_text(data = tukey_annotation, size = 5, inherit.aes = F,
            aes(label = groups, x = as.factor(Sample_standard), y = y))
p_fig4E

#### Figure 4F: Growth recovery ####

## Load the growth curve data
plate.ind <- 'Data/Growth_recovery_variants/Growth_curves_2025_description.xlsx'
file.od1 <- 'Data/Growth_recovery_variants/Growth_curves_data_26_03_25.xlsx'
file.od2 <- 'Data/Growth_recovery_variants/Growth_curves_data_01_04_25.xlsx'
file.od3 <- 'Data/Growth_recovery_variants/Growth_curves_data_04_04_25.xlsx'
file.od4 <- 'Data/Growth_recovery_variants/Growth_curves_data_09_04_25.xlsx'

## Define function to read plate data
# function to process plates ----------------------------------------------
read.my.gc <- function(file, plate.index, time_trim){
  pl <- read_excel(file,sheet = 1, skip = 27, n_max = 90, col_names = T)
  ind <- read_excel(plate.index, sheet = 1, col_names = T) 
  colnames(ind) <- c('Experiment', 'Well', 'Strain', 'State', 'Arabinose', 'TMP', 'Replicate')
  
  time <- seq(0.25,0.25*(nrow(pl)), 0.25)
  pl$Time <- time
  
  # Pivot the table to a longer format
  pl_longer <- pl %>% 
    pivot_longer(cols = colnames(pl)[3:length(colnames(pl))], 
                 names_to = 'Well', values_to = 'OD')
  
  pl_gc <- pl_longer %>% select(Time, Well, OD)
  
  pl_gc_wide <- pl_gc %>% 
    pivot_wider(names_from = 'Well', values_from = 'OD') 
  
  data.od <- SummarizeGrowthByPlate(pl_gc_wide, t_trim = time_trim) %>%
    mutate(Experiment = basename(file))
  
  data.od_final <- inner_join(data.od, ind, 
                              by = c('Experiment' = 'Experiment', 'sample' = 'Well'))
  
  return(data.od_final)
  
}

## Run the growthcurver analyses
t_trim <- 13
data.od1 <- read.my.gc(file.od1, plate.ind, t_trim) %>%
  mutate(experiment = basename(file.od1))
data.od2 <- read.my.gc(file.od2, plate.ind, t_trim) %>%
  mutate(experiment = basename(file.od2))
data.od3 <- read.my.gc(file.od3, plate.ind, t_trim) %>%
  mutate(experiment = basename(file.od3))
data.od4 <- read.my.gc(file.od4, plate.ind, t_trim) %>%
  mutate(experiment = basename(file.od4))

## Concatenate the datasets
data.od.all <- bind_rows(data.od1, data.od2, data.od3, data.od4)

## Summarize all curves with the area under the curve
data.od.summary <- data.od.all %>%
  group_by(TMP, Arabinose, sample, Strain, State, Replicate, experiment) %>%
  summarise(auc = mean(auc_e))

## Look at the areas under the curve for all controls
data_wt_noTMP <- data.od.summary %>%
  filter(Strain == 'WT', TMP == 'No')

## Check the blank as well
data_blank <- data.od.summary %>%
  filter(Strain == 'Blank')

## Show side by side the AUC for each sample with and without TMP
data.od.summary_noTMP <- data.od.summary %>% ungroup() %>%
  filter(Strain == 'WT', TMP == 'No') %>%
  group_by(TMP, Strain, State, Arabinose, experiment) %>%
  summarise(mean_auc_noTMP = mean(auc), 
            max_auc_noTMP = max(auc))

data_od_complete <- inner_join(data.od.summary, data.od.summary_noTMP, 
                               by = c('State' = 'State', 
                                      'experiment' = 'experiment', 
                                      'Arabinose' = 'Arabinose'))

data_od_complete %<>% mutate(growth_recovery = 100*(auc/mean_auc_noTMP))

## Standardize the names
data_od_complete %<>% mutate(
  Strain_standard = case_when(
    State == 'Single copy' ~ str_replace(string = Strain.x, pattern = '[\\-]', replacement = '+'),
    Strain.x == 'Blank' ~ Strain.x, 
    State == 'Duplication with WT' ~ str_c(
      str_replace(string = Strain.x, pattern = '[\\-]', replacement = '+'), 
      '/WT', sep = '')
  )
)

## Comparison for the two-copy systems at 0.0015% arabinose
data_fig_0.0015 <- data_od_complete %>% 
  filter(TMP.x == 'Yes', Arabinose == 0.0015, State == 'Duplication with WT', 
         Strain_standard %in% c('Blank', 'WT/WT', 'Q67C/WT', 'S59Y/WT', 'Q67C+S59Y/WT')) %>%
  mutate(Strain_standard = factor(Strain_standard,
                                  levels = c('Blank', 'WT/WT',
                                             'Q67C/WT', 'S59Y/WT',  
                                             'Q67C+S59Y/WT'))
  )

## Add the data for the WT/WT without TMP as a reference
data_ref <- data.od.summary %>%
  filter(TMP == 'No', Arabinose == 0.0015, State == 'Duplication with WT', Strain == 'WT')

data_ref_complete <- inner_join(data_ref, data.od.summary_noTMP, 
                               by = c('State' = 'State', 
                                      'experiment' = 'experiment', 
                                      'Arabinose' = 'Arabinose'))

data_ref_complete %<>% mutate(growth_recovery = 100*(auc/mean_auc_noTMP)) %>%
  mutate(Strain_standard = 'WT/WT (no TMP)')

## Add to full table
data_fig_0.0015 <- bind_rows(data_fig_0.0015, data_ref_complete) %>%
  mutate(Strain_standard = factor(Strain_standard,
                                  levels = c('Blank','WT/WT (no TMP)', 'WT/WT',
                                             'Q67C/WT', 'S59Y/WT',  
                                             'Q67C+S59Y/WT'))
  )

## Use an ANOVA to compare across groups ##
m <- aov(growth_recovery~Strain_standard, data=data_fig_0.0015)
anova_test <- anova(m)
write.table(anova_test, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Supp_tables/TableS7_ANOVAFig4F_growthRecovery.tsv')

tukey_test <- HSD.test(m, trt = 'Strain_standard')
tukey_test2 <- TukeyHSD(m)
write.table(tukey_test2$Strain_standard, 
            file = 'Supp_tables//TableS7_TukeyFig4F_growthRecovery.tsv',
            row.names = T, col.names = T, sep = '\t', append = F)

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = rep(110, n = 6))
tukey_annotation$Strain_standard <- rownames(tukey_annotation)

p_fig4F <- data_fig_0.0015 %>%
  ggplot(aes(x = Strain_standard, y = growth_recovery)) +
  geom_jitter(width = 0.2) +
  xlab('') + 
  ylab('Growth recovery (%)') +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(data = tukey_annotation, size = 5, inherit.aes = F,
            aes(label = groups, x = as.factor(Strain_standard), y = y))
p_fig4F

#### Put the four panels together ####

p_fig4 <- plot_grid(
  p_fig4A, p_fig4B,
  p_fig4E, p_fig4F, nrow = 2,
  labels = c('A', 'B', 'E', 'F'),
  label_size = panel_label_size, label_fontface = 'bold', 
  rel_widths = c(1, 1))

ggsave(p_fig4, device = cairo_pdf, dpi = 300, width = 22, height = 18, units = 'cm', 
       filename = 'Figures/Main_figures/Fig4A-B_E-F.pdf')
## Panels C and D were added with Inkscape

#### Fig. S5 ####

data_fig_0.0015_sc <- data_od_complete %>% 
  filter(TMP.x == 'Yes', Arabinose == 0.0015, State == 'Single copy', 
         Strain_standard %in% c('Blank', 'WT', 'Q67C', 'S59Y', 'Q67C+S59Y')) %>%
  mutate(Strain_standard = factor(Strain_standard,
                                  levels = c('Blank', 'WT',
                                             'Q67C', 'S59Y',  
                                             'Q67C+S59Y'))
  )

## Add the data for the WT/WT without TMP as a reference
data_ref_sc <- data.od.summary %>%
  filter(TMP == 'No', Arabinose == 0.0015, State == 'Single copy', Strain == 'WT')

data_ref_complete_sc <- inner_join(data_ref_sc, data.od.summary_noTMP, 
                                by = c('State' = 'State', 
                                       'experiment' = 'experiment', 
                                       'Arabinose' = 'Arabinose'))

data_ref_complete_sc %<>% mutate(growth_recovery = 100*(auc/mean_auc_noTMP)) %>%
  mutate(Strain_standard = 'WT (no TMP)')

## Add to full table
data_fig_0.0015_sc <- bind_rows(data_fig_0.0015_sc, data_ref_complete_sc) %>%
  mutate(Strain_standard = factor(Strain_standard,
                                  levels = c('Blank','WT (no TMP)', 'WT',
                                             'Q67C', 'S59Y',  
                                             'Q67C+S59Y'))
  )


## Use an ANOVA to compare across groups ##

m <- aov(growth_recovery~Strain_standard, data=data_fig_0.0015_sc)
anova_test <- anova(m)
write.table(anova_test, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Supp_tables/TableS5_ANOVA_growthRecovery_singleCopy.tsv')

tukey_test <- HSD.test(m, trt = 'Strain_standard')
tukey_test2 <- TukeyHSD(m)
write.table(tukey_test2$Strain_standard, 
            # file = 'Figures/2025-04-09_final_experiments/Tukey_test_ara0.0015_t13_singleCopy.tsv', 
            file = 'Supp_tables/TableS5_Tukey_test_growthRecovery_singleCopy.tsv', 
            row.names = T, col.names = T, sep = '\t', append = F)

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = rep(115, 6))
tukey_annotation$Strain_standard <- rownames(tukey_annotation)

p_figS5 <- data_fig_0.0015_sc %>%
  ggplot(aes(x = Strain_standard, y = growth_recovery)) +
  geom_jitter(width = 0.2) +
  xlab('') + 
  ylab('Growth recovery (%)') +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(data = tukey_annotation, size = 5, inherit.aes = F,
            aes(label = groups, x = as.factor(Strain_standard), y = y))
p_figS5
ggsave(p_figS5, device = cairo_pdf, width = 16, height = 10, units = 'cm', dpi = 300, 
       filename = 'Figures/Supp_figures/FigS5.pdf')

