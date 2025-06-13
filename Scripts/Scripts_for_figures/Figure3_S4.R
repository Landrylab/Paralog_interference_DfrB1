##################################################################
####                Figures 3, S4, S5                         ####
##################################################################

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
setwd('/path/to/Paralog_interference_DfrB1')

#### Panels A-C from the negative interference script ####

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

## Save this data frame
write.table(dms_data_final_both, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/dms_data_final_both.tsv')

dms_data_final_both_long <- dms_data_final_both %>% ungroup() %>% 
  group_by(Mutation) %>%
  pivot_longer(cols = c('Selection_coefficient_Duplicated', 'Selection_coefficient_Singleton'), 
               names_to = 'Condition', names_prefix = 'Selection_coefficient_',
               values_to = 'Selection_coefficient') %>%
  mutate(Condition = factor(Condition, levels = c('Singleton', 'Duplicated')))

## Add the distribution of stop codons as a point with errorbars
df_stop_info <- data.frame(Condition = c('Singleton', 'Duplicated'), 
                           mean_stop = c(mean_stop_EMPTY, mean_stop_CDS),
                           stop_high = c(stop_EMPTY_high, stop_CDS_high),
                           stop_low = c(stop_EMPTY_low, stop_CDS_low))

p_fig3A_boxplot <- dms_data_final_both_long %>% 
  ## Add a label to make interfering mutations blue
  mutate(color_bool = ifelse(and(Selection_coefficient < stop_CDS_low, Condition == 'Duplicated'),
                             'Negative interference', 'Other'
                             )
         ) %>%
  ggplot(aes(x = Condition, y = Selection_coefficient)) +
  geom_line(aes(group = Mutation), alpha = 0.1) +
  geom_jitter(width = 0.025, aes(colour = color_bool)) +
  scale_colour_manual(values = c('#1a53ff', 'black')) +
  geom_violin(width = 0.7, alpha = 0.3, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha= 0.3) +
  xlab('') + ylab('Selection coefficient') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(data = df_stop_info, colour = 'red', size = 3, aes(y = mean_stop)) +
  geom_errorbar(data = df_stop_info, colour = 'red', inherit.aes = F,
                aes(x = Condition, ymax = stop_high, ymin = stop_low), 
                width = 0.1, linewidth = 1.2
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = 'none', 
    legend.justification = 'center'
  ) +
  stat_compare_means( 
                     method = 'wilcox.test',
                     paired = T, aes(group = Mutation),
                     comparisons = list(c('Singleton', 'Duplicated')),
                     size = 4
  ) +
  ylim(-1, 0.7) +
  labs(colour = '') +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))
p_fig3A_boxplot

## Check which mutants in single copy were below the lower end of the stop codon distribution
single_copy_DN <- dms_data_final_both_long %>% 
  filter(Condition == 'Singleton', Selection_coefficient <= stop_EMPTY_low)
# Mutation Condition   Selection_coefficient
# W38M     Single copy                -0.946
# W38Q     Single copy                -0.911

## Same for the CDS data
neg_dom_check <- dms_data_final_both_long %>% 
  filter(Condition == 'Duplicated', Selection_coefficient <= stop_CDS_low)
## 63 mutations

#### Figure 3C: Enrichment of negative interference mutations ####

## Filter to get a count of missense mutations
missense_mut <- dms_data_final %>%
  separate(Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  filter(Mutant_Residue != '*', WT_Residue != Mutant_Residue)
# 1310 missense mutations

## Relative enrichment of negative interference mutations in the duplicated background ##

dom_eff_candidates <- dms_data_final %>%
  separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), sep = c(1, -1)) %>%
  ## Filter by mean of selection coefficient minus 2.5 sd
  filter(Selection_coefficient < (mean_stop-cutoff_sd*sd_stop), Mutant_Residue != '*')

write.table(dom_eff_candidates, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Supp_tables/TableS4_DominantNegativeCandidates.tsv')

dom_eff_candidates_regions <- left_join(x = dom_eff_candidates %>% mutate(Position = as.numeric(Position)),
                                        y = data_annotation_new, 
                                        by = c('Position' = 'Position'), relationship = 'many-to-many')

## Alternatively, filter by both the magnitude of the fitness effect and the significance
dom_eff_candidates_regions <- df_volcano_plot_regions %>% 
  ## Filter based on the mean of stop codons - 2.5*(sd of stop codons)
  filter(Selection_coefficient < (mean_stop-2.5*sd_stop), p_adj <= 0.05)

dom_eff_candidates_regions_summary <- dom_eff_candidates_regions %>% ungroup() %>%
  mutate(Site = factor(Site, levels = unique(data_annotation_new$Site))) %>%
  group_by(Site) %>%
  summarise(count_mut = n()) %>%
  complete(Site, fill = list(count_mut = 0))

## Check the number of positions in each region, excluding positions 1-9 because they are not included in the DMS
positions_per_region <- data_annotation_new %>% ungroup() %>%
  mutate(Site = factor(Site, levels = unique(data_annotation_new$Site))) %>%
  filter(Position >= 10) %>%
  group_by(Site) %>%
  summarise(position_count = n()) %>%
  complete(Site, fill = list(count_mut = 0))

positions_per_region %<>% mutate(available_missense_mut = position_count*19)

test_positions <- data_annotation_new %>% ungroup() %>%
  filter(Position >= 10)
length(table(unique(test_positions$Position)))
# We do have 69 positions

## Merge the number of total positions and the number of dominant negative candidates
dom_eff_abundance <- inner_join(x = dom_eff_candidates_regions_summary, 
                                y = positions_per_region, 
                                by = c('Site' = 'Site')
)

dom_eff_abundance_long <- dom_eff_abundance %>% 
  select(-position_count) %>%
  pivot_longer(cols = c('count_mut', 'available_missense_mut'), names_to = 'Count_type', values_to = 'Count_val')

## Add an estimation of observed vs expected ##
total_candidates <- nrow(dom_eff_candidates)
total_missense <- length(table(unique(test_positions$Position))) * 19

dom_eff_abundance %<>% 
  mutate(expected_dom_eff = round(available_missense_mut*(total_candidates / total_missense), 2))

dom_eff_abundance_long <- dom_eff_abundance %>% 
  select(-position_count) %>%
  pivot_longer(cols = c('count_mut', 'expected_dom_eff'), 
               names_to = 'Count_type', values_to = 'Count_val')

dom_eff_abundance_long <- dom_eff_abundance %>% 
  select(-position_count) %>%
  pivot_longer(cols = c('count_mut', 'expected_dom_eff', 'available_missense_mut'), 
               names_to = 'Count_type', values_to = 'Count_val')

## Another alternative, this time as a relative enrichment (observed / expected) ##
dom_eff_abundance %<>% 
  mutate(relative_enrichment = count_mut / expected_dom_eff, 
         Background = 'Duplicated') %>%
  rowwise() %>%
  mutate(
    label_enrichment = str_c(toString(round(relative_enrichment, 2)), '\n', 
                             '(', toString(count_mut), '/', toString(round(expected_dom_eff, 2)), ')')
  )

max_rel_enrichment <- max(dom_eff_abundance$relative_enrichment)

## A variable for the position of labels on Fig 3C
rel_enrichment_margin <- 3

## Refactor according to the enrichment value
dom_eff_abundance_order <- dom_eff_abundance  %>% arrange(relative_enrichment)

dom_eff_abundance %<>%
  mutate(Site = factor(Site, levels = dom_eff_abundance_order$Site))

## Draw figure 3C 
p_fig3C <- dom_eff_abundance %>%
  mutate(Site = case_when(
    Site == 'Disordered region' ~ 'Disordered\nregion',
    Site == 'Dimerization interface' ~ 'Dimerization\ninterface',
    Site == 'Binding only NADPH' ~ 'Binding only\nNADPH',
    Site == 'Buried residues' ~ 'Buried\nresidues',
    Site == 'Tetramerization interface' ~ 'Tetramerization\ninterface',
    Site == 'NADPH binding' ~ 'NADPH\nbinding',
    Site == 'Cat. residues/DHF binding' ~ 'Cat. residues/\nDHF binding', 
    Site == 'Unannotated' ~ 'Unannotated'
  )
  ) %>%
  filter(Site != 'Binding only\nNADPH') %>%
  mutate(Site = factor(Site, 
                       levels = c('Disordered\nregion', 'Unannotated', 'Dimerization\ninterface',
                                  'Buried\nresidues', 'Tetramerization\ninterface', 
                                  'NADPH\nbinding', 'Cat. residues/\nDHF binding')
  )
  ) %>%
  ggplot(aes(x = Site, y = relative_enrichment)) +
  geom_bar(stat = 'identity', aes(fill = Site)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_label(aes(label = label_enrichment, x = Site, y = relative_enrichment + 1),
             vjust = 0.5, size = 3.5
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = 'none'
  ) +
  xlab('') +
  ylab('Relative enrichment\n(neg. int. mutations)') +
  ylim(0, max_rel_enrichment + rel_enrichment_margin) +
  scale_fill_manual(values = c('#000000', # Disordered region
                               '#999999', # Unannotated
                               '#8803fc', # Dimerization interface
                               '#ff4d4d', # Buried residues
                               '#ff6600', # Tetramerization interface
                               '#56B4E9', # NADPH binding
                               '#E69F00' # Cat residues/DHF binding
  )
  )
p_fig3C
## Cartoon structures were added manually with inkscape

## Prepare the annotation ofor residues belonging to each region
ha2 <- HeatmapAnnotation(
  `Disordered region` = data_annotation_dfrb1$Disordered_region,
  `Dimerization interface` = data_annotation_dfrb1$`A,C`,
  `Buried residues` = data_annotation_dfrb1$Buried,
  `Tetramerization interface` = data_annotation_dfrb1$`A,D`,
  `NADPH binding` = data_annotation_dfrb1$NADPH,
  `Cat. residues/DHF binding` = data_annotation_dfrb1$Cat_residues,
  show_annotation_name = T,
  annotation_name_gp = gpar(fontface = 'bold', fontsize = 12),
  simple_anno_size = unit(0.75, 'cm'),
  annotation_name_side = 'left',
  show_legend = FALSE,
  col = list(
    `Disordered region` = colorRamp2(c(0, 1), c("white", "#999999")),
    `Dimerization interface` = colorRamp2(c(0, 1), c("white", "#8803fc")),
    `Buried residues` = colorRamp2(c(0, 1), c("white", "#ff4d4d")),
    `Tetramerization interface` = colorRamp2(c(0, 1), c("white", "#ff6600")),
    `NADPH binding` = colorRamp2(c(0, 1), c("white", "#56B4E9")),
    `Cat. residues/DHF binding` = colorRamp2(c(0, 1), c("white", "#E69F00"))
  ),
  gp = gpar(col = "black")
)

## Draw the annotation as a panel
m = matrix(seq(from = 77, to = 1, by = -1), 1)

column_labels = c('', '', '', '5',
                  '', '', '', '', '10',
                  '', '', '', '', '15',
                  '', '', '', '', '20',
                  '', '', '', '', '25',
                  '', '', '', '', '30',
                  '', '', '', '', '35',
                  '', '', '', '', '40',
                  '', '', '', '', '45',
                  '', '', '', '', '50',
                  '', '', '', '', '55',
                  '', '', '', '', '60',
                  '', '', '', '', '65',
                  '', '', '', '', '70',
                  '', '', '', '', '75',
                  '', '', ''
)

p_annotation_heatmap <- draw_CHeatmap(Heatmap(m, bottom_annotation = ha2, column_labels = column_labels))
## The annotation was taken from this figure and copied to figure 3C with Inkscape

#### Figure 3B: FoldX effects of negative interference mutations ####

## Load the list of dominant negative mutations
dom_neg_mut <- dom_eff_candidates

## Read the FoldX mutational effects (taken from Cisneros, et al., 2023. Sci. Adv.)
foldx_effects <- read_delim('Data/Mutational_effects/sciadv.add9109_table_s3.tsv', delim = '\t')

foldx_effects %>% filter(Position == 59, WT_Residue == 'S', Residue == 'Y') %>%
  select(Mean_ddG_int_HM_A_D)
## S59Y strongly destabilizes the tetramerization interface (ddG = 29)

## Load the table with all the selection coefficients
all_dms_cds <- df_volcano_plot_regions

#### Merge the effects ####
foldx_effects_summary <- foldx_effects %>% ungroup() %>% 
  filter(Residue != '*', Position >= 22) %>%
  group_by(Position, WT_Residue, Residue) %>%
  summarise(Mean_ddG_stab_HET = mean(as.numeric(Mean_ddG_stab_HET)), 
            Mean_ddG_int_HM_A_C = mean(as.numeric(Mean_ddG_int_HM_A_C)), 
            Mean_ddG_int_HM_A_D = mean(as.numeric(Mean_ddG_int_HM_A_D)), 
            
            Mean_ddG_int_HET_A_C = mean(as.numeric(Mean_ddG_int_HET_A_C)),
            Mean_ddG_int_HET_A_D = mean(as.numeric(Mean_ddG_int_HET_A_D))
  )

## Merge the dataframes
all_dms_foldx <- inner_join(x = all_dms_cds %>% select(Position, WT_Residue, Mutant_Residue, 
                                                       Selection_coefficient, pval, p_adj, Site), 
                            y = foldx_effects_summary, 
                            by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                   'Mutant_Residue' = 'Residue'))

#### Load the Rosetta data for binding to NADPH or DHF ####
rosetta_data_nadph <- read_delim('Data/Mutational_effects/all_results_DMS_NADPH_formatted.tsv', 
                             delim = '\t')
rosetta_data_dhf_nadph <- read_delim('Data/Mutational_effects/all_results_DMS_NADPH_DHF_formatted.tsv', 
                                 delim = '\t')

## Concatenate the Rosetta results
rosetta_data_tmp <- left_join(x = rosetta_data_dhf_nadph %>%
                                mutate(total_score_dhf = total_score) %>%
                                select(mutation_label, state, energy_unit, total_score_dhf) %>%
                                filter(state == 'ddg') %>%
                                group_by(mutation_label, state, energy_unit) %>%
                                summarise(total_score_dhf = mean(total_score_dhf)), 
                              y = rosetta_data_nadph %>%
                                mutate(total_score_nadph = total_score) %>%
                                select(mutation_label, state, energy_unit, total_score_nadph) %>%
                                filter(state == 'ddg') %>%
                                group_by(mutation_label, state, energy_unit) %>%
                                summarise(total_score_nadph = mean(total_score_nadph)), 
                              by = c('mutation_label' = 'mutation_label', 'state' = 'state', 
                                     'energy_unit' = 'energy_unit')
)

## Distinguish among configurations of corresponding homo- or heterotetramers
rosetta_data_final <- rosetta_data_tmp %>%
  separate(mutation_label, into = c('mutA', 'mutB', 'mutC', 'mutD'), sep = '_', remove = F) %>%
  mutate(het_type = case_when(
    mutA != mutB & mutB == mutC & mutC == mutD & mutA == 'WT' ~ 'het_mut3_WT1', 
    mutA != mutB & mutB == mutC & mutC == mutD & mutA != 'WT' ~ 'het_mut1_WT3',
    mutA != mutB & mutA == mutD & mutB == mutC & mutA != 'WT' ~ 'het_mut2_AD', # Tetramerization interface
    mutA != mutB & mutA == mutC & mutB == mutD & mutA != 'WT' ~ 'het_mut2_AC', # Dimerization interface
    mutA == mutB & mutB != mutC & mutC == mutD & mutA != 'WT' ~ 'het_mut2_AB', # Crossed
    mutA != 'WT' & mutB != 'WT' & mutC != 'WT' & mutD != 'WT' ~ 'hm_mut'
  )
  ) %>%
  # Prepare columns so that we can merge by the mutation
  mutate(Mutation = case_when(
    mutA != 'WT' ~ mutA, 
    mutB != 'WT' ~ mutB,
    mutC != 'WT' ~ mutC,
    mutD != 'WT' ~ mutD
  )) %>%
  separate(Mutation, into =  c('WT_Residue', 'Position', 'Residue'), sep = c(1, -1)) %>%
  mutate(Position = as.numeric(Position))

## Filter effects on binding when the mutation is present on all four chains
rosetta_data_final_hm <- rosetta_data_final %>% ungroup() %>% filter(het_type == 'hm_mut') %>%
  mutate(Mean_ddG_NADPH = total_score_nadph, 
         Mean_ddG_DHF = total_score_dhf) %>%
  select(WT_Residue, Position, Residue, Mean_ddG_NADPH, Mean_ddG_DHF)

all_dms_foldx <- left_join(x = all_dms_foldx, 
                           y = rosetta_data_final_hm, 
                           by = c('WT_Residue' = 'WT_Residue', 
                                  'Position' = 'Position',
                                  'Mutant_Residue' = 'Residue')
                           )

## Indicate the dominant negative mutations in this dataset
dom_neg_mut %<>% rowwise() %>%
  mutate(mut_str = str_c(WT_Residue, Position, Mutant_Residue, sep = ''))

all_dms_foldx_DN <- all_dms_foldx %>% rowwise() %>%
  mutate(mut_str = str_c(WT_Residue, Position, Mutant_Residue, sep = '')) %>%
  mutate(DN_bool = ifelse(mut_str %in% dom_neg_mut$mut_str, 'Neg. interference', 'Other'))

#### Let's compare the distributions of effects of dominant negative mutations versus everything else ####

all_dms_foldx_DN_long <- all_dms_foldx_DN %>% ungroup() %>% 
  select(-Mean_ddG_int_HET_A_C, -Mean_ddG_int_HET_A_D) %>%
  group_by(mut_str, DN_bool, Site) %>%
  pivot_longer(cols = c('Mean_ddG_stab_HET', 'Mean_ddG_int_HM_A_C', 'Mean_ddG_int_HM_A_D', 
                        'Mean_ddG_NADPH', 'Mean_ddG_DHF'), 
               names_to = 'ddG_type', values_to = 'ddG_value') %>%
  mutate(ddG_type = case_when(
    ddG_type == 'Mean_ddG_stab_HET' ~ 'Folding', 
    ddG_type == 'Mean_ddG_int_HM_A_C' ~ 'Dimerization',
    ddG_type == 'Mean_ddG_int_HM_A_D' ~ 'Tetramerization',
    ddG_type == 'Mean_ddG_NADPH' ~ 'NADPH binding',
    ddG_type == 'Mean_ddG_DHF' ~ 'DHF binding'
  )
  ) %>% 
  mutate(ddG_type = factor(ddG_type, levels = c('Folding', 'Dimerization', 'Tetramerization', 
                                                'NADPH binding', 'DHF binding')))

## Make sure mutations are not repeated when they appear in multiple regions
all_dms_foldx_DN_long_unique <- all_dms_foldx_DN_long %>% ungroup() %>%
  group_by(Position, WT_Residue, Mutant_Residue, DN_bool, ddG_type) %>%
  summarise(ddG_value = mean(ddG_value)) %>%
  ungroup() %>%
  mutate(ddG_type = factor(ddG_type)) %>%
  group_by(ddG_type)

## Organize the data to estimate p-values
comps <- compare_means(ddG_value~DN_bool, 
                       data = all_dms_foldx_DN_long_unique,
                       paired = F, 
                       method = 'wilcox.test',
                       group.by = 'ddG_type') %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p > 0.01, str_c('p = ', round(as.numeric(p), 2), sep = ''),
                                  sprintf("p = %2.1e", as.numeric(p)))
  ), 
  y_pos = c(11, 11, 11, 11, 11)
  ) %>%
  mutate(ddG_type = factor(ddG_type,
                           levels = c('Folding', 'Dimerization', 'Tetramerization', 
                                      'NADPH binding', 'DHF binding')))

mut_counts <- all_dms_foldx_DN_long_unique %>% ungroup() %>% 
  select(Position, Mutant_Residue, WT_Residue, DN_bool) %>%
  unique() %>% group_by(DN_bool) %>%
  summarise(count = n())
# 63 negative interference mutations
# 1076 other mutations

## Draw the figure
p_fig3B <- all_dms_foldx_DN_long_unique %>% ungroup() %>%
  mutate(ddG_type = factor(ddG_type,
                           levels = c('Folding', 'Dimerization', 'Tetramerization', 
                                      'NADPH binding', 'DHF binding')), 
         DN_bool = case_when(
           DN_bool == 'Neg. interference' ~ 'Neg. interference\n(n = 63)',
           DN_bool == 'Other' ~ 'Other (n = 1076)'
         )) %>%
  ggplot(aes(x = DN_bool, y = ddG_value)) +
  facet_wrap(~factor(ddG_type, levels = c('Folding', 'Dimerization', 'Tetramerization', 
                                          'NADPH binding', 'DHF binding')),
             ncol = 5) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.3, aes(fill = DN_bool)) +
  geom_violin(width = 0.7, alpha = 0.3, scale = 'width', aes(fill = DN_bool)) +
  scale_fill_manual(values = c('#1a53ff', 'black')) +
  scale_colour_manual(values = c('#1a53ff', 'black')) +
  ylim(-2, 12) +
  theme(legend.position = 'none', legend.justification = 'center') +
  labs(fill = '', colour = '', 
       y = expression(bold('\u0394\u0394G'))) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('') +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = 'none'
  ) +
  geom_signif(data = as.data.frame(comps) %>% ungroup() %>%
                mutate(ddG_type = as.character(ddG_type), 
                       group1 = as.character(group1), 
                       group2 = as.character(group2)
                       ) %>%
                mutate(group1 = case_when(
                  group1 == 'Neg. interference' ~ 'Neg. interference\n(n = 63)',
                  group1 == 'Other' ~ 'Other (n = 1076)'
                ), 
                group2 = case_when(
                  group2 == 'Neg. interference' ~ 'Neg. interference\n(n = 63)',
                  group2 == 'Other' ~ 'Other (n = 1076)'
                )),
              aes(xmin = group1, xmax = group2,
                  annotations=p.format,
                  y_position = y_pos
                  ),
              manual = TRUE, textsize = 4, 
              tip_length = 1)
p_fig3B
       
## Draw figure 3
p_fig3_top <- plot_grid(p_fig3A_boxplot, p_fig3B, nrow = 1,
                        labels = c('A', 'B'),
                        label_size = panel_label_size, label_fontface = 'bold',
                        rel_widths = c(0.35, 1)
)

p_fig3 <- plot_grid(p_fig3_top, p_fig3C, nrow = 2,
                    labels = c('', 'C'),
                    label_size = panel_label_size, label_fontface = 'bold',
                    rel_heights = c(1, 1))

ggsave(p_fig3, device = cairo_pdf,
       width = 26, height = 20,
       dpi = 300, units = 'cm',
       filename = 'Figures/Main_figures/Fig3_noAnnotations.pdf'
       
)
## Annotations and cartoons were added manually with Inkscape

#### Prepare supplementary table 2 with the DMS and biophysical effects ####
table_s2 <- left_join(dms_data_final_both %>%
                        separate(col = Mutation, into = c('WT_Residue', 'Position', 'Mutant_Residue'), 
                                  sep = c(1, -1)) %>%
                        mutate(Position = as.numeric(Position)), 
                      all_dms_foldx %>% 
                        select(Position, WT_Residue, Mutant_Residue, 
                               Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D, 
                               Mean_ddG_NADPH, Mean_ddG_DHF) %>%
                        unique(), 
                      by = c('WT_Residue' = 'WT_Residue', 'Position' = 'Position', 
                             'Mutant_Residue' = 'Mutant_Residue')) %>%
  arrange(Position, Mutant_Residue)

colnames(table_s2) <- c("WT_Residue", "Position", "Mutant_Residue", 
                        "Selection_coefficient_Duplicated", "Selection_coefficient_Singleton", 
                        "Foldx_ddG_folding", "Foldx_ddG_dimerization", "Foldx_ddG_tetramerization", 
                        "FlexddG_ddG_NADPH", "FlexddG_ddG_DHF")

write.table(table_s2, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Supp_tables/TableS2_DMS_s_biophysical_effects.tsv')

#### Fig. S4: Comparison of effects to other mutations from the same region ####

mutation_counts <- all_dms_foldx_DN_long %>% ungroup() %>%
  filter(Site %in% c('Cat. residues/DHF binding', 'NADPH binding', 'Tetramerization interface')) %>%
  group_by(Site, DN_bool, ddG_type) %>%
  summarise(count = n()) %>% rowwise() %>%
  mutate(label = str_c('n = ', count, sep = ''))

p_figS4 <- all_dms_foldx_DN_long %>%
  filter(Site %in% c('Cat. residues/DHF binding', 'NADPH binding', 'Tetramerization interface')) %>%
  ggplot(aes(x = DN_bool, y = ddG_value)) +
  facet_grid(ddG_type~Site) +
  geom_jitter(width = 0.15, aes(colour = DN_bool), alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.3, aes(fill = DN_bool)) +
  geom_violin(width = 0.7, alpha = 0.3, scale = 'width', aes(fill = DN_bool)) +
  scale_fill_manual(values = c('#1a53ff', 'black')) +
  scale_colour_manual(values = c('#1a53ff', 'black')) +
  ylim(-5, 10) + 
  theme(legend.position = 'none', legend.justification = 'center', 
        strip.text = element_text(size = 10)) +
  labs(fill = '', colour = '', 
       y = expression(bold('\u0394\u0394G'))) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('') +
  stat_compare_means(paired = F, 
                     method = 'wilcox.test',
                     comparisons = list(c('Neg. interference', 'Other')
                     ), label.y = 7.5, tip.length = 0.005, 
                     aes(label = str_c('p = ', ..p.format.., sep = ''))
  ) +
  geom_text(data = mutation_counts, y = -4, aes(label = label))
p_figS4
ggsave(p_foldx_regions, device = cairo_pdf, width = 20, height = 24, units = 'cm', dpi = 300,  
       filename = 'Figures/Supp_figures/FigS4_ddG_sameRegions.pdf')

