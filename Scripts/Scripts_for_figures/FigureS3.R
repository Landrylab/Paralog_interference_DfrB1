###########################################################
####                    FigureS3                       ####
###########################################################

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

## Set the path to the main directory with the results
# setwd('/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Angel/PhD_projects/Yacine/Manuscript/Github/Paralog_interference_DfrB1/')
setwd('/path/to/Paralog_interference_DfrB1')

## Load the cytometry data
metadata <- read_excel('Data/Flow_cytometry/Cytometry_2025_description.xlsx', 
                       sheet = 'Sample_description', skip = 0)

colnames(metadata) <- c('Experiment', 'Well', 'Strain', 'Arabinose', 'Replicate')

path_files <- 'Data/Flow_cytometry/2025-03-04_at_11-36-00am/'
list_files <- list.files(path_files)

all_data_expression <- c()

for(infile in list_files){
  # Load the data
  new_data <- read_delim(file.path(path_files, infile), delim = ',', col_names = T) %>%
    select(TIME, 'GRN-B-HLin', 'FSC-HLin', 'SSC-HLin', 'FSC-HLog', 'SSC-HLog')
  
  # Add the name of the file as an ID
  new_data %<>% mutate(ID = substr(x = infile, start = 1, stop = (nchar(infile) - 4)))
  
  all_data_expression <- bind_rows(all_data_expression, new_data)
}

# Add the annotation data to identify each sample
all_data_expression %<>% separate(col = ID, into = c('tmp', 'Well'), sep = 'am.')

## Remove the extra zeroes from the well names in the cytometry data
all_data_expression %<>%
  mutate(Well_new = str_replace(string = Well, pattern = '([A-Z])0([0-9])', replacement = '\\1\\2'))

all_data_expression <- inner_join(x = all_data_expression, 
                                  y = metadata,
                                  by = c('Well_new' = 'Well'))


# Add the logarithms of the green fluorescence
all_data_expression %<>% mutate(GRNBHLog = ifelse(log10(`GRN-B-HLin`) < 0, 0, log10(`GRN-B-HLin`)),
                                `FSC-HLog` = ifelse(log10(`FSC-HLin`) < 0, 0, log10(`FSC-HLin`)),
                                `SSC-HLog` = ifelse(log10(`SSC-HLin`) < 0, 0, log10(`SSC-HLin`))
)

all_data_processed_expression <- all_data_expression %>% rowwise() %>%
  mutate(circle_test = (`FSC-HLog` - 1.25)^2 + (`SSC-HLog` - 1.25)^2) %>%
  filter(circle_test < 1) %>%
  mutate(
    Arabinose = str_c(toString(Arabinose), '% arabinose', sep = '')
  )

## Calculating the mean fluorescence in each replicate and the number of cells for each
all_data_processed_summary_expression <- all_data_processed_expression %>%
  ungroup() %>%
  group_by(Strain, Arabinose, Well_new) %>%
  summarise(
    mean_fluo = mean(GRNBHLog),
    sem_fluo = sd(GRNBHLog) / sqrt(n()),
    num_cells = n()
  )
## All wells had between 4691 and 5070 cells, pretty close to the target of 5000

all_data_processed_expression %<>%
  mutate(Strain = factor(Strain, 
                         levels = c('WT (neg ctl)', 'WT-sfGFP', 
                                    'E2R-sfGFP',  'F18W-sfGFP',
                                    'N15D-sfGFP', 'S20M-sfGFP')
                         ), 
         Arabinose = factor(Arabinose, levels = c('0% arabinose', '0.2% arabinose'))
         )
  

## Draw a figure with the full distributions
wt_test <- all_data_processed_expression %>% filter(Strain == 'WT-sfGFP')
median_wt <- all_data_processed_expression %>% ungroup() %>%
  group_by(Arabinose) %>%
  summarise(med_wt = median(GRNBHLog))

## Statistical test (0% arabinose)
m <- aov(GRNBHLog~Strain, data=all_data_processed_expression %>% filter(Arabinose == '0% arabinose'))
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'Strain')
tukey_test2 <- TukeyHSD(m)

tukey_annotation_ara0 <- as.data.frame(tukey_test$groups) %>%
  mutate(y = rep(3.5, 6))
tukey_annotation_ara0$Strain <- rownames(tukey_annotation_ara0)

needed_vals <- rownames(tukey_test2$Strain)
tukey_table_ara0 <- tukey_test2$Strain %>% as.data.frame() %>%
  mutate(arabinose = '0% arabinose', 
         contrast = needed_vals)
tukey_table_ara0


## Statistical test (0.2% arabinose)
m <- aov(GRNBHLog~Strain, data=all_data_processed_expression %>% filter(Arabinose == '0.2% arabinose'))
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'Strain')
tukey_test2 <- TukeyHSD(m)

tukey_annotation_ara0.2 <- as.data.frame(tukey_test$groups) %>%
  mutate(y = rep(3.5, 6))
tukey_annotation_ara0.2$Strain <- rownames(tukey_annotation_ara0.2)

needed_vals <- rownames(tukey_test2$Strain)
tukey_table_ara0.2 <- tukey_test2$Strain %>% as.data.frame() %>%
  mutate(arabinose = '0.2% arabinose', 
         contrast = needed_vals)
tukey_table_ara0.2

tukey_table_final <- bind_rows(tukey_table_ara0 %>%
                                 select(arabinose, contrast, diff, lwr, upr, `p adj`), 
                               tukey_table_ara0.2 %>%
                                 select(arabinose, contrast, diff, lwr, upr, `p adj`)
                               )


## Reorganize columns
tukey_table_ara0.2 %<>% select(arabinose, contrast, diff, lwr, upr, `p adj`)

write.table(tukey_table_final, 
            file = 'Supp_tables/TableS3_Tukey_test_protAbundance.tsv', 
            row.names = F, col.names = T, sep = '\t', append = F)

## Annotation for both concentrations
tukey_annotation <- bind_rows(tukey_annotation_ara0 %>% mutate(Arabinose = '0% arabinose'), 
                              tukey_annotation_ara0.2 %>% mutate(Arabinose = '0.2% arabinose')
                              ) %>%
  mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose', '0.2% arabinose')))

## Annotation only for 0.2% arabinose
tukey_annotation <- tukey_annotation_ara0.2 %>% mutate(Arabinose = '0.2% arabinose') %>%
  mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose', '0.2% arabinose')))

p_figS3 <- all_data_processed_expression %>%
  ggplot(aes(x = Strain, y = GRNBHLog, fill = Strain)) + 
  facet_wrap(~Arabinose) +
  geom_boxplot(aes(group = Well_new), width = 0.7, alpha = 0.3, outlier.shape = NA) +
  geom_violin(aes(group = Well_new), width = 0.7, alpha = 0.3, scale = 'width') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = 'none') +
  geom_hline(data = median_wt, aes(yintercept = med_wt), linetype = 'dashed') +
  ylab('Fluorescence (log10)') +
  geom_text(data = tukey_annotation, size = 5, inherit.aes = F,
            aes(label = groups, x = as.factor(Strain), y = y))
p_figS3
ggsave(p_figS3, device = cairo_pdf, width = 16, height = 10, dpi = 300, units = 'cm', 
       filename = 'Manuscript/Figures/Supp_figures/FigS3.pdf')

## Run a separate ANOVA showing the effect of arabinose and variants ##
m <- aov(GRNBHLog~Strain+Arabinose, data=all_data_processed_expression)
anova_test <- anova(m)
write.table(anova_test, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Supp_tables/TableS3_ANOVA_protAbundance.tsv')


