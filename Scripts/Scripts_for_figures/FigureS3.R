#############################################################
####        Figure S3: Growth recovery plasmids          ####
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
# setwd('/path/to/Paralog_interference_DfrB1')
setwd('/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Angel/PhD_projects/Yacine/Manuscript/Github/Paralog_interference_DfrB1/')

## Load data
data_curves <- read_delim('Data/Growthcurves_WT_plasmids/Growth_curves_WT_plasmids.tsv', delim = '\t')

p_curves <- data_curves %>% filter(tmp == 10) %>%
  mutate(conc_ara = factor(conc_ara, 
                           levels = c('0% arabinose', '0.01% arabinose', '0.025% arabinose', 
                                      '0.05% arabinose', '0.2% arabinose'))) %>%
  ggplot(aes(x = time, y = OD)) +
  geom_line(aes(colour = sample, group = interaction(sample, REP, conc_ara))) +
  facet_wrap(~conc_ara, nrow = 1) +
  xlab('Time') + ylab('OD') +
  theme(legend.position = 'top', 
        legend.justification = 'center') +
  labs(colour = '')
p_curves

## Show the comparisons of AUC ##
p_auc <- data_curves %>% filter(tmp == 10) %>%
  mutate(conc_ara = factor(conc_ara, 
                           levels = c('0% arabinose', '0.01% arabinose', '0.025% arabinose', 
                                      '0.05% arabinose', '0.2% arabinose'))) %>%
  select(conc_ara, sample, REP, auc_e) %>% unique() %>%
  ggplot(aes(x = sample, y = auc_e, colour = sample)) +
  geom_point() +
  facet_wrap(~conc_ara, nrow = 1) +
  xlab('') + ylab('AUC') +
  stat_compare_means(paired = F, method = 't.test', 
                     comparisons = list(c('WT (Amp plasmid)', 'WT (Kan plasmid)')
                                        )
                     ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = 'none') +
  ylim(0, 20)
p_auc

p_figS3 <- plot_grid(p_curves, p_auc, nrow = 2, ncol = 1, labels = c('A', 'B'), 
                     label_size = panel_label_size, label_fontface = 'bold')

ggsave(p_figS3, device = cairo_pdf, width = 24, height = 20, dpi = 300, units = 'cm', 
       filename = 'Figures/Supp_figures/FigS3.pdf')

