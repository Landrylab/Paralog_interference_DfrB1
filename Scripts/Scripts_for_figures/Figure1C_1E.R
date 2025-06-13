#############################################################
####               Figure 1C,1E                          ####
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
setwd('/path/to/Paralog_interference_DfrB1')

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
p_fig1C <- df_final %>%
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
p_fig1C


#### Alternatively, use curves to illustrate this scenario as in panel C ####

# Set the mu_ratio to 2
mu_f = 1
mu_e = 1
mu_a = 1
mu_ratio = (mu_a + mu_e) / mu_f

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
p_fig1D_curves <- df_final_not_neutral %>%
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
p_fig1D_curves
## The rest of the figure was assembled with Inkscape
