library(ggplot2)
WIN_SIZE <- 1000000
MIN_COV <- 1

# Disable scientific notation
options(scipen=999)


# CCS ---------------------------------------------------------------------

# Mapped with pbmm2
## Import
ccs.pbmm2.cov <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/pbmm2/safoPUVx_001-21.ccs.pbmm2.cov_win1Mb.regions.bed.gz", 
                            header=FALSE,
                            col.names = c('CHROM', 'POS', 'END', 'COV')
)
## Add mid position
ccs.pbmm2.cov$midBIN <- ccs.pbmm2.cov$POS + WIN_SIZE/2

## Remove unplaced contigs
ccs.pbmm2.cov <- subset(ccs.pbmm2.cov, grepl(x = ccs.pbmm2.cov$CHROM, pattern = 'CM'))

## Simplify chr names
ccs.pbmm2.cov$CHROM_NUM <- substr(ccs.pbmm2.cov$CHROM, 6, 8)

## Extract low cov points
ccs_pbmm2_lowcov <- subset(ccs.pbmm2.cov, COV < MIN_COV)

## Plot
ggplot(data = ccs.pbmm2.cov) +
  facet_grid(. ~ CHROM_NUM, space = 'free_x', 
            scales = 'free' ) +
  geom_point(aes(x = midBIN, y = COV, color = (as.numeric(CHROM_NUM) %% 2 == 0)),
             size = 0.02) +
  geom_point(data = ccs_pbmm2_lowcov, aes(x = midBIN, y = COV), color = 'red',
             size = 0.04)+
  geom_hline(yintercept = MIN_COV, color = 'red', linewidth = 0.1) +
  
  ylim(0, 15) +
  
  theme(
    # Panels and background
    panel.spacing.x = unit(0.1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'grey50', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 4,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    plot.title = element_text(size = 7)
  ) + 
  
 # scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  scale_color_manual(values = c('black', 'grey60')) + 
  #scale_color_manual(values = rep(c('black', 'grey60'), 
  #                                length(unique(ccs.pbmm2.cov$CHROM)))) +
  labs(x = 'Position along each chromosome',
       y = 'Per-window coverage',
       title = 'CCS mapped with pbmm2')


# Mapped with winnowmap
## Import
ccs.winnow.cov <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/winnowmap/safoPUVx_001-21.ccs.winnow.cov_win1Mb.regions.bed.gz", 
                            header=FALSE,
                            col.names = c('CHROM', 'POS', 'END', 'COV')
)
## Add mid position
ccs.winnow.cov$midBIN <- ccs.winnow.cov$POS + WIN_SIZE/2

## Remove unplaced contigs
ccs.winnow.cov <- subset(ccs.winnow.cov, grepl(x = ccs.winnow.cov$CHROM, pattern = 'CM'))

## Simplify chr names
ccs.winnow.cov$CHROM_NUM <- substr(ccs.winnow.cov$CHROM, 6, 8)

## Extract low cov points
ccs_winnow_lowcov <- subset(ccs.winnow.cov, COV < MIN_COV)

ggplot(data = ccs.winnow.cov) +
  facet_grid(. ~ CHROM_NUM, space = 'free_x', 
             scales = 'free' ) +
  geom_point(aes(x = midBIN, y = COV, color = (as.numeric(CHROM_NUM) %% 2 == 0)),
             size = 0.02) +
  geom_point(data = ccs_winnow_lowcov, aes(x = midBIN, y = COV), color = 'red',
             size = 0.04)+
  geom_hline(yintercept = MIN_COV, color = 'red', linewidth = 0.1) +
  
  ylim(0, 15) +
  
  theme(
    # Panels and background
    panel.spacing.x = unit(0.1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'grey50', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 4,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    plot.title = element_text(size = 7)
  ) + 
  
  # scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  scale_color_manual(values = c('black', 'grey60')) + 
  #scale_color_manual(values = rep(c('black', 'grey60'), 
  #                                length(unique(ccs.winnow.cov$CHROM)))) +
  labs(x = 'Position along each chromosome',
       y = 'Per-window coverage',
       title = 'CCS mapped with winnow')


# Both 
ccs.pbmm2.cov$aligner <- 'pbmm2'
ccs.winnow.cov$aligner <- 'winnowmap'

ccs_both <- rbind(ccs.pbmm2.cov, ccs.winnow.cov)

## Extract low cov points
ccs_both_lowcov <- subset(ccs_both, COV < MIN_COV)

ggplot(data = ccs_both) +
  facet_grid(aligner ~ CHROM_NUM, space = 'free_x', 
             scales = 'free' ) +
  geom_point(aes(x = midBIN, y = COV, color = (as.numeric(CHROM_NUM) %% 2 == 0)),
             size = 0.02) +
  geom_point(data = ccs_both_lowcov, aes(x = midBIN, y = COV), color = 'red',
             size = 0.04)+
  geom_hline(yintercept = MIN_COV, color = 'red', linewidth = 0.1) +
  
  ylim(0, 10) +
  
  theme(
    # Panels and background
    panel.spacing.x = unit(0.1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'grey50', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 5,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    plot.title = element_text(size = 7)
  ) + 
  
  # scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  scale_color_manual(values = c('black', 'grey60')) + 
  #scale_color_manual(values = rep(c('black', 'grey60'), 
  #                                length(unique(ccs.winnow.cov$CHROM)))) +
  labs(x = 'Position along each chromosome',
       y = paste0('Coverage by ', WIN_SIZE/1000000, '-Mb window'),
       title = 'CCS mapped')

ggsave(filename = paste0("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/ccs_pbmm2_winnow_win", WIN_SIZE/1000000, 'Mb.png'),
       width = 4000,
       height = 2800,
       units = 'px',
       dpi = 700
)

# SUBREADS ---------------------------------------------------------------------

# Mapped with pbmm2
## Import
subr.pbmm2.cov <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/pbmm2/safoPUVx_001-21.subreads.pbmm2.cov_win1Mb.regions.bed.gz", 
                            header=FALSE,
                            col.names = c('CHROM', 'POS', 'END', 'COV')
)
## Add mid position
subr.pbmm2.cov$midBIN <- subr.pbmm2.cov$POS + WIN_SIZE/2

## Remove unplaced contigs
subr.pbmm2.cov <- subset(subr.pbmm2.cov, grepl(x = subr.pbmm2.cov$CHROM, pattern = 'CM'))

## Simplify chr names
subr.pbmm2.cov$CHROM_NUM <- substr(subr.pbmm2.cov$CHROM, 6, 8)

## Extract low cov points
subr_pbmm2_lowcov <- subset(subr.pbmm2.cov, COV < MIN_COV)

## Plot
ggplot(data = subr.pbmm2.cov) +
  facet_grid(. ~ CHROM_NUM, space = 'free_x', 
             scales = 'free' ) +
  geom_point(aes(x = midBIN, y = COV, color = (as.numeric(CHROM_NUM) %% 2 == 0)),
             size = 0.02) +
  geom_point(data = subr_pbmm2_lowcov, aes(x = midBIN, y = COV), color = 'red',
             size = 0.04)+
  geom_hline(yintercept = MIN_COV, color = 'red', linewidth = 0.1) +
  
  #ylim(0, 15) +
  
  theme(
    # Panels and background
    panel.spacing.x = unit(0.1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'grey50', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 4,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    plot.title = element_text(size = 7)
  ) + 
  
  # scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  scale_color_manual(values = c('black', 'grey60')) + 
  #scale_color_manual(values = rep(c('black', 'grey60'), 
  #                                length(unique(subr.pbmm2.cov$CHROM)))) +
  labs(x = 'Position along each chromosome',
       y = 'Per-window coverage',
       title = 'Subreads mapped with pbmm2')


# Mapped with winnowmap
## Import
subr.winnow.cov <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/winnowmap/safoPUVx_001-21.subreads.winnow.cov_win1Mb.regions.bed.gz", 
                             header=FALSE,
                             col.names = c('CHROM', 'POS', 'END', 'COV')
)
## Add mid position
subr.winnow.cov$midBIN <- subr.winnow.cov$POS + WIN_SIZE/2

## Remove unplaced contigs
subr.winnow.cov <- subset(subr.winnow.cov, grepl(x = subr.winnow.cov$CHROM, pattern = 'CM'))

## Simplify chr names
subr.winnow.cov$CHROM_NUM <- substr(subr.winnow.cov$CHROM, 6, 8)

## Extract low cov points
subr_winnow_lowcov <- subset(subr.winnow.cov, COV < MIN_COV)

ggplot(data = subr.winnow.cov) +
  facet_grid(. ~ CHROM_NUM, space = 'free_x', 
             scales = 'free' ) +
  geom_point(aes(x = midBIN, y = COV, color = (as.numeric(CHROM_NUM) %% 2 == 0)),
             size = 0.02) +
  geom_point(data = subr_winnow_lowcov, aes(x = midBIN, y = COV), color = 'red',
             size = 0.04)+
  geom_hline(yintercept = MIN_COV, color = 'red', linewidth = 0.1) +
  
  #ylim(0, 15) +
  
  theme(
    # Panels and background
    panel.spacing.x = unit(0.1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'grey50', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 5,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    plot.title = element_text(size = 7)
  ) + 
  
  # scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  scale_color_manual(values = c('black', 'grey60')) + 
  #scale_color_manual(values = rep(c('black', 'grey60'), 
  #                                length(unique(subr.winnow.cov$CHROM)))) +
  labs(x = 'Position along each chromosome',
       y = 'Per-window coverage',
       title = 'Subreads mapped with winnow')


# Both 
subr.pbmm2.cov$aligner <- 'pbmm2'
subr.winnow.cov$aligner <- 'winnowmap'

subr_both <- rbind(subr.pbmm2.cov, subr.winnow.cov)

## Extract low cov points
subr_both_lowcov <- subset(subr_both, COV < MIN_COV)

ggplot(data = subr_both) +
  facet_grid(aligner ~ CHROM_NUM, space = 'free_x', 
             scales = 'free' ) +
  geom_point(aes(x = midBIN, y = COV, color = (as.numeric(CHROM_NUM) %% 2 == 0)),
             size = 0.02) +
  geom_point(data = subr_both_lowcov, aes(x = midBIN, y = COV), color = 'red',
             size = 0.04)+
  geom_hline(yintercept = MIN_COV, color = 'red', linewidth = 0.1) +
  
  ylim(0, 100) +
  
  theme(
    # Panels and background
    panel.spacing.x = unit(0.1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'grey50', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 5,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    plot.title = element_text(size = 7)
  ) + 
  
  # scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  scale_color_manual(values = c('black', 'grey60')) + 
  #scale_color_manual(values = rep(c('black', 'grey60'), 
  #                                length(unique(subr.winnow.cov$CHROM)))) +
  labs(x = 'Position along each chromosome',
       y = paste0('Coverage by ', WIN_SIZE/1000000, '-Mb window'),
       title = 'Subreads mapped')


ggsave(filename = paste0("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/subreads_pbmm2_winnow_win", WIN_SIZE/1000000, 'Mb.png'),
              width = 4000,
              height = 2800,
              units = 'px',
              dpi = 700
       )



# SUBREADS trimmed --------------------------------------------------------


