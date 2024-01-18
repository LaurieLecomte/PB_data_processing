ccs.pbmm2.cov <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/coverage/pbmm2/safoPUVx_001-21.ccs.pbmm2.cov_win1Mb.regions.bed.gz", 
                            header=FALSE,
                            col.names = c('CHROM', 'POS', 'END', 'COV')
)
WIN_SIZE <- 1000000
ccs.pbmm2.cov$midBIN <- ccs.pbmm2.cov$POS + WIN_SIZE/2


library(ggplot2)

ggplot(data = ccs.pbmm2.cov) +
  facet_grid(. ~ CHROM, 
             scales = 'free', space = 'free_x') +
  geom_point(aes(x = midBIN, y = COV, color = CHROM),
             size = 0.04) +
  theme(
    # Panels and background
    panel.spacing.x = unit(0.6, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
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
    axis.line.x = element_line(linewidth = 0.08)
    
  ) + 
  
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('black', 'grey60'), 
                                  length(unique(ccs.pbmm2.cov$CHROM))/2)) +
  labs(x = 'Position along each chromosome',
       y = 'Per-window coverage',
       title = 'CCS mapped with pbmm2')


# Save to external file
#ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/combined_per_win_density.png',
#       width = 2800,
#       height = 3100,
#       units = 'px',
#       dpi = 700,
#        device = 'pdf'
#)
