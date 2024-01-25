# Compare pbsv SVs called from winnowmap and pbmm2 alignments



# 1. Import and format ----------------------------------------------------
ccs.pbmm2 <- read.delim("pbsv_test/safoPUVx_001-21.ccs.pbmm2.table", header=FALSE,
                        col.names = c('CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'GT', 'DP'))
ccs.winnow <- read.delim("pbsv_test/safoPUVx_001-21.ccs.winnow.table", header=FALSE,
                        col.names = c('CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'GT', 'DP'))

# Convert SVLEN to absolute size, required for merging DELs later
ccs.pbmm2$SVLEN <- as.numeric(ccs.pbmm2$SVLEN)
ccs.pbmm2$absSVLEN <- as.numeric(ccs.pbmm2$SVLEN)

ccs.winnow$SVLEN <- as.numeric(ccs.winnow$SVLEN)
ccs.winnow$absSVLEN <- as.numeric(ccs.winnow$SVLEN)

ccs.pbmm2_dups <- ccs.pbmm2[duplicated(ccs.pbmm2), ]
ccs.winnow_dups <- ccs.winnow[duplicated(ccs.winnow), ]

# Remove contigs
ccs.pbmm2 <- subset(ccs.pbmm2, grepl(ccs.pbmm2$CHROM, pattern = 'CM') )
ccs.winnow <- subset(ccs.winnow, grepl(ccs.winnow$CHROM, pattern = 'CM') )


# 2. Compare depth of SV calls --------------------------------------------
ccs.pbmm2$DP <- as.numeric(ccs.pbmm2$DP)
ccs.winnow$DP <- as.numeric(ccs.winnow$DP)

ccs.pbmm2_DP1 <- subset(ccs.pbmm2, DP > 1)
ccs.pbmm2_DP2 <- subset(ccs.pbmm2, DP > 2)
ccs.pbmm2_DP4 <- subset(ccs.pbmm2, DP >= 4)

ccs.winnow_DP1 <- subset(ccs.winnow, DP > 1)
ccs.winnow_DP2 <- subset(ccs.winnow, DP > 2)
ccs.winnow_DP4 <- subset(ccs.winnow, DP >= 4)


library(ggplot2)
ccs.pbmm2$mapper <- 'pbmm2'
ccs.winnow$mapper <- 'winnowmap'

ccs_both <- rbind(ccs.pbmm2, ccs.winnow)

ggplot(data = ccs_both) + 
  #facet_wrap(~mapper + SVTYPE) +
  facet_grid(SVTYPE ~ mapper, scales = 'free') +
  geom_histogram(aes(DP), binwidth = 1) + 
  xlim(c(0, 30))


# 3. Compare types --------------------------------------------------------
table(ccs_both$SVTYPE, ccs_both$mapper)
## less INS from winnowmap


# 4. Compare SV size ------------------------------------------------------
## Convert candidate SVs length to num and bins
SVLEN_breaks <- c(-Inf, 50, 100, 250, 500, 1000, 2500, 5000, 10000, Inf)
SVLEN_names <- c('[0-50[',
                 '[50-100[',
                 '[100-250[',
                 '[250-500[',
                 '[500-1000[',
                 '[1000-2500[',
                 '[2500-5000[',
                 '[5000-10000[',
                 '[10000+')

ccs_both$SVLEN_bin <-
  cut(abs(ccs_both$SVLEN), breaks = SVLEN_breaks, labels = SVLEN_names, right = FALSE)

# SVTYPE and SVLEN before filtering
ggplot(data = ccs_both) + 
  facet_grid(SVTYPE ~ mapper, scales = 'free') +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE), color = 'black', linewidth = 0.1) +
  theme(
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1)
  )

# SVTYPE and SVLEN after filtering
ccs_both_DP4 <- subset(ccs_both, DP >= 4)
ggplot(data = ccs_both_DP4) + 
  facet_grid(SVTYPE ~ mapper, scales = 'free') +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE), color = 'black', linewidth = 0.1) +
  theme(
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1)
  )

ggplot(data = ccs_both_DP4) + 
  #facet_wrap(~mapper + SVTYPE) +
  facet_grid(SVTYPE ~ mapper, scales = 'free') +
  geom_histogram(aes(DP), binwidth = 1) + 
  xlim(c(0, 30))

table(ccs_both_DP4$SVTYPE, ccs_both_DP4$mapper)

# 2. Match successfully genotyped SVs with a known putative SV ------------
# We match based on POS, END, SVLEN and SVTYPE
# allowing a 5 bp window around each value
OFFSET <- 10

# Add window around POS, END and SVLEN of each putative SV
ccs.pbmm2$POS_min <- ccs.pbmm2$POS - OFFSET
ccs.pbmm2$POS_max <- ccs.pbmm2$POS + OFFSET

ccs.pbmm2$absSVLEN_min <- (ccs.pbmm2$absSVLEN) - OFFSET
ccs.pbmm2$absSVLEN_max <- (ccs.pbmm2$absSVLEN) + OFFSET

ccs.pbmm2$END_min <- ccs.pbmm2$END - OFFSET
ccs.pbmm2$END_max <- ccs.pbmm2$END + OFFSET


# Reconvert to data table for merging
library(data.table)

ccs.pbmm2 <- as.data.table(ccs.pbmm2)
ccs.winnow <- as.data.table(ccs.winnow)

# Merge
merged_ccs_pbmm2_winnow <- 
  ccs.winnow[ccs.pbmm2, 
       on = .(CHROM, SVTYPE, 
              POS <= POS_max, POS >= POS_min, # allow a window on POS (VARx <|>|>=|=< VARi)
              END <= END_max, END >= END_min, # allow a window on END (VARx <|>|>=|=< VARi)
              absSVLEN <= absSVLEN_max, absSVLEN >= absSVLEN_min # allow a window on SVLEN (VARx <|>|>=|=< VARi)
       ),
       .(x.CHROM, i.POS, x.POS, # variables I need in merged output (i = pbmm2) (x = winnow)
         i.SVLEN, x.SVLEN,
         i.SVTYPE, x.SVTYPE,
         i.GT, x.GT, i.DP, x.DP),]

shared_pbmm2_winnow_win5 <- subset(merged_ccs_pbmm2_winnow, !is.na(x.CHROM))

