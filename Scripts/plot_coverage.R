library(tidyverse)

cov <- read_tsv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/Coverage/B.1.1.7_HUG_P1.cov", col_names = c("Chrom", "Position", "Coverage"))
vcf <- read_tsv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/VariantAnnotation/B.1.1.7_HUG_P1.variants.annot.tab", col_names =T)


scaleFactor <- max(cov$Coverage)/max(vcf$AF)

ggplot(aes(Position, Coverage), data=cov, color="red") +
  geom_line() + geom_point(aes(POS, AF * scaleFactor), data=vcf, color="red") +
  scale_y_continuous(name="Coverage", sec.axis=sec_axis(~./scaleFactor, name="AF")) + ggtitle("HUG_P1")

####### P2

cov <- read_tsv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/Coverage/B.1.1.7_LS_P2.cov", col_names = c("Chrom", "Position", "Coverage"))
vcf <- read_tsv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/VariantAnnotation/B.1.1.7_LS_P2.variants.annot.tab", col_names =T)


scaleFactor <- max(cov$Coverage)/max(vcf$AF)

ggplot(aes(Position, Coverage), data=cov, color="red") +
  geom_line() + geom_point(aes(POS, AF * scaleFactor), data=vcf, color="red") +
  scale_y_continuous(name="Coverage", sec.axis=sec_axis(~./scaleFactor, name="AF")) + ggtitle("P2")

####### P3
cov <- read_tsv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/Coverage/B.1.1.7_LS_P3.cov", col_names = c("Chrom", "Position", "Coverage"))
vcf <- read_tsv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/VariantAnnotation/B.1.1.7_LS_P3.variants.annot.tab", col_names =T)


scaleFactor <- max(cov$Coverage)/max(vcf$AF)

ggplot(aes(Position, Coverage), data=cov, color="red") +
  geom_line() + geom_point(aes(POS, AF * scaleFactor), data=vcf, color="red") +
  scale_y_continuous(name="Coverage", sec.axis=sec_axis(~./scaleFactor, name="AF")) + ggtitle("P3")
